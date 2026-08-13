"""Parallel simulation and reduction for the peptide-aware clustering null."""

from __future__ import annotations

import argparse
import math
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

from ..analysis.clustering_breadth import wilson_interval
from ..analysis.detection_aware_clustering import RADII
from ..core.config import load_config
from ..core.io import write_csv, write_json
from ..core.randomness import stable_seed
from .sharding import Shard, validate_complete_partition, validate_receipts, write_receipt

PRIMARY_RADIUS = 10


def parse_positions(value: object) -> np.ndarray:
    """Parse the compact, reviewable semicolon representation used by the cache."""
    token = str(value)
    if not token:
        return np.empty(0, dtype=np.int64)
    return np.fromiter((int(item) for item in token.split(";")), dtype=np.int64)


def pair_counts(positions: np.ndarray, radii: tuple[int, ...] = RADII) -> np.ndarray:
    """Count unordered pairs within every cumulative radius."""
    positions = np.sort(np.asarray(positions, dtype=np.int64))
    if len(positions) < 2:
        return np.zeros(len(radii), dtype=np.int64)
    distances = positions[np.newaxis, :] - positions[:, np.newaxis]
    upper = np.sort(distances[np.triu_indices(len(positions), 1)])
    return np.searchsorted(upper, np.asarray(radii), side="right").astype(np.int64)


def sample_catalogue(
    strata: list[tuple[np.ndarray, int]], rng: np.random.Generator
) -> np.ndarray:
    """Draw one collision-free catalogue while preserving every stratum count."""
    sampled = [rng.choice(eligible, size=count, replace=False) for eligible, count in strata]
    output = np.sort(np.concatenate(sampled).astype(np.int64, copy=False))
    if len(np.unique(output)) != len(output):
        raise AssertionError("exposure strata overlap; simulated sites are not unique")
    return output


def simulate_shard(
    proteins_path: Path,
    strata_path: Path,
    shard_index: int,
    shard_count: int,
    draws: int,
    seed: int,
    output: Path,
    receipt: Path,
) -> None:
    """Simulate the proteins assigned to one deterministic map shard."""
    shard = Shard(shard_index, shard_count)
    proteins = pd.read_csv(proteins_path).sort_values("accession", kind="stable")
    strata_frame = pd.read_csv(strata_path)
    strata_by_protein = {
        str(accession): [
            (parse_positions(row.eligible_positions), int(row.n_observed))
            for row in group.itertuples(index=False)
        ]
        for accession, group in strata_frame.groupby("accession", sort=False)
    }

    rows: list[dict[str, object]] = []
    for protein_index, protein in enumerate(proteins.itertuples(index=False)):
        if not shard.owns(protein_index):
            continue
        accession = str(protein.accession)
        observed = parse_positions(protein.observed_positions)
        strata = strata_by_protein[accession]
        if sum(count for _, count in strata) != len(observed):
            raise AssertionError(f"site-count mismatch for {accession}")
        possible_pairs = math.comb(len(observed), 2)
        observed_counts = pair_counts(observed)
        null_counts = np.empty((draws, len(RADII)), dtype=np.int32)
        rng = np.random.default_rng(stable_seed(seed, "peptide_exposure", accession))
        for draw in range(draws):
            simulated = sample_catalogue(strata, rng)
            if len(simulated) != len(observed):
                raise AssertionError(f"simulation dropped sites for {accession}")
            null_counts[draw] = pair_counts(simulated)

        for radius_index, radius in enumerate(RADII):
            null_fraction = null_counts[:, radius_index] / possible_pairs
            observed_fraction = observed_counts[radius_index] / possible_pairs
            null_mean = float(null_fraction.mean())
            effect = observed_fraction - null_mean
            if abs(effect) < 1e-12:
                effect = 0.0
            rows.append(
                {
                    "analysis_key": f"{accession}|{radius}",
                    "accession": accession,
                    "radius": radius,
                    "n_sites": len(observed),
                    "n_exposure_strata": int(protein.n_exposure_strata),
                    "n_variable_sites": int(protein.n_variable_sites),
                    "n_mapped_peptide_sites": int(protein.n_mapped_peptide_sites),
                    "n_tryptic_fallback_sites": int(protein.n_tryptic_fallback_sites),
                    "observed_close_pairs": int(observed_counts[radius_index]),
                    "possible_pairs": possible_pairs,
                    "observed_close_pair_fraction": observed_fraction,
                    "null_mean_close_pair_fraction": null_mean,
                    "null_sd_close_pair_fraction": float(null_fraction.std(ddof=1)),
                    "effect": effect,
                    "fold": observed_fraction / null_mean if null_mean else math.inf,
                    "p_value": (
                        1
                        + int(
                            (
                                null_counts[:, radius_index] >= observed_counts[radius_index]
                            ).sum()
                        )
                    )
                    / (draws + 1),
                    "simulation_draws": draws,
                }
            )

    frame = pd.DataFrame(rows).sort_values("analysis_key", kind="stable")
    write_csv(output, frame)
    write_receipt(
        receipt,
        shard=shard,
        outputs=[output],
        records=len(frame),
        metadata={"draws_per_protein": draws, "seed": seed},
    )


def bootstrap_fold_interval(
    observed: np.ndarray,
    expected: np.ndarray,
    draws: int,
    seed: int,
) -> tuple[float, float]:
    """Whole-protein bootstrap interval for a ratio of equal-protein means."""
    rng = np.random.default_rng(seed)
    values = np.empty(draws, dtype=float)
    complete = 0
    while complete < draws:
        size = min(1_000, draws - complete)
        indices = rng.integers(0, len(observed), size=(size, len(observed)))
        values[complete : complete + size] = observed[indices].mean(axis=1) / expected[
            indices
        ].mean(axis=1)
        complete += size
    return float(np.quantile(values, 0.025)), float(np.quantile(values, 0.975))


def reduce(
    config_path: Path,
    prepared_summary_path: Path,
    shard_paths: list[Path],
    receipt_paths: list[Path],
    output_dir: Path,
) -> None:
    """Validate every shard and summarize the protein population."""
    config = load_config(config_path)
    validate_receipts(receipt_paths, len(shard_paths))
    per_protein = validate_complete_partition(
        [pd.read_csv(path) for path in shard_paths], "analysis_key"
    )
    per_protein.loc[per_protein.effect.abs() < 1e-12, "effect"] = 0.0
    if len(per_protein) % len(RADII):
        raise ValueError("detection-aware output does not cover every radius")

    base_seed = int(config.values["randomness"]["analysis_base_seed"])
    bootstrap_draws = int(config.values["randomness"]["bootstrap_draws"])
    rows: list[dict[str, object]] = []
    for radius, all_proteins in per_protein.groupby("radius", sort=True):
        # A conditioned null with zero variance cannot test spatial clustering:
        # its exposure strata force the same close-pair count in every draw.
        # Retain those proteins in the audit output, but do not dilute the
        # estimand with observations that contain no experimental information.
        informative = all_proteins.null_sd_close_pair_fraction.gt(1e-12)
        group = all_proteins.loc[informative].copy()
        if group.empty:
            raise ValueError(f"no informative conditioned nulls at radius {radius}")
        observed = group.observed_close_pair_fraction.to_numpy(float)
        expected = group.null_mean_close_pair_fraction.to_numpy(float)
        effects = group.effect.to_numpy(float)
        successes = int((effects > 1e-12).sum())
        fraction_low, fraction_high = wilson_interval(successes, len(group))
        mean_effect = float(effects.mean())
        sample_se = float(effects.std(ddof=1) / math.sqrt(len(effects)))
        critical = float(stats.t.ppf(0.975, len(effects) - 1))
        fold_low, fold_high = bootstrap_fold_interval(
            observed,
            expected,
            bootstrap_draws,
            stable_seed(base_seed, "detection_aware_fold", radius),
        )
        rows.append(
            {
                "radius": int(radius),
                "n_proteins": len(group),
                "n_all_proteins": len(all_proteins),
                "n_degenerate_null_proteins": int((~informative).sum()),
                "n_sites": int(group.loc[group.radius.eq(radius), "n_sites"].sum()),
                "n_positive_excess": successes,
                "positive_excess_fraction": successes / len(group),
                "positive_excess_fraction_ci_low": fraction_low,
                "positive_excess_fraction_ci_high": fraction_high,
                "mean_observed_close_pair_fraction": float(observed.mean()),
                "mean_null_close_pair_fraction": float(expected.mean()),
                "group_fold": float(observed.mean() / expected.mean()),
                "group_fold_ci_low": fold_low,
                "group_fold_ci_high": fold_high,
                "mean_effect": mean_effect,
                "mean_effect_ci_low": mean_effect - critical * sample_se,
                "mean_effect_ci_high": mean_effect + critical * sample_se,
                "population_t_p_greater": float(
                    stats.t.sf(mean_effect / sample_se, len(effects) - 1)
                ),
                "all_protein_group_fold": float(
                    all_proteins.observed_close_pair_fraction.mean()
                    / all_proteins.null_mean_close_pair_fraction.mean()
                ),
                "resampling_unit": "protein",
            }
        )
    per_scale = pd.DataFrame(rows).sort_values("radius")
    primary = per_scale.loc[per_scale.radius.eq(PRIMARY_RADIUS)].iloc[0]
    prepared = __import__("json").loads(prepared_summary_path.read_text(encoding="utf-8"))
    write_csv(output_dir / "per_protein.csv", per_protein.drop(columns="analysis_key"))
    write_csv(output_dir / "per_scale.csv", per_scale)
    write_json(
        output_dir / "summary.json",
        {
            "primary_radius": PRIMARY_RADIUS,
            "primary": primary.to_dict(),
            "analysis_population": prepared["analysis_population"],
            "site_count_preserved": prepared["site_count_preserved"],
            "conditioning": prepared["stratum_definition"],
            "informative_population": (
                "proteins whose conditioned null has nonzero close-pair variance"
            ),
            "interpretation": (
                "Post hoc detection-bias sensitivity analysis; the primary endpoint "
                "was fixed at 10 residues before running this implementation."
            ),
        },
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    simulate_parser = commands.add_parser("simulate")
    simulate_parser.add_argument("--proteins", type=Path, required=True)
    simulate_parser.add_argument("--strata", type=Path, required=True)
    simulate_parser.add_argument("--shard", type=int, required=True)
    simulate_parser.add_argument("--shards", type=int, required=True)
    simulate_parser.add_argument("--draws", type=int, required=True)
    simulate_parser.add_argument("--seed", type=int, required=True)
    simulate_parser.add_argument("--output", type=Path, required=True)
    simulate_parser.add_argument("--receipt", type=Path, required=True)

    reduce_parser = commands.add_parser("reduce")
    reduce_parser.add_argument("--config", type=Path, required=True)
    reduce_parser.add_argument("--prepared-summary", type=Path, required=True)
    reduce_parser.add_argument("--shards", type=Path, nargs="+", required=True)
    reduce_parser.add_argument("--receipts", type=Path, nargs="+", required=True)
    reduce_parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()

    if args.command == "simulate":
        simulate_shard(
            args.proteins,
            args.strata,
            args.shard,
            args.shards,
            args.draws,
            args.seed,
            args.output,
            args.receipt,
        )
    else:
        reduce(
            args.config,
            args.prepared_summary,
            args.shards,
            args.receipts,
            args.output_dir.resolve(),
        )


if __name__ == "__main__":
    main()
