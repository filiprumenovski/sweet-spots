"""Self-clustering after restricting sites to independently corroborated evidence.

The analysis asks whether the main result persists when a site must have been
reported by at least two publications or on at least two distinct, exactly
mapped peptide sequences. These restrictions trade coverage for evidential
independence and are reported as post hoc robustness analyses.
"""

from __future__ import annotations

import argparse
import math
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd

from ..core.config import load_config
from ..core.io import write_csv, write_json
from ..core.randomness import stable_seed
from .clustering_breadth import wilson_interval
from .detection_aware_clustering import RADII, load_sequences, mapped_spans
from .self_clustering import all_eligible, close_pair_graph

PRIMARY_RADIUS = 10
MINIMUM_SITES = 3


def evidence_sets(
    atlas_path: Path, sequences: dict[str, str]
) -> tuple[dict[str, dict[str, np.ndarray]], dict[str, int]]:
    """Construct strict, multi-publication and multi-peptide site maps."""
    atlas = pd.read_csv(
        atlas_path,
        usecols=[
            "species",
            "accession",
            "position_in_protein",
            "site_residue",
            "peptide_seq",
            "pmid",
        ],
        low_memory=False,
        encoding="utf-8",
        encoding_errors="replace",
    )
    atlas = atlas.loc[atlas.species.eq("human") & atlas.site_residue.isin(["S", "T"])].copy()
    atlas["position"] = pd.to_numeric(atlas.position_in_protein, errors="coerce")
    atlas = atlas.dropna(subset=["position"])

    publications: dict[tuple[str, int], set[str]] = defaultdict(set)
    peptides: dict[tuple[str, int], set[str]] = defaultdict(set)
    valid_sites: set[tuple[str, int]] = set()
    invalid_rows = 0
    for accession, raw_position, residue, peptide, pmid in atlas[
        ["accession", "position", "site_residue", "peptide_seq", "pmid"]
    ].itertuples(index=False, name=None):
        accession = str(accession)
        position = int(raw_position)
        sequence = sequences.get(accession)
        if (
            sequence is None
            or not 1 <= position <= len(sequence)
            or sequence[position - 1] != residue
        ):
            invalid_rows += 1
            continue
        key = accession, position
        valid_sites.add(key)
        if pmid is not None and not pd.isna(pmid) and str(pmid).strip():
            publications[key].add(str(pmid).strip())
        if mapped_spans(sequence, peptide, position):
            peptides[key].add(str(peptide).strip().upper())

    definitions = {
        "all_strict_sites": valid_sites,
        "two_independent_publications": {
            key for key in valid_sites if len(publications[key]) >= 2
        },
        "two_distinct_mapped_peptides": {key for key in valid_sites if len(peptides[key]) >= 2},
    }
    output: dict[str, dict[str, np.ndarray]] = {}
    for name, sites in definitions.items():
        by_protein: dict[str, list[int]] = defaultdict(list)
        for accession, position in sorted(sites):
            by_protein[accession].append(position)
        output[name] = {
            accession: np.asarray(positions, dtype=np.int64)
            for accession, positions in by_protein.items()
        }
    return output, {
        "valid_unique_sites": len(valid_sites),
        "invalid_rows": invalid_rows,
        "multi_publication_sites": len(definitions["two_independent_publications"]),
        "multi_peptide_sites": len(definitions["two_distinct_mapped_peptides"]),
    }


def protein_table(
    sequences: dict[str, str], site_maps: dict[str, dict[str, np.ndarray]]
) -> pd.DataFrame:
    """Calculate the exact residue-opportunity expectation for each restriction."""
    rows: list[dict[str, object]] = []
    for restriction, by_protein in site_maps.items():
        for accession, observed in sorted(by_protein.items()):
            if len(observed) < MINIMUM_SITES:
                continue
            eligible = all_eligible(sequences[accession], "oglcnac")
            possible_pairs = math.comb(len(observed), 2)
            eligible_pairs = math.comb(len(eligible), 2)
            for radius in RADII:
                observed_edges, _ = close_pair_graph(observed, radius)
                eligible_edges, _ = close_pair_graph(eligible, radius)
                observed_fraction = observed_edges / possible_pairs
                null_fraction = eligible_edges / eligible_pairs
                rows.append(
                    {
                        "restriction": restriction,
                        "accession": accession,
                        "radius": radius,
                        "n_sites": len(observed),
                        "n_eligible": len(eligible),
                        "observed_close_pair_fraction": observed_fraction,
                        "null_close_pair_fraction": null_fraction,
                        "effect": observed_fraction - null_fraction,
                    }
                )
    return pd.DataFrame(rows).sort_values(["restriction", "accession", "radius"])


def fold_bootstrap(
    observed: np.ndarray,
    expected: np.ndarray,
    draws: int,
    seed: int,
) -> tuple[float, float]:
    """Whole-protein bootstrap interval for the equal-protein group fold."""
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


def summarize(per_protein: pd.DataFrame, draws: int, seed: int) -> pd.DataFrame:
    """Summarize every evidence restriction and radius at protein level."""
    rows: list[dict[str, object]] = []
    for (restriction, radius), group in per_protein.groupby(
        ["restriction", "radius"], sort=True
    ):
        observed = group.observed_close_pair_fraction.to_numpy(float)
        expected = group.null_close_pair_fraction.to_numpy(float)
        successes = int(group.effect.gt(0).sum())
        prevalence_low, prevalence_high = wilson_interval(successes, len(group))
        fold_low, fold_high = fold_bootstrap(
            observed,
            expected,
            draws,
            stable_seed(seed, "evidence_restriction", restriction, radius),
        )
        rows.append(
            {
                "restriction": restriction,
                "radius": int(radius),
                "n_proteins": len(group),
                "n_sites": int(group.n_sites.sum()),
                "n_positive_excess": successes,
                "positive_excess_fraction": successes / len(group),
                "positive_excess_fraction_ci_low": prevalence_low,
                "positive_excess_fraction_ci_high": prevalence_high,
                "mean_observed_close_pair_fraction": float(observed.mean()),
                "mean_null_close_pair_fraction": float(expected.mean()),
                "group_fold": float(observed.mean() / expected.mean()),
                "group_fold_ci_low": fold_low,
                "group_fold_ci_high": fold_high,
                "mean_effect": float((observed - expected).mean()),
                "resampling_unit": "protein",
            }
        )
    return pd.DataFrame(rows).sort_values(["restriction", "radius"])


def analyze(config_path: Path, output_dir: Path) -> None:
    """Run the complete independent-evidence robustness analysis."""
    config = load_config(config_path)
    sequences = load_sequences(config.source_root / "data/interim/fasta_human.parquet")
    maps, audit = evidence_sets(
        config.source_root / "analysis/revalidation/data/atlas_unambiguous.csv",
        sequences,
    )
    proteins = protein_table(sequences, maps)
    draws = int(config.values["randomness"]["bootstrap_draws"])
    base_seed = int(config.values["randomness"]["manuscript_base_seed"])
    scales = summarize(proteins, draws, base_seed)
    primary = scales.loc[scales.radius.eq(PRIMARY_RADIUS)].set_index("restriction")
    write_csv(output_dir / "per_protein.csv", proteins)
    write_csv(output_dir / "per_scale.csv", scales)
    write_json(
        output_dir / "summary.json",
        {
            "primary_radius": PRIMARY_RADIUS,
            "audit": audit,
            "all_strict_sites": primary.loc["all_strict_sites"].to_dict(),
            "two_independent_publications": primary.loc[
                "two_independent_publications"
            ].to_dict(),
            "two_distinct_mapped_peptides": primary.loc[
                "two_distinct_mapped_peptides"
            ].to_dict(),
            "interpretation": (
                "Post hoc independent-evidence restrictions; fewer sites and proteins "
                "change the estimand and preclude direct fold comparisons as endpoints."
            ),
        },
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    analyze(args.config, args.output_dir.resolve())


if __name__ == "__main__":
    main()
