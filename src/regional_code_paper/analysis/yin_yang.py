"""Opportunity-conditioned O-GlcNAc/phosphorylation comparisons.

The historical yin-yang question is represented by two distinct estimands:

* a local same-site test, conditioned on nearby non-O-GlcNAc S/T opportunity;
* a nested protein-level null showing how the raw overlap changes after
  conditioning on residue chemistry and disorder.

The module also implements the matched kinome comparison in Figure 3C/S2.
"""

from __future__ import annotations

import argparse
from collections import defaultdict
from pathlib import Path

import kinase_library as kl
import numpy as np
import pandas as pd

from ..core.config import load_config
from ..core.io import write_json
from ..core.randomness import stable_seed
from .regional_code import Inputs, load_inputs, local_covariates


def local_opportunity_table(inputs: Inputs) -> pd.DataFrame:
    rows = []
    for accession in sorted(set(inputs.oglcnac) & set(inputs.phospho)):
        sequence = inputs.sequences[accession]
        st_positions = np.asarray(
            [index for index, residue in enumerate(sequence, start=1) if residue in "ST"],
            dtype=int,
        )
        if len(st_positions) < 2:
            continue
        residues = np.asarray([sequence[position - 1] for position in st_positions])
        is_oglcnac = np.asarray(
            [position in inputs.oglcnac[accession] for position in st_positions]
        )
        is_phospho = np.asarray(
            [position in inputs.phospho[accession] for position in st_positions]
        )
        for window in (10, 25, 50):
            left = np.searchsorted(st_positions, st_positions - window, side="left")
            right = np.searchsorted(st_positions, st_positions + window, side="right")
            for site_index in np.flatnonzero(is_oglcnac):
                neighbourhood = slice(left[site_index], right[site_index])
                non_oglcnac = ~is_oglcnac[neighbourhood]
                for residue_matched in (True, False):
                    controls = (
                        non_oglcnac & (residues[neighbourhood] == residues[site_index])
                        if residue_matched
                        else non_oglcnac
                    )
                    n_controls = int(controls.sum())
                    if n_controls == 0:
                        continue
                    position = int(st_positions[site_index])
                    rows.append(
                        {
                            "accession": accession,
                            "position": position,
                            "window": window,
                            "residue_matched": residue_matched,
                            "observed": int(is_phospho[site_index]),
                            "expected": float(is_phospho[neighbourhood][controls].mean()),
                            "n_controls": n_controls,
                            "regional": (accession, position) in inputs.region_sites,
                        }
                    )
    return pd.DataFrame(rows)


def ratio_interval(
    frame: pd.DataFrame, *, expected: str, draws: int, seed: int
) -> dict[str, object]:
    aggregate = frame.groupby("accession")[["observed", expected]].sum()
    observed = aggregate.observed.to_numpy(float)
    null = aggregate[expected].to_numpy(float)
    if len(observed) < 3 or null.sum() == 0:
        raise ValueError("Ratio requires at least three proteins and a positive expectation")
    rng = np.random.default_rng(seed)
    estimates = np.empty(draws, dtype=float)
    for draw in range(draws):
        indices = rng.integers(0, len(observed), len(observed))
        estimates[draw] = observed[indices].sum() / null[indices].sum()
    low, high = np.quantile(estimates, [0.025, 0.975])
    return {
        "ratio": float(observed.sum() / null.sum()),
        "ci": [float(low), float(high)],
        "n_proteins": len(observed),
        "observed": int(observed.sum()),
        "expected": float(null.sum()),
    }


def nested_null_table(inputs: Inputs) -> tuple[pd.DataFrame, np.ndarray]:
    scope = sorted(set(inputs.oglcnac) & set(inputs.phospho))
    disorder_values = [
        inputs.disorder[(accession, position)]
        for accession in scope
        for position, residue in enumerate(inputs.sequences[accession], start=1)
        if residue in "ST" and (accession, position) in inputs.disorder
    ]
    cutpoints = np.quantile(np.asarray(disorder_values), [0.2, 0.4, 0.6, 0.8])

    def stratum(value: float) -> int:
        return int(np.searchsorted(cutpoints, value, side="right"))

    rows = []
    for accession in scope:
        sequence = inputs.sequences[accession]
        pool = [
            position
            for position, residue in enumerate(sequence, start=1)
            if residue in "ST" and (accession, position) in inputs.disorder
        ]
        if not pool:
            continue
        pool_set = set(pool)
        phospho = sorted(inputs.phospho[accession] & pool_set)
        if not phospho:
            continue
        pool_strata = np.asarray([stratum(inputs.disorder[(accession, p)]) for p in pool])
        phospho_strata = np.asarray([stratum(inputs.disorder[(accession, p)]) for p in phospho])
        definitions = (
            ("all", inputs.oglcnac[accession]),
            (
                "regional",
                {
                    position
                    for position in inputs.oglcnac[accession]
                    if (accession, position) in inputs.region_sites
                },
            ),
        )
        for arm, positions in definitions:
            oglcnac = sorted(positions & pool_set)
            if not oglcnac:
                continue
            oglcnac_strata = np.asarray(
                [stratum(inputs.disorder[(accession, p)]) for p in oglcnac]
            )
            expected_disorder = 0.0
            for bin_index in range(5):
                pool_size = int((pool_strata == bin_index).sum())
                if pool_size:
                    expected_disorder += (
                        float((phospho_strata == bin_index).sum())
                        * float((oglcnac_strata == bin_index).sum())
                        / pool_size
                    )
            rows.append(
                {
                    "accession": accession,
                    "arm": arm,
                    "observed": len(set(oglcnac) & set(phospho)),
                    "expected_length": len(phospho) * len(oglcnac) / len(sequence),
                    "expected_st": len(phospho) * len(oglcnac) / len(pool),
                    "expected_disorder": expected_disorder,
                    "n_oglcnac": len(oglcnac),
                    "n_phospho": len(phospho),
                    "n_st": len(pool),
                    "length": len(sequence),
                }
            )
    return pd.DataFrame(rows), cutpoints


def peptide(sequence: str, position: int, radius: int = 7) -> str:
    """Return a fixed-width kinase-library peptide, padded with underscores."""
    return "".join(
        sequence[index] if 0 <= index < len(sequence) else "_"
        for index in range(position - 1 - radius, position + radius)
    )


def matched_kinase_pairs(inputs: Inputs) -> list[tuple[str, int, int, str, str]]:
    pairs = []
    regional: dict[str, set[int]] = defaultdict(set)
    for accession, position in inputs.region_sites:
        regional[accession].add(position)
    for accession in sorted(regional):
        sequence = inputs.sequences.get(accession)
        if not sequence:
            continue
        controls = [
            position
            for position, residue in enumerate(sequence, start=1)
            if residue in "ST"
            and position not in inputs.oglcnac.get(accession, set())
            and position not in inputs.phospho.get(accession, set())
        ]
        if not controls:
            continue
        control_features = np.asarray(
            [local_covariates(inputs, accession, position, radius=12) for position in controls]
        )
        # Match S/T density over the manuscript's 51-residue window, not ±12.
        control_features[:, 1] = [
            local_covariates(inputs, accession, position, radius=25)[1] for position in controls
        ]
        scale = np.nanstd(control_features, axis=0) + 1e-9
        used = np.zeros(len(controls), dtype=bool)
        for position in sorted(regional[accession]):
            feature = np.asarray(
                [
                    local_covariates(inputs, accession, position, radius=12)[0],
                    local_covariates(inputs, accession, position, radius=25)[1],
                ]
            )
            if not np.isfinite(feature).all():
                continue
            distance = np.nansum(((control_features - feature) / scale) ** 2, axis=1)
            distance[used] = np.inf
            match = int(np.argmin(distance))
            if not np.isfinite(distance[match]):
                continue
            used[match] = True
            control = controls[match]
            pairs.append(
                (
                    accession,
                    position,
                    control,
                    peptide(sequence, position),
                    peptide(sequence, control),
                )
            )
    return pairs


def score_kinome(
    pairs: list[tuple[str, int, int, str, str]],
) -> tuple[pd.DataFrame, pd.DataFrame]:
    information = kl.get_kinome_info()
    serine_threonine = set(information.loc[information.KL_LIBRARY.eq("ser_thr"), "MATRIX_NAME"])
    family = dict(zip(information.MATRIX_NAME, information.FAMILY, strict=True))
    site_rows = []
    kinase_cluster: pd.Series[float] | None = None
    kinase_control: pd.Series[float] | None = None
    n_scored = 0
    for (
        accession,
        cluster_position,
        control_position,
        cluster_peptide,
        control_peptide,
    ) in pairs:
        try:
            cluster = kl.Substrate(cluster_peptide).percentile()
            control = kl.Substrate(control_peptide).percentile()
        except Exception:
            continue
        kinases = [
            kinase
            for kinase in cluster.index
            if kinase in serine_threonine and kinase in control.index
        ]
        cluster = cluster[kinases].astype(float)
        control = control[kinases].astype(float)
        kinase_cluster = (
            cluster.copy()
            if kinase_cluster is None
            else kinase_cluster.add(cluster, fill_value=0.0)
        )
        kinase_control = (
            control.copy()
            if kinase_control is None
            else kinase_control.add(control, fill_value=0.0)
        )
        site_rows.append(
            {
                "accession": accession,
                "cluster_position": cluster_position,
                "control_position": control_position,
                "cluster_n90": int((cluster > 90).sum()),
                "control_n90": int((control > 90).sum()),
                "cluster_max": float(cluster.max()),
                "control_max": float(control.max()),
            }
        )
        n_scored += 1
    if kinase_cluster is None or kinase_control is None:
        raise RuntimeError("No kinase peptides could be scored")
    kinases = pd.DataFrame(
        {
            "kinase": kinase_cluster.index,
            "cluster_mean_percentile": (kinase_cluster / n_scored).to_numpy(),
            "control_mean_percentile": (kinase_control / n_scored).to_numpy(),
        }
    )
    kinases["delta"] = kinases.cluster_mean_percentile - kinases.control_mean_percentile
    kinases["family"] = [family.get(name, "Other") for name in kinases.kinase]
    kinases["n_pairs"] = n_scored
    return pd.DataFrame(site_rows), kinases.sort_values("delta").reset_index(drop=True)


def kinase_protein_summary(sites: pd.DataFrame) -> pd.DataFrame:
    frame = sites.groupby("accession", as_index=False).agg(
        cluster_n90=("cluster_n90", "mean"),
        control_n90=("control_n90", "mean"),
        cluster_max=("cluster_max", "mean"),
        control_max=("control_max", "mean"),
    )
    frame["delta_n90"] = frame.cluster_n90 - frame.control_n90
    frame["delta_max"] = frame.cluster_max - frame.control_max
    return frame


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--regions", type=Path, required=True)
    parser.add_argument("--sites", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    config = load_config(args.config)
    output = args.output_dir.resolve()
    output.mkdir(parents=True, exist_ok=True)
    base_seed = int(config.values["randomness"]["manuscript_base_seed"])
    draws = int(config.values["randomness"]["figure_bootstrap_draws"])
    inputs = load_inputs(args.config, args.regions, args.sites)

    local = local_opportunity_table(inputs)
    nested, cutpoints = nested_null_table(inputs)
    kinase_sites, kinase_deltas = score_kinome(matched_kinase_pairs(inputs))
    kinase_proteins = kinase_protein_summary(kinase_sites)
    local.to_csv(output / "local_opportunity.csv", index=False)
    nested.to_csv(output / "nested_nulls.csv", index=False)
    kinase_deltas.to_csv(output / "kinase_site_scores.csv", index=False)
    kinase_proteins.to_csv(output / "kinase_summary.csv", index=False)

    local_results = {}
    for arm, selected in (("all", local), ("regional", local.loc[local.regional])):
        for window in (10, 25, 50):
            for residue_matched in (True, False):
                subset = selected.loc[
                    (selected.window == window) & (selected.residue_matched == residue_matched)
                ]
                key = f"{arm}_w{window}_{'residue_matched' if residue_matched else 'all_st'}"
                local_results[key] = ratio_interval(
                    subset,
                    expected="expected",
                    draws=draws,
                    seed=stable_seed(base_seed, "local", arm, window, residue_matched),
                )
    nested_results = {}
    for arm in ("all", "regional"):
        selected = nested.loc[nested.arm.eq(arm)]
        nested_results[arm] = {
            name: ratio_interval(
                selected,
                expected=column,
                draws=draws,
                seed=stable_seed(base_seed, "nested", arm, name),
            )
            for name, column in (
                ("N0_length", "expected_length"),
                ("N1_st", "expected_st"),
                ("N2_disorder", "expected_disorder"),
            )
        }
    delta = kinase_proteins.delta_n90.to_numpy(float)
    rng = np.random.default_rng(stable_seed(base_seed, "kinase", "bootstrap"))
    bootstrap = np.empty(20_000, dtype=float)
    for draw in range(len(bootstrap)):
        bootstrap[draw] = delta[rng.integers(0, len(delta), len(delta))].mean()
    write_json(
        output / "summary.json",
        {
            "local_opportunity": local_results,
            "nested_nulls": nested_results,
            "disorder_cutpoints": cutpoints.tolist(),
            "kinome": {
                "n_site_pairs": len(kinase_sites),
                "n_proteins": len(kinase_proteins),
                "n_kinases": len(kinase_deltas),
                "cluster_n90": float(kinase_proteins.cluster_n90.mean()),
                "control_n90": float(kinase_proteins.control_n90.mean()),
                "delta_n90": float(delta.mean()),
                "delta_n90_ci": np.quantile(bootstrap, [0.025, 0.975]).tolist(),
            },
        },
    )


if __name__ == "__main__":
    main()
