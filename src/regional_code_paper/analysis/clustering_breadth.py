"""Quantify how broadly self-clustering is distributed across proteins.

This analysis replaces thresholded per-protein Monte Carlo and pseudo-Ripley
counts with estimands that directly answer the biological question. It uses the
exact, fixed-count opportunity expectation already calculated for each protein:

* the fraction of proteins whose observed close-pair fraction exceeds its null
  expectation;
* the fraction of matched proteins in which O-GlcNAc excess exceeds the
  corresponding phosphorylation excess; and
* the fraction satisfying those inequalities at every evaluated radius.

Intervals describe the protein population. They are not site-level intervals
and do not turn the sensitivity radii into independent hypothesis tests.
"""

from __future__ import annotations

import argparse
import math
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

from ..core.io import write_csv, write_json

PRIMARY_RADIUS = 10
PRIMARY_UNIVERSE = "all_eligible_residues"
SENSITIVITY_UNIVERSE = "observed_residue_union"
FOCAL_PTM = "oglcnac"
CONTROL_PTM = "phospho"

REQUIRED_COLUMNS = {
    "universe",
    "accession",
    "ptm",
    "radius",
    "observed_close_pair_fraction",
    "null_close_pair_fraction",
    "effect",
}


def wilson_interval(
    successes: int, total: int, confidence: float = 0.95
) -> tuple[float, float]:
    """Wilson score interval for a protein-level fraction."""
    if not 0 <= successes <= total or total <= 0:
        raise ValueError("Wilson interval requires 0 <= successes <= total and total > 0")
    z = float(stats.norm.ppf(0.5 + confidence / 2))
    proportion = successes / total
    denominator = 1 + z * z / total
    center = (proportion + z * z / (2 * total)) / denominator
    half_width = (
        z
        * math.sqrt(proportion * (1 - proportion) / total + z * z / (4 * total * total))
        / denominator
    )
    return center - half_width, center + half_width


def _fraction_record(successes: int, total: int) -> dict[str, float | int]:
    low, high = wilson_interval(successes, total)
    return {
        "n_proteins": total,
        "n_successes": successes,
        "fraction": successes / total,
        "fraction_ci_low": low,
        "fraction_ci_high": high,
    }


def per_scale_breadth(per_protein: pd.DataFrame) -> pd.DataFrame:
    """Summarize positive clustering excess separately at every radius."""
    rows: list[dict[str, object]] = []
    grouped = per_protein.groupby(["universe", "ptm", "radius"], sort=True)
    for (universe, ptm, radius), group in grouped:
        successes = int(group.effect.gt(0).sum())
        record = _fraction_record(successes, len(group))
        observed = float(group.observed_close_pair_fraction.mean())
        expected = float(group.null_close_pair_fraction.mean())
        rows.append(
            {
                "universe": universe,
                "ptm": ptm,
                "radius": int(radius),
                **record,
                "mean_observed_close_pair_fraction": observed,
                "mean_null_close_pair_fraction": expected,
                "group_fold": observed / expected,
                "median_excess": float(group.effect.median()),
                "mean_excess": float(group.effect.mean()),
            }
        )
    return pd.DataFrame(rows).sort_values(["universe", "ptm", "radius"])


def matched_dominance(per_protein: pd.DataFrame) -> pd.DataFrame:
    """Summarize O-GlcNAc dominance within proteins carrying both PTMs."""
    rows: list[dict[str, object]] = []
    for (universe, radius), group in per_protein.groupby(["universe", "radius"], sort=True):
        wide = group.pivot(index="accession", columns="ptm", values="effect")
        matched = wide.dropna(subset=[FOCAL_PTM, CONTROL_PTM])
        difference = matched[FOCAL_PTM] - matched[CONTROL_PTM]
        successes = int(difference.gt(0).sum())
        rows.append(
            {
                "universe": universe,
                "radius": int(radius),
                "contrast": "oglcnac_excess_greater_than_phospho",
                **_fraction_record(successes, len(difference)),
                "median_paired_excess_difference": float(difference.median()),
                "mean_paired_excess_difference": float(difference.mean()),
            }
        )
    return pd.DataFrame(rows).sort_values(["universe", "radius"])


def multiscale_breadth(per_protein: pd.DataFrame) -> pd.DataFrame:
    """Count proteins meeting the same directional criterion at every radius."""
    rows: list[dict[str, object]] = []
    radii = sorted(int(value) for value in per_protein.radius.unique())
    for universe, universe_frame in per_protein.groupby("universe", sort=True):
        focal = universe_frame.loc[universe_frame.ptm.eq(FOCAL_PTM)].pivot(
            index="accession", columns="radius", values="effect"
        )
        focal = focal.dropna(subset=radii)
        successes = int(focal[radii].gt(0).all(axis=1).sum())
        rows.append(
            {
                "universe": universe,
                "criterion": "positive_oglcnac_excess_at_every_radius",
                "radii": ";".join(map(str, radii)),
                **_fraction_record(successes, len(focal)),
            }
        )

        paired = universe_frame.loc[universe_frame.ptm.isin([FOCAL_PTM, CONTROL_PTM])].pivot(
            index="accession", columns=["ptm", "radius"], values="effect"
        )
        required = [(ptm, radius) for ptm in (FOCAL_PTM, CONTROL_PTM) for radius in radii]
        paired = paired.dropna(subset=required)
        dominates = np.column_stack(
            [
                paired[(FOCAL_PTM, radius)].to_numpy()
                > paired[(CONTROL_PTM, radius)].to_numpy()
                for radius in radii
            ]
        ).all(axis=1)
        rows.append(
            {
                "universe": universe,
                "criterion": "oglcnac_excess_greater_than_phospho_at_every_radius",
                "radii": ";".join(map(str, radii)),
                **_fraction_record(int(dominates.sum()), len(paired)),
            }
        )
    return pd.DataFrame(rows).sort_values(["universe", "criterion"])


def analyze(input_path: Path, output_dir: Path) -> None:
    """Build all breadth products from the canonical per-protein table."""
    per_protein = pd.read_csv(input_path)
    missing = REQUIRED_COLUMNS - set(per_protein.columns)
    if missing:
        raise ValueError(f"per-protein input is missing columns: {sorted(missing)}")
    if per_protein.duplicated(["universe", "accession", "ptm", "radius"]).any():
        raise ValueError("per-protein input contains duplicate analysis keys")

    scale = per_scale_breadth(per_protein)
    dominance = matched_dominance(per_protein)
    multiscale = multiscale_breadth(per_protein)
    primary = scale.loc[
        scale.universe.eq(PRIMARY_UNIVERSE)
        & scale.ptm.eq(FOCAL_PTM)
        & scale.radius.eq(PRIMARY_RADIUS)
    ].iloc[0]
    sensitivity = scale.loc[
        scale.universe.eq(SENSITIVITY_UNIVERSE)
        & scale.ptm.eq(FOCAL_PTM)
        & scale.radius.eq(PRIMARY_RADIUS)
    ].iloc[0]
    primary_dominance = dominance.loc[
        dominance.universe.eq(PRIMARY_UNIVERSE) & dominance.radius.eq(PRIMARY_RADIUS)
    ].iloc[0]
    all_scales = multiscale.loc[
        multiscale.universe.eq(PRIMARY_UNIVERSE)
        & multiscale.criterion.eq("positive_oglcnac_excess_at_every_radius")
    ].iloc[0]

    write_csv(output_dir / "per_scale.csv", scale)
    write_csv(output_dir / "matched_dominance.csv", dominance)
    write_csv(output_dir / "multiscale.csv", multiscale)
    write_json(
        output_dir / "summary.json",
        {
            "estimand": (
                "protein-population prevalence of positive excess over the exact "
                "fixed-count residue-opportunity expectation"
            ),
            "primary_radius": PRIMARY_RADIUS,
            "primary": primary.to_dict(),
            "observed_residue_sensitivity": sensitivity.to_dict(),
            "matched_oglcnac_dominance": primary_dominance.to_dict(),
            "all_scales": all_scales.to_dict(),
            "interpretation": (
                "Descriptive protein-population breadth; sensitivity radii are not "
                "independent tests."
            ),
        },
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--per-protein", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    analyze(args.per_protein.resolve(), args.output_dir.resolve())


if __name__ == "__main__":
    main()
