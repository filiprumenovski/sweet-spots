import pandas as pd

from regional_code_paper.analysis.clustering_breadth import (
    matched_dominance,
    multiscale_breadth,
    per_scale_breadth,
    wilson_interval,
)


def example() -> pd.DataFrame:
    rows = []
    effects = {
        "A": {"oglcnac": (0.30, 0.35), "phospho": (0.10, 0.12)},
        "B": {"oglcnac": (0.20, 0.25), "phospho": (0.21, 0.26)},
        "C": {"oglcnac": (-0.10, 0.05), "phospho": (-0.20, -0.10)},
    }
    for accession, ptms in effects.items():
        for ptm, values in ptms.items():
            for radius, effect in zip((5, 10), values, strict=True):
                rows.append(
                    {
                        "universe": "all_eligible_residues",
                        "accession": accession,
                        "ptm": ptm,
                        "radius": radius,
                        "observed_close_pair_fraction": 0.5 + effect,
                        "null_close_pair_fraction": 0.5,
                        "effect": effect,
                    }
                )
    return pd.DataFrame(rows)


def test_wilson_interval_contains_observed_fraction() -> None:
    low, high = wilson_interval(8, 10)
    assert low < 0.8 < high


def test_breadth_and_dominance_use_proteins_as_rows() -> None:
    scale = per_scale_breadth(example())
    primary = scale.loc[(scale.ptm == "oglcnac") & (scale.radius == 10)].iloc[0]
    assert primary.n_proteins == 3
    assert primary.n_successes == 3

    dominance = matched_dominance(example())
    radius_10 = dominance.loc[dominance.radius == 10].iloc[0]
    assert radius_10.n_proteins == 3
    assert radius_10.n_successes == 2


def test_multiscale_requires_every_radius() -> None:
    result = multiscale_breadth(example())
    focal = result.loc[result.criterion == "positive_oglcnac_excess_at_every_radius"].iloc[0]
    assert focal.n_successes == 2
