import pandas as pd

from regional_code_paper.analysis.region_definition_sensitivity import (
    construct_definition,
    definition_id,
    interval_overlap_fraction,
)


def test_definition_identifier_is_explicit() -> None:
    assert definition_id(5, 3, 10) == "core5_min3_final10"


def test_final_gap_merges_strict_core_chains() -> None:
    regions, sites = construct_definition(
        {"P1": {100, 103, 108, 116, 120, 124}},
        core_gap=5,
        minimum_sites=3,
        final_gap=10,
    )
    assert regions[["start", "end", "span", "valence"]].to_dict("records") == [
        {"start": 100, "end": 124, "span": 25, "valence": 6}
    ]
    assert sites.position.tolist() == [100, 103, 108, 116, 120, 124]


def test_interval_overlap_counts_any_coordinate_intersection() -> None:
    primary = pd.DataFrame(
        {"accession": ["P1", "P2"], "start": [10, 10], "end": [20, 20]}
    )
    variant = pd.DataFrame({"accession": ["P1"], "start": [20], "end": [25]})
    assert interval_overlap_fraction(primary, variant) == 0.5
