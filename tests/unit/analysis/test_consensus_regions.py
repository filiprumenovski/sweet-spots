from regional_code_paper.analysis.consensus_regions import (
    consensus_components,
    segment_positions,
)


def test_segmentation_splits_only_above_the_gap() -> None:
    assert segment_positions([1, 4, 9, 20, 22, 24], maximum_gap=5, minimum_sites=3) == [
        (1, 4, 9),
        (20, 22, 24),
    ]


def test_consensus_requires_support_at_every_gap() -> None:
    positions = [1, 4, 9, 20, 22, 24]
    assert consensus_components(
        positions,
        gaps=(3, 5),
        final_gap=5,
        minimum_sites=3,
    ) == [(20, 22, 24)]


def test_nested_gap_intersection_equals_strictest_core_call() -> None:
    positions = [1, 4, 9, 17, 21, 25, 40, 47, 54]
    strict_core = segment_positions(positions, maximum_gap=5, minimum_sites=3)
    consensus = consensus_components(
        positions,
        gaps=(5, 8, 10, 12, 15),
        final_gap=5,
        minimum_sites=3,
    )
    assert consensus == strict_core
