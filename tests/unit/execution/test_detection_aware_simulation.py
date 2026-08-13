import numpy as np

from regional_code_paper.execution.detection_aware_clustering import (
    pair_counts,
    sample_catalogue,
)


def test_sampling_preserves_stratum_and_total_counts() -> None:
    strata = [
        (np.array([1, 3, 5]), 2),
        (np.array([10, 12]), 1),
    ]
    sampled = sample_catalogue(strata, np.random.default_rng(7))
    assert len(sampled) == 3
    assert len(np.unique(sampled)) == 3
    assert np.isin(sampled, [1, 3, 5, 10, 12]).all()


def test_pair_counts_are_cumulative() -> None:
    counts = pair_counts(np.array([1, 3, 10]), radii=(2, 7, 10))
    assert counts.tolist() == [1, 2, 3]
