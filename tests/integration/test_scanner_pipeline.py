"""Synthetic prepare → nested folds → receipts → reduction, without provider data."""

import json

import numpy as np
import pandas as pd
import pytest
from sklearn.ensemble import HistGradientBoostingClassifier
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import GroupKFold
from threadpoolctl import threadpool_limits

from regional_code_paper.execution import scanner as execution
from regional_code_paper.models import scanner


def test_synthetic_scanner_serial_and_sharded_agree(tmp_path, monkeypatch):
    rng = np.random.default_rng(7)
    groups = np.repeat([f"P{i:02}" for i in range(20)], 12)
    labels = np.tile([0, 1], len(groups) // 2)
    legacy = rng.normal(size=(len(groups), 3)) + labels[:, None]
    enhanced = np.column_stack([legacy, rng.normal(size=(len(groups), 3))])
    centers = np.tile(np.arange(12) * 10 + 11, 20)
    spans = {group: [(21, 51)] for group in np.unique(groups)}
    families = {"composition": [0, 1], "geometry": [2, 3], "disorder": [4, 5]}

    def splits(groups, n_splits=5):
        for train, test in GroupKFold(n_splits).split(groups, groups=groups):
            assert not set(groups[train]) & set(groups[test])
            yield train, test

    def estimator():
        return HistGradientBoostingClassifier(max_iter=4, min_samples_leaf=3, random_state=42)

    # Synthetic cohorts have their own test-only split design, never a production fallback.
    monkeypatch.setattr(scanner, "protein_group_folds", splits)
    monkeypatch.setattr(execution, "protein_group_folds", splits)
    monkeypatch.setattr(scanner, "baseline_estimator", estimator)
    monkeypatch.setattr(scanner, "make_region_classifier", estimator)
    monkeypatch.setattr(execution, "make_region_classifier", estimator)
    monkeypatch.setattr(
        execution,
        "load_dataset",
        lambda *_: (legacy, enhanced, labels, groups, groups, centers, spans, families),
    )
    cache, metadata = tmp_path / "features.npz", tmp_path / "features.json"
    execution.prepare(tmp_path / "config.yaml", tmp_path / "regions.csv", cache, metadata)
    shards, receipts = [], []
    with threadpool_limits(limits=1):
        scores, folds, calls, fold_table = scanner.evaluate_nested(
            legacy, enhanced, labels, groups
        )
        for fold in range(5):
            shards.append(tmp_path / f"fold{fold}.csv")
            receipts.append(tmp_path / f"fold{fold}.json")
            execution.run_fold(cache, metadata, fold, shards[-1], receipts[-1])
        # Held-out labels cannot affect selection, calibrated cutoffs or predictions.
        changed = labels.copy()
        _, test = next(iter(splits(groups)))
        changed[test] = 1 - changed[test]
        before, info = scanner.fit_outer_fold(legacy, enhanced, labels, groups, 0)
        after, changed_info = scanner.fit_outer_fold(legacy, enhanced, changed, groups, 0)
        pd.testing.assert_frame_equal(before, after)
        assert info["thresholds"] == changed_info["thresholds"]
        assert info["selected"] == changed_info["selected"]
    output = tmp_path / "results"
    execution.reduce_folds(cache, metadata, shards[::-1], receipts[::-1], output)
    predictions = pd.read_csv(output / "outer_predictions.csv")
    np.testing.assert_array_equal(predictions.fold, folds)
    for name in scores:
        np.testing.assert_allclose(predictions[f"{name}_score"], scores[name], atol=1e-15)
    distributed_folds = pd.read_csv(output / "outer_folds.csv")
    identity = ["fold", "n_test", "n_test_positive", "selected", "inner_average_precision"]
    pd.testing.assert_frame_equal(
        distributed_folds[identity], fold_table[identity], check_dtype=False
    )
    # Metrics use the published CSV boundary; decimal parsing can split floating-point ties.
    for fold in range(5):
        subset = predictions.loc[predictions.fold == fold]
        assert distributed_folds.loc[fold, "baseline_auroc"] == pytest.approx(
            roc_auc_score(subset.label, subset.baseline_score)
        )
    assert json.loads((output / "summary.json").read_text())["cohort"]["tiles"] == len(labels)
    with pytest.raises(ValueError, match="five outer folds"):
        execution.reduce_folds(cache, metadata, shards[:-1], receipts[:-1], output)
