import json

import numpy as np
import pytest

from regional_code_paper.core.validation import FOLD_MANIFEST, protein_group_folds


def test_archived_scanner_folds_preserve_proteins_and_survive_row_reordering():
    partitions = json.loads(FOLD_MANIFEST.read_text())["partitions"]
    cohort = next(p["cohort"] for p in partitions.values() if "scanner.outer" in p["names"])
    groups = np.repeat([row[0] for row in cohort], [row[1] for row in cohort])
    expected = np.repeat([row[2] for row in cohort], [row[1] for row in cohort])
    order = np.random.default_rng(17).permutation(len(groups))
    groups, expected = groups[order], expected[order]
    seen = np.zeros(len(groups), dtype=int)
    for fold, (train, test) in enumerate(protein_group_folds(groups)):
        assert not set(groups[train]) & set(groups[test])
        assert np.all(expected[test] == fold)
        seen[test] += 1
    assert np.all(seen == 1)
    known = {protein: fold for protein, _, fold in cohort}
    assert {p: known[p] for p in ["Q8WZ42", "Q8WXI7", "Q8NF91", "Q5VST9", "Q03001"]} == {
        "Q8WZ42": 0,
        "Q8WXI7": 1,
        "Q8NF91": 2,
        "Q5VST9": 3,
        "Q03001": 4,
    }
    with pytest.raises(ValueError, match="changed cohort"):
        list(protein_group_folds(groups[:-1]))


def test_unknown_cohort_is_not_silently_resplit():
    with pytest.raises(ValueError, match="No archived validation partition"):
        list(protein_group_folds(np.array(["new-protein"])))


def test_nested_scanner_partitions_exclude_each_outer_test_set():
    partitions = list(json.loads(FOLD_MANIFEST.read_text())["partitions"].values())
    outer = next(p["cohort"] for p in partitions if "scanner.outer" in p["names"])
    for fold in range(5):
        inner = next(p for p in partitions if f"scanner.inner.{fold}" in p["names"])
        assert inner["n_splits"] == 3
        assert {(p, count) for p, count, _ in inner["cohort"]} == {
            (p, count) for p, count, assigned in outer if assigned != fold
        }
