"""Frozen protein-grouped folds, independent of CPU-specific NumPy tie sorting."""

from __future__ import annotations

import hashlib
import json
from functools import cache
from pathlib import Path

import numpy as np

FOLD_MANIFEST = Path(__file__).with_name("validation_folds.json")


@cache
def _partitions() -> dict:
    return json.loads(FOLD_MANIFEST.read_text(encoding="utf-8"))["partitions"]


def protein_group_folds(groups: np.ndarray, n_splits: int = 5):
    """Yield archived splits; reject changes to the frozen cohort or tile counts.

    GroupKFold sorts group sizes with an unstable, CPU-dispatched argsort.
    Explicit memberships preserve its validated AVX-512 partitions on all CPUs.
    Labels and features are never used to select a partition.
    """
    unique, inverse, counts = np.unique(groups, return_inverse=True, return_counts=True)
    cohort = list(zip(unique.tolist(), counts.tolist(), strict=True))
    encoded = json.dumps([n_splits, cohort], separators=(",", ":")).encode()
    key = hashlib.sha256(encoded).hexdigest()
    try:
        partition = _partitions()[key]
    except KeyError as error:
        raise ValueError(
            "No archived validation partition for this protein cohort and tile counts; "
            "a changed cohort requires an explicit new validation design."
        ) from error
    group_folds = np.asarray([row[2] for row in partition["cohort"]], dtype=int)
    if set(group_folds) != set(range(n_splits)) or len(group_folds) != len(unique):
        raise ValueError("Invalid archived protein-grouped validation partition")
    folds = group_folds[inverse]
    for fold in range(n_splits):
        yield np.flatnonzero(folds != fold), np.flatnonzero(folds == fold)
