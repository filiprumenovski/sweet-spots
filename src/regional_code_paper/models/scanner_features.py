"""Leakage-resistant feature extraction for the O-GlcNAc region scanner.

This module deliberately contains no O-GlcNAc coordinates.  Feature extraction
accepts sequence, center, and (optionally) independently predicted disorder only.
The feature-family slices keep the scientific decomposition explicit:

* ``composition``: non-S/T chemistry, normalized within the non-S/T residues;
* ``geometry``: S/T density and S/T-mask arrangement only;
* ``disorder``: predicted disorder only;
* ``full``: all of the above plus reduced-alphabet arrangement features.
"""

from __future__ import annotations

import math
from collections import Counter
from collections.abc import Mapping
from dataclasses import dataclass

import numpy as np
from sklearn.ensemble import HistGradientBoostingClassifier

AMINO_ACIDS = "ACDEFGHIKLMNPQRSTVWY"
NON_ST_AMINO_ACIDS = tuple(aa for aa in AMINO_ACIDS if aa not in "ST")
WINDOW_RADII = (10, 20, 40)
KD = dict(
    zip(
        "ARNDCQEGHILKMFPSTWYV",
        [
            1.8,
            -4.5,
            -3.5,
            -3.5,
            2.5,
            -3.5,
            -3.5,
            -0.4,
            -3.2,
            4.5,
            3.8,
            -3.9,
            1.9,
            2.8,
            -1.6,
            -0.8,
            -0.7,
            -0.9,
            -1.3,
            4.2,
        ],
        strict=True,
    )
)
REDUCED_CLASS = {}
for _group, _residues in {
    "H": "AVLIM",
    "R": "FWY",
    "P": "P",
    "G": "G",
    "O": "ST",
    "+": "KRH",
    "-": "DE",
    "N": "NQC",
}.items():
    for _aa in _residues:
        REDUCED_CLASS[_aa] = _group
REDUCED_CLASSES = "HRPGO+-N"
REDUCED_PAIRS = tuple(a + b for a in REDUCED_CLASSES for b in REDUCED_CLASSES)


@dataclass(frozen=True)
class FeatureVector:
    values: np.ndarray
    names: tuple[str, ...]
    families: Mapping[str, tuple[int, ...]]


def make_region_classifier() -> HistGradientBoostingClassifier:
    """Return the validated full regional classifier.

    Average-precision early stopping is intentional for the ~1%-prevalence
    task.  Class weighting is omitted: it did not improve inner grouped-CV AP
    consistently and produces scores on a different scale.
    """

    return HistGradientBoostingClassifier(
        max_iter=600,
        max_depth=4,
        learning_rate=0.03,
        l2_regularization=1.0,
        early_stopping=True,
        scoring="average_precision",
        validation_fraction=0.1,
        n_iter_no_change=30,
        random_state=42,
    )


def _window(sequence: str, center: int, radius: int) -> str:
    """Return the historical 1-indexed-center window used by P33."""

    return sequence[max(0, center - 1 - radius) : min(len(sequence), center + radius)]


def _scd(segment: str) -> float:
    charge = np.array([1 if aa in "KR" else (-1 if aa in "DE" else 0) for aa in segment])
    charged = np.flatnonzero(charge)
    if len(charged) < 2:
        return 0.0
    total = 0.0
    for left in range(len(charged)):
        for right in range(left + 1, len(charged)):
            i, j = charged[left], charged[right]
            total += charge[i] * charge[j] * math.sqrt(j - i)
    return total / len(segment)


def _charge_blockiness(segment: str, block_size: int = 5) -> float:
    charge = np.array([1 if aa in "KR" else (-1 if aa in "DE" else 0) for aa in segment])
    if charge.size < 2 * block_size or np.abs(charge).sum() == 0:
        return 0.0
    block_means = [
        charge[start : start + block_size].mean()
        for start in range(0, len(charge) - block_size + 1, block_size)
    ]
    return float(np.var(block_means))


def _reduced_pair_frequencies(segment: str, gap: int = 2) -> list[float]:
    counts: Counter[str] = Counter()
    total = 0
    for i in range(len(segment) - gap - 1):
        left = REDUCED_CLASS.get(segment[i])
        right = REDUCED_CLASS.get(segment[i + gap + 1])
        if left is not None and right is not None:
            counts[left + right] += 1
            total += 1
    return [counts[pair] / total if total else 0.0 for pair in REDUCED_PAIRS]


def _st_geometry(segment: str) -> list[float]:
    """S/T-mask features; these cannot encode the identity of non-S/T residues."""

    length = len(segment)
    mask = np.fromiter((aa in "ST" for aa in segment), dtype=float, count=length)
    positions = np.flatnonzero(mask)
    gaps = np.diff(positions)
    if len(gaps):
        gap_mean = float(gaps.mean())
        gap_cv = float(gaps.std() / (gap_mean + 1e-9))
        adjacent = float(np.mean(gaps == 1))
        close = float(np.mean(gaps <= 3))
    else:
        gap_mean, gap_cv, adjacent, close = float(length), 1.0, 0.0, 0.0
    autocorrelation = [
        float(np.mean(mask[:-lag] * mask[lag:])) if length > lag else 0.0
        for lag in (1, 2, 3, 5, 8, 10)
    ]
    return [float(mask.mean()), gap_mean, gap_cv, adjacent, close, *autocorrelation]


def extract_region_features(
    sequence: str,
    center: int,
    disorder_by_position: Mapping[int, float] | None = None,
) -> FeatureVector | None:
    """Extract model inputs without consulting PTM sites or labels.

    ``center`` and disorder keys are 1-indexed.  Composition features are
    conditional on a residue being non-S/T, preventing the composition-only
    family from reconstructing S/T density by summing residue fractions.
    """

    values: list[float] = []
    names: list[str] = []
    family_indices: dict[str, list[int]] = {
        "composition": [],
        "geometry": [],
        "disorder": [],
        "arrangement": [],
    }

    def add(family: str, name: str, value: float) -> None:
        family_indices[family].append(len(values))
        names.append(name)
        values.append(float(value))

    for radius in WINDOW_RADII:
        segment = _window(sequence, center, radius)
        if len(segment) < 10:
            return None
        counts = Counter(segment)
        non_st_total = sum(counts[aa] for aa in NON_ST_AMINO_ACIDS)
        denominator = max(non_st_total, 1)
        for aa in NON_ST_AMINO_ACIDS:
            add("composition", f"w{radius}_nonst_{aa}", counts[aa] / denominator)
        non_st_hydropathy = sum(KD[aa] * counts[aa] for aa in NON_ST_AMINO_ACIDS)
        add("composition", f"w{radius}_nonst_hydropathy", non_st_hydropathy / denominator)

        geometry = _st_geometry(segment)
        geometry_names = (
            "st_fraction",
            "st_gap_mean",
            "st_gap_cv",
            "st_adjacent_fraction",
            "st_close_fraction",
            "st_autocorr_1",
            "st_autocorr_2",
            "st_autocorr_3",
            "st_autocorr_5",
            "st_autocorr_8",
            "st_autocorr_10",
        )
        for name, value in zip(geometry_names, geometry, strict=True):
            add("geometry", f"w{radius}_{name}", value)

    arrangement_segment = _window(sequence, center, 25)
    add("arrangement", "scd_w25", _scd(arrangement_segment))
    add("arrangement", "charge_blockiness_w25", _charge_blockiness(arrangement_segment))
    pair_segment = _window(sequence, center, 15)
    for pair, value in zip(REDUCED_PAIRS, _reduced_pair_frequencies(pair_segment), strict=True):
        add("arrangement", f"gap2_pair_{pair}", value)

    disorder_values = []
    if disorder_by_position is not None:
        disorder_values = [
            disorder_by_position[position]
            for position in range(center - 15, center + 16)
            if position in disorder_by_position
        ]
    add(
        "disorder", "disorder_mean_w15", np.mean(disorder_values) if disorder_values else np.nan
    )

    composition = tuple(family_indices["composition"])
    geometry = tuple(family_indices["geometry"])
    disorder = tuple(family_indices["disorder"])
    arrangement = tuple(family_indices["arrangement"])
    families = {
        "composition": composition,
        "geometry": geometry,
        "disorder": disorder,
        "full": composition + geometry + disorder + arrangement,
    }
    return FeatureVector(np.asarray(values, dtype=float), tuple(names), families)


def precision_recall_at_predictions(
    y: np.ndarray, predicted: np.ndarray
) -> tuple[float, float]:
    true_positive = int(np.sum((y == 1) & predicted))
    precision = true_positive / max(int(predicted.sum()), 1)
    recall = true_positive / max(int(y.sum()), 1)
    return precision, recall
