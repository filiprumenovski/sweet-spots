"""Nested, protein-grouped validation of the enhanced regional scanner.

Candidate selection and operating-threshold calibration occur only inside each
outer-training partition. The outer test fold is used once, for evaluation.
"""

from __future__ import annotations

import argparse
import json
import math
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.base import clone
from sklearn.ensemble import HistGradientBoostingClassifier
from sklearn.metrics import average_precision_score, roc_auc_score
from sklearn.model_selection import GroupKFold

from ..core.config import load_config
from ..core.io import write_json
from .scanner_features import (
    REDUCED_CLASS,
    REDUCED_PAIRS,
    extract_region_features,
    make_region_classifier,
    precision_recall_at_predictions,
)

QUANTILES = (0.90, 0.95, 0.98, 0.99, 0.995)


def legacy_composition(segment: str) -> list[float]:
    length = len(segment)
    counts = Counter(segment)

    def fraction(residues: str) -> float:
        return sum(counts[residue] for residue in residues) / length

    hydropathy = dict(
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
    return [
        fraction("ST"),
        fraction("DEFWY"),
        (counts["K"] + counts["R"] - counts["D"] - counts["E"]) / length,
        sum(hydropathy.get(residue, 0.0) * count for residue, count in counts.items()) / length,
        counts["P"] / length,
        fraction("QN"),
    ]


def legacy_scd(segment: str) -> float:
    charge = np.asarray(
        [1 if residue in "KR" else (-1 if residue in "DE" else 0) for residue in segment]
    )
    charged = np.flatnonzero(charge)
    if len(charged) < 2:
        return 0.0
    total = 0.0
    for left in range(len(charged)):
        for right in range(left + 1, len(charged)):
            i, j = charged[left], charged[right]
            total += charge[i] * charge[j] * math.sqrt(j - i)
    return total / len(segment)


def legacy_pairs(segment: str, gap: int = 2) -> list[float]:
    counts = Counter()
    total = 0
    for index in range(len(segment) - gap - 1):
        left = REDUCED_CLASS.get(segment[index])
        right = REDUCED_CLASS.get(segment[index + gap + 1])
        if left is not None and right is not None:
            counts[left + right] += 1
            total += 1
    return [counts[pair] / total if total else 0.0 for pair in REDUCED_PAIRS]


def legacy_features(
    sequence: str, center: int, disorder: dict[int, float]
) -> list[float] | None:
    values = []
    for radius in (10, 20, 40):
        segment = sequence[max(0, center - 1 - radius) : min(len(sequence), center + radius)]
        if len(segment) < 10:
            return None
        values.extend(legacy_composition(segment))
    values.append(legacy_scd(sequence[max(0, center - 26) : min(len(sequence), center + 25)]))
    disorder_window = [
        disorder[position]
        for position in range(center - 15, center + 16)
        if position in disorder
    ]
    values.append(float(np.mean(disorder_window)) if disorder_window else np.nan)
    values.extend(legacy_pairs(sequence[max(0, center - 16) : min(len(sequence), center + 15)]))
    return values


def load_dataset(config_path: Path, regions_path: Path):
    config = load_config(config_path)
    root = config.source_root
    fasta = pd.read_parquet(
        root / "data/interim/fasta_all.parquet",
        columns=["accession", "sequence", "taxon_id", "is_canonical"],
    )
    fasta = fasta.loc[(fasta.taxon_id == 9606) & fasta.is_canonical]
    fasta = fasta.drop_duplicates("accession")
    sequences = dict(zip(fasta.accession.astype(str), fasta.sequence.astype(str), strict=True))
    disorder_frame = pd.read_parquet(
        root / "data/interim/iupred_residue_scores.parquet",
        columns=["accession", "position", "disorder_score"],
    )
    disorder: dict[str, dict[int, float]] = defaultdict(dict)
    for accession, position, score in disorder_frame.itertuples(index=False, name=None):
        disorder[str(accession)][int(position)] = float(score)
    regions = pd.read_csv(regions_path)
    spans: dict[str, list[tuple[int, int]]] = defaultdict(list)
    for accession, start, end in regions[["accession", "start", "end"]].itertuples(
        index=False, name=None
    ):
        spans[str(accession)].append((int(start), int(end)))

    # The scanner's estimand is the complete validated atlas universe, not only
    # proteins that contain a positive region.  Keeping the two collections
    # separate is important: ``spans`` defines labels, whereas ``site_proteins``
    # defines the population over which false positives are measured.
    atlas = pd.read_csv(
        root / "analysis/revalidation/data/atlas_unambiguous.csv",
        usecols=["species", "accession", "position_in_protein", "site_residue"],
        dtype=str,
        encoding_errors="replace",
    )
    atlas = atlas.loc[atlas.species.eq("human") & atlas.site_residue.isin(["S", "T"])].copy()
    atlas["position"] = pd.to_numeric(atlas.position_in_protein, errors="coerce")
    site_proteins: set[str] = set()
    for accession, position, residue in (
        atlas[["accession", "position", "site_residue"]]
        .dropna(subset=["position"])
        .drop_duplicates()
        .itertuples(index=False, name=None)
    ):
        sequence = sequences.get(str(accession))
        position = int(position)
        if (
            sequence is not None
            and 1 <= position <= len(sequence)
            and sequence[position - 1] == residue
        ):
            site_proteins.add(str(accession))

    legacy, enhanced, labels, groups, accessions, centers = [], [], [], [], [], []
    family_indices = None
    # Sorting is essential. HistGradientBoosting's internal validation split is
    # row-order-sensitive; unordered table iteration makes repeated runs drift.
    for accession in sorted(site_proteins):
        sequence = sequences[accession]
        for center in range(11, len(sequence) - 10, 10):
            old = legacy_features(sequence, center, disorder[accession])
            new = extract_region_features(sequence, center, disorder[accession])
            if old is None or new is None:
                continue
            family_indices = new.families if family_indices is None else family_indices
            legacy.append(old)
            enhanced.append(new.values)
            labels.append(
                int(any(start - 15 <= center <= end + 15 for start, end in spans[accession]))
            )
            groups.append(accession)
            accessions.append(accession)
            centers.append(center)
    if family_indices is None:
        raise RuntimeError("No scanner features were extracted")
    return (
        np.asarray(legacy, dtype=float),
        np.asarray(enhanced, dtype=float),
        np.asarray(labels, dtype=int),
        np.asarray(groups),
        np.asarray(accessions),
        np.asarray(centers, dtype=int),
        dict(spans),
        family_indices,
    )


def baseline_estimator() -> HistGradientBoostingClassifier:
    return HistGradientBoostingClassifier(
        max_iter=600,
        max_depth=4,
        learning_rate=0.03,
        l2_regularization=1.0,
        random_state=42,
    )


def inner_oof(estimator, matrix, labels, groups) -> np.ndarray:
    predictions = np.full(len(labels), np.nan)
    for train, validation in GroupKFold(3).split(matrix, labels, groups):
        fitted = clone(estimator).fit(matrix[train], labels[train])
        predictions[validation] = fitted.predict_proba(matrix[validation])[:, 1]
    return predictions


def evaluate_nested(legacy, enhanced, labels, groups):
    candidates = {
        "legacy": (legacy, baseline_estimator()),
        "enhanced": (enhanced, make_region_classifier()),
        "enhanced_balanced": (
            enhanced,
            make_region_classifier().set_params(class_weight="balanced"),
        ),
    }
    scores = {
        name: np.full(len(labels), np.nan) for name in ("baseline", "enhanced", "selected")
    }
    fold_assignment = np.full(len(labels), -1, dtype=int)
    predicted = {
        name: {quantile: np.zeros(len(labels), dtype=bool) for quantile in QUANTILES}
        for name in scores
    }
    fold_rows = []
    for fold_index, (outer_train, outer_test) in enumerate(
        GroupKFold(5).split(legacy, labels, groups)
    ):
        inner_predictions = {}
        inner_ap = {}
        for name, (matrix, estimator) in candidates.items():
            values = inner_oof(
                estimator,
                matrix[outer_train],
                labels[outer_train],
                groups[outer_train],
            )
            inner_predictions[name] = values
            inner_ap[name] = float(average_precision_score(labels[outer_train], values))
        selected_name = max(inner_ap, key=inner_ap.get)
        selected_matrix, selected_estimator = candidates[selected_name]
        fitted = {
            "baseline": baseline_estimator().fit(legacy[outer_train], labels[outer_train]),
            "enhanced": make_region_classifier().fit(
                enhanced[outer_train], labels[outer_train]
            ),
            "selected": clone(selected_estimator).fit(
                selected_matrix[outer_train], labels[outer_train]
            ),
        }
        score_matrices = {
            "baseline": legacy,
            "enhanced": enhanced,
            "selected": selected_matrix,
        }
        candidate_names = {
            "baseline": "legacy",
            "enhanced": "enhanced",
            "selected": selected_name,
        }
        for name in scores:
            scores[name][outer_test] = fitted[name].predict_proba(
                score_matrices[name][outer_test]
            )[:, 1]
            calibration = inner_predictions[candidate_names[name]][labels[outer_train] == 0]
            for quantile in QUANTILES:
                threshold = float(np.quantile(calibration, quantile))
                predicted[name][quantile][outer_test] = scores[name][outer_test] >= threshold
        fold_assignment[outer_test] = fold_index
        fold_rows.append(
            {
                "fold": fold_index,
                "n_test": len(outer_test),
                "n_test_positive": int(labels[outer_test].sum()),
                "selected": selected_name,
                "inner_average_precision": json.dumps(inner_ap, sort_keys=True),
                "baseline_auroc": roc_auc_score(
                    labels[outer_test], scores["baseline"][outer_test]
                ),
                "baseline_average_precision": average_precision_score(
                    labels[outer_test], scores["baseline"][outer_test]
                ),
                "enhanced_auroc": roc_auc_score(
                    labels[outer_test], scores["enhanced"][outer_test]
                ),
                "enhanced_average_precision": average_precision_score(
                    labels[outer_test], scores["enhanced"][outer_test]
                ),
            }
        )
    return scores, fold_assignment, predicted, pd.DataFrame(fold_rows)


def cluster_recall(predicted, accessions, centers, spans) -> dict[str, object]:
    centers_by_accession: dict[str, list[int]] = defaultdict(list)
    for is_predicted, accession, center in zip(predicted, accessions, centers, strict=True):
        if is_predicted:
            centers_by_accession[accession].append(int(center))
    hit = 0
    total = 0
    for accession, regions in spans.items():
        for start, end in regions:
            total += 1
            hit += int(
                any(
                    start - 15 <= center <= end + 15
                    for center in centers_by_accession[accession]
                )
            )
    return {"hit": hit, "total": total, "recall": hit / total}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--regions", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    output = args.output_dir.resolve()
    output.mkdir(parents=True, exist_ok=True)
    legacy, enhanced, labels, groups, accessions, centers, spans, families = load_dataset(
        args.config, args.regions
    )
    scores, folds, predicted, fold_table = evaluate_nested(legacy, enhanced, labels, groups)
    prediction_table = pd.DataFrame(
        {
            "accession": accessions,
            "center": centers,
            "label": labels,
            "fold": folds,
            "baseline_score": scores["baseline"],
            "enhanced_score": scores["enhanced"],
            "selected_score": scores["selected"],
        }
    )
    operating_rows = []
    operating_summary = {}
    for model in scores:
        operating_summary[model] = {}
        for quantile in QUANTILES:
            precision, recall = precision_recall_at_predictions(
                labels, predicted[model][quantile]
            )
            region_recall = cluster_recall(
                predicted[model][quantile], accessions, centers, spans
            )
            row = {
                "model": model,
                "negative_quantile": quantile,
                "precision": precision,
                "tile_recall": recall,
                "n_predicted": int(predicted[model][quantile].sum()),
                "region_hit": region_recall["hit"],
                "region_total": region_recall["total"],
                "region_recall": region_recall["recall"],
            }
            operating_rows.append(row)
            operating_summary[model][str(quantile)] = row
    family_results = {}
    splitter = GroupKFold(5)
    for family in ("composition", "geometry", "disorder"):
        matrix = enhanced[:, np.asarray(families[family])]
        family_scores = np.full(len(labels), np.nan)
        for train, test in splitter.split(matrix, labels, groups):
            fitted = make_region_classifier().fit(matrix[train], labels[train])
            family_scores[test] = fitted.predict_proba(matrix[test])[:, 1]
        family_results[family] = {
            "n_features": matrix.shape[1],
            "auroc": float(roc_auc_score(labels, family_scores)),
            "average_precision": float(average_precision_score(labels, family_scores)),
        }
    prediction_table.to_csv(output / "outer_predictions.csv", index=False)
    fold_table.to_csv(output / "outer_folds.csv", index=False)
    pd.DataFrame(operating_rows).to_csv(output / "operating_points.csv", index=False)
    write_json(
        output / "summary.json",
        {
            "cohort": {
                "tiles": len(labels),
                "positive_tiles": int(labels.sum()),
                "prevalence": float(labels.mean()),
                "proteins": len(np.unique(groups)),
                "regions": sum(map(len, spans.values())),
            },
            "baseline": {
                "auroc": float(roc_auc_score(labels, scores["baseline"])),
                "average_precision": float(average_precision_score(labels, scores["baseline"])),
            },
            "enhanced": {
                "auroc": float(roc_auc_score(labels, scores["enhanced"])),
                "average_precision": float(average_precision_score(labels, scores["enhanced"])),
            },
            "selected": {
                "auroc": float(roc_auc_score(labels, scores["selected"])),
                "average_precision": float(average_precision_score(labels, scores["selected"])),
            },
            "operating_points": operating_summary,
            "feature_families": family_results,
            "outer_selection": fold_table.selected.tolist(),
        },
    )


if __name__ == "__main__":
    main()
