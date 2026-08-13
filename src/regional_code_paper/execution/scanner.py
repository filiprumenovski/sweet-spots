"""Restartable map/reduce implementation of nested scanner validation.

``prepare`` extracts the expensive sequence features once. ``fold`` owns one
outer GroupKFold partition, including all inner model selection and threshold
calibration. ``reduce`` refuses incomplete or duplicated row partitions before
publishing the manuscript-facing tables.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.base import clone
from sklearn.metrics import average_precision_score, roc_auc_score
from sklearn.model_selection import GroupKFold

from ..core.io import write_csv, write_json, write_npz
from ..models.scanner import (
    QUANTILES,
    baseline_estimator,
    cluster_recall,
    inner_oof,
    load_dataset,
)
from ..models.scanner_features import make_region_classifier, precision_recall_at_predictions
from .sharding import Shard, validate_complete_partition, validate_receipts, write_receipt


def prepare(config: Path, regions: Path, cache: Path, metadata: Path) -> None:
    legacy, enhanced, labels, groups, accessions, centers, spans, families = load_dataset(
        config, regions
    )
    write_npz(
        cache,
        legacy=legacy,
        enhanced=enhanced,
        labels=labels,
        groups=groups,
        accessions=accessions,
        centers=centers,
    )
    write_json(
        metadata,
        {
            "schema_version": 1,
            "rows": len(labels),
            "positive_rows": int(labels.sum()),
            "groups": len(np.unique(groups)),
            "families": families,
            "spans": {key: value for key, value in sorted(spans.items())},
        },
    )


def run_fold(cache: Path, metadata: Path, fold: int, output: Path, receipt: Path) -> None:
    shard = Shard(fold, 5)
    with np.load(cache) as data:
        legacy = data["legacy"]
        enhanced = data["enhanced"]
        labels = data["labels"]
        groups = data["groups"]
    document = json.loads(metadata.read_text(encoding="utf-8"))
    families = document["families"]
    outer_train, outer_test = list(GroupKFold(5).split(legacy, labels, groups))[fold]
    candidates = {
        "legacy": (legacy, baseline_estimator()),
        "enhanced": (enhanced, make_region_classifier()),
        "enhanced_balanced": (
            enhanced,
            make_region_classifier().set_params(class_weight="balanced"),
        ),
    }
    inner_predictions: dict[str, np.ndarray] = {}
    inner_ap: dict[str, float] = {}
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
    estimators = {
        "baseline": baseline_estimator().fit(legacy[outer_train], labels[outer_train]),
        "enhanced": make_region_classifier().fit(enhanced[outer_train], labels[outer_train]),
        "selected": clone(selected_estimator).fit(
            selected_matrix[outer_train], labels[outer_train]
        ),
    }
    matrices = {"baseline": legacy, "enhanced": enhanced, "selected": selected_matrix}
    candidate_names = {
        "baseline": "legacy",
        "enhanced": "enhanced",
        "selected": selected_name,
    }
    rows = pd.DataFrame({"row_id": outer_test, "fold": fold})
    thresholds: dict[str, dict[str, float]] = {}
    for name, estimator in estimators.items():
        scores = estimator.predict_proba(matrices[name][outer_test])[:, 1]
        rows[f"{name}_score"] = scores
        calibration = inner_predictions[candidate_names[name]][labels[outer_train] == 0]
        thresholds[name] = {}
        for quantile in QUANTILES:
            threshold = float(np.quantile(calibration, quantile))
            thresholds[name][str(quantile)] = threshold
            rows[f"{name}_q{quantile:g}"] = scores >= threshold
    for family in ("composition", "geometry", "disorder"):
        columns = np.asarray(families[family], dtype=int)
        estimator = make_region_classifier().fit(
            enhanced[outer_train][:, columns], labels[outer_train]
        )
        rows[f"family_{family}_score"] = estimator.predict_proba(
            enhanced[outer_test][:, columns]
        )[:, 1]
    rows = rows.sort_values("row_id", kind="stable")
    write_csv(output, rows)
    write_receipt(
        receipt,
        shard=shard,
        outputs=[output],
        records=len(rows),
        metadata={
            "selected": selected_name,
            "inner_average_precision": inner_ap,
            "thresholds": thresholds,
            "test_positive": int(labels[outer_test].sum()),
        },
    )


def reduce_folds(
    cache: Path,
    metadata: Path,
    shards: list[Path],
    receipts: list[Path],
    output_dir: Path,
) -> None:
    if len(shards) != 5 or len(receipts) != 5:
        raise ValueError("scanner reduction requires exactly five outer folds")
    receipt_documents = validate_receipts(receipts, 5)
    combined = validate_complete_partition([pd.read_csv(path) for path in shards], "row_id")
    with np.load(cache) as data:
        labels = data["labels"]
        groups = data["groups"]
        accessions = data["accessions"]
        centers = data["centers"]
    if not np.array_equal(combined.row_id.to_numpy(int), np.arange(len(labels))):
        raise ValueError("scanner folds do not cover every prepared row exactly once")
    document = json.loads(metadata.read_text(encoding="utf-8"))
    spans = {
        accession: [tuple(interval) for interval in intervals]
        for accession, intervals in document["spans"].items()
    }
    output_dir.mkdir(parents=True, exist_ok=True)
    prediction_table = pd.DataFrame(
        {
            "accession": accessions,
            "center": centers,
            "label": labels,
            "fold": combined.fold.to_numpy(int),
            "baseline_score": combined.baseline_score,
            "enhanced_score": combined.enhanced_score,
            "selected_score": combined.selected_score,
        }
    )
    fold_rows = []
    for receipt in sorted(receipt_documents, key=lambda item: item["shard"]):
        fold = int(receipt["shard"])
        subset = combined.loc[combined.fold.eq(fold)]
        meta = receipt["metadata"]
        fold_rows.append(
            {
                "fold": fold,
                "n_test": len(subset),
                "n_test_positive": int(labels[subset.row_id.to_numpy(int)].sum()),
                "selected": meta["selected"],
                "inner_average_precision": json.dumps(
                    meta["inner_average_precision"], sort_keys=True
                ),
                "baseline_auroc": roc_auc_score(
                    labels[subset.row_id.to_numpy(int)], subset.baseline_score
                ),
                "baseline_average_precision": average_precision_score(
                    labels[subset.row_id.to_numpy(int)], subset.baseline_score
                ),
                "enhanced_auroc": roc_auc_score(
                    labels[subset.row_id.to_numpy(int)], subset.enhanced_score
                ),
                "enhanced_average_precision": average_precision_score(
                    labels[subset.row_id.to_numpy(int)], subset.enhanced_score
                ),
            }
        )
    fold_table = pd.DataFrame(fold_rows)
    operating_rows: list[dict[str, object]] = []
    operating_summary: dict[str, dict[str, object]] = {}
    for model in ("baseline", "enhanced", "selected"):
        operating_summary[model] = {}
        for quantile in QUANTILES:
            predicted = combined[f"{model}_q{quantile:g}"].to_numpy(bool)
            precision, recall = precision_recall_at_predictions(labels, predicted)
            region_recall = cluster_recall(predicted, accessions, centers, spans)
            row = {
                "model": model,
                "negative_quantile": quantile,
                "precision": precision,
                "tile_recall": recall,
                "n_predicted": int(predicted.sum()),
                "region_hit": region_recall["hit"],
                "region_total": region_recall["total"],
                "region_recall": region_recall["recall"],
            }
            operating_rows.append(row)
            operating_summary[model][str(quantile)] = row
    family_results = {}
    for family in ("composition", "geometry", "disorder"):
        scores = combined[f"family_{family}_score"].to_numpy(float)
        family_results[family] = {
            "n_features": len(document["families"][family]),
            "auroc": float(roc_auc_score(labels, scores)),
            "average_precision": float(average_precision_score(labels, scores)),
        }
    write_csv(output_dir / "outer_predictions.csv", prediction_table)
    write_csv(output_dir / "outer_folds.csv", fold_table)
    write_csv(output_dir / "operating_points.csv", pd.DataFrame(operating_rows))
    score = {
        model: combined[f"{model}_score"].to_numpy(float)
        for model in ("baseline", "enhanced", "selected")
    }
    write_json(
        output_dir / "summary.json",
        {
            "cohort": {
                "tiles": len(labels),
                "positive_tiles": int(labels.sum()),
                "prevalence": float(labels.mean()),
                "proteins": len(np.unique(groups)),
                "regions": sum(map(len, spans.values())),
            },
            **{
                model: {
                    "auroc": float(roc_auc_score(labels, score[model])),
                    "average_precision": float(average_precision_score(labels, score[model])),
                }
                for model in score
            },
            "operating_points": operating_summary,
            "feature_families": family_results,
            "outer_selection": fold_table.selected.tolist(),
        },
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)
    prepare_parser = subparsers.add_parser("prepare")
    prepare_parser.add_argument("--config", type=Path, required=True)
    prepare_parser.add_argument("--regions", type=Path, required=True)
    prepare_parser.add_argument("--cache", type=Path, required=True)
    prepare_parser.add_argument("--metadata", type=Path, required=True)
    fold_parser = subparsers.add_parser("fold")
    fold_parser.add_argument("--cache", type=Path, required=True)
    fold_parser.add_argument("--metadata", type=Path, required=True)
    fold_parser.add_argument("--fold", type=int, required=True)
    fold_parser.add_argument("--output", type=Path, required=True)
    fold_parser.add_argument("--receipt", type=Path, required=True)
    reduce_parser = subparsers.add_parser("reduce")
    reduce_parser.add_argument("--cache", type=Path, required=True)
    reduce_parser.add_argument("--metadata", type=Path, required=True)
    reduce_parser.add_argument("--shards", type=Path, nargs="+", required=True)
    reduce_parser.add_argument("--receipts", type=Path, nargs="+", required=True)
    reduce_parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    if args.command == "prepare":
        prepare(args.config, args.regions, args.cache, args.metadata)
    elif args.command == "fold":
        run_fold(args.cache, args.metadata, args.fold, args.output, args.receipt)
    else:
        reduce_folds(args.cache, args.metadata, args.shards, args.receipts, args.output_dir)


if __name__ == "__main__":
    main()
