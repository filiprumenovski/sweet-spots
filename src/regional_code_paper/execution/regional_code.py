"""Deterministic map/reduce execution for the regional-code analyses."""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.ensemble import HistGradientBoostingClassifier
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import GroupKFold

from ..analysis.regional_code import (
    MODEL_FEATURES,
    acceptor_models,
    bootstrap_interval,
    build_tiles,
    load_inputs,
    per_protein_auc,
    phosphosite_distance_profile,
    proline_profiles,
)
from ..core.config import load_config
from ..core.io import write_csv, write_json, write_npz
from ..core.randomness import stable_seed
from .sharding import Shard, validate_complete_partition, validate_receipts, write_receipt


def prepare(
    config_path: Path,
    regions_path: Path,
    sites_path: Path,
    output_dir: Path,
) -> None:
    config = load_config(config_path)
    inputs = load_inputs(config_path, regions_path, sites_path)
    base_seed = int(config.values["randomness"]["manuscript_base_seed"])
    tiles = build_tiles(inputs).reset_index(drop=True)
    tiles.insert(0, "row_id", np.arange(len(tiles)))
    feature_names = sorted(
        {column for columns in MODEL_FEATURES.values() for column in columns}
    )
    write_npz(
        output_dir / "tile_features.npz",
        matrix=tiles[feature_names].to_numpy(float),
        feature_names=np.asarray(feature_names),
        labels=tiles.label.to_numpy(int),
        groups=tiles.accession.to_numpy(str),
        row_ids=tiles.row_id.to_numpy(int),
    )
    acceptors, counts = acceptor_models(inputs)
    proline = proline_profiles(inputs)
    distance_matches, distance = phosphosite_distance_profile(
        inputs, draws=5_000, base_seed=base_seed
    )
    write_csv(output_dir / "tile_features.csv", tiles)
    write_csv(output_dir / "acceptor_model.csv", acceptors)
    write_csv(output_dir / "proline_profiles.csv", proline)
    write_csv(output_dir / "distance_profile.csv", distance)
    write_csv(output_dir / "distance_matches.csv", distance_matches)
    write_json(
        output_dir / "prepared.json",
        {
            "schema_version": 1,
            "tiles": len(tiles),
            "acceptor_model": counts,
            "distance": {
                "n_pairs": len(distance_matches),
                "n_proteins": distance_matches.accession.nunique(),
            },
        },
    )


def run_fold(
    cache_path: Path,
    fold: int,
    seed: int,
    output: Path,
    receipt: Path,
) -> None:
    shard = Shard(fold, 5)
    with np.load(cache_path) as data:
        matrix = data["matrix"]
        feature_names = data["feature_names"].astype(str).tolist()
        labels = data["labels"]
        groups = data["groups"]
        row_ids = data["row_ids"]
    column_index = {name: index for index, name in enumerate(feature_names)}
    train, test = list(GroupKFold(5).split(matrix, labels, groups))[fold]
    rows = pd.DataFrame({"row_id": row_ids[test], "fold": fold})
    for name, columns in MODEL_FEATURES.items():
        selected = [column_index[column] for column in columns]
        model = HistGradientBoostingClassifier(
            max_iter=300,
            max_depth=3,
            learning_rate=0.04,
            l2_regularization=1.0,
            random_state=seed,
        )
        model.fit(matrix[train][:, selected], labels[train])
        rows[f"{name}_score"] = model.predict_proba(matrix[test][:, selected])[:, 1]
    rows = rows.sort_values("row_id", kind="stable")
    write_csv(output, rows)
    write_receipt(receipt, shard=shard, outputs=[output], records=len(rows))


def reduce(
    config_path: Path,
    prepared_dir: Path,
    fold_paths: list[Path],
    fold_receipts: list[Path],
    output_dir: Path,
) -> None:
    config = load_config(config_path)
    base_seed = int(config.values["randomness"]["manuscript_base_seed"])
    figure_draws = int(config.values["randomness"]["figure_bootstrap_draws"])
    if len(fold_paths) != 5 or len(fold_receipts) != 5:
        raise ValueError("regional-code prediction requires exactly five fold shards")
    validate_receipts(fold_receipts, 5)
    tiles = pd.read_csv(prepared_dir / "tile_features.csv")
    predictions = validate_complete_partition(
        [pd.read_csv(path) for path in fold_paths], "row_id"
    )
    if not np.array_equal(predictions.row_id.to_numpy(int), tiles.row_id.to_numpy(int)):
        raise ValueError("regional-code folds do not cover all prepared tiles")
    for name in MODEL_FEATURES:
        tiles[f"{name}_score"] = predictions[f"{name}_score"].to_numpy(float)
    proteins = per_protein_auc(tiles)
    prepared = __import__("json").loads(
        (prepared_dir / "prepared.json").read_text(encoding="utf-8")
    )
    labels = tiles.label.to_numpy(int)
    matched = (tiles.disorder >= tiles.disorder.median()) & (
        tiles.st_fraction >= tiles.st_fraction.median()
    )
    within_values = proteins.within_protein_auroc.to_numpy(float)
    within_ci = bootstrap_interval(
        within_values,
        draws=figure_draws,
        seed=stable_seed(base_seed, "within_protein_auroc"),
    )
    output_dir.mkdir(parents=True, exist_ok=True)
    write_csv(output_dir / "held_out_tiles.csv", tiles.drop(columns=["row_id"]))
    write_csv(output_dir / "per_protein.csv", proteins)
    for name in (
        "acceptor_model.csv",
        "proline_profiles.csv",
        "distance_profile.csv",
        "distance_matches.csv",
    ):
        write_csv(output_dir / name, pd.read_csv(prepared_dir / name))
    scores = {name: tiles[f"{name}_score"].to_numpy(float) for name in MODEL_FEATURES}
    write_json(
        output_dir / "summary.json",
        {
            "tiles": {
                "n_tiles": len(tiles),
                "n_positive": int(labels.sum()),
                "n_proteins": tiles.accession.nunique(),
                "global_auroc": {
                    name: float(roc_auc_score(labels, values))
                    for name, values in scores.items()
                },
                "matched": {
                    "n_tiles": int(matched.sum()),
                    "composition_auroc": float(
                        roc_auc_score(labels[matched], scores["composition"][matched])
                    ),
                    "geometry_auroc": float(
                        roc_auc_score(labels[matched], scores["geometry"][matched])
                    ),
                },
            },
            "within_protein": {
                "mean_auroc": float(within_values.mean()),
                "ci": within_ci,
                "n_proteins": len(proteins),
            },
            "acceptor_model": prepared["acceptor_model"],
            "distance": prepared["distance"],
        },
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    prepare_parser = commands.add_parser("prepare")
    for item in (prepare_parser,):
        item.add_argument("--config", type=Path, required=True)
        item.add_argument("--regions", type=Path, required=True)
        item.add_argument("--sites", type=Path, required=True)
    prepare_parser.add_argument("--output-dir", type=Path, required=True)
    fold_parser = commands.add_parser("fold")
    fold_parser.add_argument("--cache", type=Path, required=True)
    fold_parser.add_argument("--fold", type=int, required=True)
    fold_parser.add_argument("--seed", type=int, required=True)
    fold_parser.add_argument("--output", type=Path, required=True)
    fold_parser.add_argument("--receipt", type=Path, required=True)
    reduce_parser = commands.add_parser("reduce")
    reduce_parser.add_argument("--config", type=Path, required=True)
    reduce_parser.add_argument("--prepared-dir", type=Path, required=True)
    reduce_parser.add_argument("--folds", type=Path, nargs="+", required=True)
    reduce_parser.add_argument("--fold-receipts", type=Path, nargs="+", required=True)
    reduce_parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    if args.command == "prepare":
        prepare(args.config, args.regions, args.sites, args.output_dir)
    elif args.command == "fold":
        run_fold(args.cache, args.fold, args.seed, args.output, args.receipt)
    else:
        reduce(
            args.config,
            args.prepared_dir,
            args.folds,
            args.fold_receipts,
            args.output_dir,
        )


if __name__ == "__main__":
    main()
