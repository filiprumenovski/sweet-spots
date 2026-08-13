"""Batch kinase scoring without changing the yin-yang estimands."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import kinase_library as kl
import numpy as np
import pandas as pd

from ..analysis.regional_code import load_inputs
from ..analysis.yin_yang import (
    kinase_protein_summary,
    local_opportunity_table,
    matched_kinase_pairs,
    nested_null_table,
    ratio_interval,
)
from ..core.config import load_config
from ..core.io import write_csv, write_json, write_npz
from ..core.randomness import stable_seed
from .sharding import Shard, validate_complete_partition, validate_receipts, write_receipt

PAIR_COLUMNS = (
    "accession",
    "cluster_position",
    "control_position",
    "cluster_peptide",
    "control_peptide",
)


def prepare(config: Path, regions: Path, sites: Path, output_dir: Path) -> None:
    inputs = load_inputs(config, regions, sites)
    local = local_opportunity_table(inputs)
    nested, cutpoints = nested_null_table(inputs)
    pairs = pd.DataFrame(matched_kinase_pairs(inputs), columns=PAIR_COLUMNS)
    pairs.insert(0, "pair_id", np.arange(len(pairs)))
    write_csv(output_dir / "local_opportunity.csv", local)
    write_csv(output_dir / "nested_nulls.csv", nested)
    write_csv(output_dir / "kinase_pairs.csv", pairs)
    write_json(
        output_dir / "prepared.json",
        {
            "schema_version": 1,
            "pairs": len(pairs),
            "disorder_cutpoints": cutpoints.tolist(),
        },
    )


def run_kinase_shard(
    pairs_path: Path,
    shard_index: int,
    shard_count: int,
    sites_output: Path,
    sums_output: Path,
    receipt: Path,
) -> None:
    shard = Shard(shard_index, shard_count)
    pairs = pd.read_csv(pairs_path)
    pairs = pairs.loc[pairs.pair_id.map(lambda value: shard.owns(int(value)))]
    information = kl.get_kinome_info()
    serine_threonine = set(information.loc[information.KL_LIBRARY.eq("ser_thr"), "MATRIX_NAME"])
    family = dict(zip(information.MATRIX_NAME, information.FAMILY, strict=True))
    site_rows: list[dict[str, object]] = []
    score_pair_ids: list[int] = []
    score_kinases: list[str] = []
    score_cluster: list[float] = []
    score_control: list[float] = []
    score_families: list[str] = []
    for row in pairs.itertuples(index=False):
        try:
            cluster = kl.Substrate(row.cluster_peptide).percentile()
            control = kl.Substrate(row.control_peptide).percentile()
        except Exception:
            continue
        kinases = [
            kinase
            for kinase in cluster.index
            if kinase in serine_threonine and kinase in control.index
        ]
        cluster = cluster[kinases].astype(float)
        control = control[kinases].astype(float)
        for kinase in kinases:
            score_pair_ids.append(int(row.pair_id))
            score_kinases.append(kinase)
            score_cluster.append(float(cluster[kinase]))
            score_control.append(float(control[kinase]))
            score_families.append(family.get(kinase, "Other"))
        site_rows.append(
            {
                "pair_id": row.pair_id,
                "accession": row.accession,
                "cluster_position": row.cluster_position,
                "control_position": row.control_position,
                "cluster_n90": int((cluster > 90).sum()),
                "control_n90": int((control > 90).sum()),
                "cluster_max": float(cluster.max()),
                "control_max": float(control.max()),
            }
        )
    site_frame = pd.DataFrame(site_rows)
    write_csv(sites_output, site_frame)
    write_npz(
        sums_output,
        pair_ids=np.asarray(score_pair_ids, dtype=np.int64),
        kinases=np.asarray(score_kinases),
        cluster=np.asarray(score_cluster, dtype=np.float64),
        control=np.asarray(score_control, dtype=np.float64),
        families=np.asarray(score_families),
    )
    write_receipt(
        receipt,
        shard=shard,
        outputs=[sites_output, sums_output],
        records=len(site_frame),
        metadata={"assigned_pairs": len(pairs), "successfully_scored_pairs": len(site_frame)},
    )


def reduce(
    config_path: Path,
    prepared_dir: Path,
    site_shards: list[Path],
    sum_shards: list[Path],
    receipts: list[Path],
    output_dir: Path,
) -> None:
    if not (len(site_shards) == len(sum_shards) == len(receipts)):
        raise ValueError("kinase site, sum, and receipt shard counts must match")
    validate_receipts(receipts, len(receipts))
    kinase_sites = validate_complete_partition(
        [pd.read_csv(path) for path in site_shards], "pair_id"
    )
    arrays = []
    for path in sum_shards:
        with np.load(path) as data:
            arrays.append({name: data[name] for name in data.files})
    pair_ids = np.concatenate([item["pair_ids"] for item in arrays])
    kinases = np.concatenate([item["kinases"] for item in arrays])
    cluster = np.concatenate([item["cluster"] for item in arrays])
    control = np.concatenate([item["control"] for item in arrays])
    families = np.concatenate([item["families"] for item in arrays])
    order = np.argsort(pair_ids, kind="stable")
    accumulators: dict[str, list[object]] = {}
    for index in order:
        kinase = str(kinases[index])
        aggregate = accumulators.setdefault(kinase, [0.0, 0.0, 0, str(families[index])])
        # Reconstruct the original pair-order sum. Shard count then affects
        # scheduling and file layout, never floating-point reduction order.
        aggregate[0] = float(aggregate[0]) + float(cluster[index])
        aggregate[1] = float(aggregate[1]) + float(control[index])
        aggregate[2] = int(aggregate[2]) + 1
    kinase_deltas = pd.DataFrame(
        [
            {
                "kinase": kinase,
                "cluster_mean_percentile": float(values[0]) / int(values[2]),
                "control_mean_percentile": float(values[1]) / int(values[2]),
                "family": values[3],
                "n_pairs": int(values[2]),
            }
            for kinase, values in accumulators.items()
        ]
    )
    kinase_deltas["delta"] = (
        kinase_deltas.cluster_mean_percentile - kinase_deltas.control_mean_percentile
    )
    kinase_deltas = kinase_deltas[
        [
            "kinase",
            "cluster_mean_percentile",
            "control_mean_percentile",
            "delta",
            "family",
            "n_pairs",
        ]
    ].sort_values("delta", kind="stable")
    kinase_proteins = kinase_protein_summary(kinase_sites)
    config = load_config(config_path)
    base_seed = int(config.values["randomness"]["analysis_base_seed"])
    draws = int(config.values["randomness"]["figure_bootstrap_draws"])
    local = pd.read_csv(prepared_dir / "local_opportunity.csv")
    nested = pd.read_csv(prepared_dir / "nested_nulls.csv")
    prepared = json.loads((prepared_dir / "prepared.json").read_text(encoding="utf-8"))
    local_results = {}
    for arm, selected in (("all", local), ("regional", local.loc[local.regional])):
        for window in (10, 25, 50):
            for residue_matched in (True, False):
                subset = selected.loc[
                    (selected.window == window) & (selected.residue_matched == residue_matched)
                ]
                key = f"{arm}_w{window}_{'residue_matched' if residue_matched else 'all_st'}"
                local_results[key] = ratio_interval(
                    subset,
                    expected="expected",
                    draws=draws,
                    seed=stable_seed(base_seed, "local", arm, window, residue_matched),
                )
    nested_results = {}
    for arm in ("all", "regional"):
        selected = nested.loc[nested.arm.eq(arm)]
        nested_results[arm] = {
            name: ratio_interval(
                selected,
                expected=column,
                draws=draws,
                seed=stable_seed(base_seed, "nested", arm, name),
            )
            for name, column in (
                ("N0_length", "expected_length"),
                ("N1_st", "expected_st"),
                ("N2_disorder", "expected_disorder"),
            )
        }
    delta = kinase_proteins.delta_n90.to_numpy(float)
    rng = np.random.default_rng(stable_seed(base_seed, "kinase", "bootstrap"))
    bootstrap = np.empty(20_000, dtype=float)
    for draw in range(len(bootstrap)):
        bootstrap[draw] = delta[rng.integers(0, len(delta), len(delta))].mean()
    output_dir.mkdir(parents=True, exist_ok=True)
    write_csv(output_dir / "local_opportunity.csv", local)
    write_csv(output_dir / "nested_nulls.csv", nested)
    write_csv(output_dir / "kinase_site_scores.csv", kinase_deltas)
    write_csv(output_dir / "kinase_summary.csv", kinase_proteins)
    write_json(
        output_dir / "summary.json",
        {
            "local_opportunity": local_results,
            "nested_nulls": nested_results,
            "disorder_cutpoints": prepared["disorder_cutpoints"],
            "kinome": {
                "n_site_pairs": len(kinase_sites),
                "n_proteins": len(kinase_proteins),
                "n_kinases": len(kinase_deltas),
                "cluster_n90": float(kinase_proteins.cluster_n90.mean()),
                "control_n90": float(kinase_proteins.control_n90.mean()),
                "delta_n90": float(delta.mean()),
                "delta_n90_ci": np.quantile(bootstrap, [0.025, 0.975]).tolist(),
            },
        },
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    prepare_parser = commands.add_parser("prepare")
    prepare_parser.add_argument("--config", type=Path, required=True)
    prepare_parser.add_argument("--regions", type=Path, required=True)
    prepare_parser.add_argument("--sites", type=Path, required=True)
    prepare_parser.add_argument("--output-dir", type=Path, required=True)
    shard_parser = commands.add_parser("kinase-shard")
    shard_parser.add_argument("--pairs", type=Path, required=True)
    shard_parser.add_argument("--shard", type=int, required=True)
    shard_parser.add_argument("--shards", type=int, required=True)
    shard_parser.add_argument("--sites-output", type=Path, required=True)
    shard_parser.add_argument("--sums-output", type=Path, required=True)
    shard_parser.add_argument("--receipt", type=Path, required=True)
    reduce_parser = commands.add_parser("reduce")
    reduce_parser.add_argument("--config", type=Path, required=True)
    reduce_parser.add_argument("--prepared-dir", type=Path, required=True)
    reduce_parser.add_argument("--site-shards", type=Path, nargs="+", required=True)
    reduce_parser.add_argument("--sum-shards", type=Path, nargs="+", required=True)
    reduce_parser.add_argument("--receipts", type=Path, nargs="+", required=True)
    reduce_parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    if args.command == "prepare":
        prepare(args.config, args.regions, args.sites, args.output_dir)
    elif args.command == "kinase-shard":
        run_kinase_shard(
            args.pairs,
            args.shard,
            args.shards,
            args.sites_output,
            args.sums_output,
            args.receipt,
        )
    else:
        reduce(
            args.config,
            args.prepared_dir,
            args.site_shards,
            args.sum_shards,
            args.receipts,
            args.output_dir,
        )


if __name__ == "__main__":
    main()
