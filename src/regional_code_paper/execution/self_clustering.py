"""PTM-parallel execution and deterministic reduction of self-clustering."""

from __future__ import annotations

import argparse
import platform
from datetime import UTC, datetime
from pathlib import Path

import numpy as np
import pandas as pd
import scipy

from ..analysis.self_clustering import (
    BASE_SEED,
    MIN_SITES,
    PRIMARY_RADIUS,
    PTMS,
    RADII,
    direct_contrasts,
    load_inputs,
    matched_fold_ratios,
    protein_rows,
    self_test,
    sha256,
    summarize,
)
from ..core.config import load_config
from ..core.io import write_csv, write_json
from .sharding import Shard, validate_complete_partition, validate_receipts, write_receipt


def run_ptm(config_path: Path, ptm: str, output: Path, receipt: Path) -> None:
    if ptm not in PTMS:
        raise ValueError(f"unknown PTM {ptm!r}; expected one of {PTMS}")
    config = load_config(config_path)
    self_test()
    sequences, sites, _ = load_inputs(config.source_root)
    frame = protein_rows(sequences, sites, ptms=(ptm,))
    row_ids = [
        f"{row.universe}|{row.ptm}|{row.accession}|{row.radius}" for row in frame.itertuples()
    ]
    frame.insert(0, "row_id", row_ids)
    write_csv(output, frame)
    index = PTMS.index(ptm)
    write_receipt(
        receipt,
        shard=Shard(index, len(PTMS)),
        outputs=[output],
        records=len(frame),
        metadata={"ptm": ptm},
    )


def reduce(config_path: Path, shards: list[Path], receipts: list[Path], output: Path) -> None:
    if len(shards) != len(PTMS) or len(receipts) != len(PTMS):
        raise ValueError("self-clustering reduction requires one shard per PTM")
    validate_receipts(receipts, len(PTMS))
    config = load_config(config_path)
    permutations = int(config.values["randomness"]["bootstrap_draws"])
    configured_seed = int(config.values["randomness"]["archived_clustering_seed"])
    if configured_seed != BASE_SEED:
        raise ValueError(f"archived clustering seed must remain {BASE_SEED}")
    frames = [pd.read_csv(path) for path in shards]
    # Validate with an order-neutral table, then reconstruct the monolithic
    # loop order so floating-point reductions remain bit-for-bit comparable.
    validate_complete_partition(frames, "row_id")
    by_ptm = {str(frame.ptm.iloc[0]): frame for frame in frames}
    per_protein = pd.concat(
        [
            by_ptm[ptm].loc[by_ptm[ptm].universe.eq(universe)]
            for universe in ("all_eligible_residues", "observed_residue_union")
            for ptm in PTMS
        ],
        ignore_index=True,
    ).drop(columns="row_id")
    summary = summarize(per_protein)
    contrasts = direct_contrasts(per_protein, permutations)
    fold_ratios = matched_fold_ratios(per_protein, permutations)
    _, _, provenance = load_inputs(config.source_root)
    script_path = Path(__file__).resolve()
    provenance.update(
        {
            "script": str(script_path.relative_to(config.project_root)),
            "script_sha256": sha256(script_path),
            "generated_utc": datetime.now(UTC).isoformat(),
            "config": {
                "ptms": PTMS,
                "radii": RADII,
                "primary_radius": PRIMARY_RADIUS,
                "minimum_sites_per_protein": MIN_SITES,
                "permutations": permutations,
                "base_seed": BASE_SEED,
                "weighting": "equal_protein",
                "execution": "one deterministic map job per PTM; validated reduce",
            },
            "environment": {
                "python": platform.python_version(),
                "platform": platform.platform(),
                "numpy": np.__version__,
                "pandas": pd.__version__,
                "scipy": scipy.__version__,
            },
        }
    )
    output.mkdir(parents=True, exist_ok=True)
    write_csv(output / "per_protein.csv", per_protein)
    write_csv(output / "summary.csv", summary)
    write_csv(output / "contrasts.csv", contrasts)
    write_csv(output / "fold_ratios.csv", fold_ratios)
    write_json(output / "provenance.json", provenance)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    ptm_parser = commands.add_parser("ptm")
    ptm_parser.add_argument("--config", type=Path, required=True)
    ptm_parser.add_argument("--ptm", required=True)
    ptm_parser.add_argument("--output", type=Path, required=True)
    ptm_parser.add_argument("--receipt", type=Path, required=True)
    reduce_parser = commands.add_parser("reduce")
    reduce_parser.add_argument("--config", type=Path, required=True)
    reduce_parser.add_argument("--shards", type=Path, nargs="+", required=True)
    reduce_parser.add_argument("--receipts", type=Path, nargs="+", required=True)
    reduce_parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    if args.command == "ptm":
        run_ptm(args.config, args.ptm, args.output, args.receipt)
    else:
        reduce(args.config, args.shards, args.receipts, args.output_dir)


if __name__ == "__main__":
    main()
