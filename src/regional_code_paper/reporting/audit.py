"""Audit numerical reproduction and figure completeness."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import pandas as pd

from ..core.config import load_config
from ..core.io import sha256_file, write_json

EXPECTED = {
    "object.regions": 824,
    "object.region_sites": 4196,
    "clustering.r10": 4.073557885388373,
    "regional_code.full.auroc": 0.8788975288086618,
    "regional_code.within_protein.auroc": 0.8743789354729757,
    "regional_code.composition_within_protein.auroc": 0.792152662273436,
    "kinome.delta_n90": -5.291756173032191,
    "evolution.worm.gap": 0.401840490797546,
    "transfer.human_to_target.rice": 0.8065612089584693,
    "fg_nup.recovery_fraction": 0.8,
}


def numerical_checks(metrics: pd.DataFrame) -> list[dict[str, object]]:
    indexed = metrics.set_index("claim_id")
    checks = []
    for claim, expected in EXPECTED.items():
        observed = float(indexed.at[claim, "value"]) if claim in indexed.index else math.nan
        tolerance = max(1e-10, abs(expected) * 1e-9)
        checks.append(
            {
                "claim_id": claim,
                "expected": expected,
                "observed": observed,
                "absolute_error": abs(observed - expected),
                "passed": math.isfinite(observed) and abs(observed - expected) <= tolerance,
            }
        )
    return checks


def figure_checks(manifest_path: Path, project_root: Path) -> list[dict[str, object]]:
    with manifest_path.open(encoding="utf-8") as handle:
        manifest = json.load(handle)
    checks = []
    for item in manifest["figures"]:
        for file in item["files"]:
            recorded_path = Path(file["path"])
            path = (
                recorded_path if recorded_path.is_absolute() else project_root / recorded_path
            )
            checks.append(
                {
                    "figure": item["figure"],
                    "path": str(recorded_path),
                    "exists": path.exists(),
                    "nonempty": path.exists() and path.stat().st_size > 1_000,
                    "hash_matches_manifest": path.exists()
                    and sha256_file(path) == file["sha256"],
                }
            )
    return checks


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--metrics", type=Path, required=True)
    parser.add_argument("--figure-manifest", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    # The CSV is the Python boundary; its Parquet sibling is created solely by
    # the DuckDB CLI rule and is the distribution artifact.
    metrics_csv = args.metrics.with_suffix(".csv")
    metrics = pd.read_csv(metrics_csv)
    numbers = numerical_checks(metrics)
    config = load_config(args.config)
    figures = figure_checks(args.figure_manifest, config.project_root)
    computational_pass = all(item["passed"] for item in numbers) and all(
        item["exists"] and item["nonempty"] and item["hash_matches_manifest"]
        for item in figures
    )
    write_json(
        args.output.resolve(),
        {
            "status": "reproduced" if computational_pass else "failed",
            "computational_pass": computational_pass,
            "numerical_checks": numbers,
            "figure_checks": figures,
        },
    )
    if not computational_pass:
        raise SystemExit("One or more computational reproduction checks failed")


if __name__ == "__main__":
    main()
