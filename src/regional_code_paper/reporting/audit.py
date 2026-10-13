"""Compare the complete scientific record with the immutable archived reference."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import pandas as pd

from ..core.config import load_config
from ..core.io import sha256_file, write_json
from .input_manifest import verify_manifest

FIGURES = [f"Figure_{name}" for name in ("1", "2", "3", "4", "5", "S1", "S2", "S3", "S4")]


def table_checks(observed: pd.DataFrame, expected: pd.DataFrame, key: str) -> list[dict]:
    """Compare every cell, rejecting missing/extra/duplicate rows and columns."""
    checks = []
    inventory = (
        set(observed.columns) == set(expected.columns)
        and key in observed
        and not observed[key].duplicated().any()
        and set(observed[key]) == set(expected[key])
    )
    checks.append({"check": "table inventory", "passed": inventory})
    if not inventory:
        return checks
    observed = observed.set_index(key)
    expected = expected.set_index(key)
    for identity, row in expected.iterrows():
        for column, value in row.items():
            actual = observed.at[identity, column]
            if pd.isna(value):
                passed = bool(pd.isna(actual))
            elif pd.api.types.is_numeric_dtype(expected[column]):
                passed = bool(
                    pd.notna(actual)
                    and math.isclose(float(actual), float(value), rel_tol=1e-9, abs_tol=1e-10)
                )
            else:
                passed = actual == value
            checks.append(
                {
                    "row": str(identity),
                    "column": column,
                    "expected": None
                    if pd.isna(value)
                    else value.item()
                    if hasattr(value, "item")
                    else value,
                    "observed": None
                    if pd.isna(actual)
                    else actual.item()
                    if hasattr(actual, "item")
                    else actual,
                    "passed": bool(passed),
                }
            )
    return checks


def figure_checks(manifest_path: Path, project_root: Path) -> list[dict]:
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    expected = {
        f"results/figures/{name}.{extension}"
        for name in FIGURES
        for extension in ("pdf", "png")
    }
    files = [file for item in manifest["figures"] for file in item["files"]]
    checks = [
        {
            "check": "figure inventory",
            "passed": len(files) == len(expected)
            and {f["path"] for f in files} == expected
            and sorted(item["figure"] for item in manifest["figures"]) == sorted(FIGURES),
        }
    ]
    for file in files:
        path = project_root / file["path"]
        checks.append(
            {
                "path": file["path"],
                "passed": file["path"] in expected
                and path.is_file()
                and path.stat().st_size > 1000
                and sha256_file(path) == file["sha256"],
            }
        )
    return checks


def build_report(project_root: Path, results: Path, reference: Path) -> dict:
    checks = {}
    tables = {
        "claims": ("tables/manuscript_result_ledger.csv", "claim_id"),
        "catalogue_sensitivity": (
            "analysis/region_definition_sensitivity/catalogue_grid.csv",
            "definition_id",
        ),
        "model_sensitivity": (
            "analysis/region_definition_sensitivity/model_grid.csv",
            "definition_id",
        ),
    }
    for name, (relative, key) in tables.items():
        checks[name] = table_checks(
            pd.read_csv(results / relative), pd.read_csv(reference / Path(relative).name), key
        )
    fingerprints = json.loads((reference / "predictions.json").read_text())
    checks["predictions"] = [
        {
            "path": path,
            "expected": digest,
            "observed": sha256_file(results / path),
            "passed": sha256_file(results / path) == digest,
        }
        for path, digest in fingerprints.items()
    ]
    try:
        verify_manifest(
            json.loads((results / "provenance/input_manifest.json").read_text()),
            json.loads((reference / "input_manifest.json").read_text()),
        )
        checks["inputs"] = [{"passed": True}]
    except ValueError as error:
        checks["inputs"] = [{"passed": False, "error": str(error)}]
    checks["figures"] = figure_checks(results / "figures/figure_manifest.json", project_root)
    passed = all(check["passed"] for group in checks.values() for check in group)
    return {
        "status": "reproduced" if passed else "failed",
        "computational_pass": passed,
        "reference_commit": "c03b09bcb170767f8319bf08df18e1bedbed979e",
        "checks": checks,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    config = load_config(args.config)
    try:
        report = build_report(
            config.project_root, config.results_root, config.project_root / "reference"
        )
    except (OSError, ValueError, KeyError, TypeError) as error:
        report = {"status": "failed", "computational_pass": False, "error": str(error)}
    # The downstream gate fails; this successfully written report remains inspectable.
    write_json(args.output, report)


if __name__ == "__main__":
    main()
