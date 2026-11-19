"""Compare the complete scientific record with the immutable archived reference."""

from __future__ import annotations

import argparse
import gzip
import hashlib
import io
import json
import math
from pathlib import Path

import numpy as np
import pandas as pd

from ..core.config import load_config
from ..core.io import sha256_file, write_json
from .input_manifest import verify_manifest

FIGURES = [f"Figure_{name}" for name in ("1", "2", "3", "4", "5", "S1", "S2", "S3", "S4")]


def json_cell(value):
    if pd.isna(value):
        return None
    value = value.item() if hasattr(value, "item") else value
    return str(value) if isinstance(value, float) and not math.isfinite(value) else value


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
                    "expected": json_cell(value),
                    "observed": json_cell(actual),
                    "passed": bool(passed),
                }
            )
    return checks


def archived_table(reference: Path, relative: str, digest: str) -> pd.DataFrame | None:
    """Load the compressed archived copy of ``results/<relative>``.

    The decompressed bytes must carry the archived SHA-256, so the copy is provably the
    original x86 output rather than a regenerated file.
    """
    path = reference / "tables" / f"{relative}.gz"
    if not path.is_file():
        return None
    data = gzip.decompress(path.read_bytes())
    if hashlib.sha256(data).hexdigest() != digest:
        raise ValueError(f"{path} does not decompress to the archived {relative}")
    return pd.read_csv(io.BytesIO(data))


def frame_equivalence(observed: pd.DataFrame, expected: pd.DataFrame) -> tuple[bool, float]:
    """Same rows, columns and text; every number within the audit tolerances.

    Floating-point results can differ in the last bit between CPU architectures
    (for example x86 AVX2 and ARM64), so byte identity is recorded but not required.
    """
    if list(observed.columns) != list(expected.columns) or len(observed) != len(expected):
        return False, math.inf
    worst = 0.0
    for column in expected.columns:
        left, right = observed[column], expected[column]
        if pd.api.types.is_numeric_dtype(left) and pd.api.types.is_numeric_dtype(right):
            x, y = left.to_numpy(float), right.to_numpy(float)
            if not np.allclose(x, y, rtol=1e-9, atol=1e-10, equal_nan=True):
                return False, math.inf
            finite = np.isfinite(x) & np.isfinite(y)
            if finite.any():
                worst = max(worst, float(np.abs(x[finite] - y[finite]).max()))
        elif not left.equals(right):
            return False, math.inf
    return True, worst


def file_check(results: Path, reference: Path, relative: str, digest: str) -> dict:
    """Byte identity if possible, otherwise complete values against the archived copy."""
    observed = sha256_file(results / relative)
    check = {"path": relative, "expected": digest, "observed": observed,
             "byte_identical": observed == digest}
    if observed == digest:
        return {**check, "method": "SHA-256", "passed": True}
    archived = archived_table(reference, relative, digest)
    if archived is None:
        return {**check, "method": "SHA-256", "passed": False}
    passed, worst = frame_equivalence(pd.read_csv(results / relative), archived)
    return {**check, "method": "complete values (rtol=1e-9, atol=1e-10)",
            "max_abs_diff": worst if math.isfinite(worst) else None, "passed": passed}


def figure_checks(manifest_path: Path, project_root: Path) -> list[dict]:
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    expected = {
        str(manifest_path.parent.relative_to(project_root) / f"{name}.{extension}")
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
    from .editorial import editorial_checks

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
        file_check(results, reference, path, digest) for path, digest in fingerprints.items()
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
    checks["editorial"] = editorial_checks(project_root, results)
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
