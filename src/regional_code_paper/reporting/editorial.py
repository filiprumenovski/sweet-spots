"""Verify retained editorial artwork and its links to reproduced source products."""

from __future__ import annotations

import json
import math
from pathlib import Path

import pandas as pd

from ..core.io import sha256_file
from .audit import FIGURES, table_checks


def equal_values(left, right) -> bool:
    """Compare complete JSON records with the scientific audit's tolerances."""
    if isinstance(left, dict) and isinstance(right, dict):
        return left.keys() == right.keys() and all(
            equal_values(left[k], right[k]) for k in left
        )
    if isinstance(left, list) and isinstance(right, list):
        return len(left) == len(right) and all(
            equal_values(a, b) for a, b in zip(left, right, strict=True)
        )
    if isinstance(left, bool) or isinstance(right, bool):
        return type(left) is type(right) and left == right
    if isinstance(left, int | float) and isinstance(right, int | float):
        return math.isclose(left, right, rel_tol=1e-9, abs_tol=1e-10)
    return left == right


def editorial_checks(project: Path, results: Path) -> list[dict]:
    directory = project / "assets/review_figures"
    manifest = json.loads((directory / "manifest.json").read_text())
    expected_panels = dict(zip(FIGURES, [6, 5, 3, 6, 4, 5, 3, 5, 4], strict=True))
    figures = manifest["figures"]
    expected_assets = {
        f"assets/review_figures/{figure}.{extension}"
        for figure in FIGURES
        for extension in ("pdf", "png", "svg")
    }
    checks = [
        {
            "check": "editorial inventory",
            "passed": sorted(f["figure"] for f in figures) == sorted(FIGURES)
            and expected_assets <= manifest["files"].keys()
            and all(
                [p["panel"] for p in f["panels"]]
                == list("ABCDEF"[: expected_panels[f["figure"]]])
                for f in figures
            )
            and all(p["sources"] for f in figures for p in f["panels"])
            and {s for f in figures for p in f["panels"] for s in p["sources"]}
            <= manifest["sources"].keys(),
        }
    ]
    for path, digest in manifest["files"].items():
        file = project / path
        checks.append(
            {
                "path": path,
                "method": "asset SHA-256",
                "passed": file.is_file() and sha256_file(file) == digest,
            }
        )
    for path, digest in manifest["sources"].items():
        file = (
            results / path.removeprefix("results/")
            if path.startswith("results/")
            else project / path
        )
        same = file.is_file() and sha256_file(file) == digest
        passed, method = same, "source SHA-256"
        reference = directory / "reference" / path.removeprefix("results/analysis/")
        if not same and file.is_file() and reference.is_file():
            method = "complete source values (rtol=1e-9, atol=1e-10)"
            if file.suffix == ".json":
                passed = equal_values(
                    json.loads(file.read_text()), json.loads(reference.read_text())
                )
            elif file.suffix == ".csv":
                actual, expected = (
                    pd.read_csv(file).reset_index(),
                    pd.read_csv(reference).reset_index(),
                )
                passed = all(c["passed"] for c in table_checks(actual, expected, "index"))
        checks.append(
            {"path": path, "method": method, "passed": passed, "byte_identical": same}
        )
    return checks
