"""Resolve the complete production DAG without accessing provider data."""

import json
import shutil
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]


def test_full_workflow_resolves_with_declared_inputs(tmp_path):
    for folder in ("workflow", "config"):
        shutil.copytree(ROOT / folder, tmp_path / folder)
    for name in ("src", "reference", "assets", "pyproject.toml", "uv.lock", "Makefile"):
        (tmp_path / name).symlink_to(ROOT / name)
    manifest = json.loads((ROOT / "reference/input_manifest.json").read_text())
    for record in manifest["files"].values():
        path = tmp_path / "inputs" / record["path"]
        path.parent.mkdir(parents=True, exist_ok=True)
        path.touch()
    for record in manifest["directories"].values():
        path = tmp_path / "inputs" / record["path"] / record["pattern"].replace("*", "fixture")
        path.parent.mkdir(parents=True, exist_ok=True)
        path.touch()
    # Empty files are for graph resolution only; the production preflight would reject them.
    result = subprocess.run(
        [
            sys.executable,
            "-m",
            "snakemake",
            "-s",
            "workflow/Snakefile",
            "--dry-run",
            "--cores",
            "2",
            "all",
        ],
        cwd=tmp_path,
        text=True,
        capture_output=True,
    )
    assert result.returncode == 0, result.stdout + result.stderr
    for rule in ("scanner_fold", "audit_report", "audit", "reviewer_report"):
        assert rule in result.stdout
    assert not (tmp_path / "results").exists()
