"""Typed loading and path resolution for the paper workflow configuration."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any

import yaml


@dataclass(frozen=True)
class PaperConfig:
    """Resolved workflow configuration."""

    project_root: Path
    source_root: Path
    results_root: Path
    values: dict[str, Any]


def load_config(path: Path) -> PaperConfig:
    """Load YAML once and resolve every path against the project directory."""
    path = path.resolve()
    values = yaml.safe_load(path.read_text(encoding="utf-8"))
    project_root = path.parent.parent

    def resolve(value: str) -> Path:
        candidate = Path(value)
        return (
            candidate.resolve()
            if candidate.is_absolute()
            else (project_root / candidate).resolve()
        )

    return PaperConfig(
        project_root=project_root,
        source_root=resolve(values["source_root"]),
        results_root=resolve(values["results_root"]),
        values=values,
    )
