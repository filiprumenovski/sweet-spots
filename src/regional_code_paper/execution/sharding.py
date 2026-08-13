"""Deterministic map/reduce primitives shared by expensive analyses."""

from __future__ import annotations

import json
from collections.abc import Iterable, Sequence
from dataclasses import dataclass
from pathlib import Path

import pandas as pd

from ..core.io import sha256_file, write_json


@dataclass(frozen=True)
class Shard:
    """One zero-indexed member of a fixed-size partition."""

    index: int
    count: int

    def __post_init__(self) -> None:
        if self.count < 1:
            raise ValueError("shard count must be positive")
        if not 0 <= self.index < self.count:
            raise ValueError(f"shard index {self.index} is outside [0, {self.count})")

    def owns(self, global_index: int) -> bool:
        """Assign sorted records round-robin without Python's salted hash."""
        return global_index % self.count == self.index


def partition_sorted(values: Iterable[str], shard: Shard) -> list[str]:
    """Return the stable, disjoint subset owned by ``shard``."""
    return [value for index, value in enumerate(sorted(values)) if shard.owns(index)]


def validate_complete_partition(frames: Sequence[pd.DataFrame], key: str) -> pd.DataFrame:
    """Concatenate shards and reject absent, null, or duplicate record keys."""
    if not frames:
        raise ValueError("at least one shard is required")
    columns = list(frames[0].columns)
    for index, frame in enumerate(frames):
        if list(frame.columns) != columns:
            raise ValueError(f"shard {index} has a different schema")
        if key not in frame:
            raise ValueError(f"partition key {key!r} is absent from shard {index}")
        if frame[key].isna().any():
            raise ValueError(f"shard {index} contains null {key!r} values")
    combined = pd.concat(frames, ignore_index=True)
    duplicates = combined[key].duplicated(keep=False)
    if duplicates.any():
        examples = combined.loc[duplicates, key].astype(str).head(5).tolist()
        raise ValueError(f"duplicate partition keys: {examples}")
    return combined.sort_values(key, kind="stable").reset_index(drop=True)


def write_receipt(
    path: Path,
    *,
    shard: Shard,
    outputs: Sequence[Path],
    records: int,
    metadata: dict[str, object] | None = None,
) -> None:
    """Write a machine-verifiable completion receipt after all shard outputs."""
    missing = [str(output) for output in outputs if not output.is_file()]
    if missing:
        raise FileNotFoundError(f"cannot receipt missing outputs: {missing}")
    payload: dict[str, object] = {
        "schema_version": 1,
        "shard": shard.index,
        "shards": shard.count,
        "records": records,
        "outputs": {
            str(output): {"bytes": output.stat().st_size, "sha256": sha256_file(output)}
            for output in outputs
        },
    }
    if metadata:
        payload["metadata"] = metadata
    write_json(path, payload)


def validate_receipts(paths: Sequence[Path], expected_shards: int) -> list[dict]:
    """Verify shard coverage and every byte-count/digest recorded by workers."""
    documents = [json.loads(path.read_text(encoding="utf-8")) for path in paths]
    if len(documents) != expected_shards:
        raise ValueError(f"expected {expected_shards} receipts, received {len(documents)}")
    observed = sorted(int(document["shard"]) for document in documents)
    if observed != list(range(expected_shards)):
        raise ValueError(f"receipts cover shards {observed}, not 0..{expected_shards - 1}")
    for document in documents:
        if int(document["shards"]) != expected_shards:
            raise ValueError("receipt was created under a different shard count")
        for name, claimed in document["outputs"].items():
            output = Path(name)
            if not output.is_file():
                raise FileNotFoundError(f"receipted output is missing: {output}")
            if output.stat().st_size != int(claimed["bytes"]):
                raise ValueError(f"byte-count mismatch for {output}")
            if sha256_file(output) != claimed["sha256"]:
                raise ValueError(f"SHA-256 mismatch for {output}")
    return documents
