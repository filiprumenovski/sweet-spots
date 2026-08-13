from __future__ import annotations

import pandas as pd
import pytest

from regional_code_paper.execution.sharding import (
    Shard,
    partition_sorted,
    validate_complete_partition,
    validate_receipts,
    write_receipt,
)


def test_partition_is_complete_disjoint_and_order_independent() -> None:
    values = ["delta", "alpha", "charlie", "bravo"]
    parts = [partition_sorted(reversed(values), Shard(index, 3)) for index in range(3)]
    flattened = [value for part in parts for value in part]
    assert sorted(flattened) == sorted(values)
    assert len(flattened) == len(set(flattened))


def test_reducer_rejects_duplicate_keys() -> None:
    shards = [pd.DataFrame({"row_id": [0]}), pd.DataFrame({"row_id": [0]})]
    with pytest.raises(ValueError, match="duplicate partition keys"):
        validate_complete_partition(shards, "row_id")


def test_invalid_shard_fails_early() -> None:
    with pytest.raises(ValueError, match="outside"):
        Shard(4, 4)


def test_receipt_detects_post_commit_tampering(tmp_path) -> None:
    output = tmp_path / "shard.csv"
    output.write_text("row_id\n1\n", encoding="utf-8")
    receipt = tmp_path / "receipt.json"
    write_receipt(receipt, shard=Shard(0, 1), outputs=[output], records=1)
    validate_receipts([receipt], 1)
    output.write_text("row_id\n2\n", encoding="utf-8")
    with pytest.raises(ValueError, match="SHA-256 mismatch"):
        validate_receipts([receipt], 1)
