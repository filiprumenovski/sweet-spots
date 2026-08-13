"""Platform-independent random seed derivation."""

from __future__ import annotations

import hashlib


def stable_seed(base_seed: int, *tokens: object) -> int:
    """Derive a NumPy-compatible seed for one named analysis."""
    label = "|".join(str(token) for token in tokens).encode("utf-8")
    offset = int(hashlib.sha256(label).hexdigest()[:8], 16)
    return (base_seed + offset) % (2**32)
