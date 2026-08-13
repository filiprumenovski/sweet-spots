"""Record the immutable analysis-ready inputs used by the paper workflow."""

from __future__ import annotations

import argparse
import platform
from pathlib import Path

import numpy
import pandas
import polars
import scipy
import sklearn

from ..core.config import load_config
from ..core.io import sha256_file, write_json


def directory_digest(path: Path, pattern: str) -> tuple[str, int, int]:
    """Hash relative names and file hashes for a deterministic directory digest."""
    import hashlib

    digest = hashlib.sha256()
    files = sorted(item for item in path.glob(pattern) if item.is_file())
    for item in files:
        record = f"{item.relative_to(path)}\t{sha256_file(item)}\n".encode()
        digest.update(record)
    return digest.hexdigest(), len(files), sum(item.stat().st_size for item in files)


def build_manifest(config_path: Path) -> dict[str, object]:
    config = load_config(config_path)
    source = config.source_root
    files = {
        "canonical_human_sequences": source / "data/interim/fasta_human.parquet",
        "canonical_multispecies_sequences": source / "data/interim/fasta_all.parquet",
        "residue_disorder": source / "data/interim/iupred_residue_scores.parquet",
        "oglcnac_atlas_strict": source / "analysis/revalidation/data/atlas_unambiguous.csv",
        "unified_ptm_catalogue": source / "data/processed/landscape/ptm_unified.parquet",
        "reference_consensus_regions": (
            source / "data/processed/regions/oglcnac_consensus_regions.parquet"
        ),
        "rice_sequences": (source / "data/external/multispecies_oglcnac/rice_sequences.fasta"),
        "arabidopsis_sequences": (
            source / "data/external/multispecies_oglcnac/arabidopsis_sequences.fasta"
        ),
        "drosophila_sequences": (
            source / "data/external/multispecies_oglcnac/drosophila_sequences.fasta"
        ),
        "celegans_sequences": (
            source / "data/external/multispecies_oglcnac/celegans_sequences.fasta"
        ),
        "ogt_human_sequence": (
            source / "analysis/revalidation/data/ogt_orthologs/O15294.fasta"
        ),
        "ogt_ortholog_alignment": (
            source / "analysis/revalidation/data/ogt_orthologs/orthologs_aln.fasta"
        ),
    }
    missing = [str(path) for path in files.values() if not path.is_file()]
    if missing:
        raise FileNotFoundError("Missing required paper input(s):\n" + "\n".join(missing))

    msa_digest, msa_files, msa_bytes = directory_digest(source / "data/interim/msa", "*.afa")
    structure_digest, structure_files, structure_bytes = directory_digest(
        source / "analysis/revalidation/data/ogt_structures", "*.cif"
    )
    if not msa_files:
        raise FileNotFoundError("No human ortholog alignments matched data/interim/msa/*.afa")
    if not structure_files:
        raise FileNotFoundError(
            "No OGT structures matched analysis/revalidation/data/ogt_structures/*.cif"
        )
    return {
        "contract": "analysis-ready inputs; no network access is permitted by the DAG",
        "files": {
            name: {
                "path": str(path.relative_to(source)),
                "bytes": path.stat().st_size,
                "sha256": sha256_file(path),
            }
            for name, path in files.items()
        },
        "directories": {
            "human_ortholog_alignments": {
                "path": "data/interim/msa",
                "pattern": "*.afa",
                "files": msa_files,
                "bytes": msa_bytes,
                "tree_sha256": msa_digest,
            },
            "ogt_structures": {
                "path": "analysis/revalidation/data/ogt_structures",
                "pattern": "*.cif",
                "files": structure_files,
                "bytes": structure_bytes,
                "tree_sha256": structure_digest,
            },
        },
        "runtime": {
            "python": platform.python_version(),
            "numpy": numpy.__version__,
            "pandas": pandas.__version__,
            "polars": polars.__version__,
            "scipy": scipy.__version__,
            "scikit_learn": sklearn.__version__,
        },
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    write_json(args.output, build_manifest(args.config))


if __name__ == "__main__":
    main()
