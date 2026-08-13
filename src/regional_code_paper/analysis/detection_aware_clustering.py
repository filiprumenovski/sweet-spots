"""Peptide-exposure-conditioned null for O-GlcNAc self-clustering.

Mass-spectrometry site catalogues do not expose every serine and threonine
equally. This module constructs exact, non-overlapping opportunity strata from
the observed peptide intervals. Two residues share a stratum only when they
have the same residue identity and are covered by the same set of observed
peptide intervals. Sampling within those strata therefore preserves, for every
protein, the site count, serine/threonine composition, and peptide-detection
opportunity without collisions or dropped sites.
"""

from __future__ import annotations

import argparse
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd

from ..core.config import load_config
from ..core.io import write_csv, write_json

MINIMUM_SITES = 3
RADII = (5, 10, 15, 20, 25, 30)


def tryptic_span(sequence: str, position: int) -> tuple[int, int]:
    """Return the one-based theoretical tryptic interval containing a position."""
    if not 1 <= position <= len(sequence):
        raise ValueError("position is outside the sequence")
    start = 1
    for index, residue in enumerate(sequence, start=1):
        next_residue = sequence[index] if index < len(sequence) else None
        cleaves = residue in "KR" and next_residue != "P"
        if cleaves:
            if start <= position <= index:
                return start, index
            start = index + 1
    return start, len(sequence)


def mapped_spans(sequence: str, peptide: object, position: int) -> set[tuple[int, int]]:
    """Find every exact peptide occurrence that contains the reported site."""
    if peptide is None or pd.isna(peptide):
        return set()
    token = str(peptide).strip().upper()
    if not token or any(residue not in "ACDEFGHIKLMNPQRSTVWY" for residue in token):
        return set()
    output: set[tuple[int, int]] = set()
    offset = sequence.find(token)
    while offset >= 0:
        span = offset + 1, offset + len(token)
        if span[0] <= position <= span[1]:
            output.add(span)
        offset = sequence.find(token, offset + 1)
    return output


def load_sequences(path: Path) -> dict[str, str]:
    frame = pd.read_parquet(
        path,
        columns=["accession", "sequence", "taxon_id", "is_canonical"],
    )
    frame = frame.loc[(frame.taxon_id == 9606) & frame.is_canonical]
    frame = frame.drop_duplicates("accession")
    return dict(zip(frame.accession.astype(str), frame.sequence.astype(str), strict=True))


def validated_site_spans(
    atlas_path: Path, sequences: dict[str, str]
) -> tuple[
    dict[tuple[str, int], set[tuple[int, int]]],
    set[tuple[str, int]],
    dict[str, int],
]:
    """Validate atlas rows and collect mapped peptide intervals per unique site."""
    atlas = pd.read_csv(
        atlas_path,
        usecols=[
            "species",
            "accession",
            "position_in_protein",
            "site_residue",
            "peptide_seq",
        ],
        low_memory=False,
        encoding="utf-8",
        encoding_errors="replace",
    )
    atlas = atlas.loc[atlas.species.eq("human") & atlas.site_residue.isin(["S", "T"])].copy()
    atlas["position"] = pd.to_numeric(atlas.position_in_protein, errors="coerce")
    atlas = atlas.dropna(subset=["position"])

    spans: dict[tuple[str, int], set[tuple[int, int]]] = defaultdict(set)
    valid_sites: set[tuple[str, int]] = set()
    audit = {
        "input_rows": len(atlas),
        "valid_rows": 0,
        "invalid_rows": 0,
        "valid_unique_sites": 0,
        "sites_with_mapped_peptide": 0,
        "sites_with_tryptic_fallback": 0,
    }
    for accession, raw_position, residue, peptide in atlas[
        ["accession", "position", "site_residue", "peptide_seq"]
    ].itertuples(index=False, name=None):
        sequence = sequences.get(str(accession))
        position = int(raw_position)
        if (
            sequence is None
            or not 1 <= position <= len(sequence)
            or sequence[position - 1] != residue
        ):
            audit["invalid_rows"] += 1
            continue
        audit["valid_rows"] += 1
        key = str(accession), position
        valid_sites.add(key)
        spans[key].update(mapped_spans(sequence, peptide, position))

    mapped_sites = {key for key in valid_sites if spans[key]}
    for accession, position in sorted(valid_sites):
        key = accession, position
        if spans[key]:
            audit["sites_with_mapped_peptide"] += 1
        else:
            spans[key].add(tryptic_span(sequences[accession], position))
            audit["sites_with_tryptic_fallback"] += 1
    audit["valid_unique_sites"] = len(valid_sites)
    return dict(spans), mapped_sites, audit


def exposure_tables(
    sequences: dict[str, str],
    site_spans: dict[tuple[str, int], set[tuple[int, int]]],
    mapped_sites: set[tuple[str, int]],
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Build disjoint residue-by-peptide-coverage strata for each protein."""
    sites_by_protein: dict[str, list[int]] = defaultdict(list)
    for accession, position in site_spans:
        sites_by_protein[accession].append(position)

    protein_rows: list[dict[str, object]] = []
    stratum_rows: list[dict[str, object]] = []
    for accession, observed_list in sorted(sites_by_protein.items()):
        observed = np.asarray(sorted(set(observed_list)), dtype=np.int64)
        if len(observed) < MINIMUM_SITES:
            continue
        sequence = sequences[accession]
        intervals = sorted(
            {span for position in observed for span in site_spans[(accession, int(position))]}
        )

        eligible_by_key: dict[tuple[str, tuple[int, ...]], list[int]] = defaultdict(list)
        for position, residue in enumerate(sequence, start=1):
            if residue not in "ST":
                continue
            signature = tuple(
                index
                for index, (start, end) in enumerate(intervals)
                if start <= position <= end
            )
            if signature:
                eligible_by_key[(residue, signature)].append(position)

        observed_by_key: dict[tuple[str, tuple[int, ...]], list[int]] = defaultdict(list)
        for position in observed:
            residue = sequence[int(position) - 1]
            signature = tuple(
                index
                for index, (start, end) in enumerate(intervals)
                if start <= position <= end
            )
            if not signature:
                raise AssertionError(f"observed site lacks an exposure interval: {accession}")
            observed_by_key[(residue, signature)].append(int(position))

        variable_sites = 0
        for stratum_index, (key, observed_positions) in enumerate(
            sorted(observed_by_key.items(), key=lambda item: item[0])
        ):
            residue, signature = key
            eligible = sorted(eligible_by_key[key])
            if len(eligible) < len(observed_positions):
                raise AssertionError(f"invalid exposure stratum: {accession}, {key}")
            if len(eligible) > len(observed_positions):
                variable_sites += len(observed_positions)
            signature_label = "|".join(
                f"{intervals[index][0]}-{intervals[index][1]}" for index in signature
            )
            stratum_rows.append(
                {
                    "accession": accession,
                    "stratum": stratum_index,
                    "residue": residue,
                    "peptide_coverage_signature": signature_label,
                    "n_observed": len(observed_positions),
                    "n_eligible": len(eligible),
                    "eligible_positions": ";".join(map(str, eligible)),
                }
            )

        mapped_count = sum((accession, int(position)) in mapped_sites for position in observed)
        protein_rows.append(
            {
                "accession": accession,
                "n_sites": len(observed),
                "observed_positions": ";".join(map(str, observed)),
                "n_exposure_strata": len(observed_by_key),
                "n_variable_sites": variable_sites,
                "n_mapped_peptide_sites": mapped_count,
                "n_tryptic_fallback_sites": len(observed) - mapped_count,
            }
        )

    proteins = pd.DataFrame(protein_rows).sort_values("accession").reset_index(drop=True)
    strata = pd.DataFrame(stratum_rows).sort_values(["accession", "stratum"])
    if int(proteins.n_sites.sum()) != int(strata.n_observed.sum()):
        raise AssertionError("exposure strata do not preserve the observed site count")
    return proteins, strata.reset_index(drop=True)


def prepare(config_path: Path, output_dir: Path) -> None:
    """Materialize the reviewable detection-opportunity cache."""
    config = load_config(config_path)
    fasta_path = config.source_root / "data/interim/fasta_human.parquet"
    atlas_path = config.source_root / "analysis/revalidation/data/atlas_unambiguous.csv"
    sequences = load_sequences(fasta_path)
    site_spans, mapped_sites, audit = validated_site_spans(atlas_path, sequences)
    proteins, strata = exposure_tables(sequences, site_spans, mapped_sites)
    write_csv(output_dir / "proteins.csv", proteins)
    write_csv(output_dir / "strata.csv", strata)
    write_json(
        output_dir / "summary.json",
        {
            "audit": audit,
            "analysis_population": {
                "proteins": len(proteins),
                "sites": int(proteins.n_sites.sum()),
                "exposure_strata": len(strata),
                "variable_sites": int(proteins.n_variable_sites.sum()),
                "mapped_peptide_sites": int(proteins.n_mapped_peptide_sites.sum()),
                "tryptic_fallback_sites": int(proteins.n_tryptic_fallback_sites.sum()),
            },
            "stratum_definition": (
                "same protein, residue identity, and exact observed peptide-coverage signature"
            ),
            "site_count_preserved": int(proteins.n_sites.sum()) == int(strata.n_observed.sum()),
        },
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    prepare(args.config, args.output_dir.resolve())


if __name__ == "__main__":
    main()
