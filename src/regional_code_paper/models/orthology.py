"""Alignment and multispecies helpers shared by evolutionary analyses."""

from __future__ import annotations

from collections import Counter, defaultdict
from collections.abc import Iterable
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd

from ..analysis.consensus_regions import consensus_components

TAXON_TO_SPECIES = {
    "10090": "mouse",
    "10116": "rat",
    "7955": "zebrafish",
    "7227": "fly",
    "6239": "worm",
}
SPECIES_TO_TAXON = {value: int(key) for key, value in TAXON_TO_SPECIES.items()}
AMINO_ACIDS = "ACDEFGHIKLMNPQRSTVWY"


@dataclass(frozen=True)
class Alignment:
    human: str
    human_position_to_column: dict[int, int]
    orthologs: dict[str, tuple[str, str]]


@dataclass(frozen=True)
class EvolutionInputs:
    sequences_by_taxon: dict[tuple[int, str], str]
    genes_by_taxon: dict[tuple[int, str], str]
    human_sequences: dict[str, str]
    species_sites: dict[str, dict[str, set[int]]]
    spans: dict[str, list[tuple[int, int]]]
    human_region_sites: dict[str, set[int]]
    alignments: dict[str, Alignment]
    disorder: dict[tuple[str, int], float]
    external_sequences: dict[str, str]


def parse_fasta(path: Path) -> dict[str, str]:
    records: dict[str, str] = {}
    header: str | None = None
    sequence: list[str] = []
    for raw_line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        line = raw_line.strip()
        if line.startswith(">"):
            if header is not None:
                records[header] = "".join(sequence)
            header = line[1:].split()[0]
            sequence = []
        else:
            sequence.append(line)
    if header is not None:
        records[header] = "".join(sequence)
    return records


def position_to_column(aligned_sequence: str) -> dict[int, int]:
    output: dict[int, int] = {}
    position = 0
    for column, residue in enumerate(aligned_sequence):
        if residue != "-":
            position += 1
            output[position] = column
    return output


def load_alignments(msa_directory: Path, accessions: Iterable[str]) -> dict[str, Alignment]:
    """Load the exact accession-named MSA products used by the archived analysis.

    The source directory also contains hashed cache variants. The published run
    selected accession-named files only; reproducing that selection is therefore
    part of the input contract, not an arbitrary glob.
    """
    output = {}
    for accession in sorted(set(accessions)):
        path = msa_directory / f"{accession}.afa"
        if not path.is_file():
            continue
        records = parse_fasta(path)
        human = records.get("human") or records.get(accession)
        if not human:
            continue
        orthologs = {}
        for header, sequence in records.items():
            taxon, separator, ortholog_accession = header.partition("_")
            if separator and taxon in TAXON_TO_SPECIES:
                orthologs[taxon] = (ortholog_accession, sequence)
        if orthologs:
            output[accession] = Alignment(human, position_to_column(human), orthologs)
    return output


def load_evolution_inputs(
    source_root: Path, regions_path: Path, sites_path: Path
) -> EvolutionInputs:
    fasta = pd.read_parquet(
        source_root / "data/interim/fasta_all.parquet",
        columns=["taxon_id", "accession", "sequence", "gene", "is_canonical"],
    )
    fasta = fasta.loc[fasta.is_canonical].drop_duplicates(["taxon_id", "accession"])
    sequences_by_taxon = {
        (int(taxon), str(accession)): str(sequence)
        for taxon, accession, sequence in fasta[
            ["taxon_id", "accession", "sequence"]
        ].itertuples(index=False, name=None)
    }
    genes_by_taxon = {
        (int(taxon), str(accession)): "" if pd.isna(gene) else str(gene).upper()
        for taxon, accession, gene in fasta[["taxon_id", "accession", "gene"]].itertuples(
            index=False, name=None
        )
    }
    human_sequences = {
        accession: sequence
        for (taxon, accession), sequence in sequences_by_taxon.items()
        if taxon == 9606
    }

    regions = pd.read_csv(regions_path)
    spans: dict[str, list[tuple[int, int]]] = defaultdict(list)
    for accession, start, end in regions[["accession", "start", "end"]].itertuples(
        index=False, name=None
    ):
        spans[str(accession)].append((int(start), int(end)))
    region_sites = pd.read_csv(sites_path, usecols=["accession", "position"])
    human_region_sites: dict[str, set[int]] = defaultdict(set)
    for accession, position in region_sites.itertuples(index=False, name=None):
        human_region_sites[str(accession)].add(int(position))

    atlas = pd.read_csv(
        source_root / "analysis/revalidation/data/atlas_unambiguous.csv",
        usecols=["species", "accession", "position_in_protein", "site_residue"],
        low_memory=False,
        encoding_errors="replace",
    )
    atlas = atlas.loc[atlas.site_residue.isin(["S", "T"])].copy()
    atlas["position"] = pd.to_numeric(atlas.position_in_protein, errors="coerce")
    atlas = atlas.dropna(subset=["position"])
    species_sites: dict[str, dict[str, set[int]]] = defaultdict(lambda: defaultdict(set))
    for species, accession, position in (
        atlas[["species", "accession", "position"]]
        .drop_duplicates()
        .itertuples(index=False, name=None)
    ):
        species_sites[str(species)][str(accession)].add(int(position))

    disorder_frame = pd.read_parquet(
        source_root / "data/interim/iupred_residue_scores.parquet",
        columns=["accession", "position", "disorder_score"],
    )
    disorder = {
        (str(accession), int(position)): float(score)
        for accession, position, score in disorder_frame.itertuples(index=False, name=None)
    }
    external_sequences: dict[str, str] = {}
    external_directory = source_root / "data/external/multispecies_oglcnac"
    for name in (
        "rice_sequences.fasta",
        "arabidopsis_sequences.fasta",
        "drosophila_sequences.fasta",
        "celegans_sequences.fasta",
    ):
        external_sequences.update(parse_fasta(external_directory / name))

    return EvolutionInputs(
        sequences_by_taxon=sequences_by_taxon,
        genes_by_taxon=genes_by_taxon,
        human_sequences=human_sequences,
        species_sites={species: dict(values) for species, values in species_sites.items()},
        spans=dict(spans),
        human_region_sites=dict(human_region_sites),
        alignments=load_alignments(source_root / "data/interim/msa", spans),
        disorder=disorder,
        external_sequences=external_sequences,
    )


def species_consensus_spans(positions: Iterable[int]) -> list[tuple[int, int]]:
    components = consensus_components(
        positions,
        gaps=(5, 8, 10, 12, 15),
        final_gap=10,
        minimum_sites=3,
    )
    return [(component[0], component[-1]) for component in components]


def sequence_features(sequence: str, center: int) -> list[float]:
    """Alignment-free 24-feature representation at radii 10, 20 and 40."""
    output = []
    for radius in (10, 20, 40):
        segment = sequence[max(0, center - 1 - radius) : min(len(sequence), center + radius)]
        length = max(1, len(segment))
        counts = Counter(segment)
        output.extend(
            [
                (counts["S"] + counts["T"]) / length,
                (counts["D"] + counts["E"]) / length,
                (counts["F"] + counts["W"] + counts["Y"]) / length,
                counts["A"] / length,
                counts["P"] / length,
                (counts["Q"] + counts["N"]) / length,
                (counts["K"] + counts["R"] - counts["D"] - counts["E"]) / length,
                sum(counts[residue] for residue in "PESQKRGDN") / length,
            ]
        )
    return output


def composition_vector(residues: Iterable[str]) -> np.ndarray:
    counts = Counter(residues)
    return np.asarray([counts[residue] for residue in AMINO_ACIDS], dtype=float)


def raw_pearson(left: np.ndarray, right: np.ndarray) -> float:
    if left.std() == 0 or right.std() == 0:
        return np.nan
    return float(np.corrcoef(left, right)[0, 1])


def cosine(left: np.ndarray, right: np.ndarray) -> float:
    if left.sum() == 0 or right.sum() == 0:
        return np.nan
    p, q = left / left.sum(), right / right.sum()
    denominator = np.linalg.norm(p) * np.linalg.norm(q)
    return float(p @ q / denominator) if denominator else np.nan


def one_minus_js(left: np.ndarray, right: np.ndarray) -> float:
    if left.sum() == 0 or right.sum() == 0:
        return np.nan
    p, q = left / left.sum(), right / right.sum()
    midpoint = 0.5 * (p + q)

    def divergence(values: np.ndarray, reference: np.ndarray) -> float:
        nonzero = values > 0
        return float(np.sum(values[nonzero] * np.log2(values[nonzero] / reference[nonzero])))

    return float(1.0 - 0.5 * (divergence(p, midpoint) + divergence(q, midpoint)))


def position_correlation(left: np.ndarray, right: np.ndarray) -> float | None:
    if left.std() == 0 or right.std() == 0:
        return None
    return float(np.corrcoef(left, right)[0, 1])
