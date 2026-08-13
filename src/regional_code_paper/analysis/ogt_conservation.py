"""Reconstruct OGT's peptide-reading surface and score its conservation.

No network access is used. Peptide contacts are recomputed from six cached OGT
complexes at a 5 Å heavy-atom cutoff. Surface controls are derived from monomeric
solvent-accessible surface area in 3PE3, and the 7YEA dimer interface is mapped
independently at the same cutoff.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from Bio.Align import PairwiseAligner, substitution_matrices
from Bio.PDB import MMCIFParser, ShrakeRupley

from ..core.config import load_config
from ..core.io import write_json
from ..models.orthology import parse_fasta, position_to_column

COMPLEXES = {
    "3PE4": ("A", "B"),
    "5HGV": ("A", "B"),
    "4N39": ("A", "B"),
    "4N3A": ("A", "B"),
    "4N3B": ("A", "B"),
    "4N3C": ("A", "B"),
}
CATALYTIC = {508, 568, 852}
MAXIMUM_ASA = {
    "A": 129,
    "R": 274,
    "N": 195,
    "D": 193,
    "C": 167,
    "Q": 225,
    "E": 223,
    "G": 104,
    "H": 224,
    "I": 197,
    "L": 201,
    "K": 236,
    "M": 224,
    "F": 240,
    "P": 159,
    "S": 155,
    "T": 172,
    "W": 285,
    "Y": 263,
    "V": 174,
}
THREE_TO_ONE = {
    "ALA": "A",
    "ARG": "R",
    "ASN": "N",
    "ASP": "D",
    "CYS": "C",
    "GLN": "Q",
    "GLU": "E",
    "GLY": "G",
    "HIS": "H",
    "ILE": "I",
    "LEU": "L",
    "LYS": "K",
    "MET": "M",
    "PHE": "F",
    "PRO": "P",
    "SER": "S",
    "THR": "T",
    "TRP": "W",
    "TYR": "Y",
    "VAL": "V",
}


def amino_acid_residues(chain) -> list:
    return [
        residue for residue in chain if residue.id[0] == " " and residue.resname in THREE_TO_ONE
    ]


def sequence_and_ids(chain) -> tuple[str, list[int]]:
    residues = amino_acid_residues(chain)
    return "".join(THREE_TO_ONE[residue.resname] for residue in residues), [
        residue.id[1] for residue in residues
    ]


def structure_to_uniprot(chain, human_sequence: str) -> dict[int, int]:
    sequence, residue_ids = sequence_and_ids(chain)
    aligner = PairwiseAligner()
    aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
    aligner.open_gap_score = -10.0
    aligner.extend_gap_score = -0.5
    aligner.mode = "global"
    alignment = aligner.align(human_sequence, sequence)[0]
    mapping = {}
    # `aligned` gives corresponding half-open blocks without parsing formatted text.
    human_blocks, structure_blocks = alignment.aligned
    for (human_start, human_end), (structure_start, structure_end) in zip(
        human_blocks, structure_blocks, strict=True
    ):
        for offset in range(min(human_end - human_start, structure_end - structure_start)):
            mapping[residue_ids[structure_start + offset]] = human_start + offset + 1
    return mapping


def minimum_distance(residue, atoms: np.ndarray) -> float:
    coordinates = np.asarray([atom.coord for atom in residue if atom.element != "H"])
    if not len(coordinates) or not len(atoms):
        return np.inf
    return float(
        np.sqrt(((coordinates[:, None, :] - atoms[None, :, :]) ** 2).sum(axis=2)).min()
    )


def peptide_contacts(structures: Path, human_sequence: str) -> set[int]:
    parser = MMCIFParser(QUIET=True)
    contacts = set()
    for pdb_id, (ogt_chain, peptide_chain) in COMPLEXES.items():
        model = parser.get_structure(pdb_id, structures / f"{pdb_id}.cif")[0]
        ogt = model[ogt_chain]
        peptide_atoms = np.asarray(
            [
                atom.coord
                for residue in amino_acid_residues(model[peptide_chain])
                for atom in residue
                if atom.element != "H"
            ]
        )
        mapping = structure_to_uniprot(ogt, human_sequence)
        for residue in amino_acid_residues(ogt):
            if minimum_distance(residue, peptide_atoms) <= 5.0 and residue.id[1] in mapping:
                contacts.add(mapping[residue.id[1]])
    return contacts


def surface_categories(structures: Path, human_sequence: str, reading_surface: set[int]):
    parser = MMCIFParser(QUIET=True)
    model = parser.get_structure("3PE3", structures / "3PE3.cif")[0]
    chain = model["A"]
    mapping = structure_to_uniprot(chain, human_sequence)
    ShrakeRupley(n_points=200).compute(chain, level="R")
    surface, buried = set(), set()
    for residue in amino_acid_residues(chain):
        position = mapping.get(residue.id[1])
        if position is None:
            continue
        residue_type = THREE_TO_ONE[residue.resname]
        relative_sasa = float(residue.sasa) / MAXIMUM_ASA[residue_type]
        (surface if relative_sasa >= 0.20 else buried).add(position)
    return surface - reading_surface, buried - reading_surface


def dimer_interface(structures: Path, human_sequence: str) -> set[int]:
    parser = MMCIFParser(QUIET=True)
    model = parser.get_structure("7YEA", structures / "7YEA.cif")[0]
    large_chains = sorted(
        [chain for chain in model if len(amino_acid_residues(chain)) >= 300],
        key=lambda chain: chain.id,
    )
    if len(large_chains) < 2:
        raise RuntimeError("7YEA does not contain two OGT-sized chains")
    first, second = large_chains[:2]
    mapping = structure_to_uniprot(first, human_sequence)
    second_atoms = np.asarray(
        [
            atom.coord
            for residue in amino_acid_residues(second)
            for atom in residue
            if atom.element != "H"
        ]
    )
    return {
        mapping[residue.id[1]]
        for residue in amino_acid_residues(first)
        if residue.id[1] in mapping and minimum_distance(residue, second_atoms) <= 5.0
    }


def conservation_rows(
    alignment: dict[str, str], categories: dict[str, set[int]]
) -> pd.DataFrame:
    human = alignment["human"]
    position_columns = position_to_column(human)
    mouse = alignment["mouse"]
    rows = []
    for category, positions in categories.items():
        for position in sorted(positions):
            column = position_columns.get(position)
            if column is None or column >= len(mouse):
                continue
            rows.append(
                {
                    "category": category,
                    "position": position,
                    "human_residue": human[column],
                    "mouse_residue": mouse[column],
                    "identical": human[column] == mouse[column],
                }
            )
    return pd.DataFrame(rows)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    config = load_config(args.config)
    output = args.output_dir.resolve()
    output.mkdir(parents=True, exist_ok=True)
    source = config.source_root / "analysis/revalidation/data"
    structures = source / "ogt_structures"
    orthologs = source / "ogt_orthologs"
    human_sequence = next(iter(parse_fasta(orthologs / "O15294.fasta").values()))
    alignment = parse_fasta(orthologs / "orthologs_aln.fasta")

    contacts = peptide_contacts(structures, human_sequence)
    ladder = {
        position
        for position in contacts
        if position <= 475 and human_sequence[position - 1] == "N"
    }
    reading_surface = contacts | CATALYTIC
    other_surface, buried = surface_categories(structures, human_sequence, reading_surface)
    interface = dimer_interface(structures, human_sequence)
    categories = {
        "peptide_contact_channel": contacts,
        "asn_ladder": ladder,
        "catalytic_triad": CATALYTIC,
        "reading_surface": reading_surface,
        "other_surface": other_surface,
        "buried_core": buried,
    }
    residues = conservation_rows(alignment, categories)
    summary = (
        residues.groupby("category", as_index=False)
        .agg(n=("position", "size"), identity=("identical", "mean"))
        .sort_values("category")
    )
    residues.to_csv(output / "residue_conservation.csv", index=False)
    summary.to_csv(output / "category_summary.csv", index=False)
    write_json(
        output / "summary.json",
        {
            "peptide_contacts": len(contacts),
            "asn_ladder": [int(position) for position in sorted(ladder)],
            "reading_surface": len(reading_surface),
            "dimer_interface": len(interface),
            "dimer_reading_surface_overlap": len(interface & reading_surface),
            "human_mouse_categories": [
                {
                    "category": str(row.category),
                    "n": int(row.n),
                    "identity": float(row.identity),
                }
                for row in summary.itertuples(index=False)
            ],
        },
    )


if __name__ == "__main__":
    main()
