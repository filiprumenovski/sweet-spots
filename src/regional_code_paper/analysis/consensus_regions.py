"""Construct the operational O-GlcNAc region object used throughout the paper.

The definition has two stages:

1. Validate human serine/threonine sites against canonical sequences.
2. Identify strict core sites from components of at least three sites at gap 5.
3. Re-segment those core sites at gap 10 to obtain the reported regions.

For provenance, step 2 is implemented as the intersection of qualifying
components at gaps 5, 8, 10, 12, and 15. These scalar thresholds are nested, so
their intersection is mathematically equivalent to the strictest, gap-5 pass.
The gap-10 pass changes final grouping, not core membership.

All coordinates are one-based and inclusive. CSV is used at the Python boundary;
the workflow performs any Parquet materialisation with DuckDB.
"""

from __future__ import annotations

import argparse
from collections import defaultdict
from collections.abc import Iterable
from dataclasses import dataclass
from pathlib import Path

import pandas as pd
from scipy.optimize import linear_sum_assignment

from ..core.config import load_config
from ..core.io import ensure_parent, sha256_file, write_json

DEFINITION_VERSION = "og_consensus_gaps5-8-10-12-15_core-gap10_min3_strict-v1"


@dataclass(frozen=True)
class SiteEvidence:
    accession: str
    position: int
    residue: str
    pmids: frozenset[str]


def segment_positions(
    positions: Iterable[int], *, maximum_gap: int, minimum_sites: int
) -> list[tuple[int, ...]]:
    """Return single-linkage components meeting the minimum-site threshold."""
    ordered = sorted(set(int(position) for position in positions))
    if not ordered:
        return []
    components: list[tuple[int, ...]] = []
    current = [ordered[0]]
    for position in ordered[1:]:
        if position - current[-1] > maximum_gap:
            if len(current) >= minimum_sites:
                components.append(tuple(current))
            current = []
        current.append(position)
    if len(current) >= minimum_sites:
        components.append(tuple(current))
    return components


def consensus_components(
    positions: Iterable[int], *, gaps: tuple[int, ...], final_gap: int, minimum_sites: int
) -> list[tuple[int, ...]]:
    """Intersect nested core calls, then group the surviving sites at the final gap.

    With a common ``minimum_sites`` threshold, membership is monotone in the
    maximum gap. The intersection over ``gaps`` therefore equals the call at
    ``min(gaps)``. The explicit loop is retained because it mirrors the archived
    definition and records the complete configured construction.
    """
    votes: dict[int, int] = defaultdict(int)
    ordered = sorted(set(positions))
    for gap in gaps:
        for component in segment_positions(
            ordered, maximum_gap=gap, minimum_sites=minimum_sites
        ):
            for position in component:
                votes[position] += 1
    core = [position for position, vote in votes.items() if vote == len(gaps)]
    return segment_positions(core, maximum_gap=final_gap, minimum_sites=minimum_sites)


def load_sequences(path: Path) -> dict[str, str]:
    frame = pd.read_parquet(path, columns=["accession", "sequence", "taxon_id", "is_canonical"])
    frame = frame.loc[(frame.taxon_id == 9606) & frame.is_canonical]
    frame = frame.drop_duplicates("accession")
    return dict(zip(frame.accession.astype(str), frame.sequence.astype(str), strict=True))


def load_site_evidence(
    path: Path, sequences: dict[str, str]
) -> tuple[dict[str, set[int]], dict[tuple[str, int], SiteEvidence], dict[str, int]]:
    frame = pd.read_csv(
        path,
        usecols=["species", "accession", "position_in_protein", "site_residue", "pmid"],
        low_memory=False,
        encoding_errors="replace",
    )
    frame = frame.loc[frame.species.eq("human") & frame.site_residue.isin(["S", "T"])]
    frame["position"] = pd.to_numeric(frame.position_in_protein, errors="coerce")
    frame = frame.dropna(subset=["position"])

    pmids: dict[tuple[str, int], set[str]] = defaultdict(set)
    residues: dict[tuple[str, int], str] = {}
    audit = defaultdict(int)
    audit["input_rows"] = len(frame)
    for accession, position, residue, pmid in frame[
        ["accession", "position", "site_residue", "pmid"]
    ].itertuples(index=False, name=None):
        accession = str(accession)
        position = int(position)
        sequence = sequences.get(accession)
        if sequence is None:
            audit["missing_sequence"] += 1
            continue
        if not 1 <= position <= len(sequence):
            audit["out_of_range"] += 1
            continue
        if sequence[position - 1] != residue:
            audit["residue_mismatch"] += 1
            continue
        key = (accession, position)
        residues[key] = residue
        if pd.notna(pmid) and str(pmid).strip():
            pmids[key].add(str(int(pmid)) if isinstance(pmid, float) else str(pmid))
        else:
            pmids[key]

    evidence = {
        key: SiteEvidence(key[0], key[1], residues[key], frozenset(values))
        for key, values in pmids.items()
    }
    sites: dict[str, set[int]] = defaultdict(set)
    for accession, position in evidence:
        sites[accession].add(position)
    audit["validated_sites"] = len(evidence)
    audit["validated_proteins"] = len(sites)
    return dict(sites), evidence, dict(audit)


def maximum_publication_matching(
    component: tuple[int, ...], accession: str, evidence: dict[tuple[str, int], SiteEvidence]
) -> int:
    """Maximum sites that can be assigned to distinct supporting publications."""
    publications = sorted(
        set().union(*(evidence[(accession, position)].pmids for position in component))
    )
    if not publications:
        return 0
    cost = [
        [0.0 if pmid in evidence[(accession, position)].pmids else 1.0 for pmid in publications]
        for position in component
    ]
    rows, columns = linear_sum_assignment(cost)
    return sum(cost[row][column] == 0 for row, column in zip(rows, columns, strict=True))


def survives_each_publication_removal(
    component: tuple[int, ...],
    accession: str,
    evidence: dict[tuple[str, int], SiteEvidence],
    *,
    final_gap: int,
    minimum_sites: int,
) -> bool:
    publications = set().union(
        *(evidence[(accession, position)].pmids for position in component)
    )
    if not publications:
        return False
    for removed in publications:
        surviving = [
            position
            for position in component
            if evidence[(accession, position)].pmids - {removed}
        ]
        if not segment_positions(surviving, maximum_gap=final_gap, minimum_sites=minimum_sites):
            return False
    return True


def build_regions(
    sequences: dict[str, str],
    site_map: dict[str, set[int]],
    evidence: dict[tuple[str, int], SiteEvidence],
    *,
    gaps: tuple[int, ...],
    final_gap: int,
    minimum_sites: int,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    region_rows: list[dict[str, object]] = []
    site_rows: list[dict[str, object]] = []
    for accession in sorted(site_map):
        components = consensus_components(
            site_map[accession], gaps=gaps, final_gap=final_gap, minimum_sites=minimum_sites
        )
        for region_index, component in enumerate(components, start=1):
            start, end = component[0], component[-1]
            region_id = f"{DEFINITION_VERSION}:{accession}:{start}-{end}"
            publications = sorted(
                set().union(*(evidence[(accession, position)].pmids for position in component))
            )
            region_rows.append(
                {
                    "region_id": region_id,
                    "accession": accession,
                    "region_index": region_index,
                    "start": start,
                    "end": end,
                    "span": end - start + 1,
                    "valence": len(component),
                    "positions": ";".join(map(str, component)),
                    "sequence": sequences[accession][start - 1 : end],
                    "n_contributing_pmids": len(publications),
                    "pmids": ";".join(publications),
                    "max_sites_distinct_pmids": maximum_publication_matching(
                        component, accession, evidence
                    ),
                    "survives_every_single_pmid_removal": survives_each_publication_removal(
                        component,
                        accession,
                        evidence,
                        final_gap=final_gap,
                        minimum_sites=minimum_sites,
                    ),
                    "definition_version": DEFINITION_VERSION,
                }
            )
            for position in component:
                item = evidence[(accession, position)]
                site_rows.append(
                    {
                        "region_id": region_id,
                        "accession": accession,
                        "position": position,
                        "residue": item.residue,
                        "n_pmids": len(item.pmids),
                        "pmids": ";".join(sorted(item.pmids)),
                        "definition_version": DEFINITION_VERSION,
                    }
                )
    return pd.DataFrame(region_rows), pd.DataFrame(site_rows)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--regions", type=Path, required=True)
    parser.add_argument("--sites", type=Path, required=True)
    parser.add_argument("--summary", type=Path, required=True)
    args = parser.parse_args()

    config = load_config(args.config)
    settings = config.values["analysis"]
    gaps = tuple(int(value) for value in settings["consensus_gaps"])
    final_gap = int(settings["final_gap"])
    minimum_sites = int(settings["minimum_region_sites"])
    fasta = config.source_root / "data/interim/fasta_human.parquet"
    atlas = config.source_root / "analysis/revalidation/data/atlas_unambiguous.csv"
    sequences = load_sequences(fasta)
    site_map, evidence, audit = load_site_evidence(atlas, sequences)
    regions, sites = build_regions(
        sequences,
        site_map,
        evidence,
        gaps=gaps,
        final_gap=final_gap,
        minimum_sites=minimum_sites,
    )
    ensure_parent([args.regions, args.sites, args.summary])
    regions.to_csv(args.regions, index=False)
    sites.to_csv(args.sites, index=False)

    reference = config.source_root / "data/processed/regions/oglcnac_consensus_regions.parquet"
    summary = {
        "definition": {
            "gaps": gaps,
            "final_gap": final_gap,
            "minimum_sites": minimum_sites,
            "coordinate_system": "one-based inclusive",
        },
        "audit": audit,
        "results": {
            "regions": len(regions),
            "region_sites": len(sites),
            "region_bearing_proteins": regions.accession.nunique(),
            "atlas_proteins": len(site_map),
            "regions_with_two_pmids": int((regions.n_contributing_pmids >= 2).sum()),
            "regions_with_three_pmids": int((regions.n_contributing_pmids >= 3).sum()),
            "regions_surviving_each_publication_removal": int(
                regions.survives_every_single_pmid_removal.sum()
            ),
        },
        "inputs": {"fasta_sha256": sha256_file(fasta), "atlas_sha256": sha256_file(atlas)},
        "reference_region_object_sha256": sha256_file(reference),
    }
    write_json(args.summary, summary)


if __name__ == "__main__":
    main()
