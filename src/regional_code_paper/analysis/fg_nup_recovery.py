"""Test whether the regional object recovers canonical FG nucleoporins.

This module deliberately excludes exploratory ontology scans.  The manuscript
claims one prespecified biological recovery: enrichment of ten canonical
FG/FxFG nucleoporins.  We report that test against the atlas universe, a
covariate-matched universe, and a covariate-adjusted logistic model.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy.stats import hypergeom

from ..core.config import load_config
from ..core.io import write_json
from ..core.randomness import stable_seed

FG_NUPS = frozenset(
    {
        "NUP62",
        "NUP54",
        "NUP58",
        "NUP98",
        "NUP153",
        "NUP214",
        "NUP42",
        "NUP50",
        "RANBP2",
        "POM121",
    }
)
FG_NUP_COMPLEX = frozenset(
    {
        "NUP62",
        "NUP54",
        "NUP58",
        "NUP98",
        "NUP153",
        "NUP214",
        "NUP42",
        "NUP50",
        "NUP35",
        "RANBP2",
        "POM121",
        "NUP88",
        "NUP93",
        "NUP205",
        "NUP188",
        "NUP155",
        "NUP160",
        "NUP133",
        "NUP107",
        "NUP85",
        "NUP43",
        "NUP37",
        "SEH1L",
        "SEC13",
        "TPR",
        "AHCTF1",
        "NDC1",
        "POM121C",
    }
)


def protein_table(config_path: Path, regions_path: Path) -> pd.DataFrame:
    """Return one row per validated human O-GlcNAc protein."""
    config = load_config(config_path)
    root = config.source_root
    fasta = pd.read_parquet(
        root / "data/interim/fasta_all.parquet",
        columns=["accession", "gene", "sequence", "taxon_id", "is_canonical"],
    )
    fasta = fasta.loc[(fasta.taxon_id == 9606) & fasta.is_canonical]
    fasta = fasta.drop_duplicates("accession").set_index("accession")
    atlas = pd.read_csv(
        root / "analysis/revalidation/data/atlas_unambiguous.csv",
        usecols=["species", "accession", "position_in_protein", "site_residue"],
        dtype=str,
        encoding_errors="replace",
    )
    atlas = atlas.loc[atlas.species.eq("human") & atlas.site_residue.isin(["S", "T"])]
    atlas["position"] = pd.to_numeric(atlas.position_in_protein, errors="coerce")
    counts: dict[str, int] = {}
    for accession, group in atlas.dropna(subset=["position"]).groupby("accession", sort=True):
        if accession not in fasta.index:
            continue
        sequence = str(fasta.at[accession, "sequence"])
        positions = {
            int(position)
            for position, residue in group[["position", "site_residue"]].itertuples(
                index=False, name=None
            )
            if 1 <= int(position) <= len(sequence) and sequence[int(position) - 1] == residue
        }
        if positions:
            counts[str(accession)] = len(positions)
    disorder = pd.read_parquet(
        root / "data/interim/iupred_residue_scores.parquet",
        columns=["accession", "disorder_score"],
    )
    disorder = disorder.assign(disordered=disorder.disorder_score.gt(0.5))
    disorder_fraction = disorder.groupby("accession").disordered.mean()
    region_accessions = set(
        pd.read_csv(regions_path, usecols=["accession"]).accession.astype(str)
    )

    rows = []
    fallback_disorder = float(disorder_fraction.reindex(counts).median())
    for accession, n_sites in sorted(counts.items()):
        sequence = str(fasta.at[accession, "sequence"])
        gene = str(fasta.at[accession, "gene"])
        rows.append(
            {
                "accession": accession,
                "gene": gene,
                "has_region": accession in region_accessions,
                "is_fg_nup": gene in FG_NUPS,
                "is_fg_nup_complex": gene in FG_NUP_COMPLEX,
                "length": len(sequence),
                "disorder_fraction": float(disorder_fraction.get(accession, fallback_disorder)),
                "st_fraction": (sequence.count("S") + sequence.count("T")) / len(sequence),
                "n_sites": n_sites,
            }
        )
    return pd.DataFrame(rows)


def matched_background(table: pd.DataFrame, *, seed: int) -> tuple[pd.DataFrame, float]:
    """Greedily match each region protein to three unused non-region proteins."""
    features = table[["disorder_fraction", "st_fraction", "length", "n_sites"]].copy()
    features[["length", "n_sites"]] = np.log10(features[["length", "n_sites"]] + 1)
    z = (features - features.mean()) / features.std(ddof=0).replace(0, 1)
    positive = table.index[table.has_region].tolist()
    donors = set(table.index[~table.has_region].tolist())
    rng = np.random.default_rng(seed)
    # A seeded jitter makes exact-distance ties explicit and reproducible.
    jitter = dict(zip(table.index, rng.uniform(0, 1e-12, len(table)), strict=True))
    controls: list[int] = []
    for index in sorted(positive, key=lambda item: (-z.at[item, "disorder_fraction"], item)):
        candidates = []
        for donor in donors:
            distance = np.abs(z.loc[donor] - z.loc[index])
            if bool((distance <= 0.5).all()):
                candidates.append((float(np.square(distance).sum()), jitter[donor], donor))
        for _, _, donor in sorted(candidates)[:3]:
            controls.append(donor)
            donors.remove(donor)
    matched = table.loc[sorted(set(positive + controls))].copy()
    fg_total = int(matched.is_fg_nup_complex.sum())
    fg_positive = int((matched.is_fg_nup_complex & matched.has_region).sum())
    p_value = float(
        hypergeom.sf(fg_positive - 1, len(matched), fg_total, int(matched.has_region.sum()))
    )
    return matched, p_value


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--regions", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    output = args.output_dir.resolve()
    output.mkdir(parents=True, exist_ok=True)
    config = load_config(args.config)
    table = protein_table(args.config, args.regions)
    matched, matched_p = matched_background(
        table,
        seed=stable_seed(int(config.values["randomness"]["manuscript_base_seed"]), "fg_nup"),
    )
    design = pd.DataFrame(
        {
            "intercept": 1.0,
            "has_region": table.has_region.astype(float),
            "disorder_fraction": table.disorder_fraction,
            "st_fraction": table.st_fraction,
            "log10_length": np.log10(table.length),
        }
    )
    model = sm.Logit(table.is_fg_nup_complex.astype(int), design).fit(disp=False)
    present = table.loc[table.is_fg_nup]
    recovered = present.loc[present.has_region]
    hypergeometric_p = float(
        hypergeom.sf(
            len(recovered) - 1,
            len(table),
            len(present),
            int(table.has_region.sum()),
        )
    )
    table.to_csv(output / "proteins.csv", index=False)
    matched.to_csv(output / "matched_universe.csv", index=False)
    write_json(
        output / "summary.json",
        {
            "fg_nups_present": sorted(present.gene.tolist()),
            "fg_nups_recovered": sorted(recovered.gene.tolist()),
            "n_present": len(present),
            "n_recovered": len(recovered),
            "recovery_fraction": len(recovered) / len(present),
            "background_region_fraction": float(table.has_region.mean()),
            "hypergeometric_p": hypergeometric_p,
            "matched_p": matched_p,
            "adjusted_odds_ratio": float(np.exp(model.params["has_region"])),
            "adjusted_p": float(model.pvalues["has_region"]),
            "universe_proteins": len(table),
            "matched_universe_proteins": len(matched),
        },
    )


if __name__ == "__main__":
    main()
