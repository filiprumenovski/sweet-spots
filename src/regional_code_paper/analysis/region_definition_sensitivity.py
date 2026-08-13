"""Post hoc sensitivity analysis for the operational region definition.

The catalogue grid varies strict-core gap, minimum component size, and final
grouping gap. Predictive models are refit for every strict-core definition at
the declared final grouping gap. Every prediction remains protein-held-out.
"""

from __future__ import annotations

import argparse
import itertools
import json
import tempfile
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.ensemble import HistGradientBoostingClassifier
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import GroupKFold

from ..core.config import load_config
from ..core.io import write_csv, write_json
from ..core.randomness import stable_seed
from .consensus_regions import load_sequences, load_site_evidence, segment_positions
from .regional_code import COMPOSITION_FEATURES, bootstrap_interval, build_tiles, load_inputs

MODEL_FEATURES = {
    "composition": COMPOSITION_FEATURES,
    "covariates": ("disorder", "st_fraction"),
    "adjusted_composition": ("disorder", "st_fraction", *COMPOSITION_FEATURES),
}


def definition_id(core_gap: int, minimum_sites: int, final_gap: int) -> str:
    """Return a sortable identifier for one operational definition."""
    return f"core{core_gap}_min{minimum_sites}_final{final_gap}"


def construct_definition(
    site_map: dict[str, set[int]],
    *,
    core_gap: int,
    minimum_sites: int,
    final_gap: int,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Build minimal region and core-site tables for one sensitivity definition."""
    region_rows: list[dict[str, object]] = []
    site_rows: list[dict[str, object]] = []
    identifier = definition_id(core_gap, minimum_sites, final_gap)
    for accession in sorted(site_map):
        strict_components = segment_positions(
            site_map[accession], maximum_gap=core_gap, minimum_sites=minimum_sites
        )
        core_sites = [position for component in strict_components for position in component]
        final_components = segment_positions(
            core_sites, maximum_gap=final_gap, minimum_sites=minimum_sites
        )
        for region_index, component in enumerate(final_components, start=1):
            start, end = component[0], component[-1]
            region_id = f"{identifier}:{accession}:{start}-{end}"
            region_rows.append(
                {
                    "definition_id": identifier,
                    "region_id": region_id,
                    "accession": accession,
                    "region_index": region_index,
                    "start": start,
                    "end": end,
                    "span": end - start + 1,
                    "valence": len(component),
                    "positions": ";".join(map(str, component)),
                }
            )
            site_rows.extend(
                {
                    "definition_id": identifier,
                    "region_id": region_id,
                    "accession": accession,
                    "position": position,
                }
                for position in component
            )
    return pd.DataFrame(region_rows), pd.DataFrame(site_rows)


def site_keys(sites: pd.DataFrame) -> set[tuple[str, int]]:
    """Return accession-position keys from a core-site table."""
    return {
        (str(accession), int(position))
        for accession, position in sites[["accession", "position"]].itertuples(
            index=False, name=None
        )
    }


def interval_overlap_fraction(primary: pd.DataFrame, variant: pd.DataFrame) -> float:
    """Fraction of primary intervals touched by at least one variant interval."""
    by_accession: dict[str, list[tuple[int, int]]] = defaultdict(list)
    for accession, start, end in variant[["accession", "start", "end"]].itertuples(
        index=False, name=None
    ):
        by_accession[str(accession)].append((int(start), int(end)))
    recovered = 0
    for accession, start, end in primary[["accession", "start", "end"]].itertuples(
        index=False, name=None
    ):
        recovered += any(
            other_start <= int(end) and other_end >= int(start)
            for other_start, other_end in by_accession.get(str(accession), [])
        )
    return recovered / len(primary)


def load_validated_sites(config_path: Path):
    """Load the same sequence-validated atlas used by the primary definition."""
    config = load_config(config_path)
    fasta = config.source_root / "data/interim/fasta_human.parquet"
    atlas = config.source_root / "analysis/revalidation/data/atlas_unambiguous.csv"
    sequences = load_sequences(fasta)
    site_map, evidence, audit = load_site_evidence(atlas, sequences)
    return config, sequences, site_map, evidence, audit


def build_catalogue_grid(config_path: Path, output: Path) -> None:
    """Evaluate catalogue size and overlap for the complete parameter grid."""
    config, _, site_map, _, audit = load_validated_sites(config_path)
    settings = config.values["analysis"]
    sensitivity = settings["region_sensitivity"]
    primary_key = (
        min(int(value) for value in settings["consensus_gaps"]),
        int(settings["minimum_region_sites"]),
        int(settings["final_gap"]),
    )
    definitions: dict[tuple[int, int, int], tuple[pd.DataFrame, pd.DataFrame]] = {}
    for core_gap, minimum_sites, final_gap in itertools.product(
        sensitivity["core_gaps"],
        sensitivity["minimum_sites"],
        sensitivity["final_gaps"],
    ):
        key = (int(core_gap), int(minimum_sites), int(final_gap))
        definitions[key] = construct_definition(
            site_map,
            core_gap=key[0],
            minimum_sites=key[1],
            final_gap=key[2],
        )
    primary_regions, primary_sites = definitions[primary_key]
    primary_keys = site_keys(primary_sites)
    rows = []
    for (core_gap, minimum_sites, final_gap), (regions, sites) in definitions.items():
        variant_keys = site_keys(sites)
        intersection = len(primary_keys & variant_keys)
        union = len(primary_keys | variant_keys)
        rows.append(
            {
                "definition_id": definition_id(core_gap, minimum_sites, final_gap),
                "core_gap": core_gap,
                "minimum_sites": minimum_sites,
                "final_gap": final_gap,
                "is_primary": (core_gap, minimum_sites, final_gap) == primary_key,
                "regions": len(regions),
                "core_sites": len(sites),
                "region_bearing_proteins": regions.accession.nunique(),
                "median_span": float(regions.span.median()),
                "maximum_span": int(regions.span.max()),
                "fraction_validated_sites": len(sites) / int(audit["validated_sites"]),
                "primary_site_recall": intersection / len(primary_keys),
                "variant_site_precision": intersection / len(variant_keys),
                "site_jaccard": intersection / union,
                "primary_interval_overlap": interval_overlap_fraction(
                    primary_regions, regions
                ),
            }
        )
    write_csv(output, pd.DataFrame(rows).sort_values(
        ["minimum_sites", "core_gap", "final_gap"], kind="stable"
    ))


def held_out_scores(tiles: pd.DataFrame, seed: int) -> pd.DataFrame:
    """Fit the sensitivity models with five protein-grouped outer folds."""
    labels = tiles.label.to_numpy(int)
    groups = tiles.accession.to_numpy(str)
    scored = tiles[["accession", "center", "label", "disorder", "st_fraction"]].copy()
    for model_name, features in MODEL_FEATURES.items():
        scores = np.empty(len(tiles), dtype=float)
        for train, test in GroupKFold(5).split(tiles, labels, groups):
            model = HistGradientBoostingClassifier(
                max_iter=300,
                max_depth=3,
                learning_rate=0.04,
                l2_regularization=1.0,
                random_state=seed,
            )
            model.fit(tiles.iloc[train][list(features)], labels[train])
            scores[test] = model.predict_proba(tiles.iloc[test][list(features)])[:, 1]
        scored[f"{model_name}_score"] = scores
    return scored


def macro_within_auc(scored: pd.DataFrame, score_column: str) -> tuple[float, int]:
    """Macro-average held-out AUROC across informative proteins."""
    values = []
    for _, group in scored.groupby("accession", sort=True):
        if group.label.nunique() == 2:
            values.append(roc_auc_score(group.label, group[score_column]))
    return float(np.mean(values)), len(values)


def within_auc_by_protein(scored: pd.DataFrame, score_column: str) -> pd.Series:
    """Return held-out AUROC indexed by informative protein."""
    values = {}
    for accession, group in scored.groupby("accession", sort=True):
        if group.label.nunique() == 2:
            values[str(accession)] = roc_auc_score(group.label, group[score_column])
    return pd.Series(values, dtype=float)


def fit_definition_model(
    config_path: Path,
    *,
    core_gap: int,
    minimum_sites: int,
    output: Path,
) -> None:
    """Refit regional prediction for one core definition at the primary final gap."""
    config, _, site_map, _, _ = load_validated_sites(config_path)
    final_gap = int(config.values["analysis"]["final_gap"])
    regions, sites = construct_definition(
        site_map,
        core_gap=core_gap,
        minimum_sites=minimum_sites,
        final_gap=final_gap,
    )
    with tempfile.TemporaryDirectory(prefix="region-sensitivity-") as temporary:
        temporary_path = Path(temporary)
        regions_path = temporary_path / "regions.csv"
        sites_path = temporary_path / "sites.csv"
        regions.to_csv(regions_path, index=False)
        sites.to_csv(sites_path, index=False)
        inputs = load_inputs(config_path, regions_path, sites_path)
        tiles = build_tiles(inputs)
    seed = int(config.values["randomness"]["analysis_base_seed"])
    bootstrap_draws = int(config.values["randomness"]["figure_bootstrap_draws"])
    scored = held_out_scores(tiles, seed)
    labels = scored.label.to_numpy(int)
    matched = (scored.disorder >= scored.disorder.median()) & (
        scored.st_fraction >= scored.st_fraction.median()
    )
    row: dict[str, object] = {
        "definition_id": definition_id(core_gap, minimum_sites, final_gap),
        "core_gap": core_gap,
        "minimum_sites": minimum_sites,
        "final_gap": final_gap,
        "is_primary": core_gap == 5 and minimum_sites == 3,
        "regions": len(regions),
        "core_sites": len(sites),
        "region_bearing_proteins": regions.accession.nunique(),
        "tiles": len(scored),
        "positive_tiles": int(labels.sum()),
    }
    for model_name in MODEL_FEATURES:
        score_column = f"{model_name}_score"
        row[f"{model_name}_global_auroc"] = roc_auc_score(labels, scored[score_column])
        within_values = within_auc_by_protein(scored, score_column)
        row[f"{model_name}_within_protein_auroc"] = float(within_values.mean())
        interval = bootstrap_interval(
            within_values.to_numpy(float),
            draws=bootstrap_draws,
            seed=stable_seed(
                seed,
                f"region_definition_sensitivity:{core_gap}:{minimum_sites}:{model_name}",
            ),
        )
        row[f"{model_name}_within_protein_ci_low"] = interval[0]
        row[f"{model_name}_within_protein_ci_high"] = interval[1]
        row["within_protein_count"] = len(within_values)
    row["composition_matched_auroc"] = roc_auc_score(
        labels[matched], scored.loc[matched, "composition_score"]
    )
    row["adjusted_composition_increment"] = (
        float(row["adjusted_composition_within_protein_auroc"])
        - float(row["covariates_within_protein_auroc"])
    )
    adjusted_by_protein = within_auc_by_protein(
        scored, "adjusted_composition_score"
    )
    covariates_by_protein = within_auc_by_protein(scored, "covariates_score")
    paired_increment = adjusted_by_protein - covariates_by_protein
    increment_interval = bootstrap_interval(
        paired_increment.to_numpy(float),
        draws=bootstrap_draws,
        seed=stable_seed(
            seed, f"region_definition_sensitivity:{core_gap}:{minimum_sites}"
        ),
    )
    row["adjusted_composition_increment_ci_low"] = increment_interval[0]
    row["adjusted_composition_increment_ci_high"] = increment_interval[1]
    row["adjusted_composition_increment_positive_fraction"] = float(
        (paired_increment > 0).mean()
    )
    write_json(output, row)


def reduce_models(inputs: list[Path], output: Path) -> None:
    """Combine one JSON result per strict-core definition."""
    rows = [json.loads(path.read_text(encoding="utf-8")) for path in inputs]
    write_csv(
        output,
        pd.DataFrame(rows).sort_values(["minimum_sites", "core_gap"], kind="stable"),
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    catalogue = commands.add_parser("catalogue")
    catalogue.add_argument("--config", type=Path, required=True)
    catalogue.add_argument("--output", type=Path, required=True)
    model = commands.add_parser("model")
    model.add_argument("--config", type=Path, required=True)
    model.add_argument("--core-gap", type=int, required=True)
    model.add_argument("--minimum-sites", type=int, required=True)
    model.add_argument("--output", type=Path, required=True)
    reduce_parser = commands.add_parser("reduce")
    reduce_parser.add_argument("--inputs", type=Path, nargs="+", required=True)
    reduce_parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.command == "catalogue":
        build_catalogue_grid(args.config, args.output)
    elif args.command == "model":
        fit_definition_model(
            args.config,
            core_gap=args.core_gap,
            minimum_sites=args.minimum_sites,
            output=args.output,
        )
    else:
        reduce_models(args.inputs, args.output)


if __name__ == "__main__":
    main()
