"""Evolutionary conservation, transfer, and null-calibrated regional grammar."""

from __future__ import annotations

import argparse
from collections import defaultdict
from collections.abc import Callable
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.ensemble import HistGradientBoostingClassifier
from sklearn.metrics import roc_auc_score

from ..core.config import load_config
from ..core.io import write_json
from ..core.randomness import stable_seed
from ..models.orthology import (
    TAXON_TO_SPECIES,
    EvolutionInputs,
    composition_vector,
    cosine,
    load_evolution_inputs,
    one_minus_js,
    position_correlation,
    raw_pearson,
    sequence_features,
    species_consensus_spans,
)

DIVERGENCE_MILLION_YEARS = {"mouse": 90, "rat": 90, "zebrafish": 450, "fly": 700, "worm": 700}


def site_conservation(inputs: EvolutionInputs) -> pd.DataFrame:
    rows = []
    for accession, sites in inputs.human_region_sites.items():
        alignment = inputs.alignments.get(accession)
        if not alignment:
            continue
        for taxon, (_, ortholog) in alignment.orthologs.items():
            species = TAXON_TO_SPECIES[taxon]
            for position in sites:
                column = alignment.human_position_to_column.get(position)
                if column is None or column >= len(ortholog):
                    continue
                neighbourhood = [
                    alignment.human_position_to_column[nearby]
                    for nearby in range(position - 5, position + 6)
                    if nearby in alignment.human_position_to_column
                    and alignment.human_position_to_column[nearby] < len(ortholog)
                ]
                rows.append(
                    {
                        "accession": accession,
                        "species": species,
                        "position": position,
                        "exact": int(ortholog[column] in "ST"),
                        "regional": int(
                            any(ortholog[index] in "ST" for index in neighbourhood)
                        ),
                    }
                )
    return pd.DataFrame(rows)


def target_sequences(inputs: EvolutionInputs, species: str, taxon: int) -> dict[str, str]:
    sequences = {
        accession: sequence
        for (candidate_taxon, accession), sequence in inputs.sequences_by_taxon.items()
        if candidate_taxon == taxon
    }
    if species in {"rice", "Arabidopsis"}:
        sequences.update(inputs.external_sequences)
    return sequences


def build_species_tiles(
    inputs: EvolutionInputs, species: str, taxon: int
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    sequences = target_sequences(inputs, species, taxon)
    if species == "human":
        spans = inputs.spans
    else:
        spans = {
            accession: species_consensus_spans(positions)
            for accession, positions in inputs.species_sites.get(species, {}).items()
        }
    features, labels, accessions, genes = [], [], [], []
    for accession in sorted(spans):
        sequence = sequences.get(accession)
        regions = spans[accession]
        if not sequence or not regions:
            continue
        gene = inputs.genes_by_taxon.get((taxon, accession), "")
        for center in range(11, len(sequence) - 10, 10):
            features.append(sequence_features(sequence, center))
            labels.append(int(any(start - 15 <= center <= end + 15 for start, end in regions)))
            accessions.append(accession)
            genes.append(gene)
    return (
        np.asarray(features, dtype=float),
        np.asarray(labels, dtype=int),
        np.asarray(accessions),
        np.asarray(genes),
    )


def scanner() -> HistGradientBoostingClassifier:
    return HistGradientBoostingClassifier(
        max_iter=500,
        max_depth=4,
        learning_rate=0.03,
        l2_regularization=1.0,
        random_state=42,
    )


def transfer_predictions(inputs: EvolutionInputs) -> pd.DataFrame:
    human = build_species_tiles(inputs, "human", 9606)
    human_features, human_labels, human_accessions, human_genes = human
    model = scanner().fit(human_features, human_labels)
    training_genes = set(human_genes[human_genes != ""])
    rows = []
    target_data = {}
    for species, taxon in (("mouse", 10090), ("Arabidopsis", 3702), ("rice", 39947)):
        features, labels, accessions, genes = build_species_tiles(inputs, species, taxon)
        target_data[species] = (features, labels, accessions, genes)
        scores = model.predict_proba(features)[:, 1]
        keep = np.ones(len(labels), dtype=bool)
        if species == "mouse":
            keep = np.asarray([bool(gene) and gene not in training_genes for gene in genes])
        for accession, label, score in zip(
            accessions[keep], labels[keep], scores[keep], strict=True
        ):
            rows.append(
                {
                    "direction": "human_to_target",
                    "species": species,
                    "accession": accession,
                    "label": int(label),
                    "score": float(score),
                    "gene_disjoint": species == "mouse",
                }
            )
    rice_features, rice_labels, _, _ = target_data["rice"]
    rice_model = scanner().fit(rice_features, rice_labels)
    rice_to_human = rice_model.predict_proba(human_features)[:, 1]
    for accession, label, score in zip(
        human_accessions, human_labels, rice_to_human, strict=True
    ):
        rows.append(
            {
                "direction": "rice_to_human",
                "species": "human",
                "accession": accession,
                "label": int(label),
                "score": float(score),
                "gene_disjoint": False,
            }
        )
    return pd.DataFrame(rows)


def bootstrap_by_protein(
    frame: pd.DataFrame,
    statistic: Callable[[pd.DataFrame], float],
    *,
    draws: int,
    seed: int,
) -> tuple[float, float]:
    accessions = frame.accession.unique()
    groups = {accession: frame.loc[frame.accession.eq(accession)] for accession in accessions}
    rng = np.random.default_rng(seed)
    estimates = []
    for _ in range(draws):
        selection = rng.choice(accessions, len(accessions), replace=True)
        sample = pd.concat([groups[accession] for accession in selection], ignore_index=True)
        estimates.append(statistic(sample))
    low, high = np.quantile(estimates, [0.025, 0.975])
    return float(low), float(high)


def local_disorder(inputs: EvolutionInputs, accession: str, start: int, end: int) -> float:
    values = [inputs.disorder.get((accession, position)) for position in range(start, end + 1)]
    observed = [value for value in values if value is not None]
    return float(np.mean(observed)) if observed else np.nan


def matched_control_windows(
    inputs: EvolutionInputs, accession: str, start: int, end: int
) -> list[tuple[int, int]]:
    sequence = inputs.human_sequences[accession]
    length = end - start + 1
    region_sequence = sequence[start - 1 : end]
    region_st = (region_sequence.count("S") + region_sequence.count("T")) / length
    region_disorder = local_disorder(inputs, accession, start, end)
    controls = []
    for control_start in range(1, len(sequence) - length + 2, max(1, length // 2)):
        control_end = control_start + length - 1
        if any(
            not (control_end < other_start - length or control_start > other_end + length)
            for other_start, other_end in inputs.spans[accession]
        ):
            continue
        control_sequence = sequence[control_start - 1 : control_end]
        control_st = (control_sequence.count("S") + control_sequence.count("T")) / length
        control_disorder = local_disorder(inputs, accession, control_start, control_end)
        if (
            np.isfinite(control_disorder)
            and abs(control_st - region_st) <= 0.10
            and abs(control_disorder - region_disorder) <= 0.15
        ):
            controls.append((control_start, control_end))
    return controls


def aligned_segment(alignment, ortholog: str, start: int, end: int) -> str:
    return "".join(
        ortholog[alignment.human_position_to_column[position]]
        for position in range(start, end + 1)
        if position in alignment.human_position_to_column
        and alignment.human_position_to_column[position] < len(ortholog)
        and ortholog[alignment.human_position_to_column[position]] != "-"
    )


def directional_substitutions(inputs: EvolutionInputs) -> pd.DataFrame:
    sticker = set("DEFWY")
    rung_taxa = {
        "mammal": {"10090", "10116"},
        "vertebrate": {"10090", "10116", "7955"},
        "all": set(TAXON_TO_SPECIES),
    }
    rows = []
    for rung, taxa in rung_taxa.items():
        for accession, regions in inputs.spans.items():
            alignment = inputs.alignments.get(accession)
            if not alignment:
                continue
            orthologs = [
                sequence
                for taxon, (_, sequence) in alignment.orthologs.items()
                if taxon in taxa
            ]
            if not orthologs:
                continue
            for start, end in regions:
                controls = matched_control_windows(inputs, accession, start, end)
                if not controls:
                    continue
                control = controls[len(controls) // 2]
                for arm, (left, right) in (("region", (start, end)), ("control", control)):
                    for ortholog in orthologs:
                        for position in range(left, right + 1):
                            column = alignment.human_position_to_column.get(position)
                            if (
                                column is None
                                or column >= len(ortholog)
                                or ortholog[column] == "-"
                            ):
                                continue
                            if "-" in ortholog[max(0, column - 3) : column + 4]:
                                continue
                            human_residue = inputs.human_sequences[accession][position - 1]
                            ortholog_residue = ortholog[column]
                            if (
                                human_residue != ortholog_residue
                                and human_residue not in sticker
                            ):
                                rows.append(
                                    {
                                        "accession": accession,
                                        "rung": rung,
                                        "arm": arm,
                                        "gain": int(ortholog_residue in sticker),
                                    }
                                )
    return pd.DataFrame(rows)


def composition_position_nulls(inputs: EvolutionInputs, seed: int) -> pd.DataFrame:
    rng = np.random.default_rng(seed)
    sequence_pool: dict[str, list[str]] = defaultdict(list)
    for alignment in inputs.alignments.values():
        for taxon, (_, ortholog) in alignment.orthologs.items():
            ungapped = ortholog.replace("-", "")
            if len(ungapped) >= 60:
                sequence_pool[taxon].append(ungapped)

    rows = []
    for accession, regions in inputs.spans.items():
        alignment = inputs.alignments.get(accession)
        if not alignment:
            continue
        for taxon, (ortholog_accession, ortholog) in alignment.orthologs.items():
            species = TAXON_TO_SPECIES[taxon]
            sites = inputs.species_sites.get(species, {}).get(ortholog_accession, set())
            if not sites:
                continue
            ortholog_position_by_column = {}
            position = 0
            for column, residue in enumerate(ortholog):
                if residue != "-":
                    position += 1
                    ortholog_position_by_column[column] = position
            for start, end in regions:
                columns = [
                    alignment.human_position_to_column[p]
                    for p in range(start, end + 1)
                    if p in alignment.human_position_to_column
                    and alignment.human_position_to_column[p] < len(ortholog)
                    and ortholog[alignment.human_position_to_column[p]] != "-"
                ]
                if len(columns) < 8:
                    continue
                human_positions = [
                    p
                    for p in range(start, end + 1)
                    if p in alignment.human_position_to_column
                    and alignment.human_position_to_column[p] in columns
                ]
                human_sites = np.asarray(
                    [int(p in inputs.human_region_sites[accession]) for p in human_positions]
                )
                ortholog_sites = np.asarray(
                    [int(ortholog_position_by_column[column] in sites) for column in columns]
                )
                length = min(len(human_sites), len(ortholog_sites))
                human_sites, ortholog_sites = human_sites[:length], ortholog_sites[:length]
                if length < 8:
                    continue
                human_vector = composition_vector(
                    inputs.human_sequences[accession][start - 1 : end]
                )
                true_vector = composition_vector(ortholog[column] for column in columns)
                true_start = min(ortholog_position_by_column[column] for column in columns)
                true_end = max(ortholog_position_by_column[column] for column in columns)
                ungapped = ortholog.replace("-", "")
                null_a_vector = None
                null_a_sites = None
                if len(ungapped) > length + 20:
                    for _ in range(40):
                        null_start = int(rng.integers(1, len(ungapped) - length + 1))
                        if (
                            null_start + length - 1 < true_start - 10
                            or null_start > true_end + 10
                        ):
                            null_a_vector = composition_vector(
                                ungapped[null_start - 1 : null_start + length - 1]
                            )
                            null_a_sites = np.asarray(
                                [int(null_start + offset in sites) for offset in range(length)]
                            )
                            break
                null_b_vector = None
                candidates = sequence_pool.get(taxon, [])
                if candidates:
                    for _ in range(20):
                        candidate = candidates[int(rng.integers(0, len(candidates)))]
                        if len(candidate) > length + 5:
                            null_start = int(rng.integers(1, len(candidate) - length + 1))
                            null_b_vector = composition_vector(
                                candidate[null_start - 1 : null_start + length - 1]
                            )
                            break
                record = {
                    "accession": accession,
                    "species": species,
                    "n_columns": length,
                    "true_pearson": raw_pearson(human_vector, true_vector),
                    "true_cosine": cosine(human_vector, true_vector),
                    "true_js": one_minus_js(human_vector, true_vector),
                    "true_position": position_correlation(human_sites, ortholog_sites),
                }
                if null_a_vector is not None and null_a_sites is not None:
                    record.update(
                        {
                            "null_a_pearson": raw_pearson(human_vector, null_a_vector),
                            "null_a_cosine": cosine(human_vector, null_a_vector),
                            "null_a_js": one_minus_js(human_vector, null_a_vector),
                            "null_a_position": position_correlation(
                                human_sites, null_a_sites[:length]
                            ),
                        }
                    )
                if null_b_vector is not None:
                    record.update(
                        {
                            "null_b_pearson": raw_pearson(human_vector, null_b_vector),
                            "null_b_cosine": cosine(human_vector, null_b_vector),
                            "null_b_js": one_minus_js(human_vector, null_b_vector),
                        }
                    )
                rows.append(record)
    return pd.DataFrame(rows)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--regions", type=Path, required=True)
    parser.add_argument("--sites", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    config = load_config(args.config)
    output = args.output_dir.resolve()
    output.mkdir(parents=True, exist_ok=True)
    base_seed = int(config.values["randomness"]["analysis_base_seed"])
    inputs = load_evolution_inputs(config.source_root, args.regions, args.sites)

    conservation = site_conservation(inputs)
    transfer = transfer_predictions(inputs)
    substitutions = directional_substitutions(inputs)
    nulls = composition_position_nulls(inputs, base_seed)
    conservation.to_csv(output / "site_conservation.csv", index=False)
    transfer.to_csv(output / "transfer.csv", index=False)
    substitutions.to_csv(output / "substitutions.csv", index=False)
    nulls.to_csv(output / "composition_position_nulls.csv", index=False)

    conservation_results = []
    for species in ("mouse", "rat", "zebrafish", "fly", "worm"):
        group = conservation.loc[conservation.species.eq(species)]
        if group.empty:
            continue
        gap = group.regional - group.exact
        low, high = bootstrap_by_protein(
            group.assign(gap=gap),
            lambda frame: float(frame.gap.mean()),
            draws=10_000,
            seed=stable_seed(base_seed, "conservation", species),
        )
        conservation_results.append(
            {
                "species": species,
                "divergence_million_years": DIVERGENCE_MILLION_YEARS[species],
                "exact": float(group.exact.mean()),
                "regional": float(group.regional.mean()),
                "gap": float(gap.mean()),
                "gap_ci": [low, high],
                "n_sites": len(group),
                "n_proteins": group.accession.nunique(),
            }
        )
    transfer_results = []
    for (direction, species), group in transfer.groupby(["direction", "species"], sort=True):
        auc = float(roc_auc_score(group.label, group.score))
        low, high = bootstrap_by_protein(
            group,
            lambda frame: float(roc_auc_score(frame.label, frame.score)),
            draws=2_000,
            seed=stable_seed(base_seed, "transfer", direction, species),
        )
        transfer_results.append(
            {
                "direction": direction,
                "species": species,
                "auroc": auc,
                "ci": [low, high],
                "n_tiles": len(group),
                "n_proteins": group.accession.nunique(),
            }
        )

    null_results = {}
    for metric in ("pearson", "cosine", "js"):
        valid = nulls[["accession", f"true_{metric}", f"null_a_{metric}"]].dropna()
        valid = valid.assign(excess=valid[f"true_{metric}"] - valid[f"null_a_{metric}"])
        low, high = bootstrap_by_protein(
            valid,
            lambda frame: float(frame.excess.mean()),
            draws=10_000,
            seed=stable_seed(base_seed, "composition_null", metric),
        )
        null_results[metric] = {
            "true": float(valid[f"true_{metric}"].mean()),
            "null_a": float(valid[f"null_a_{metric}"].mean()),
            "excess": float(valid.excess.mean()),
            "excess_ci": [low, high],
        }
    position = nulls[["accession", "true_position", "null_a_position"]].fillna(0.0)
    position = position.assign(excess=position.true_position - position.null_a_position)
    low, high = bootstrap_by_protein(
        position,
        lambda frame: float(frame.excess.mean()),
        draws=10_000,
        seed=stable_seed(base_seed, "position_null"),
    )
    null_results["position"] = {
        "true": float(position.true_position.mean()),
        "null_a": float(position.null_a_position.mean()),
        "excess": float(position.excess.mean()),
        "excess_ci": [low, high],
    }
    write_json(
        output / "summary.json",
        {
            "alignments_used": len(inputs.alignments),
            "site_conservation": conservation_results,
            "transfer": transfer_results,
            "composition_position_nulls": null_results,
            "substitution_events": len(substitutions),
        },
    )


if __name__ == "__main__":
    main()
