"""Regional-code prediction, acceptor comparison, and local sequence summaries.

This module implements the analyses behind Figure 2, Figure 3B and Figure 3F.
Every cross-validation split is grouped by accession, and every interval
resamples whole proteins.
"""

from __future__ import annotations

import math
from collections import defaultdict
from collections.abc import Callable, Iterable
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.api as sm
from sklearn.ensemble import HistGradientBoostingClassifier
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import GroupKFold

from ..core.config import load_config
from ..core.randomness import stable_seed

COMPOSITION_FEATURES = (
    "acidic",
    "alanine",
    "aromatic",
    "proline",
    "qn",
    "charge",
    "hydrophobic",
)
GEOMETRY_FEATURES = ("geo_st", "geo_gapmean", "geo_gapsd", "geo_close")
MODEL_FEATURES = {
    "disorder": ("disorder",),
    "st_only": ("st_fraction",),
    "covariates": ("disorder", "st_fraction"),
    "composition": COMPOSITION_FEATURES,
    "geometry": GEOMETRY_FEATURES,
    "full": ("disorder", "st_fraction", *COMPOSITION_FEATURES, *GEOMETRY_FEATURES),
}


@dataclass(frozen=True)
class Inputs:
    sequences: dict[str, str]
    disorder: dict[tuple[str, int], float]
    oglcnac: dict[str, set[int]]
    phospho: dict[str, set[int]]
    spans: dict[str, list[tuple[int, int]]]
    region_sites: set[tuple[str, int]]


def _site_map(frame: pd.DataFrame) -> dict[str, set[int]]:
    output: dict[str, set[int]] = defaultdict(set)
    for accession, position in frame[["accession", "position"]].itertuples(
        index=False, name=None
    ):
        output[str(accession)].add(int(position))
    return dict(output)


def load_inputs(config_path: Path, regions_path: Path, sites_path: Path) -> Inputs:
    config = load_config(config_path)
    root = config.source_root
    fasta = pd.read_parquet(
        root / "data/interim/fasta_human.parquet",
        columns=["accession", "sequence", "taxon_id", "is_canonical"],
    )
    fasta = fasta.loc[(fasta.taxon_id == 9606) & fasta.is_canonical]
    fasta = fasta.drop_duplicates("accession")
    sequences = dict(zip(fasta.accession.astype(str), fasta.sequence.astype(str), strict=True))

    disorder_frame = pd.read_parquet(
        root / "data/interim/iupred_residue_scores.parquet",
        columns=["accession", "position", "disorder_score"],
    )
    disorder = {
        (str(accession), int(position)): float(score)
        for accession, position, score in disorder_frame.itertuples(index=False, name=None)
    }

    atlas = pd.read_csv(
        root / "analysis/revalidation/data/atlas_unambiguous.csv",
        usecols=["species", "accession", "position_in_protein", "site_residue"],
        low_memory=False,
        encoding_errors="replace",
    )
    atlas = atlas.loc[atlas.species.eq("human") & atlas.site_residue.isin(["S", "T"])]
    atlas["position"] = pd.to_numeric(atlas.position_in_protein, errors="coerce")
    atlas = atlas.dropna(subset=["position"])
    valid_oglcnac = []
    for accession, position, residue in (
        atlas[["accession", "position", "site_residue"]]
        .drop_duplicates()
        .itertuples(index=False, name=None)
    ):
        accession, position = str(accession), int(position)
        sequence = sequences.get(accession)
        if sequence and 1 <= position <= len(sequence) and sequence[position - 1] == residue:
            valid_oglcnac.append((accession, position))
    oglcnac = _site_map(pd.DataFrame(valid_oglcnac, columns=["accession", "position"]))

    unified = pd.read_parquet(
        root / "data/processed/landscape/ptm_unified.parquet",
        columns=["accession", "position", "residue", "mod_type", "taxon_id", "ambiguous"],
    )
    unified = unified.loc[
        (unified.taxon_id == 9606)
        & ~unified.ambiguous
        & unified.mod_type.eq("phospho")
        & unified.residue.isin(["S", "T"])
    ]
    valid_phospho = []
    for accession, position, residue in (
        unified[["accession", "position", "residue"]]
        .drop_duplicates()
        .itertuples(index=False, name=None)
    ):
        accession, position = str(accession), int(position)
        sequence = sequences.get(accession)
        if sequence and 1 <= position <= len(sequence) and sequence[position - 1] == residue:
            valid_phospho.append((accession, position))
    phospho = _site_map(pd.DataFrame(valid_phospho, columns=["accession", "position"]))

    regions = pd.read_csv(regions_path)
    spans: dict[str, list[tuple[int, int]]] = defaultdict(list)
    for accession, start, end in regions[["accession", "start", "end"]].itertuples(
        index=False, name=None
    ):
        spans[str(accession)].append((int(start), int(end)))
    region_sites_frame = pd.read_csv(sites_path, usecols=["accession", "position"])
    region_sites = {
        (str(accession), int(position))
        for accession, position in region_sites_frame.itertuples(index=False, name=None)
    }
    return Inputs(sequences, disorder, oglcnac, phospho, dict(spans), region_sites)


def local_covariates(
    inputs: Inputs, accession: str, position: int, radius: int = 25
) -> tuple[float, float]:
    """Mean disorder and S/T fraction in a clipped, residue-centred window."""
    sequence = inputs.sequences[accession]
    segment = sequence[max(0, position - 1 - radius) : min(len(sequence), position + radius)]
    values = [
        inputs.disorder.get((accession, index))
        for index in range(max(1, position - radius), min(len(sequence), position + radius) + 1)
    ]
    observed = [value for value in values if value is not None]
    disorder = float(np.mean(observed)) if observed else math.nan
    st_fraction = (segment.count("S") + segment.count("T")) / len(segment)
    return disorder, st_fraction


def sequence_features(segment: str) -> tuple[list[float], list[float]]:
    """Return deliberately disjoint composition and S/T-geometry feature families."""
    length = max(1, len(segment))
    fraction = {residue: segment.count(residue) / length for residue in set(segment)}

    def f(residue: str) -> float:
        return fraction.get(residue, 0.0)

    composition = [
        f("D") + f("E"),
        f("A"),
        f("F") + f("W") + f("Y"),
        f("P"),
        f("Q") + f("N"),
        f("K") + f("R") - f("D") - f("E"),
        sum(f(residue) for residue in "AVILM"),
    ]
    st_positions = np.asarray(
        [index for index, residue in enumerate(segment) if residue in "ST"], dtype=int
    )
    gaps = np.diff(st_positions)
    geometry = [
        len(st_positions) / length,
        float(gaps.mean()) if len(gaps) else float(length),
        float(gaps.std()) if len(gaps) > 1 else 0.0,
        float((gaps <= 5).mean()) if len(gaps) else 0.0,
    ]
    return composition, geometry


def build_tiles(inputs: Inputs) -> pd.DataFrame:
    columns = [
        "accession",
        "center",
        "label",
        "disorder",
        "st_fraction",
        *COMPOSITION_FEATURES,
        *GEOMETRY_FEATURES,
    ]
    rows: list[list[object]] = []
    for accession in sorted(inputs.spans):
        sequence = inputs.sequences.get(accession)
        if not sequence:
            continue
        for center in range(26, len(sequence) - 24, 10):
            segment = sequence[center - 26 : center + 25]
            composition, geometry = sequence_features(segment)
            disorder, st_fraction = local_covariates(inputs, accession, center)
            if not np.isfinite(disorder):
                continue
            label = int(
                any(start - 15 <= center <= end + 15 for start, end in inputs.spans[accession])
            )
            rows.append(
                [accession, center, label, disorder, st_fraction, *composition, *geometry]
            )
    return pd.DataFrame(rows, columns=columns)


def grouped_predictions(tiles: pd.DataFrame, seed: int) -> dict[str, np.ndarray]:
    labels = tiles.label.to_numpy(int)
    groups = tiles.accession.to_numpy(str)
    output: dict[str, np.ndarray] = {}
    for name, columns in MODEL_FEATURES.items():
        scores = np.zeros(len(tiles), dtype=float)
        for train, test in GroupKFold(5).split(tiles, labels, groups):
            model = HistGradientBoostingClassifier(
                max_iter=300,
                max_depth=3,
                learning_rate=0.04,
                l2_regularization=1.0,
                random_state=seed,
            )
            model.fit(tiles.iloc[train][list(columns)], labels[train])
            scores[test] = model.predict_proba(tiles.iloc[test][list(columns)])[:, 1]
        output[name] = scores
    return output


def bootstrap_interval(
    values: np.ndarray,
    statistic: Callable[[np.ndarray], float] = np.mean,
    *,
    draws: int,
    seed: int,
) -> tuple[float, float]:
    rng = np.random.default_rng(seed)
    estimates = np.empty(draws, dtype=float)
    for draw in range(draws):
        estimates[draw] = statistic(values[rng.integers(0, len(values), len(values))])
    low, high = np.quantile(estimates, [0.025, 0.975])
    return float(low), float(high)


def per_protein_auc(tiles: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for accession, group in tiles.groupby("accession", sort=True):
        if group.label.nunique() == 2:
            rows.append(
                {
                    "accession": accession,
                    "within_protein_auroc": roc_auc_score(group.label, group.full_score),
                }
            )
    return pd.DataFrame(rows)


def acceptor_models(inputs: Inputs) -> tuple[pd.DataFrame, dict[str, object]]:
    rows = []
    for accession, positions in inputs.oglcnac.items():
        sequence = inputs.sequences[accession]
        for position in positions:
            segment = sequence[max(0, position - 26) : min(len(sequence), position + 25)]
            composition, _ = sequence_features(segment)
            disorder, st_fraction = local_covariates(inputs, accession, position)
            if np.isfinite(disorder):
                rows.append(
                    [
                        accession,
                        position,
                        int((accession, position) in inputs.region_sites),
                        disorder,
                        st_fraction,
                        composition[0],
                        composition[1],
                        composition[2],
                    ]
                )
    sites = pd.DataFrame(
        rows,
        columns=[
            "accession",
            "position",
            "cluster",
            "disorder",
            "st_fraction",
            "acidic",
            "alanine",
            "aromatic",
        ],
    )
    estimates = []
    for outcome in ("acidic", "alanine", "aromatic"):
        design = sm.add_constant(sites[["cluster", "disorder", "st_fraction"]])
        fit = sm.OLS(sites[outcome], design).fit(
            cov_type="cluster", cov_kwds={"groups": sites.accession}
        )
        low, high = fit.conf_int().loc["cluster"]
        estimates.append(
            {
                "outcome": outcome,
                "coefficient": float(fit.params["cluster"]),
                "ci_low": float(low),
                "ci_high": float(high),
                "p_value": float(fit.pvalues["cluster"]),
                "n_sites": len(sites),
                "n_proteins": sites.accession.nunique(),
            }
        )
    return pd.DataFrame(estimates), {
        "n_sites": len(sites),
        "n_proteins": sites.accession.nunique(),
    }


def proline_profiles(inputs: Inputs) -> pd.DataFrame:
    offsets = tuple(range(-5, 6))
    flank_offsets = {-4, -3, -2, 2, 3, 4}
    rows: list[dict[str, object]] = []
    definitions = (
        ("phospho", inputs.phospho, None),
        ("oglcnac", inputs.oglcnac, None),
        ("cluster", inputs.oglcnac, inputs.region_sites),
    )
    for label, site_map, restriction in definitions:
        for accession, positions in site_map.items():
            sequence = inputs.sequences.get(accession, "")
            values = {offset: [] for offset in offsets}
            for position in positions:
                if restriction is not None and (accession, position) not in restriction:
                    continue
                for offset in offsets:
                    index = position - 1 + offset
                    if 0 <= index < len(sequence):
                        values[offset].append(int(sequence[index] == "P"))
            if not values[1]:
                continue
            row: dict[str, object] = {"accession": accession, "set": label}
            row.update(
                {
                    f"offset_{offset:+d}": float(np.mean(values[offset]))
                    if values[offset]
                    else math.nan
                    for offset in offsets
                }
            )
            row["plus_one"] = float(np.mean(values[1]))
            row["flank"] = float(
                np.mean([value for offset in flank_offsets for value in values[offset]])
            )
            rows.append(row)
    return pd.DataFrame(rows)


def distance_to_spans(position: int, spans: Iterable[tuple[int, int]]) -> int:
    return min(
        0 if start <= position <= end else min(abs(position - start), abs(position - end))
        for start, end in spans
    )


def phosphosite_distance_profile(
    inputs: Inputs, *, draws: int, base_seed: int
) -> tuple[pd.DataFrame, pd.DataFrame]:
    rows = []
    for accession in sorted(set(inputs.spans) & set(inputs.phospho)):
        sequence = inputs.sequences[accession]
        phosphosites = sorted(inputs.phospho[accession])
        controls = [
            index
            for index, residue in enumerate(sequence, start=1)
            if residue in "ST" and index not in inputs.phospho[accession]
        ]
        if not controls:
            continue
        control_features = np.asarray(
            [local_covariates(inputs, accession, position) for position in controls]
        )
        scale = np.nanstd(control_features, axis=0) + 1e-9
        used = np.zeros(len(controls), dtype=bool)
        for phosphosite in phosphosites:
            feature = np.asarray(local_covariates(inputs, accession, phosphosite))
            distance = np.nansum(((control_features - feature) / scale) ** 2, axis=1)
            distance[used] = np.inf
            match = int(np.argmin(distance))
            if not np.isfinite(distance[match]):
                continue
            used[match] = True
            control = controls[match]
            rows.append(
                {
                    "accession": accession,
                    "phosphosite": phosphosite,
                    "control": control,
                    "phosphosite_distance": distance_to_spans(
                        phosphosite, inputs.spans[accession]
                    ),
                    "control_distance": distance_to_spans(control, inputs.spans[accession]),
                    "phosphosite_disorder": feature[0],
                    "phosphosite_st_fraction": feature[1],
                    "control_disorder": control_features[match, 0],
                    "control_st_fraction": control_features[match, 1],
                }
            )
    matches = pd.DataFrame(rows)
    profile = []
    for radius in range(61):
        differences = (
            matches.assign(
                difference=(matches.phosphosite_distance <= radius).astype(float)
                - (matches.control_distance <= radius).astype(float)
            )
            .groupby("accession")
            .difference.mean()
            .to_numpy()
        )
        low, high = bootstrap_interval(
            differences,
            draws=draws,
            seed=stable_seed(base_seed, "distance", radius),
        )
        profile.append(
            {
                "radius": radius,
                "phosphosite_cdf": float((matches.phosphosite_distance <= radius).mean()),
                "control_cdf": float((matches.control_distance <= radius).mean()),
                "difference": float(differences.mean()),
                "ci_low": low,
                "ci_high": high,
                "n_proteins": len(differences),
            }
        )
    return matches, pd.DataFrame(profile)
