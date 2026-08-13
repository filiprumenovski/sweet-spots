"""Protein-balanced test of within-PTM spatial self-clustering.

This is the complete implementation behind Figures 1B, 1C and 3A. The null is
analytic: for each protein, the observed fraction of close site pairs is compared
with the exact expectation under fixed-count sampling from chemically eligible
residues. Whole proteins, never sites, are the resampling unit.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import platform
from datetime import UTC, datetime
from pathlib import Path

import numpy as np
import pandas as pd
import scipy
from scipy import stats

from ..core.config import load_config

PTMS = ("oglcnac", "phospho", "acetylation", "ubiquitination")
RADII = (5, 10, 15, 20, 25, 30)
PRIMARY_RADIUS = 10
MIN_SITES = 3
BASE_SEED = 20260814


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def falling_ratio(k: int, n: int, order: int) -> float:
    """Return (k)_order / (n)_order, including the k < order case."""
    if k < order:
        return 0.0
    value = 1.0
    for offset in range(order):
        value *= (k - offset) / (n - offset)
    return value


def close_pair_graph(positions: np.ndarray, radius: int) -> tuple[int, int]:
    """Number of close edges and unordered edge pairs sharing a vertex."""
    positions = np.asarray(positions, dtype=np.int64)
    left = np.searchsorted(positions, positions - radius, side="left")
    right = np.searchsorted(positions, positions + radius, side="right")
    degrees = right - left - 1
    edges = int(degrees.sum() // 2)
    adjacent_edge_pairs = int(np.sum(degrees * (degrees - 1) // 2))
    return edges, adjacent_edge_pairs


def null_pair_moments(
    n_eligible: int,
    n_sites: int,
    eligible_edges: int,
    adjacent_edge_pairs: int,
) -> tuple[float, float]:
    """Exact mean/variance of selected close edges under fixed-size sampling."""
    p2 = falling_ratio(n_sites, n_eligible, 2)
    p3 = falling_ratio(n_sites, n_eligible, 3)
    p4 = falling_ratio(n_sites, n_eligible, 4)
    all_edge_pairs = eligible_edges * (eligible_edges - 1) // 2
    disjoint_edge_pairs = all_edge_pairs - adjacent_edge_pairs
    mean = eligible_edges * p2
    second = eligible_edges * p2 + 2 * adjacent_edge_pairs * p3 + 2 * disjoint_edge_pairs * p4
    variance = max(0.0, second - mean * mean)
    return float(mean), float(variance)


def validate_site_map(
    frame: pd.DataFrame,
    sequences: dict[str, str],
    allowed_residues: set[str],
) -> tuple[dict[str, np.ndarray], dict[str, int]]:
    output: dict[str, set[int]] = {}
    audit = {
        "input_unique_sites": len(frame),
        "strict_sites": 0,
        "missing_canonical": 0,
        "out_of_range": 0,
        "residue_mismatch": 0,
        "wrong_residue": 0,
    }
    for accession, position, residue in frame[["accession", "position", "residue"]].itertuples(
        index=False, name=None
    ):
        position = int(position)
        if residue not in allowed_residues:
            audit["wrong_residue"] += 1
        elif accession not in sequences:
            audit["missing_canonical"] += 1
        elif not 1 <= position <= len(sequences[accession]):
            audit["out_of_range"] += 1
        elif sequences[accession][position - 1] != residue:
            audit["residue_mismatch"] += 1
        else:
            output.setdefault(accession, set()).add(position)
    arrays = {
        accession: np.asarray(sorted(positions), dtype=np.int64)
        for accession, positions in output.items()
    }
    audit["strict_sites"] = int(sum(len(x) for x in arrays.values()))
    audit["strict_proteins"] = len(arrays)
    return arrays, audit


def load_inputs(repo: Path) -> tuple[dict[str, str], dict[str, dict[str, np.ndarray]], dict]:
    fasta_path = repo / "data/interim/fasta_human.parquet"
    unified_path = repo / "data/processed/landscape/ptm_unified.parquet"
    atlas_path = repo / "analysis/revalidation/data/atlas_unambiguous.csv"

    fasta = pd.read_parquet(
        fasta_path,
        columns=["accession", "sequence", "taxon_id", "is_canonical"],
    )
    fasta = fasta.loc[(fasta.taxon_id == 9606) & fasta.is_canonical]
    fasta = fasta.drop_duplicates("accession")
    sequences = dict(zip(fasta.accession.astype(str), fasta.sequence.astype(str), strict=True))

    atlas = pd.read_csv(
        atlas_path,
        usecols=["accession", "position_in_protein", "site_residue", "species"],
        low_memory=False,
        encoding="utf-8",
        encoding_errors="replace",
    )
    atlas = atlas.loc[
        atlas.species.eq("human") & atlas.site_residue.isin(["S", "T"]),
        ["accession", "position_in_protein", "site_residue"],
    ].rename(columns={"position_in_protein": "position", "site_residue": "residue"})
    atlas["position"] = pd.to_numeric(atlas.position, errors="coerce")
    atlas = atlas.dropna(subset=["position"]).drop_duplicates()
    oglcnac, oglcnac_audit = validate_site_map(atlas, sequences, {"S", "T"})

    unified = pd.read_parquet(
        unified_path,
        columns=["accession", "position", "residue", "mod_type", "taxon_id", "ambiguous"],
    )
    unified = unified.loc[(unified.taxon_id == 9606) & ~unified.ambiguous]
    allowed = {
        "phospho": {"S", "T"},
        "acetylation": {"K"},
        "ubiquitination": {"K"},
    }
    sites = {"oglcnac": oglcnac}
    audits = {"oglcnac_atlas": oglcnac_audit}
    for ptm, residues in allowed.items():
        subset = unified.loc[
            unified.mod_type.eq(ptm), ["accession", "position", "residue"]
        ].drop_duplicates()
        sites[ptm], audits[ptm] = validate_site_map(subset, sequences, residues)

    provenance = {
        "inputs": {
            str(fasta_path.relative_to(repo)): sha256(fasta_path),
            str(unified_path.relative_to(repo)): sha256(unified_path),
            str(atlas_path.relative_to(repo)): sha256(atlas_path),
        },
        "audits": audits,
        "n_canonical_sequences": len(sequences),
    }
    return sequences, sites, provenance


def all_eligible(sequence: str, ptm: str) -> np.ndarray:
    residues = {"S", "T"} if ptm in {"oglcnac", "phospho"} else {"K"}
    return np.fromiter(
        (index for index, residue in enumerate(sequence, start=1) if residue in residues),
        dtype=np.int64,
    )


def observed_union_pools(
    sites: dict[str, dict[str, np.ndarray]],
) -> dict[str, dict[str, np.ndarray]]:
    pools: dict[str, dict[str, np.ndarray]] = {ptm: {} for ptm in PTMS}
    for first, second in (("oglcnac", "phospho"), ("acetylation", "ubiquitination")):
        for accession in sorted(set(sites[first]) & set(sites[second])):
            union = np.union1d(sites[first][accession], sites[second][accession])
            # Require opportunity beyond each PTM's already observed set.
            for ptm in (first, second):
                if len(union) > len(sites[ptm][accession]):
                    pools[ptm][accession] = union
    return pools


def protein_rows(
    sequences: dict[str, str],
    sites: dict[str, dict[str, np.ndarray]],
    *,
    ptms: tuple[str, ...] = PTMS,
) -> pd.DataFrame:
    union_pools = observed_union_pools(sites)
    rows = []
    for universe in ("all_eligible_residues", "observed_residue_union"):
        for ptm in ptms:
            for accession, observed in sites[ptm].items():
                if len(observed) < MIN_SITES:
                    continue
                if universe == "all_eligible_residues":
                    eligible = all_eligible(sequences[accession], ptm)
                else:
                    eligible = union_pools[ptm].get(accession)
                    if eligible is None:
                        continue
                n, k = len(eligible), len(observed)
                if n < k or k < MIN_SITES:
                    raise AssertionError((universe, ptm, accession, n, k))
                site_pairs = math.comb(k, 2)
                eligible_pairs = math.comb(n, 2)
                for radius in RADII:
                    observed_edges, _ = close_pair_graph(observed, radius)
                    eligible_edges, adjacent = close_pair_graph(eligible, radius)
                    null_edges, null_edge_variance = null_pair_moments(
                        n, k, eligible_edges, adjacent
                    )
                    observed_fraction = observed_edges / site_pairs
                    null_fraction = eligible_edges / eligible_pairs
                    # These must agree algebraically; retain the assertion as a guardrail.
                    if not np.isclose(null_edges / site_pairs, null_fraction, atol=1e-12):
                        raise AssertionError("Pair-null expectation mismatch")
                    rows.append(
                        {
                            "universe": universe,
                            "accession": accession,
                            "ptm": ptm,
                            "radius": radius,
                            "n_sites": k,
                            "n_eligible": n,
                            "observed_close_pairs": observed_edges,
                            "possible_site_pairs": site_pairs,
                            "eligible_close_pairs": eligible_edges,
                            "possible_eligible_pairs": eligible_pairs,
                            "observed_close_pair_fraction": observed_fraction,
                            "null_close_pair_fraction": null_fraction,
                            "effect": observed_fraction - null_fraction,
                            "null_variance_effect": null_edge_variance / (site_pairs**2),
                        }
                    )
    return pd.DataFrame(rows)


def holm_adjust(values: pd.Series) -> np.ndarray:
    p = values.to_numpy(dtype=float)
    order = np.argsort(p)
    ranked = p[order]
    adjusted = np.maximum.accumulate(ranked * (len(p) - np.arange(len(p))))
    output = np.empty_like(adjusted)
    output[order] = np.minimum(adjusted, 1.0)
    return output


def bh_adjust(values: pd.Series) -> np.ndarray:
    p = values.to_numpy(dtype=float)
    order = np.argsort(p)
    ranked = p[order]
    adjusted = np.minimum.accumulate((ranked * len(p) / np.arange(1, len(p) + 1))[::-1])[::-1]
    output = np.empty_like(adjusted)
    output[order] = np.minimum(adjusted, 1.0)
    return output


def summarize(per_protein: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for (universe, ptm, radius), group in per_protein.groupby(
        ["universe", "ptm", "radius"], sort=True
    ):
        effects = group.effect.to_numpy(dtype=float)
        n = len(group)
        mean_effect = float(effects.mean())
        sample_se = float(effects.std(ddof=1) / math.sqrt(n)) if n > 1 else np.nan
        t_stat = mean_effect / sample_se if sample_se > 0 else np.inf
        population_p = float(stats.t.sf(t_stat, df=n - 1)) if n > 1 else np.nan
        exact_null_se = math.sqrt(float(group.null_variance_effect.sum())) / n
        exact_z = mean_effect / exact_null_se if exact_null_se > 0 else np.inf
        exact_p = float(stats.norm.sf(exact_z))
        exact_log10_p = float(stats.norm.logsf(exact_z) / math.log(10))
        observed = float(group.observed_close_pair_fraction.mean())
        null = float(group.null_close_pair_fraction.mean())
        t_crit = float(stats.t.ppf(0.975, df=n - 1)) if n > 1 else np.nan
        rows.append(
            {
                "universe": universe,
                "ptm": ptm,
                "radius": int(radius),
                "n_proteins": n,
                "n_sites": int(group.n_sites.sum()),
                "mean_observed_close_pair_fraction": observed,
                "mean_null_close_pair_fraction": null,
                "mean_effect": mean_effect,
                "fold_of_equal_protein_means": observed / null if null else np.nan,
                "population_ci_low": mean_effect - t_crit * sample_se,
                "population_ci_high": mean_effect + t_crit * sample_se,
                "population_t_p_greater": population_p,
                "population_log10_p_greater": float(
                    stats.t.logsf(t_stat, df=n - 1) / math.log(10)
                ),
                "exact_selection_null_se": exact_null_se,
                "exact_selection_z": exact_z,
                "exact_selection_p_greater": exact_p,
                "exact_selection_log10_p_greater": exact_log10_p,
                "weighting": "equal_protein",
            }
        )
    summary = pd.DataFrame(rows).sort_values(["universe", "radius", "ptm"])
    summary["global_bh_q_population"] = bh_adjust(summary.population_t_p_greater)
    summary["global_bh_q_exact_null"] = bh_adjust(summary.exact_selection_p_greater)
    summary["primary_holm_p_population"] = np.nan
    summary["primary_holm_p_exact_null"] = np.nan
    mask = summary.universe.eq("all_eligible_residues") & summary.radius.eq(PRIMARY_RADIUS)
    summary.loc[mask, "primary_holm_p_population"] = holm_adjust(
        summary.loc[mask, "population_t_p_greater"]
    )
    summary.loc[mask, "primary_holm_p_exact_null"] = holm_adjust(
        summary.loc[mask, "exact_selection_p_greater"]
    )
    return summary.reset_index(drop=True)


def sign_flip_p(differences: np.ndarray, permutations: int, seed: int) -> float:
    rng = np.random.default_rng(seed)
    observed = abs(float(differences.mean()))
    extreme = 0
    complete = 0
    chunk = 1000
    while complete < permutations:
        size = min(chunk, permutations - complete)
        signs = rng.choice(np.array([-1.0, 1.0]), size=(size, len(differences)))
        values = np.abs((signs * differences).mean(axis=1))
        extreme += int(np.sum(values >= observed))
        complete += size
    return (extreme + 1) / (permutations + 1)


def bootstrap_ci(differences: np.ndarray, permutations: int, seed: int) -> tuple[float, float]:
    rng = np.random.default_rng(seed)
    output = np.empty(permutations, dtype=float)
    complete = 0
    chunk = 1000
    while complete < permutations:
        size = min(chunk, permutations - complete)
        indices = rng.integers(0, len(differences), size=(size, len(differences)))
        output[complete : complete + size] = differences[indices].mean(axis=1)
        complete += size
    return float(np.quantile(output, 0.025)), float(np.quantile(output, 0.975))


def stable_seed(*tokens: object) -> int:
    payload = "|".join(map(str, tokens)).encode()
    return (BASE_SEED + int(hashlib.sha256(payload).hexdigest()[:8], 16)) % (2**32)


def direct_contrasts(per_protein: pd.DataFrame, permutations: int) -> pd.DataFrame:
    rows = []
    for universe in ("all_eligible_residues", "observed_residue_union"):
        subset = per_protein.loc[per_protein.universe.eq(universe)]
        for radius in RADII:
            wide = subset.loc[subset.radius.eq(radius)].pivot(
                index="accession", columns="ptm", values="effect"
            )
            for control in ("phospho", "acetylation", "ubiquitination"):
                matched = wide.dropna(subset=["oglcnac", control])
                if matched.empty:
                    continue
                differences = (matched.oglcnac - matched[control]).to_numpy(dtype=float)
                ci_low, ci_high = bootstrap_ci(
                    differences, permutations, stable_seed(universe, radius, control, "boot")
                )
                rows.append(
                    {
                        "universe": universe,
                        "radius": radius,
                        "contrast": f"oglcnac_minus_{control}",
                        "mean_effect_difference": float(differences.mean()),
                        "bootstrap_ci_low": ci_low,
                        "bootstrap_ci_high": ci_high,
                        "sign_flip_p_two_sided": sign_flip_p(
                            differences,
                            permutations,
                            stable_seed(universe, radius, control, "sign"),
                        ),
                        "n_common_proteins": len(differences),
                        "n_permutations": permutations,
                    }
                )
    output = pd.DataFrame(rows).sort_values(["universe", "radius", "contrast"])
    output["global_bh_q"] = bh_adjust(output.sign_flip_p_two_sided)
    output["primary_holm_p"] = np.nan
    mask = output.universe.eq("all_eligible_residues") & output.radius.eq(PRIMARY_RADIUS)
    output.loc[mask, "primary_holm_p"] = holm_adjust(output.loc[mask, "sign_flip_p_two_sided"])
    return output.reset_index(drop=True)


def matched_fold_ratios(per_protein: pd.DataFrame, permutations: int) -> pd.DataFrame:
    """Ratio of equal-protein fold enrichments with paired protein bootstrap.

    The folds are recomputed inside every joint protein resample. This retains
    proteins with zero observed close pairs, for which a literal per-protein fold
    ratio would be zero or infinite and would require an arbitrary pseudocount.
    """
    rows = []
    for universe in ("all_eligible_residues", "observed_residue_union"):
        subset = per_protein.loc[per_protein.universe.eq(universe)]
        for radius in RADII:
            radius_frame = subset.loc[subset.radius.eq(radius)]
            observed = radius_frame.pivot(
                index="accession", columns="ptm", values="observed_close_pair_fraction"
            )
            expected = radius_frame.pivot(
                index="accession", columns="ptm", values="null_close_pair_fraction"
            )
            for control in ("phospho", "acetylation", "ubiquitination"):
                accessions = observed.dropna(subset=["oglcnac", control]).index
                accessions = accessions.intersection(
                    expected.dropna(subset=["oglcnac", control]).index
                )
                if accessions.empty:
                    continue
                og_obs = observed.loc[accessions, "oglcnac"].to_numpy(float)
                og_null = expected.loc[accessions, "oglcnac"].to_numpy(float)
                ct_obs = observed.loc[accessions, control].to_numpy(float)
                ct_null = expected.loc[accessions, control].to_numpy(float)
                og_fold = float(og_obs.mean() / og_null.mean())
                control_fold = float(ct_obs.mean() / ct_null.mean())
                ratio = og_fold / control_fold
                rng = np.random.default_rng(
                    stable_seed(universe, radius, control, "fold_ratio_boot")
                )
                boot = np.empty(permutations, dtype=float)
                complete = 0
                chunk = 1000
                while complete < permutations:
                    size = min(chunk, permutations - complete)
                    indices = rng.integers(0, len(accessions), size=(size, len(accessions)))
                    boot_og_fold = og_obs[indices].mean(axis=1) / og_null[indices].mean(axis=1)
                    boot_ct_fold = ct_obs[indices].mean(axis=1) / ct_null[indices].mean(axis=1)
                    boot[complete : complete + size] = boot_og_fold / boot_ct_fold
                    complete += size
                rows.append(
                    {
                        "universe": universe,
                        "radius": radius,
                        "contrast": f"oglcnac_over_{control}",
                        "oglcnac_fold": og_fold,
                        "control_fold": control_fold,
                        "fold_ratio": ratio,
                        "bootstrap_ci_low": float(np.quantile(boot, 0.025)),
                        "bootstrap_ci_high": float(np.quantile(boot, 0.975)),
                        "n_common_proteins": len(accessions),
                        "n_bootstrap": permutations,
                        "resampling_unit": "matched_protein",
                        "estimand": "ratio_of_equal_protein_mean_fold_enrichments",
                    }
                )
    return (
        pd.DataFrame(rows)
        .sort_values(["universe", "radius", "contrast"])
        .reset_index(drop=True)
    )


def self_test() -> None:
    positions = np.asarray([1, 3, 10], dtype=np.int64)
    assert close_pair_graph(positions, 2) == (1, 0)
    assert close_pair_graph(positions, 9) == (3, 3)
    mean, variance = null_pair_moments(3, 3, 1, 0)
    assert np.isclose(mean, 1.0) and np.isclose(variance, 0.0)
    mean, variance = null_pair_moments(4, 2, 2, 1)
    # Brute force: the six pairs contain two close edges, hence Bernoulli(1/3).
    assert np.isclose(mean, 1 / 3) and np.isclose(variance, 2 / 9)


def render_report(
    summary: pd.DataFrame,
    contrasts: pd.DataFrame,
    fold_ratios: pd.DataFrame,
    provenance: dict,
) -> str:
    primary = summary.loc[
        summary.universe.eq("all_eligible_residues") & summary.radius.eq(PRIMARY_RADIUS)
    ].sort_values("ptm")
    direct = contrasts.loc[
        contrasts.universe.eq("all_eligible_residues") & contrasts.radius.eq(PRIMARY_RADIUS)
    ].sort_values("contrast")
    ratios = fold_ratios.loc[
        fold_ratios.universe.eq("all_eligible_residues") & fold_ratios.radius.eq(PRIMARY_RADIUS)
    ].sort_values("contrast")
    sensitivity = summary.loc[
        summary.universe.eq("observed_residue_union") & summary.radius.eq(PRIMARY_RADIUS)
    ].sort_values("ptm")
    lines = [
        "# PTM Self-Clustering Rerun",
        "",
        f"Generated: {datetime.now(UTC).isoformat()}",
        "",
        "## Question and estimand",
        "",
        "This tests whether sites of each PTM lie close to other sites of the same PTM.",
        "For every protein with at least three sites, the endpoint is the fraction of all",
        "same-PTM site pairs separated by at most 10 residues minus its exact expectation",
        "when the same number of sites is sampled without replacement from eligible S/T",
        "or K residues in that protein. Proteins receive equal weight.",
        "",
        "## Primary result: all eligible residues",
        "",
        "| PTM | Observed pair fraction | Null | Excess | Fold | 95% population CI | "
        "Exact-null log10(p) | Population p | Holm population p | Proteins | Sites |",
        "|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for row in primary.itertuples(index=False):
        lines.append(
            f"| {row.ptm} | {row.mean_observed_close_pair_fraction:.4f} | "
            f"{row.mean_null_close_pair_fraction:.4f} | {row.mean_effect:.4f} | "
            f"{row.fold_of_equal_protein_means:.3f} | "
            f"[{row.population_ci_low:.4f}, {row.population_ci_high:.4f}] | "
            f"{row.exact_selection_log10_p_greater:.1f} | {row.population_t_p_greater:.3g} | "
            f"{row.primary_holm_p_population:.3g} | {row.n_proteins} | {row.n_sites} |"
        )
    lines.extend(
        [
            "",
            "## Direct O-GlcNAc contrasts on the same proteins",
            "",
            "The primary comparison is the fold-enrichment ratio, recomputed in",
            "20,000 joint matched-protein bootstrap samples. No sites are pooled as",
            "independent observations.",
            "",
            "| Contrast | O-GlcNAc fold | Control fold | Fold ratio | "
            "95% protein-bootstrap CI | Common proteins |",
            "|---|---:|---:|---:|---:|---:|",
        ]
    )
    for row in ratios.itertuples(index=False):
        lines.append(
            f"| {row.contrast} | {row.oglcnac_fold:.3f} | {row.control_fold:.3f} | "
            f"{row.fold_ratio:.3f} | [{row.bootstrap_ci_low:.3f}, "
            f"{row.bootstrap_ci_high:.3f}] | {row.n_common_proteins} |"
        )
    lines.extend(
        [
            "",
            "Because individual proteins can have zero observed close pairs, literal",
            "per-protein fold ratios can be infinite. Joint protein resampling estimates",
            "the requested ratio without dropping zero proteins or adding a pseudocount.",
            "The same-chemistry O-GlcNAc/phospho comparison is primary. K-based controls",
            "are secondary because K and S/T have different sequence geometry.",
            "",
            "The corresponding additive null-adjusted comparisons are:",
            "",
            "| Contrast | O-GlcNAc excess minus control excess | 95% bootstrap CI | "
            "Sign-flip p | Holm p | Common proteins |",
            "|---|---:|---:|---:|---:|---:|",
        ]
    )
    for row in direct.itertuples(index=False):
        lines.append(
            f"| {row.contrast} | {row.mean_effect_difference:.4f} | "
            f"[{row.bootstrap_ci_low:.4f}, {row.bootstrap_ci_high:.4f}] | "
            f"{row.sign_flip_p_two_sided:.3g} | {row.primary_holm_p:.3g} | "
            f"{row.n_common_proteins} |"
        )
    lines.extend(
        [
            "",
            "## Detectability-sensitive check: observed-residue union",
            "",
            "Here S/T opportunities are restricted to the union of observed O-GlcNAc and",
            "phosphosites, and K opportunities to the union of observed acetylation and",
            "ubiquitination sites. Proteins must contain both PTM types, so this is a more",
            "conservative but smaller analysis.",
            "",
            "| PTM | Observed | Null | Excess | Fold | 95% population CI | Proteins |",
            "|---|---:|---:|---:|---:|---:|---:|",
        ]
    )
    for row in sensitivity.itertuples(index=False):
        lines.append(
            f"| {row.ptm} | {row.mean_observed_close_pair_fraction:.4f} | "
            f"{row.mean_null_close_pair_fraction:.4f} | {row.mean_effect:.4f} | "
            f"{row.fold_of_equal_protein_means:.3f} | "
            f"[{row.population_ci_low:.4f}, {row.population_ci_high:.4f}] | "
            f"{row.n_proteins} |"
        )
    lines.extend(
        [
            "",
            "## Guardrails",
            "",
            "- This is self-clustering, not distance from O-GlcNAc to a different PTM.",
            "- The exact selection null conditions on protein, PTM site count, eligible",
            "  residue chemistry, and the actual spacing of eligible residues.",
            "- The population interval treats proteins as independent observational units.",
            "- Catalog ascertainment can create apparent clustering; the observed-union",
            "  sensitivity reduces but cannot eliminate study and peptide-detection bias.",
            "- Radii 5, 15, 20, 25, and 30 are sensitivity checks; 10 is primary.",
            "",
            "## Provenance",
            "",
            f"- Direct-contrast permutations: {provenance['config']['permutations']:,}",
            f"- Base seed: {BASE_SEED}",
        ]
    )
    for path, digest in provenance["inputs"].items():
        lines.append(f"- `{path}`: `{digest}`")
    lines.append(f"- Script SHA-256: `{provenance['script_sha256']}`")
    return "\n".join(lines) + "\n"


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    config = load_config(args.config)
    permutations = int(config.values["randomness"]["bootstrap_draws"])
    configured_seed = int(config.values["randomness"]["archived_clustering_seed"])
    if configured_seed != BASE_SEED:
        raise ValueError(
            f"Archived clustering seed must remain {BASE_SEED}; received {configured_seed}"
        )
    if permutations < 100:
        raise SystemExit("At least 100 direct-contrast permutations are required")
    self_test()
    repo = config.source_root
    output = args.output_dir.resolve()
    output.mkdir(parents=True, exist_ok=True)

    sequences, sites, provenance = load_inputs(repo)
    per_protein = protein_rows(sequences, sites)
    summary = summarize(per_protein)
    contrasts = direct_contrasts(per_protein, permutations)
    fold_ratios = matched_fold_ratios(per_protein, permutations)

    script_path = Path(__file__).resolve()
    provenance.update(
        {
            "script": str(script_path.relative_to(config.project_root)),
            "script_sha256": sha256(script_path),
            "generated_utc": datetime.now(UTC).isoformat(),
            "config": {
                "ptms": PTMS,
                "radii": RADII,
                "primary_radius": PRIMARY_RADIUS,
                "minimum_sites_per_protein": MIN_SITES,
                "permutations": permutations,
                "base_seed": BASE_SEED,
                "weighting": "equal_protein",
                "primary_null": (
                    "fixed-count sampling without replacement from all eligible residues"
                ),
                "sensitivity_null": (
                    "fixed-count sampling from observed same-chemistry PTM union"
                ),
            },
            "environment": {
                "python": platform.python_version(),
                "platform": platform.platform(),
                "numpy": np.__version__,
                "pandas": pd.__version__,
                "scipy": scipy.__version__,
            },
        }
    )

    # CSV is the explicit Python boundary. DuckDB performs every Parquet
    # materialisation declared by the workflow.
    per_protein.to_csv(output / "per_protein.csv", index=False)
    summary.to_csv(output / "summary.csv", index=False)
    contrasts.to_csv(output / "contrasts.csv", index=False)
    fold_ratios.to_csv(output / "fold_ratios.csv", index=False)
    (output / "provenance.json").write_text(
        json.dumps(provenance, indent=2, sort_keys=True) + "\n"
    )
    (output / "RUN_REPORT.md").write_text(
        render_report(summary, contrasts, fold_ratios, provenance)
    )
    print(f"wrote {len(per_protein):,} protein-scale rows to {output}", flush=True)


if __name__ == "__main__":
    main()
