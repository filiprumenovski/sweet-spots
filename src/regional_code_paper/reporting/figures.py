"""Render every main and supplemental figure from declared workflow products.

The plotting layer performs display transformations only.  Statistical
estimation lives in the analysis modules; this file reads their tidy tables and
draws them without synthetic curves, jittered pseudo-data, or copied point
estimates.  PDF is the archival output and PNG is a review convenience.
"""

from __future__ import annotations

import argparse
import json
from collections.abc import Callable
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.axes import Axes
from matplotlib.figure import Figure
from sklearn.metrics import roc_curve

from ..core.config import load_config
from ..core.io import sha256_file, write_json

MM = 1 / 25.4
WIDTH = 174 * MM
INK = "#20242a"
BLUE = "#21618c"
BLUE_LIGHT = "#dbeaf3"
RED = "#bb4d3e"
TEAL = "#188977"
GOLD = "#c18b2f"
GREY = "#8b9298"
LIGHT = "#e8ebed"


def configure_style() -> None:
    """Set a compact journal style without mutating scientific content."""
    plt.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": ["Arial", "Liberation Sans", "DejaVu Sans"],
            "font.size": 6.2,
            "axes.labelsize": 6.4,
            "xtick.labelsize": 5.8,
            "ytick.labelsize": 5.8,
            "legend.fontsize": 5.6,
            "axes.linewidth": 0.55,
            "lines.linewidth": 0.85,
            "savefig.dpi": 400,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )


def panel(ax: Axes, label: str) -> None:
    ax.text(-0.16, 1.08, label, transform=ax.transAxes, fontweight="bold", fontsize=8)


def null_line(ax: Axes, value: float = 1.0, *, axis: str = "y") -> None:
    draw = ax.axhline if axis == "y" else ax.axvline
    draw(value, color=GREY, linewidth=0.6, linestyle=(0, (3, 2)), zorder=0)


def clean(ax: Axes) -> None:
    ax.spines[["top", "right"]].set_visible(False)
    ax.tick_params(length=2.2, width=0.5)


def save(fig: Figure, output: Path, name: str) -> list[Path]:
    output.mkdir(parents=True, exist_ok=True)
    paths = [output / f"{name}.pdf", output / f"{name}.png"]
    for path in paths:
        fig.savefig(path, bbox_inches="tight", pad_inches=0.03)
    plt.close(fig)
    return paths


def read_json(path: Path) -> dict:
    with path.open(encoding="utf-8") as handle:
        return json.load(handle)


def figure_1(root: Path) -> Figure:
    regions = pd.read_csv(root / "analysis/consensus_regions/consensus_regions.csv")
    region_summary = read_json(root / "analysis/consensus_regions/summary.json")
    proteins = pd.read_csv(root / "analysis/self_clustering/per_protein.csv")
    ratios = pd.read_csv(root / "analysis/self_clustering/fold_ratios.csv")
    evidence = pd.read_csv(root / "analysis/evidence_restricted_clustering/per_scale.csv")
    detection = read_json(root / "analysis/detection_aware_clustering/summary.json")["primary"]
    deletion = pd.read_csv(root / "analysis/adversarial_deletion/adversarial_curve.csv")
    envelope = pd.read_csv(root / "analysis/adversarial_deletion/random_envelope.csv")

    fig, axes = plt.subplots(2, 3, figsize=(WIDTH, 104 * MM), constrained_layout=True)
    ax = axes[0, 0]
    panel(ax, "a")
    values = [
        region_summary["audit"]["validated_sites"],
        region_summary["results"]["region_sites"],
        region_summary["results"]["regions"],
    ]
    labels = [
        "validated\nO-GlcNAc sites",
        "strict core\nsites (gap 5)",
        "reported\nregions (gap 10)",
    ]
    ax.bar(range(3), values, color=[GREY, BLUE, TEAL], width=0.62)
    for index, value in enumerate(values):
        ax.text(index, value * 1.03, f"{value:,}", ha="center", fontsize=5.5)
    ax.set_xticks(range(3), labels)
    ax.set_ylabel("count")
    clean(ax)

    ax = axes[0, 1]
    panel(ax, "b")
    primary = proteins.loc[
        proteins.universe.eq("all_eligible_residues") & proteins.radius.eq(10)
    ].copy()
    primary["fold"] = primary.observed_close_pair_fraction / primary.null_close_pair_fraction
    paired = primary.pivot(index="accession", columns="ptm", values="fold").dropna(
        subset=["oglcnac", "phospho"]
    )
    visible = paired.loc[(paired.oglcnac > 0) & (paired.phospho > 0)]
    ax.scatter(visible.phospho, visible.oglcnac, s=4, color=BLUE, alpha=0.35, linewidth=0)
    bounds = (0.15, 150)
    ax.plot(bounds, bounds, color=GREY, linestyle=(0, (3, 2)), linewidth=0.6)
    ax.set(xscale="log", yscale="log", xlim=bounds, ylim=bounds)
    ax.set_xlabel("phosphorylation close-pair fold")
    ax.set_ylabel("O-GlcNAc close-pair fold")
    clean(ax)

    ax = axes[0, 2]
    panel(ax, "c")
    curve = ratios.loc[
        ratios.universe.eq("all_eligible_residues") & ratios.contrast.eq("oglcnac_over_phospho")
    ].sort_values("radius")
    ax.fill_between(
        curve.radius,
        curve.bootstrap_ci_low,
        curve.bootstrap_ci_high,
        color=BLUE_LIGHT,
    )
    ax.plot(curve.radius, curve.fold_ratio, "o-", color=BLUE, markersize=3)
    null_line(ax)
    ax.set(xlabel="radius (residues)", ylabel="O-GlcNAc / phosphorylation")
    clean(ax)

    ax = axes[1, 0]
    panel(ax, "d")
    evidence_primary = evidence.loc[evidence.radius.eq(10)].set_index("restriction")
    conditions = [
        ("all strict\nsites", evidence_primary.loc["all_strict_sites"]),
        (
            "≥2 independent\npublications",
            evidence_primary.loc["two_independent_publications"],
        ),
        (
            "≥2 mapped\npeptides",
            evidence_primary.loc["two_distinct_mapped_peptides"],
        ),
    ]
    values = [float(item.group_fold) for _, item in conditions] + [
        float(detection["group_fold"])
    ]
    low = [float(item.group_fold_ci_low) for _, item in conditions] + [
        float(detection["group_fold_ci_low"])
    ]
    high = [float(item.group_fold_ci_high) for _, item in conditions] + [
        float(detection["group_fold_ci_high"])
    ]
    positive = [float(item.positive_excess_fraction) for _, item in conditions] + [
        float(detection["positive_excess_fraction"])
    ]
    labels = [
        "all validated sites",
        "≥2 publications",
        "≥2 mapped peptides",
        "peptide-conditioned",
    ]
    y = np.arange(len(labels))[::-1]
    ax.barh(y, values, color=[BLUE, TEAL, TEAL, GOLD], height=0.58)
    ax.errorbar(
        values,
        y,
        xerr=[np.asarray(values) - low, high - np.asarray(values)],
        fmt="none",
        ecolor=INK,
        elinewidth=0.65,
        capsize=2,
    )
    for index, (value, fraction) in enumerate(zip(values, positive, strict=True)):
        ax.text(
            high[index] + 0.3,
            y[index],
            f"{value:.1f}-fold; {fraction:.0%} positive",
            ha="left",
            va="center",
            fontsize=4.7,
        )
    null_line(ax, axis="x")
    ax.set_yticks(y, labels)
    ax.set(xlabel="group close-pair fold", xlim=(0, 17.5))
    clean(ax)

    ax = axes[1, 1]
    panel(ax, "e")
    size = np.clip(regions.n_contributing_pmids.to_numpy(float), 1, 8)
    ax.scatter(regions.span, regions.valence, s=3 + 2 * size, alpha=0.35, color=BLUE)
    ax.set(xlabel="region span (residues)", ylabel="validated sites per region")
    ax.set_xscale("log")
    clean(ax)

    ax = axes[1, 2]
    panel(ax, "f")
    ax.fill_between(
        envelope.removed_pct,
        envelope.fold_q025,
        envelope.fold_q975,
        color=LIGHT,
        label="random removal, 95% envelope",
    )
    ax.plot(deletion.removed_pct, deletion.group_fold, color=BLUE, label="strongest first")
    null_line(ax)
    ax.set(xlabel="proteins removed (%)", ylabel="O-GlcNAc group fold", yscale="log")
    ax.legend(frameon=False)
    clean(ax)
    return fig


def draw_roc(ax: Axes, labels: np.ndarray, scores: np.ndarray, label: str, color: str) -> None:
    false_positive, true_positive, _ = roc_curve(labels, scores)
    ax.plot(false_positive, true_positive, label=label, color=color)


def figure_2(root: Path) -> Figure:
    tiles = pd.read_csv(root / "analysis/regional_code/held_out_tiles.csv")
    coefficients = pd.read_csv(root / "analysis/regional_code/acceptor_model.csv")
    per_protein = pd.read_csv(root / "analysis/regional_code/per_protein.csv")
    summary = read_json(root / "analysis/regional_code/summary.json")
    regions = pd.read_csv(root / "analysis/consensus_regions/consensus_regions.csv")

    fig = plt.figure(figsize=(WIDTH, 102 * MM), constrained_layout=True)
    grid = fig.add_gridspec(2, 6)
    axes = [
        fig.add_subplot(grid[0, :2]),
        fig.add_subplot(grid[0, 2:4]),
        fig.add_subplot(grid[0, 4:]),
        fig.add_subplot(grid[1, :3]),
        fig.add_subplot(grid[1, 3:]),
    ]
    for ax, label in zip(axes, "abcde", strict=True):
        panel(ax, label)

    draw_roc(axes[0], tiles.label, tiles.composition_score, "composition", BLUE)
    draw_roc(axes[0], tiles.label, tiles.geometry_score, "acceptor geometry", GOLD)
    draw_roc(axes[0], tiles.label, tiles.full_score, "full", TEAL)
    axes[0].plot([0, 1], [0, 1], color=GREY, linestyle=(0, (2, 2)), linewidth=0.6)
    axes[0].set(xlabel="false-positive rate", ylabel="true-positive rate")
    axes[0].legend(frameon=False)

    y = np.arange(len(coefficients))[::-1]
    axes[1].errorbar(
        coefficients.coefficient,
        y,
        xerr=[
            coefficients.coefficient - coefficients.ci_low,
            coefficients.ci_high - coefficients.coefficient,
        ],
        fmt="o",
        color=BLUE,
        capsize=2,
    )
    axes[1].axvline(0, color=GREY, linewidth=0.6)
    axes[1].set_yticks(y, coefficients.outcome)
    axes[1].set_xlabel("cluster coefficient")

    axes[2].hist(
        per_protein.composition_within_protein_auroc,
        bins=24,
        color=BLUE_LIGHT,
        edgecolor=BLUE,
    )
    axes[2].axvline(summary["composition_within_protein"]["mean_auroc"], color=BLUE)
    axes[2].set(xlabel="composition-only within-protein AUROC", ylabel="proteins")

    informative = (
        tiles.groupby("accession")
        .filter(lambda group: group.label.nunique() == 2)
        .groupby("accession")
        .agg(
            n=("label", "size"),
            positives=("label", "sum"),
            spread=("composition_score", "std"),
        )
    )
    selected = informative.sort_values(["positives", "spread"], ascending=False).head(6).index
    for row, accession in enumerate(selected):
        group = tiles.loc[tiles.accession.eq(accession)].sort_values("center")
        offset = row * 1.1
        axes[3].plot(
            group.center,
            group.composition_score + offset,
            color=BLUE,
            linewidth=0.7,
        )
        axes[3].scatter(
            group.loc[group.label.eq(1), "center"],
            group.loc[group.label.eq(1), "composition_score"] + offset,
            s=4,
            color=RED,
        )
        axes[3].text(group.center.min(), offset + 0.85, accession, fontsize=4.5)
    axes[3].set(xlabel="residue position", ylabel="held-out score, offset by protein")
    axes[3].set_yticks([])

    positive = tiles.loc[tiles.label.eq(1)]
    negative = tiles.loc[tiles.label.eq(0)]
    metrics = {
        "S/T/P fraction": (
            positive.st_fraction + positive.proline,
            negative.st_fraction + negative.proline,
        ),
        "aromatic fraction": (positive.aromatic, negative.aromatic),
        "mean S/T gap": (positive.geo_gapmean, negative.geo_gapmean),
    }
    x = np.arange(len(metrics))
    axes[4].bar(x - 0.18, [item[0].mean() for item in metrics.values()], 0.36, color=BLUE)
    axes[4].bar(x + 0.18, [item[1].mean() for item in metrics.values()], 0.36, color=LIGHT)
    axes[4].set_xticks(x, metrics, rotation=18, ha="right")
    axes[4].legend(["regional tiles", "same-protein background"], frameon=False)
    axes[4].text(
        0.98,
        0.96,
        f"{len(regions):,} regions",
        transform=axes[4].transAxes,
        ha="right",
        va="top",
        fontsize=5,
    )
    for ax in axes:
        clean(ax)
    return fig


def figure_3(root: Path) -> Figure:
    clustering = pd.read_csv(root / "analysis/self_clustering/summary.csv")
    profiles = pd.read_csv(root / "analysis/regional_code/proline_profiles.csv")
    kinase = pd.read_csv(root / "analysis/yin_yang/kinase_summary.csv")
    local = read_json(root / "analysis/yin_yang/summary.json")["local_opportunity"]
    nested = read_json(root / "analysis/yin_yang/summary.json")["nested_nulls"]
    distance = pd.read_csv(root / "analysis/regional_code/distance_profile.csv")

    fig, axes = plt.subplots(2, 3, figsize=(WIDTH, 104 * MM), constrained_layout=True)
    for ax, label in zip(axes.ravel(), "abcdef", strict=True):
        panel(ax, label)
    primary = clustering.loc[
        clustering.universe.eq("all_eligible_residues") & clustering.radius.eq(10)
    ]
    axes[0, 0].bar(
        primary.ptm,
        primary.fold_of_equal_protein_means,
        color=[BLUE if item == "oglcnac" else GREY for item in primary.ptm],
    )
    axes[0, 0].tick_params(axis="x", rotation=25)
    axes[0, 0].set_ylabel("group close-pair fold")
    null_line(axes[0, 0])

    offset_columns = [f"offset_{value:+d}" for value in range(-5, 6)]
    for label, group in profiles.groupby("set"):
        axes[0, 1].plot(range(-5, 6), group[offset_columns].mean(), label=label)
    axes[0, 1].set(xlabel="offset from acceptor", ylabel="proline occupancy")
    axes[0, 1].legend(frameon=False)

    axes[0, 2].plot(
        kinase.control_n90,
        kinase.cluster_n90,
        "o",
        color=BLUE,
        alpha=0.25,
        markersize=2,
    )
    limit = max(kinase.control_n90.max(), kinase.cluster_n90.max())
    axes[0, 2].plot([0, limit], [0, limit], color=GREY, linestyle=(0, (3, 2)))
    axes[0, 2].set(xlabel="matched control", ylabel="cluster acceptor")

    rows = []
    for name, item in local.items():
        arm, window, matching = (
            name.split("_w", maxsplit=1)[0],
            name.split("_w")[1].split("_")[0],
            name.endswith("residue_matched"),
        )
        rows.append((arm, int(window), matching, item))
    for index, (arm, _window, matched, item) in enumerate(sorted(rows)):
        marker = "o" if matched else "s"
        color = BLUE if arm == "all" else TEAL
        axes[1, 0].errorbar(
            index,
            item["ratio"],
            yerr=[[item["ratio"] - item["ci"][0]], [item["ci"][1] - item["ratio"]]],
            fmt=marker,
            color=color,
            capsize=1.5,
        )
    null_line(axes[1, 0])
    axes[1, 0].set(xlabel="window / matching configuration", ylabel="observed / expected")
    axes[1, 0].set_xticks([])

    names = ["N0 length", "N1 S/T", "N2 disorder"]
    for offset, arm in enumerate(("all", "regional")):
        items = list(nested[arm].values())
        axes[1, 1].plot(
            np.arange(3) + (offset - 0.5) * 0.08,
            [item["ratio"] for item in items],
            "o-",
            label=arm,
        )
    null_line(axes[1, 1])
    axes[1, 1].set_xticks(range(3), names, rotation=20, ha="right")
    axes[1, 1].set_ylabel("same-site observed / expected")
    axes[1, 1].legend(frameon=False)

    axes[1, 2].fill_between(
        distance.radius, distance.ci_low, distance.ci_high, color=BLUE_LIGHT
    )
    axes[1, 2].plot(distance.radius, distance.difference, color=BLUE)
    axes[1, 2].axhline(0, color=GREY, linewidth=0.6)
    axes[1, 2].set(
        xlabel="distance to region (residues)", ylabel="phosphosite minus control CDF"
    )
    for ax in axes.ravel():
        clean(ax)
    return fig


def figure_4(root: Path) -> Figure:
    categories = pd.read_csv(root / "analysis/ogt_conservation/category_summary.csv")
    conservation = pd.read_csv(root / "analysis/evolution/site_conservation.csv")
    transfer = pd.read_csv(root / "analysis/evolution/transfer.csv")
    nulls = pd.read_csv(root / "analysis/evolution/composition_position_nulls.csv")
    substitutions = pd.read_csv(root / "analysis/evolution/substitutions.csv")

    fig, axes = plt.subplots(2, 3, figsize=(WIDTH, 104 * MM), constrained_layout=True)
    for ax, label in zip(axes.ravel(), "abcdef", strict=True):
        panel(ax, label)
    category_labels = {
        "asn_ladder": "Asn ladder",
        "buried_core": "buried core",
        "catalytic_triad": "catalytic triad",
        "other_surface": "other surface",
        "peptide_contact_channel": "peptide contacts",
        "reading_surface": "reading surface",
    }
    categories = categories.sort_values("identity").assign(
        label=lambda frame: frame.category.map(category_labels)
    )
    axes[0, 0].barh(categories.label, categories.identity * 100, color=BLUE)
    axes[0, 0].set(xlabel="human to mouse identity (%)", xlim=(95, 100.2))

    summary = conservation.groupby("species").agg(
        exact=("exact", "mean"), regional=("regional", "mean")
    )
    order = [
        item for item in ("mouse", "rat", "zebrafish", "fly", "worm") if item in summary.index
    ]
    x = np.arange(len(order))
    axes[0, 1].plot(x, summary.loc[order, "exact"], "o-", label="exact site", color=RED)
    axes[0, 1].plot(x, summary.loc[order, "regional"], "o-", label="within ±5", color=BLUE)
    axes[0, 1].set_xticks(x, order, rotation=25, ha="right")
    axes[0, 1].set_ylabel("fraction conserved")
    axes[0, 1].legend(frameon=False)

    transfer_colors = [BLUE, RED, GOLD, TEAL]
    for color, ((_direction, species), group) in zip(
        transfer_colors,
        transfer.groupby(["direction", "species"], sort=True),
        strict=True,
    ):
        draw_roc(
            axes[0, 2],
            group.label.to_numpy(),
            group.score.to_numpy(),
            f"{species}",
            color,
        )
    axes[0, 2].plot([0, 1], [0, 1], color=GREY, linestyle=(0, (2, 2)))
    axes[0, 2].set(xlabel="false-positive rate", ylabel="true-positive rate")
    axes[0, 2].legend(frameon=False, ncol=2)

    metrics = ["pearson", "cosine", "js", "position"]
    true_values = [
        nulls[f"true_{metric}"].fillna(0).mean()
        if metric == "position"
        else nulls[f"true_{metric}"].mean()
        for metric in metrics
    ]
    null_values = [
        nulls[f"null_a_{metric}"].fillna(0).mean()
        if metric == "position"
        else nulls[f"null_a_{metric}"].mean()
        for metric in metrics
    ]
    axes[1, 0].bar(np.arange(4) - 0.18, true_values, 0.36, color=BLUE, label="orthologous")
    axes[1, 0].bar(
        np.arange(4) + 0.18, null_values, 0.36, color=LIGHT, label="within-protein null"
    )
    axes[1, 0].set_xticks(range(4), metrics, rotation=20)
    axes[1, 0].set_ylabel("similarity / correlation")
    axes[1, 0].legend(frameon=False, loc="lower left")

    excess = np.asarray(true_values) - np.asarray(null_values)
    axes[1, 1].bar(metrics, excess, color=[BLUE, BLUE, BLUE, TEAL])
    axes[1, 1].set_ylabel("excess over within-protein null")
    axes[1, 1].tick_params(axis="x", rotation=20)

    rates = substitutions.groupby(["rung", "arm"]).gain.mean().unstack()
    rates["ratio"] = rates.region / rates.control
    axes[1, 2].bar(rates.index, rates.ratio, color=BLUE)
    null_line(axes[1, 2])
    axes[1, 2].set_ylabel("acidic/aromatic introduction rate ratio")
    axes[1, 2].tick_params(axis="x", rotation=20)
    for ax in axes.ravel():
        clean(ax)
    return fig


def figure_5(root: Path) -> Figure:
    fg = read_json(root / "analysis/fg_nup_recovery/summary.json")
    folds = pd.read_csv(root / "analysis/scanner/outer_folds.csv")
    operating = pd.read_csv(root / "analysis/scanner/operating_points.csv")

    fig, axes = plt.subplots(1, 4, figsize=(WIDTH, 55 * MM), constrained_layout=True)
    for ax, label in zip(axes, "abcd", strict=True):
        panel(ax, label)
    axes[0].bar(
        [0, 1],
        [fg["recovery_fraction"] * 100, fg["background_region_fraction"] * 100],
        color=[BLUE, LIGHT],
    )
    axes[0].set_xticks([0, 1], ["FG Nups", "atlas\nbackground"])
    axes[0].set_ylabel("with a region (%)")
    axes[0].text(
        0,
        fg["recovery_fraction"] * 100 + 2,
        f"{fg['n_recovered']} / {fg['n_present']}",
        ha="center",
    )

    p_values = [fg["hypergeometric_p"], fg["matched_p"], fg["adjusted_p"]]
    axes[1].barh([2, 1, 0], -np.log10(p_values), color=[GREY, TEAL, BLUE])
    axes[1].axvline(-np.log10(0.05), color=GREY, linestyle=(0, (3, 2)))
    axes[1].set_yticks([2, 1, 0], ["hypergeometric", "matched", "adjusted"])
    axes[1].set_xlabel("-log10 p")

    for _, row in folds.iterrows():
        axes[2].plot(
            [0, 1],
            [row.baseline_average_precision, row.enhanced_average_precision],
            "o-",
            color=BLUE,
            alpha=0.7,
        )
    axes[2].set_xticks([0, 1], ["legacy", "enhanced"])
    axes[2].set_ylabel("outer-fold average precision")

    for model, group in operating.groupby("model"):
        axes[3].plot(group.tile_recall, group.precision, "o-", label=model)
    axes[3].set(xlabel="tile recall", ylabel="catalogue precision")
    axes[3].legend(frameon=False)
    for ax in axes:
        clean(ax)
    return fig


def figure_s1(root: Path) -> Figure:
    regions = pd.read_csv(root / "analysis/consensus_regions/consensus_regions.csv")
    ratios = pd.read_csv(root / "analysis/self_clustering/fold_ratios.csv")
    fig, axes = plt.subplots(1, 3, figsize=(WIDTH, 50 * MM), constrained_layout=True)
    for ax, label in zip(axes, "abc", strict=True):
        panel(ax, label)
    axes[0].hist(regions.span, bins=35, color=BLUE_LIGHT, edgecolor=BLUE)
    axes[0].set(xlabel="region length", ylabel="regions")
    axes[1].bar(
        ["≥2 PMIDs", "≥3 PMIDs", "any PMID\nremoved"],
        [
            regions.n_contributing_pmids.ge(2).mean(),
            regions.n_contributing_pmids.ge(3).mean(),
            regions.survives_every_single_pmid_removal.astype(bool).mean(),
        ],
        color=[BLUE, BLUE, TEAL],
    )
    axes[1].set_ylabel("fraction of regions")
    primary = ratios.loc[ratios.contrast.eq("oglcnac_over_phospho")]
    axes[2].errorbar(
        primary.radius,
        primary.fold_ratio,
        yerr=[
            primary.fold_ratio - primary.bootstrap_ci_low,
            primary.bootstrap_ci_high - primary.fold_ratio,
        ],
        fmt="o-",
        color=BLUE,
        capsize=2,
    )
    axes[2].set(xlabel="radius", ylabel="O-GlcNAc / phospho")
    for ax in axes:
        clean(ax)
    return fig


def figure_s2(root: Path) -> Figure:
    deltas = pd.read_csv(root / "analysis/yin_yang/kinase_site_scores.csv").sort_values("delta")
    local = pd.read_csv(root / "analysis/yin_yang/local_opportunity.csv")
    fig, axes = plt.subplots(1, 3, figsize=(WIDTH, 50 * MM), constrained_layout=True)
    for ax, label in zip(axes, "abc", strict=True):
        panel(ax, label)
    axes[0].plot(np.arange(len(deltas)), deltas.delta, color=BLUE)
    axes[0].axhline(0, color=GREY, linewidth=0.6)
    axes[0].set(xlabel="kinases, ranked", ylabel="cluster - control percentile")
    family = deltas.groupby("family").delta.mean().sort_values()
    axes[1].barh(family.index, family, color=BLUE)
    axes[1].set_xlabel("mean percentile shift")
    summary = (
        local.groupby(["window", "residue_matched", "regional"])[["observed", "expected"]]
        .sum()
        .assign(ratio=lambda frame: frame.observed / frame.expected)
        .reset_index()
    )
    for key, group in summary.groupby(["residue_matched", "regional"]):
        axes[2].plot(
            group.window, group.ratio, "o-", label=f"matched={key[0]}, regional={key[1]}"
        )
    null_line(axes[2])
    axes[2].set(xlabel="window", ylabel="observed / expected")
    axes[2].legend(frameon=False, fontsize=4.5)
    for ax in axes:
        clean(ax)
    return fig


def figure_s3(root: Path) -> Figure:
    nulls = pd.read_csv(root / "analysis/evolution/composition_position_nulls.csv")
    conservation = pd.read_csv(root / "analysis/evolution/site_conservation.csv")
    fig, axes = plt.subplots(1, 3, figsize=(WIDTH, 50 * MM), constrained_layout=True)
    for ax, label in zip(axes, "abc", strict=True):
        panel(ax, label)
    metrics = ("pearson", "cosine", "js")
    for index, metric in enumerate(metrics):
        axes[0].boxplot(
            [nulls[f"true_{metric}"].dropna(), nulls[f"null_a_{metric}"].dropna()],
            positions=[index * 3, index * 3 + 1],
            widths=0.7,
            showfliers=False,
        )
    axes[0].set_xticks([0.5, 3.5, 6.5], metrics)
    axes[0].set_ylabel("similarity")
    position = nulls.true_position
    axes[1].hist(
        position.fillna(0), bins=25, color=BLUE_LIGHT, edgecolor=BLUE, label="zero-imputed"
    )
    axes[1].hist(position.dropna(), bins=25, histtype="step", color=RED, label="non-degenerate")
    axes[1].set_xlabel("modified-position correlation")
    axes[1].legend(frameon=False)
    grouped = conservation.groupby("species").agg(
        exact=("exact", "mean"), regional=("regional", "mean")
    )
    grouped.plot.bar(ax=axes[2], color=[RED, BLUE])
    axes[2].set_ylabel("fraction conserved")
    axes[2].legend(frameon=False)
    for ax in axes:
        clean(ax)
    return fig


def figure_s4(root: Path) -> Figure:
    """Region-definition sensitivity across the complete post hoc grid."""
    catalogue = pd.read_csv(
        root / "analysis/region_definition_sensitivity/catalogue_grid.csv"
    )
    models = pd.read_csv(root / "analysis/region_definition_sensitivity/model_grid.csv")
    final_ten = catalogue.loc[catalogue.final_gap.eq(10)]
    fig, axes = plt.subplots(2, 2, figsize=(WIDTH, 90 * MM), constrained_layout=True)
    for ax, label in zip(axes.ravel(), "abcd", strict=True):
        ax.text(0, 1.04, label, transform=ax.transAxes, fontweight="bold", fontsize=8)

    colors = {3: BLUE, 4: TEAL, 5: GOLD}
    for minimum_sites, group in final_ten.groupby("minimum_sites"):
        ordered = group.sort_values("core_gap")
        axes[0, 0].plot(
            ordered.core_gap,
            ordered.core_sites,
            "o-",
            color=colors[int(minimum_sites)],
            label=f"minimum {int(minimum_sites)} sites",
        )
    axes[0, 0].set(xlabel="strict-core gap", ylabel="retained core sites")
    axes[0, 0].set_xticks(sorted(final_ten.core_gap.unique()))
    axes[0, 0].legend(frameon=False)

    jaccard = final_ten.pivot(
        index="minimum_sites", columns="core_gap", values="site_jaccard"
    ).sort_index(ascending=False)
    image = axes[0, 1].imshow(jaccard, vmin=0, vmax=1, cmap="Blues", aspect="auto")
    axes[0, 1].set_xticks(range(len(jaccard.columns)), jaccard.columns)
    axes[0, 1].set_yticks(range(len(jaccard.index)), jaccard.index)
    axes[0, 1].set(xlabel="strict-core gap", ylabel="minimum sites")
    for row in range(jaccard.shape[0]):
        for column in range(jaccard.shape[1]):
            value = jaccard.iloc[row, column]
            axes[0, 1].text(
                column,
                row,
                f"{value:.2f}",
                ha="center",
                va="center",
                color="white" if value > 0.65 else INK,
                fontsize=5,
            )
    fig.colorbar(image, ax=axes[0, 1], label="site Jaccard vs primary")

    for minimum_sites, group in models.groupby("minimum_sites"):
        ordered = group.sort_values("core_gap")
        axes[1, 0].plot(
            ordered.core_gap,
            ordered.composition_within_protein_auroc,
            "o-",
            color=colors[int(minimum_sites)],
            label=f"minimum {int(minimum_sites)} sites",
        )
    axes[1, 0].axhline(0.5, color=GREY, linestyle=(0, (3, 2)), linewidth=0.6)
    axes[1, 0].set(
        xlabel="strict-core gap",
        ylabel="composition-only within-protein AUROC",
        ylim=(0.5, 0.9),
    )
    axes[1, 0].set_xticks(sorted(models.core_gap.unique()))
    primary = models.loc[models.is_primary].iloc[0]
    axes[1, 0].scatter(
        [primary.core_gap],
        [primary.composition_within_protein_auroc],
        s=45,
        facecolors="none",
        edgecolors=RED,
        linewidths=0.9,
        zorder=4,
    )
    axes[1, 0].legend(frameon=False)

    for minimum_sites, group in models.groupby("minimum_sites"):
        ordered = group.sort_values("core_gap")
        axes[1, 1].errorbar(
            ordered.core_gap,
            ordered.adjusted_composition_increment,
            yerr=[
                ordered.adjusted_composition_increment
                - ordered.adjusted_composition_increment_ci_low,
                ordered.adjusted_composition_increment_ci_high
                - ordered.adjusted_composition_increment,
            ],
            fmt="o-",
            capsize=2,
            color=colors[int(minimum_sites)],
            label=f"minimum {int(minimum_sites)} sites",
        )
    axes[1, 1].axhline(0, color=GREY, linestyle=(0, (3, 2)), linewidth=0.6)
    axes[1, 1].set(
        xlabel="strict-core gap",
        ylabel="composition increment in within-protein AUROC",
    )
    axes[1, 1].set_xticks(sorted(models.core_gap.unique()))
    axes[1, 1].scatter(
        [primary.core_gap],
        [primary.adjusted_composition_increment],
        s=45,
        facecolors="none",
        edgecolors=RED,
        linewidths=0.9,
        zorder=4,
    )
    axes[1, 1].legend(frameon=False)

    for ax in axes.ravel():
        clean(ax)
    return fig


BUILDERS: dict[str, Callable[[Path], Figure]] = {
    "Figure_1": figure_1,
    "Figure_2": figure_2,
    "Figure_3": figure_3,
    "Figure_4": figure_4,
    "Figure_5": figure_5,
    "Figure_S1": figure_s1,
    "Figure_S2": figure_s2,
    "Figure_S3": figure_s3,
    "Figure_S4": figure_s4,
}


def render_all(results: Path, output: Path, project_root: Path) -> list[dict[str, object]]:
    manifest = []
    for name, builder in BUILDERS.items():
        paths = save(builder(results), output, name)
        manifest.append(
            {
                "figure": name,
                "files": [
                    {
                        "path": str(path.relative_to(project_root)),
                        "sha256": sha256_file(path),
                        "bytes": path.stat().st_size,
                    }
                    for path in paths
                ],
            }
        )
    return manifest


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    config = load_config(args.config)
    configure_style()
    output = args.output_dir.resolve()
    manifest = render_all(config.results_root, output, config.project_root)
    write_json(output / "figure_manifest.json", {"figures": manifest})


if __name__ == "__main__":
    main()
