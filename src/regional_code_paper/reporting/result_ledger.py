"""Assemble a tidy, machine-readable ledger of manuscript-facing results.

Each row is one estimand.  The table is intentionally boring: no display
rounding, no prose scraping, and no inferred units.  Figures and the audit both
consume this same ledger so a number cannot drift independently in three places.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

import pandas as pd

from ..core.config import load_config
from ..core.io import ensure_parent

FIELDS = (
    "claim_id",
    "figure_panel",
    "analysis",
    "estimand",
    "value",
    "ci_low",
    "ci_high",
    "unit",
    "inferential_unit",
    "source_file",
)


def read_json(path: Path) -> dict[str, Any]:
    with path.open(encoding="utf-8") as handle:
        return json.load(handle)


def source_label(path: Path, results_root: Path) -> str:
    """Return a repository-relative path for a retained result source."""
    return str(Path(results_root.name) / path.relative_to(results_root))


class Ledger:
    """Small helper that makes omissions and accidental schema drift obvious."""

    def __init__(self) -> None:
        self.rows: list[dict[str, object]] = []

    def add(
        self,
        claim_id: str,
        panel: str,
        analysis: str,
        estimand: str,
        value: float,
        *,
        interval: tuple[float, float] | list[float] | None = None,
        unit: str,
        inferential_unit: str,
        source: str,
    ) -> None:
        low, high = interval if interval is not None else (None, None)
        self.rows.append(
            dict(
                zip(
                    FIELDS,
                    (
                        claim_id,
                        panel,
                        analysis,
                        estimand,
                        value,
                        low,
                        high,
                        unit,
                        inferential_unit,
                        source,
                    ),
                    strict=True,
                )
            )
        )


def assemble(root: Path) -> pd.DataFrame:
    ledger = Ledger()
    regions_path = root / "analysis/consensus_regions/summary.json"
    regions = read_json(regions_path)
    for metric, unit in (
        ("regions", "regions"),
        ("region_sites", "sites"),
        ("region_bearing_proteins", "proteins"),
        ("atlas_proteins", "proteins"),
    ):
        ledger.add(
            f"object.{metric}",
            "1A;S1",
            "consensus_regions",
            metric,
            regions["results"][metric],
            unit=unit,
            inferential_unit="site catalogue",
            source=source_label(regions_path, root),
        )

    clustering_path = root / "analysis/self_clustering/fold_ratios.csv"
    clustering = pd.read_csv(clustering_path)
    for row in clustering.loc[
        clustering.contrast.eq("oglcnac_over_phospho")
        & clustering.universe.eq("all_eligible_residues")
    ].to_dict("records"):
        ledger.add(
            f"clustering.r{row['radius']}",
            "1B;1C",
            "self_clustering",
            "O-GlcNAc to phosphorylation close-pair fold ratio",
            row["fold_ratio"],
            interval=[row["bootstrap_ci_low"], row["bootstrap_ci_high"]],
            unit="fold",
            inferential_unit="protein",
            source=source_label(clustering_path, root),
        )

    breadth_path = root / "analysis/clustering_breadth/summary.json"
    breadth = read_json(breadth_path)
    for claim, block, estimand in (
        (
            "primary_positive_fraction",
            "primary",
            "fraction of proteins with positive O-GlcNAc clustering excess",
        ),
        (
            "observed_union_positive_fraction",
            "observed_residue_sensitivity",
            "positive-excess fraction under observed-residue opportunities",
        ),
        (
            "matched_dominance_fraction",
            "matched_oglcnac_dominance",
            "fraction of matched proteins with greater O-GlcNAc than phospho excess",
        ),
        (
            "all_scales_positive_fraction",
            "all_scales",
            "fraction with positive O-GlcNAc excess at every tested radius",
        ),
    ):
        values = breadth[block]
        ledger.add(
            f"clustering_breadth.{claim}",
            "1F",
            "clustering_breadth",
            estimand,
            values["fraction"],
            interval=[values["fraction_ci_low"], values["fraction_ci_high"]],
            unit="fraction",
            inferential_unit="protein",
            source=source_label(breadth_path, root),
        )

    detection_path = root / "analysis/detection_aware_clustering/summary.json"
    detection = read_json(detection_path)["primary"]
    ledger.add(
        "detection_aware_clustering.group_fold",
        "1D;S1",
        "detection_aware_clustering",
        "O-GlcNAc self-clustering fold under peptide-exposure-conditioned null",
        detection["group_fold"],
        interval=[detection["group_fold_ci_low"], detection["group_fold_ci_high"]],
        unit="fold",
        inferential_unit="protein",
        source=source_label(detection_path, root),
    )
    ledger.add(
        "detection_aware_clustering.positive_fraction",
        "1D;S1",
        "detection_aware_clustering",
        "fraction of proteins with positive excess under peptide-aware null",
        detection["positive_excess_fraction"],
        interval=[
            detection["positive_excess_fraction_ci_low"],
            detection["positive_excess_fraction_ci_high"],
        ],
        unit="fraction",
        inferential_unit="protein",
        source=source_label(detection_path, root),
    )

    evidence_path = root / "analysis/evidence_restricted_clustering/summary.json"
    evidence = read_json(evidence_path)
    for key, label in (
        ("two_independent_publications", "at least two independent publications"),
        ("two_distinct_mapped_peptides", "at least two distinct mapped peptides"),
    ):
        values = evidence[key]
        ledger.add(
            f"evidence_restricted_clustering.{key}.group_fold",
            "1D;S1",
            "evidence_restricted_clustering",
            f"O-GlcNAc self-clustering fold among sites supported by {label}",
            values["group_fold"],
            interval=[values["group_fold_ci_low"], values["group_fold_ci_high"]],
            unit="fold",
            inferential_unit="protein",
            source=source_label(evidence_path, root),
        )
        ledger.add(
            f"evidence_restricted_clustering.{key}.positive_fraction",
            "1D;S1",
            "evidence_restricted_clustering",
            f"positive-excess fraction among sites supported by {label}",
            values["positive_excess_fraction"],
            interval=[
                values["positive_excess_fraction_ci_low"],
                values["positive_excess_fraction_ci_high"],
            ],
            unit="fraction",
            inferential_unit="protein",
            source=source_label(evidence_path, root),
        )

    regional_path = root / "analysis/regional_code/summary.json"
    regional = read_json(regional_path)
    tile_values = {
        "full": regional["tiles"]["global_auroc"]["full"],
        "composition": regional["tiles"]["global_auroc"]["composition"],
        "geometry": regional["tiles"]["global_auroc"]["geometry"],
        "composition_matched": regional["tiles"]["matched"]["composition_auroc"],
        "geometry_matched": regional["tiles"]["matched"]["geometry_auroc"],
    }
    for key, value in tile_values.items():
        ledger.add(
            f"regional_code.{key}.auroc",
            "2A",
            "regional_code",
            f"{key} held-out AUROC",
            value,
            unit="AUROC",
            inferential_unit="protein-grouped tile",
            source=source_label(regional_path, root),
        )
    ledger.add(
        "regional_code.within_protein.auroc",
        "2C",
        "regional_code",
        "full-model mean within-protein AUROC",
        regional["within_protein"]["mean_auroc"],
        interval=regional["within_protein"]["ci"],
        unit="AUROC",
        inferential_unit="protein",
        source=source_label(regional_path, root),
    )
    ledger.add(
        "regional_code.composition_within_protein.auroc",
        "2C",
        "regional_code",
        "composition-only mean within-protein AUROC",
        regional["composition_within_protein"]["mean_auroc"],
        interval=regional["composition_within_protein"]["ci"],
        unit="AUROC",
        inferential_unit="protein",
        source=source_label(regional_path, root),
    )

    yin_path = root / "analysis/yin_yang/summary.json"
    yin = read_json(yin_path)
    for name, item in yin["local_opportunity"].items():
        ledger.add(
            f"yin_yang.local.{name}",
            "3D;S2C",
            "yin_yang",
            "same-site observed to expected ratio",
            item["ratio"],
            interval=item["ci"],
            unit="ratio",
            inferential_unit="protein",
            source=source_label(yin_path, root),
        )
    for arm, nulls in yin["nested_nulls"].items():
        for null, item in nulls.items():
            ledger.add(
                f"yin_yang.nested.{arm}.{null}",
                "3E",
                "yin_yang",
                f"{arm} same-site ratio under {null}",
                item["ratio"],
                interval=item["ci"],
                unit="ratio",
                inferential_unit="protein",
                source=source_label(yin_path, root),
            )
    for key in ("cluster_n90", "control_n90", "delta_n90"):
        ledger.add(
            f"kinome.{key}",
            "3C;S2A",
            "kinome",
            key,
            yin["kinome"][key],
            interval=yin["kinome"].get("delta_n90_ci") if key == "delta_n90" else None,
            unit="kinases above 90th percentile",
            inferential_unit="protein",
            source=source_label(yin_path, root),
        )

    evolution_path = root / "analysis/evolution/summary.json"
    evolution = read_json(evolution_path)
    for item in evolution["site_conservation"]:
        ledger.add(
            f"evolution.{item['species']}.gap",
            "4B",
            "evolution",
            "regional minus exact conservation",
            item["gap"],
            interval=item["gap_ci"],
            unit="fraction",
            inferential_unit="protein",
            source=source_label(evolution_path, root),
        )
    for item in evolution["transfer"]:
        ledger.add(
            f"transfer.{item['direction']}.{item['species']}",
            "4C",
            "evolution",
            "alignment-free transfer AUROC",
            item["auroc"],
            interval=item["ci"],
            unit="AUROC",
            inferential_unit="protein-grouped tile",
            source=source_label(evolution_path, root),
        )
    for metric, item in evolution["composition_position_nulls"].items():
        ledger.add(
            f"evolution.null_excess.{metric}",
            "4E;S3",
            "evolution",
            f"{metric} conservation excess over within-protein null",
            item["excess"],
            interval=item["excess_ci"],
            unit="correlation-scale difference",
            inferential_unit="protein",
            source=source_label(evolution_path, root),
        )

    scanner_path = root / "analysis/scanner/summary.json"
    scanner = read_json(scanner_path)
    for model in ("baseline", "enhanced", "selected"):
        for metric in ("auroc", "average_precision"):
            ledger.add(
                f"scanner.{model}.{metric}",
                "5C;5D",
                "scanner",
                f"{model} {metric}",
                scanner[model][metric],
                unit=metric,
                inferential_unit="outer protein-grouped fold",
                source=source_label(scanner_path, root),
            )

    ogt_path = root / "analysis/ogt_conservation/summary.json"
    ogt = read_json(ogt_path)
    for item in ogt["human_mouse_categories"]:
        key, value = item["category"], item["identity"]
        ledger.add(
            f"ogt.identity.{key}",
            "4A",
            "ogt_conservation",
            f"human to mouse identity: {key}",
            value,
            unit="fraction",
            inferential_unit="OGT residue",
            source=source_label(ogt_path, root),
        )

    fg_path = root / "analysis/fg_nup_recovery/summary.json"
    fg = read_json(fg_path)
    for key, unit in (
        ("recovery_fraction", "fraction"),
        ("background_region_fraction", "fraction"),
        ("hypergeometric_p", "p-value"),
        ("matched_p", "p-value"),
        ("adjusted_odds_ratio", "odds ratio"),
        ("adjusted_p", "p-value"),
    ):
        ledger.add(
            f"fg_nup.{key}",
            "5A;5B",
            "fg_nup_recovery",
            key,
            fg[key],
            unit=unit,
            inferential_unit="protein",
            source=source_label(fg_path, root),
        )

    deletion_path = root / "analysis/adversarial_deletion/summary.csv"
    deletion = pd.read_csv(deletion_path).iloc[0]
    deletion_metrics = {
        "full_fold": "full_group_fold",
        "fold_after_50_percent": "fold_after_top_50pct",
        "first_below_null_percent": "first_below_null_pct",
    }
    for key, column in deletion_metrics.items():
        ledger.add(
            f"deletion.{key}",
            "1F",
            "adversarial_deletion",
            key,
            float(deletion[column]),
            unit="fold" if "fold" in key else "percent removed",
            inferential_unit="protein",
            source=source_label(deletion_path, root),
        )

    return pd.DataFrame(ledger.rows, columns=FIELDS).sort_values("claim_id")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    config = load_config(args.config)
    output = args.output.resolve()
    ensure_parent(output)
    assemble(config.results_root).to_csv(output, index=False)


if __name__ == "__main__":
    main()
