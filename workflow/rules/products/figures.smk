"""Figure rules. Every panel reads only products declared by the analysis DAG."""

FIGURE_INPUTS = [
    rules.consensus_regions.output.regions,
    rules.consensus_regions.output.sites,
    rules.consensus_regions.output.summary,
    rules.self_clustering.output.proteins,
    rules.self_clustering.output.summary,
    rules.self_clustering.output.ratios,
    rules.clustering_breadth.output.per_scale,
    rules.clustering_breadth.output.dominance,
    rules.clustering_breadth.output.multiscale,
    rules.detection_aware_clustering.output.scales,
    rules.detection_aware_clustering.output.summary,
    rules.evidence_restricted_clustering.output.scales,
    rules.evidence_restricted_clustering.output.summary,
    rules.adversarial_deletion.output.curve,
    rules.adversarial_deletion.output.envelope,
    rules.regional_code.output.tiles,
    rules.regional_code.output.proteins,
    rules.regional_code.output.acceptors,
    rules.regional_code.output.proline,
    rules.regional_code.output.distance,
    rules.regional_code.output.summary,
    rules.yin_yang.output.local,
    rules.yin_yang.output.nested,
    rules.yin_yang.output.kinase_sites,
    rules.yin_yang.output.kinase_summary,
    rules.yin_yang.output.summary,
    rules.evolution.output.site_conservation,
    rules.evolution.output.transfer,
    rules.evolution.output.composition_nulls,
    rules.evolution.output.substitutions,
    rules.scanner.output.predictions,
    rules.scanner.output.folds,
    rules.scanner.output.operating_points,
    rules.ogt_conservation.output.categories,
    rules.fg_nup_recovery.output.proteins,
    rules.fg_nup_recovery.output.matched,
    rules.fg_nup_recovery.output.summary,
]


rule figures:
    input:
        FIGURE_INPUTS,
    output:
        pdf=expand(f"{RESULTS}/figures/{{figure}}.pdf", figure=FIGURES),
        png=expand(f"{RESULTS}/figures/{{figure}}.png", figure=FIGURES),
        manifest=f"{RESULTS}/figures/figure_manifest.json",
    log:
        f"{RESULTS}/logs/figures.log",
    resources:
        mem_mb=8000,
        runtime=120,
        disk_mb=10000,
    shell:
        "uv run --frozen python -m regional_code_paper.reporting.figures "
        "--config config/config.yaml "
        "--output-dir {RESULTS}/figures > {log:q} 2>&1"
