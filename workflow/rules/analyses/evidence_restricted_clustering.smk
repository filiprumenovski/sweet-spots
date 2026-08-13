"""Independent-publication and independent-peptide clustering restrictions."""


rule evidence_restricted_clustering:
    input:
        manifest=rules.input_manifest.output,
    output:
        proteins=f"{RESULTS}/analysis/evidence_restricted_clustering/per_protein.csv",
        scales=f"{RESULTS}/analysis/evidence_restricted_clustering/per_scale.csv",
        summary=f"{RESULTS}/analysis/evidence_restricted_clustering/summary.json",
    log:
        f"{RESULTS}/logs/evidence_restricted_clustering.log",
    resources:
        mem_mb=6000,
        runtime=60,
        disk_mb=4000,
    shell:
        "uv run --frozen python -m "
        "regional_code_paper.analysis.evidence_restricted_clustering "
        "--config config/config.yaml "
        "--output-dir {RESULTS}/analysis/evidence_restricted_clustering > {log:q} 2>&1"
