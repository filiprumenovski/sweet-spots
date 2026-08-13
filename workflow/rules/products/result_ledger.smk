"""Single review ledger assembled only after every paper analysis succeeds."""


rule result_ledger_csv:
    input:
        regions=rules.consensus_regions.output.summary,
        clustering=rules.self_clustering.output.ratios,
        clustering_breadth=rules.clustering_breadth.output.summary,
        detection_aware_clustering=rules.detection_aware_clustering.output.summary,
        evidence_restricted_clustering=rules.evidence_restricted_clustering.output.summary,
        regional_code=rules.regional_code.output.summary,
        yin_yang=rules.yin_yang.output.summary,
        evolution=rules.evolution.output.summary,
        scanner=rules.scanner.output.summary,
        ogt=rules.ogt_conservation.output.summary,
        fg_nup=rules.fg_nup_recovery.output.summary,
        deletion=rules.adversarial_deletion.output.summary,
    output:
        f"{RESULTS}/tables/manuscript_result_ledger.csv",
    log:
        f"{RESULTS}/logs/metrics.log",
    resources:
        mem_mb=2000,
        runtime=15,
        disk_mb=2000,
    shell:
        "uv run --frozen python -m regional_code_paper.reporting.result_ledger "
        "--config config/config.yaml "
        "--output {output:q} > {log:q} 2>&1"


rule result_ledger_parquet:
    input:
        rules.result_ledger_csv.output,
    output:
        f"{RESULTS}/tables/manuscript_result_ledger.parquet",
    log:
        f"{RESULTS}/logs/metrics_parquet.log",
    resources:
        mem_mb=2000,
        runtime=15,
        disk_mb=2000,
    shell:
        "mkdir -p {RESULTS}/tables && "
        "uv run --frozen duckdb -c \"COPY (SELECT * FROM read_csv_auto('{input}')) "
        "TO '{output}' (FORMAT PARQUET, COMPRESSION ZSTD);\" > {log:q} 2>&1"
