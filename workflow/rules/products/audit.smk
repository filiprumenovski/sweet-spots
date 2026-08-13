"""Fail-closed checks of the numerical results and frozen figure set."""


rule audit:
    input:
        metrics=rules.result_ledger_parquet.output,
        figure_manifest=rules.figures.output.manifest,
        figures=rules.figures.output.pdf + rules.figures.output.png,
    output:
        f"{RESULTS}/audit/reproduction_report.json",
    log:
        f"{RESULTS}/logs/audit.log",
    resources:
        mem_mb=4000,
        runtime=30,
        disk_mb=4000,
    shell:
        "uv run --frozen python -m regional_code_paper.reporting.audit "
        "--config config/config.yaml "
        "--metrics {input.metrics:q} --figure-manifest {input.figure_manifest:q} "
        "--output {output:q} > {log:q} 2>&1"
