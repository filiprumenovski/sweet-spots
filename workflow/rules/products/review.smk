"""An offline, source-linked report of every claim and figure."""

from regional_code_paper.reporting.review import source_files


rule reviewer_report:
    input:
        audit=rules.audit_report.output,
        code=[str(path) for path in source_files(Path("."))],
        sources=rules.audit_report.input,
    output:
        html=f"{RESULTS}/review/index.html",
        dot=f"{RESULTS}/review/workflow.dot",
        svg=f"{RESULTS}/review/workflow.svg",
        hashes=f"{RESULTS}/review/source_hashes.json",
    log:
        f"{RESULTS}/logs/review.log",
    resources:
        mem_mb=4000,
        runtime=15,
        disk_mb=2000,
    shell:
        "uv run --frozen python -m regional_code_paper.reporting.review "
        "--config config/config.yaml > {log:q} 2>&1"
