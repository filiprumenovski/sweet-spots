"""Conservation of OGT substrate-reading and catalytic surfaces."""


rule ogt_conservation:
    input:
        manifest=rules.input_manifest.output,
    output:
        residues=f"{RESULTS}/analysis/ogt_conservation/residue_conservation.csv",
        categories=f"{RESULTS}/analysis/ogt_conservation/category_summary.csv",
        summary=f"{RESULTS}/analysis/ogt_conservation/summary.json",
    log:
        f"{RESULTS}/logs/ogt.log",
    resources:
        mem_mb=6000,
        runtime=60,
        disk_mb=4000,
    shell:
        "uv run --frozen python -m regional_code_paper.analysis.ogt_conservation "
        "--config config/config.yaml "
        "--output-dir {RESULTS}/analysis/ogt_conservation > {log:q} 2>&1"
