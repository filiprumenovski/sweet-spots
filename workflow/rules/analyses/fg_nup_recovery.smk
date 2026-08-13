"""Recovery of FG-nucleoporins as an external biological benchmark."""


rule fg_nup_recovery:
    input:
        manifest=rules.input_manifest.output,
        regions=rules.consensus_regions.output.regions,
    output:
        proteins=f"{RESULTS}/analysis/fg_nup_recovery/proteins.csv",
        matched=f"{RESULTS}/analysis/fg_nup_recovery/matched_universe.csv",
        summary=f"{RESULTS}/analysis/fg_nup_recovery/summary.json",
    log:
        f"{RESULTS}/logs/fg_nup.log",
    resources:
        mem_mb=10000,
        runtime=120,
        disk_mb=8000,
    shell:
        "uv run --frozen python -m regional_code_paper.analysis.fg_nup_recovery "
        "--config config/config.yaml "
        "--regions {input.regions:q} --output-dir {RESULTS}/analysis/fg_nup_recovery "
        "> {log:q} 2>&1"
