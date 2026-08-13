"""Consensus-region definition shared by every downstream branch."""


rule consensus_regions:
    input:
        manifest=rules.input_manifest.output,
    output:
        regions=f"{RESULTS}/analysis/consensus_regions/consensus_regions.csv",
        sites=f"{RESULTS}/analysis/consensus_regions/consensus_region_sites.csv",
        summary=f"{RESULTS}/analysis/consensus_regions/summary.json",
    log:
        f"{RESULTS}/logs/regions.log",
    resources:
        mem_mb=6000,
        runtime=60,
        disk_mb=4000,
    shell:
        "uv run --frozen python -m regional_code_paper.analysis.consensus_regions "
        "--config config/config.yaml "
        "--regions {output.regions:q} --sites {output.sites:q} --summary {output.summary:q} "
        "> {log:q} 2>&1"
