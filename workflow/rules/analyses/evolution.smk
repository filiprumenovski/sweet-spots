"""Evolutionary conservation and alignment-free transfer analyses."""


rule evolution:
    input:
        manifest=rules.input_manifest.output,
        regions=rules.consensus_regions.output.regions,
        sites=rules.consensus_regions.output.sites,
    output:
        site_conservation=f"{RESULTS}/analysis/evolution/site_conservation.csv",
        transfer=f"{RESULTS}/analysis/evolution/transfer.csv",
        composition_nulls=f"{RESULTS}/analysis/evolution/composition_position_nulls.csv",
        substitutions=f"{RESULTS}/analysis/evolution/substitutions.csv",
        summary=f"{RESULTS}/analysis/evolution/summary.json",
    threads: 2
    log:
        f"{RESULTS}/logs/evolution.log",
    benchmark:
        f"{RESULTS}/benchmarks/evolution.tsv",
    resources:
        mem_mb=16000,
        runtime=360,
        disk_mb=12000,
    shell:
        "OMP_NUM_THREADS={threads} OPENBLAS_NUM_THREADS={threads} "
        "uv run --frozen python -m regional_code_paper.analysis.evolution "
        "--config config/config.yaml "
        "--regions {input.regions:q} --sites {input.sites:q} "
        "--output-dir {RESULTS}/analysis/evolution > {log:q} 2>&1"
