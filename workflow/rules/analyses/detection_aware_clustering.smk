"""Peptide-exposure-conditioned O-GlcNAc self-clustering sensitivity."""

DETECTION_PREP = f"{RESULTS}/work/detection_aware_clustering/prepared"


rule detection_aware_clustering_prepare:
    input:
        manifest=rules.input_manifest.output,
    output:
        proteins=f"{DETECTION_PREP}/proteins.csv",
        strata=f"{DETECTION_PREP}/strata.csv",
        summary=f"{DETECTION_PREP}/summary.json",
    log:
        f"{RESULTS}/logs/detection_aware_clustering/prepare.log",
    benchmark:
        f"{RESULTS}/benchmarks/detection_aware_clustering/prepare.tsv",
    resources:
        mem_mb=6000,
        runtime=60,
        disk_mb=4000,
    shell:
        "uv run --frozen python -m regional_code_paper.analysis.detection_aware_clustering "
        "--config config/config.yaml --output-dir {DETECTION_PREP} > {log:q} 2>&1"


rule detection_aware_clustering_shard:
    input:
        proteins=rules.detection_aware_clustering_prepare.output.proteins,
        strata=rules.detection_aware_clustering_prepare.output.strata,
    output:
        data=f"{RESULTS}/work/detection_aware_clustering/shard_{{shard}}.csv",
        receipt=(
            f"{RESULTS}/work/detection_aware_clustering/shard_{{shard}}.receipt.json"
        ),
    wildcard_constraints:
        shard="\\d+",
    log:
        f"{RESULTS}/logs/detection_aware_clustering/shard_{{shard}}.log",
    benchmark:
        f"{RESULTS}/benchmarks/detection_aware_clustering/shard_{{shard}}.tsv",
    resources:
        mem_mb=config["resources"]["simulation_shard"]["mem_mb"],
        runtime=config["resources"]["simulation_shard"]["runtime"],
        disk_mb=config["resources"]["simulation_shard"]["disk_mb"],
    shell:
        "uv run --frozen python -m regional_code_paper.execution.detection_aware_clustering "
        "simulate --proteins {input.proteins:q} --strata {input.strata:q} "
        "--shard {wildcards.shard} --shards {config[parallel][detection_null_shards]} "
        "--draws {config[randomness][detection_null_draws]} "
        "--seed {config[randomness][manuscript_base_seed]} --output {output.data:q} "
        "--receipt {output.receipt:q} > {log:q} 2>&1"


rule detection_aware_clustering:
    input:
        prepared_summary=rules.detection_aware_clustering_prepare.output.summary,
        shards=expand(
            f"{RESULTS}/work/detection_aware_clustering/shard_{{shard}}.csv",
            shard=DETECTION_SHARDS,
        ),
        receipts=expand(
            f"{RESULTS}/work/detection_aware_clustering/shard_{{shard}}.receipt.json",
            shard=DETECTION_SHARDS,
        ),
    output:
        proteins=f"{RESULTS}/analysis/detection_aware_clustering/per_protein.csv",
        scales=f"{RESULTS}/analysis/detection_aware_clustering/per_scale.csv",
        summary=f"{RESULTS}/analysis/detection_aware_clustering/summary.json",
    log:
        f"{RESULTS}/logs/detection_aware_clustering/reduce.log",
    resources:
        mem_mb=4000,
        runtime=30,
        disk_mb=4000,
    shell:
        "uv run --frozen python -m regional_code_paper.execution.detection_aware_clustering "
        "reduce --config config/config.yaml --prepared-summary {input.prepared_summary:q} "
        "--shards {input.shards:q} --receipts {input.receipts:q} "
        "--output-dir {RESULTS}/analysis/detection_aware_clustering > {log:q} 2>&1"
