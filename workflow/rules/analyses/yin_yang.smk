"""Opportunity nulls plus horizontally scalable kinase-library scoring."""

YIN_PREP = f"{RESULTS}/work/yin_yang/prepared"


rule yin_yang_prepare:
    input:
        manifest=rules.input_manifest.output,
        regions=rules.consensus_regions.output.regions,
        sites=rules.consensus_regions.output.sites,
    output:
        local=f"{YIN_PREP}/local_opportunity.csv",
        nested=f"{YIN_PREP}/nested_nulls.csv",
        pairs=f"{YIN_PREP}/kinase_pairs.csv",
        metadata=f"{YIN_PREP}/prepared.json",
    log:
        f"{RESULTS}/logs/yin_yang/prepare.log",
    resources:
        mem_mb=10000,
        runtime=120,
        disk_mb=8000,
    shell:
        "uv run --frozen python -m regional_code_paper.execution.yin_yang prepare "
        "--config config/config.yaml --regions {input.regions:q} --sites {input.sites:q} "
        "--output-dir {YIN_PREP} > {log:q} 2>&1"


rule kinase_scoring_shard:
    input:
        pairs=rules.yin_yang_prepare.output.pairs,
    output:
        sites=f"{RESULTS}/work/yin_yang/kinase/shard_{{shard}}.sites.csv",
        sums=f"{RESULTS}/work/yin_yang/kinase/shard_{{shard}}.scores.npz",
        receipt=f"{RESULTS}/work/yin_yang/kinase/shard_{{shard}}.receipt.json",
    wildcard_constraints:
        shard="\\d+",
    log:
        f"{RESULTS}/logs/yin_yang/kinase_{{shard}}.log",
    benchmark:
        f"{RESULTS}/benchmarks/yin_yang/kinase_{{shard}}.tsv",
    resources:
        mem_mb=config["resources"]["kinase_shard"]["mem_mb"],
        runtime=config["resources"]["kinase_shard"]["runtime"],
        disk_mb=config["resources"]["kinase_shard"]["disk_mb"],
    shell:
        "uv run --frozen python -m regional_code_paper.execution.yin_yang kinase-shard "
        "--pairs {input.pairs:q} --shard {wildcards.shard} "
        "--shards {config[parallel][kinase_scoring_shards]} "
        "--sites-output {output.sites:q} --sums-output {output.sums:q} "
        "--receipt {output.receipt:q} > {log:q} 2>&1"


rule yin_yang:
    input:
        prepared=rules.yin_yang_prepare.output,
        site_shards=expand(f"{RESULTS}/work/yin_yang/kinase/shard_{{shard}}.sites.csv", shard=KINASE_SHARDS),
        sum_shards=expand(f"{RESULTS}/work/yin_yang/kinase/shard_{{shard}}.scores.npz", shard=KINASE_SHARDS),
        receipts=expand(f"{RESULTS}/work/yin_yang/kinase/shard_{{shard}}.receipt.json", shard=KINASE_SHARDS),
    output:
        local=f"{RESULTS}/analysis/yin_yang/local_opportunity.csv",
        nested=f"{RESULTS}/analysis/yin_yang/nested_nulls.csv",
        kinase_sites=f"{RESULTS}/analysis/yin_yang/kinase_site_scores.csv",
        kinase_summary=f"{RESULTS}/analysis/yin_yang/kinase_summary.csv",
        summary=f"{RESULTS}/analysis/yin_yang/summary.json",
    log:
        f"{RESULTS}/logs/yin_yang/reduce.log",
    resources:
        mem_mb=8000,
        runtime=90,
        disk_mb=6000,
    shell:
        "uv run --frozen python -m regional_code_paper.execution.yin_yang reduce "
        "--config config/config.yaml --prepared-dir {YIN_PREP} "
        "--site-shards {input.site_shards:q} --sum-shards {input.sum_shards:q} "
        "--receipts {input.receipts:q} --output-dir {RESULTS}/analysis/yin_yang "
        "> {log:q} 2>&1"
