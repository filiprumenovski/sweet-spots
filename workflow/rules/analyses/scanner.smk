"""Shared feature cache and five independently restartable nested-CV folds."""

SCANNER_WORK = f"{RESULTS}/work/scanner"


rule scanner_prepare:
    input:
        manifest=rules.input_manifest.output,
        regions=rules.consensus_regions.output.regions,
    output:
        cache=f"{SCANNER_WORK}/features.npz",
        metadata=f"{SCANNER_WORK}/features.json",
    log:
        f"{RESULTS}/logs/scanner/prepare.log",
    benchmark:
        f"{RESULTS}/benchmarks/scanner/prepare.tsv",
    resources:
        mem_mb=config["resources"]["feature_extraction"]["mem_mb"],
        runtime=config["resources"]["feature_extraction"]["runtime"],
        disk_mb=config["resources"]["feature_extraction"]["disk_mb"],
    shell:
        "uv run --frozen python -m regional_code_paper.execution.scanner prepare "
        "--config config/config.yaml --regions {input.regions:q} --cache {output.cache:q} "
        "--metadata {output.metadata:q} > {log:q} 2>&1"


rule scanner_fold:
    input:
        cache=rules.scanner_prepare.output.cache,
        metadata=rules.scanner_prepare.output.metadata,
    output:
        data=f"{SCANNER_WORK}/folds/fold_{{fold}}.csv",
        receipt=f"{SCANNER_WORK}/folds/fold_{{fold}}.receipt.json",
    wildcard_constraints:
        fold="[0-4]",
    # One numerical thread per fit; the five outer folds are the parallel unit.
    threads: 1
    log:
        f"{RESULTS}/logs/scanner/fold_{{fold}}.log",
    benchmark:
        f"{RESULTS}/benchmarks/scanner/fold_{{fold}}.tsv",
    resources:
        mem_mb=config["resources"]["model_fold"]["mem_mb"],
        runtime=config["resources"]["model_fold"]["runtime"],
        disk_mb=config["resources"]["model_fold"]["disk_mb"],
    shell:
        "OMP_NUM_THREADS={threads} OPENBLAS_NUM_THREADS={threads} "
        "uv run --frozen python -m regional_code_paper.execution.scanner fold "
        "--cache {input.cache:q} --metadata {input.metadata:q} --fold {wildcards.fold} "
        "--output {output.data:q} --receipt {output.receipt:q} > {log:q} 2>&1"


rule scanner:
    input:
        cache=rules.scanner_prepare.output.cache,
        metadata=rules.scanner_prepare.output.metadata,
        shards=expand(f"{SCANNER_WORK}/folds/fold_{{fold}}.csv", fold=OUTER_FOLDS),
        receipts=expand(f"{SCANNER_WORK}/folds/fold_{{fold}}.receipt.json", fold=OUTER_FOLDS),
    output:
        predictions=f"{RESULTS}/analysis/scanner/outer_predictions.csv",
        folds=f"{RESULTS}/analysis/scanner/outer_folds.csv",
        operating_points=f"{RESULTS}/analysis/scanner/operating_points.csv",
        summary=f"{RESULTS}/analysis/scanner/summary.json",
    log:
        f"{RESULTS}/logs/scanner/reduce.log",
    resources:
        mem_mb=10000,
        runtime=60,
        disk_mb=8000,
    shell:
        "uv run --frozen python -m regional_code_paper.execution.scanner reduce "
        "--cache {input.cache:q} --metadata {input.metadata:q} "
        "--shards {input.shards:q} --receipts {input.receipts:q} "
        "--output-dir {RESULTS}/analysis/scanner > {log:q} 2>&1"
