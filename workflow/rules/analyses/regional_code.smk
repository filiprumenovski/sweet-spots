"""Feature preparation, protein-grouped CV maps, and validated reduction."""

REGIONAL_PREP = f"{RESULTS}/work/regional_code/prepared"


rule regional_code_prepare:
    input:
        manifest=rules.input_manifest.output,
        regions=rules.consensus_regions.output.regions,
        sites=rules.consensus_regions.output.sites,
    output:
        tiles=f"{REGIONAL_PREP}/tile_features.csv",
        cache=f"{REGIONAL_PREP}/tile_features.npz",
        acceptors=f"{REGIONAL_PREP}/acceptor_model.csv",
        proline=f"{REGIONAL_PREP}/proline_profiles.csv",
        distance=f"{REGIONAL_PREP}/distance_profile.csv",
        matches=f"{REGIONAL_PREP}/distance_matches.csv",
        metadata=f"{REGIONAL_PREP}/prepared.json",
    log:
        f"{RESULTS}/logs/regional_code/prepare.log",
    benchmark:
        f"{RESULTS}/benchmarks/regional_code/prepare.tsv",
    resources:
        mem_mb=config["resources"]["feature_extraction"]["mem_mb"],
        runtime=config["resources"]["feature_extraction"]["runtime"],
        disk_mb=config["resources"]["feature_extraction"]["disk_mb"],
    shell:
        "uv run --frozen python -m regional_code_paper.execution.regional_code prepare "
        "--config config/config.yaml --regions {input.regions:q} --sites {input.sites:q} "
        "--output-dir {REGIONAL_PREP} > {log:q} 2>&1"


rule regional_code_fold:
    input:
        cache=rules.regional_code_prepare.output.cache,
    output:
        data=f"{RESULTS}/work/regional_code/folds/fold_{{fold}}.csv",
        receipt=f"{RESULTS}/work/regional_code/folds/fold_{{fold}}.receipt.json",
    wildcard_constraints:
        fold="[0-4]",
    # HistGradientBoosting is stable across fold jobs but can drift with
    # OpenMP reduction order inside a fit. Parallelize folds, not one fit.
    threads: 1
    log:
        f"{RESULTS}/logs/regional_code/fold_{{fold}}.log",
    benchmark:
        f"{RESULTS}/benchmarks/regional_code/fold_{{fold}}.tsv",
    resources:
        mem_mb=config["resources"]["model_fold"]["mem_mb"],
        runtime=config["resources"]["model_fold"]["runtime"],
        disk_mb=config["resources"]["model_fold"]["disk_mb"],
    shell:
        "OMP_NUM_THREADS={threads} OPENBLAS_NUM_THREADS={threads} "
        "uv run --frozen python -m regional_code_paper.execution.regional_code fold "
        "--cache {input.cache:q} --fold {wildcards.fold} "
        "--seed {config[randomness][analysis_base_seed]} --output {output.data:q} "
        "--receipt {output.receipt:q} > {log:q} 2>&1"


rule regional_code:
    input:
        prepared=rules.regional_code_prepare.output,
        folds=expand(f"{RESULTS}/work/regional_code/folds/fold_{{fold}}.csv", fold=OUTER_FOLDS),
        fold_receipts=expand(f"{RESULTS}/work/regional_code/folds/fold_{{fold}}.receipt.json", fold=OUTER_FOLDS),
    output:
        tiles=f"{RESULTS}/analysis/regional_code/held_out_tiles.csv",
        proteins=f"{RESULTS}/analysis/regional_code/per_protein.csv",
        acceptors=f"{RESULTS}/analysis/regional_code/acceptor_model.csv",
        proline=f"{RESULTS}/analysis/regional_code/proline_profiles.csv",
        distance=f"{RESULTS}/analysis/regional_code/distance_profile.csv",
        distance_matches=f"{RESULTS}/analysis/regional_code/distance_matches.csv",
        summary=f"{RESULTS}/analysis/regional_code/summary.json",
    log:
        f"{RESULTS}/logs/regional_code/reduce.log",
    resources:
        mem_mb=10000,
        runtime=90,
        disk_mb=8000,
    shell:
        "uv run --frozen python -m regional_code_paper.execution.regional_code reduce "
        "--config config/config.yaml --prepared-dir {REGIONAL_PREP} "
        "--folds {input.folds:q} --fold-receipts {input.fold_receipts:q} "
        "--output-dir {RESULTS}/analysis/regional_code > {log:q} 2>&1"
