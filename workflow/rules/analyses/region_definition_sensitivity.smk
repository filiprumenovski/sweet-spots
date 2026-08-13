"""Post hoc robustness grid for strict-core membership and final grouping."""


REGION_SENSITIVITY = f"{RESULTS}/analysis/region_definition_sensitivity"
REGION_SENSITIVITY_WORK = f"{RESULTS}/work/region_definition_sensitivity"


rule region_definition_catalogue_sensitivity:
    input:
        manifest=rules.input_manifest.output,
        primary=rules.consensus_regions.output.regions,
    output:
        f"{REGION_SENSITIVITY}/catalogue_grid.csv",
    log:
        f"{RESULTS}/logs/region_definition_sensitivity/catalogue.log",
    resources:
        mem_mb=config["resources"]["default"]["mem_mb"],
        runtime=config["resources"]["default"]["runtime"],
        disk_mb=config["resources"]["default"]["disk_mb"],
    shell:
        "uv run --frozen python -m regional_code_paper.analysis.region_definition_sensitivity "
        "catalogue --config config/config.yaml --output {output:q} > {log:q} 2>&1"


rule region_definition_model_sensitivity:
    input:
        manifest=rules.input_manifest.output,
        primary=rules.consensus_regions.output.regions,
    output:
        f"{REGION_SENSITIVITY_WORK}/core{{core_gap}}_min{{minimum_sites}}.json",
    wildcard_constraints:
        core_gap="|".join(REGION_SENSITIVITY_CORE_GAPS),
        minimum_sites="|".join(REGION_SENSITIVITY_MINIMUM_SITES),
    threads: 1
    log:
        f"{RESULTS}/logs/region_definition_sensitivity/core{{core_gap}}_min{{minimum_sites}}.log",
    resources:
        mem_mb=config["resources"]["model_fold"]["mem_mb"],
        runtime=config["resources"]["model_fold"]["runtime"],
        disk_mb=config["resources"]["model_fold"]["disk_mb"],
    shell:
        "OMP_NUM_THREADS={threads} OPENBLAS_NUM_THREADS={threads} "
        "uv run --frozen python -m regional_code_paper.analysis.region_definition_sensitivity "
        "model --config config/config.yaml --core-gap {wildcards.core_gap} "
        "--minimum-sites {wildcards.minimum_sites} --output {output:q} > {log:q} 2>&1"


rule region_definition_model_sensitivity_reduce:
    input:
        expand(
            f"{REGION_SENSITIVITY_WORK}/core{{core_gap}}_min{{minimum_sites}}.json",
            core_gap=REGION_SENSITIVITY_CORE_GAPS,
            minimum_sites=REGION_SENSITIVITY_MINIMUM_SITES,
        ),
    output:
        f"{REGION_SENSITIVITY}/model_grid.csv",
    log:
        f"{RESULTS}/logs/region_definition_sensitivity/reduce.log",
    resources:
        mem_mb=config["resources"]["default"]["mem_mb"],
        runtime=config["resources"]["default"]["runtime"],
        disk_mb=config["resources"]["default"]["disk_mb"],
    shell:
        "uv run --frozen python -m regional_code_paper.analysis.region_definition_sensitivity "
        "reduce --inputs {input:q} --output {output:q} > {log:q} 2>&1"
