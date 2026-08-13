"""Independent PTM maps followed by protein-level statistical reduction."""


rule self_clustering_ptm:
    input:
        manifest=rules.input_manifest.output,
    output:
        data=f"{RESULTS}/work/self_clustering/{{ptm}}.csv",
        receipt=f"{RESULTS}/work/self_clustering/{{ptm}}.receipt.json",
    wildcard_constraints:
        ptm="|".join(PTMS),
    log:
        f"{RESULTS}/logs/self_clustering/{{ptm}}.log",
    benchmark:
        f"{RESULTS}/benchmarks/self_clustering/{{ptm}}.tsv",
    resources:
        mem_mb=6000,
        runtime=90,
        disk_mb=4000,
    shell:
        "uv run --frozen python -m regional_code_paper.execution.self_clustering ptm "
        "--config config/config.yaml --ptm {wildcards.ptm:q} --output {output.data:q} "
        "--receipt {output.receipt:q} > {log:q} 2>&1"


rule self_clustering:
    input:
        shards=expand(f"{RESULTS}/work/self_clustering/{{ptm}}.csv", ptm=PTMS),
        receipts=expand(f"{RESULTS}/work/self_clustering/{{ptm}}.receipt.json", ptm=PTMS),
    output:
        proteins=f"{RESULTS}/analysis/self_clustering/per_protein.csv",
        summary=f"{RESULTS}/analysis/self_clustering/summary.csv",
        contrasts=f"{RESULTS}/analysis/self_clustering/contrasts.csv",
        ratios=f"{RESULTS}/analysis/self_clustering/fold_ratios.csv",
        provenance=f"{RESULTS}/analysis/self_clustering/provenance.json",
    log:
        f"{RESULTS}/logs/self_clustering/reduce.log",
    resources:
        mem_mb=8000,
        runtime=120,
        disk_mb=6000,
    shell:
        "uv run --frozen python -m regional_code_paper.execution.self_clustering reduce "
        "--config config/config.yaml --shards {input.shards:q} --receipts {input.receipts:q} "
        "--output-dir {RESULTS}/analysis/self_clustering > {log:q} 2>&1"


rule clustering_breadth:
    input:
        proteins=rules.self_clustering.output.proteins,
    output:
        per_scale=f"{RESULTS}/analysis/clustering_breadth/per_scale.csv",
        dominance=f"{RESULTS}/analysis/clustering_breadth/matched_dominance.csv",
        multiscale=f"{RESULTS}/analysis/clustering_breadth/multiscale.csv",
        summary=f"{RESULTS}/analysis/clustering_breadth/summary.json",
    log:
        f"{RESULTS}/logs/clustering_breadth.log",
    resources:
        mem_mb=2000,
        runtime=15,
        disk_mb=2000,
    shell:
        "uv run --frozen python -m regional_code_paper.analysis.clustering_breadth "
        "--per-protein {input.proteins:q} "
        "--output-dir {RESULTS}/analysis/clustering_breadth > {log:q} 2>&1"


rule adversarial_deletion:
    input:
        proteins=rules.self_clustering.output.proteins,
    output:
        curve=f"{RESULTS}/analysis/adversarial_deletion/adversarial_curve.csv",
        envelope=f"{RESULTS}/analysis/adversarial_deletion/random_envelope.csv",
        summary=f"{RESULTS}/analysis/adversarial_deletion/summary.csv",
    params:
        sql="src/regional_code_paper/sql/adversarial_deletion.sql",
    log:
        f"{RESULTS}/logs/adversarial_deletion.log",
    resources:
        mem_mb=4000,
        runtime=30,
        disk_mb=4000,
    shell:
        "mkdir -p {RESULTS}/analysis/adversarial_deletion && "
        "PAPER_CLUSTERING_INPUT={input.proteins:q} "
        "PAPER_OUTPUT_ROOT={RESULTS}/analysis/adversarial_deletion "
        "uv run --frozen duckdb < {params.sql:q} > {log:q} 2>&1"
