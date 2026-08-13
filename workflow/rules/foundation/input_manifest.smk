"""Immutable inputs and cryptographic provenance."""


rule input_manifest:
    input:
        fasta_human=f"{SOURCE_ROOT}/data/interim/fasta_human.parquet",
        fasta_all=f"{SOURCE_ROOT}/data/interim/fasta_all.parquet",
        disorder=f"{SOURCE_ROOT}/data/interim/iupred_residue_scores.parquet",
        atlas=f"{SOURCE_ROOT}/analysis/revalidation/data/atlas_unambiguous.csv",
        ptm=f"{SOURCE_ROOT}/data/processed/landscape/ptm_unified.parquet",
        reference_regions=(
            f"{SOURCE_ROOT}/data/processed/regions/oglcnac_consensus_regions.parquet"
        ),
        human_ortholog_alignment_directory=f"{SOURCE_ROOT}/data/interim/msa",
        external_sequence_fastas=EXTERNAL_SEQUENCE_FASTAS,
        ogt_structure_directory=(
            f"{SOURCE_ROOT}/analysis/revalidation/data/ogt_structures"
        ),
        ogt_human_sequence=(
            f"{SOURCE_ROOT}/analysis/revalidation/data/ogt_orthologs/O15294.fasta"
        ),
        ogt_ortholog_alignment=(
            f"{SOURCE_ROOT}/analysis/revalidation/data/ogt_orthologs/orthologs_aln.fasta"
        ),
    output:
        f"{RESULTS}/provenance/input_manifest.json",
    log:
        f"{RESULTS}/logs/input_manifest.log",
    resources:
        mem_mb=2000,
        runtime=15,
        disk_mb=2000,
    shell:
        "uv run --frozen python -m regional_code_paper.reporting.input_manifest "
        "--config config/config.yaml "
        "--output {output:q} > {log:q} 2>&1"
