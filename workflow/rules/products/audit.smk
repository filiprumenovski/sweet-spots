"""Write durable evidence before the separate fail-closed gate."""

rule audit_report:
    input:
        ledger=rules.result_ledger_csv.output,
        parquet=rules.result_ledger_parquet.output,
        manifest=rules.input_manifest.output,
        figures=rules.figures.output,
        predictions=f"{RESULTS}/analysis/scanner/outer_predictions.csv",
        regional=f"{RESULTS}/analysis/regional_code/held_out_tiles.csv",
        catalogue=f"{RESULTS}/analysis/region_definition_sensitivity/catalogue_grid.csv",
        models=f"{RESULTS}/analysis/region_definition_sensitivity/model_grid.csv",
        reference=sorted(str(path) for path in __import__("pathlib").Path("reference").glob("*")),
    output:
        f"{RESULTS}/audit/reproduction_report.json",
    log:
        f"{RESULTS}/logs/audit.log",
    resources:
        mem_mb=4000,
        runtime=30,
        disk_mb=4000,
    shell:
        "uv run --frozen python -m regional_code_paper.reporting.audit "
        "--config config/config.yaml --output {output:q} > {log:q} 2>&1"

rule audit:
    input:
        rules.audit_report.output,
    output:
        f"{RESULTS}/audit/PASS",
    run:
        import json
        from pathlib import Path
        report = json.loads(Path(input[0]).read_text())
        if not report["computational_pass"]:
            raise ValueError(f"Reproduction failed; inspect {input[0]}")
        Path(output[0]).write_text("All frozen scientific checks passed.\n")
