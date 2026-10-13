import json
import subprocess
import sys
from pathlib import Path

from regional_code_paper.reporting.audit import main

ROOT = Path(__file__).resolve().parents[2]


def test_failed_audit_report_survives_snakemake_gate(tmp_path, monkeypatch):
    config = tmp_path / "config/config.yaml"
    config.parent.mkdir()
    config.write_text("source_root: inputs\nresults_root: results\n")
    report = tmp_path / "results/audit/reproduction_report.json"
    monkeypatch.setattr(
        sys, "argv", ["audit", "--config", str(config), "--output", str(report)]
    )
    main()
    assert json.loads(report.read_text())["computational_pass"] is False
    # Exercise the actual production gate, with only its upstream reference substituted.
    rule = (ROOT / "workflow/rules/products/audit.smk").read_text().split("rule audit:\n")[1]
    rule = rule.replace(
        "rules.audit_report.output", 'f"{RESULTS}/audit/reproduction_report.json"'
    )
    snakefile = tmp_path / "Snakefile"
    snakefile.write_text('RESULTS = "results"\nrule audit:\n' + rule)
    result = subprocess.run(
        [
            sys.executable,
            "-m",
            "snakemake",
            "-s",
            str(snakefile),
            "--cores",
            "1",
            "results/audit/PASS",
        ],
        cwd=tmp_path,
        capture_output=True,
        text=True,
    )
    assert result.returncode != 0
    assert "Reproduction failed" in result.stderr + result.stdout
    assert report.is_file()
    assert not (report.parent / "PASS").exists()
