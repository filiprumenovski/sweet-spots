import json
from pathlib import Path

import pandas as pd
import pytest

from regional_code_paper.reporting.audit import figure_checks, table_checks
from regional_code_paper.reporting.input_manifest import verify_manifest

REFERENCE = Path(__file__).resolve().parents[3] / "reference"


@pytest.mark.parametrize(
    "change", ["value", "ci_low", "missing", "extra", "duplicate", "metadata"]
)
def test_audit_catches_drift_anywhere_in_complete_ledger(change):
    expected = pd.read_csv(REFERENCE / "manuscript_result_ledger.csv")
    observed = expected.copy()
    if change in ("value", "ci_low"):
        row = observed[change].first_valid_index()
        observed.loc[row, change] += 0.001
    elif change == "missing":
        observed = observed.iloc[1:]
    elif change == "extra":
        observed.loc[0, "claim_id"] = "unexpected"
    elif change == "duplicate":
        observed = pd.concat([observed, observed.iloc[:1]])
    else:
        observed.loc[0, "inferential_unit"] = "wrong"
    assert not all(item["passed"] for item in table_checks(observed, expected, "claim_id"))


def test_complete_ledger_passes_and_is_json_serializable():
    expected = pd.read_csv(REFERENCE / "manuscript_result_ledger.csv")
    checks = table_checks(expected, expected, "claim_id")
    assert all(item["passed"] for item in checks)
    json.dumps(checks, allow_nan=False)
    assert expected.claim_id.nunique() == 82


def test_empty_figure_inventory_fails(tmp_path):
    manifest = tmp_path / "manifest.json"
    manifest.write_text('{"figures": []}')
    assert not all(item["passed"] for item in figure_checks(manifest, tmp_path))


def test_input_content_and_directory_membership_are_enforced():
    expected = json.loads((REFERENCE / "input_manifest.json").read_text())
    observed = json.loads(json.dumps(expected))
    observed["runtime"] = {"python": "another platform"}
    verify_manifest(observed, expected)
    observed["files"]["canonical_human_sequences"]["sha256"] = "changed"
    with pytest.raises(ValueError, match="files"):
        verify_manifest(observed, expected)
    observed["files"] = expected["files"]
    observed["directories"]["human_ortholog_alignments"]["files"] += 1
    with pytest.raises(ValueError, match="directories"):
        verify_manifest(observed, expected)


def test_nonfinite_observation_fails_but_report_remains_serializable():
    expected = pd.DataFrame({"claim_id": ["a"], "value": [0.8]})
    observed = expected.copy()
    observed.loc[0, "value"] = float("inf")
    checks = table_checks(observed, expected, "claim_id")
    assert not all(item["passed"] for item in checks)
    json.dumps(checks, allow_nan=False)
