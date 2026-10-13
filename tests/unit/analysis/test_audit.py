import pandas as pd

from regional_code_paper.reporting.audit import EXPECTED, numerical_checks


def test_audit_catches_scanner_drift_even_when_regional_checks_match():
    metrics = pd.DataFrame({"claim_id": list(EXPECTED), "value": list(EXPECTED.values())})
    metrics.loc[metrics.claim_id == "scanner.enhanced.average_precision", "value"] += 0.001
    failures = [item["claim_id"] for item in numerical_checks(metrics) if not item["passed"]]
    assert failures == ["scanner.enhanced.average_precision"]
