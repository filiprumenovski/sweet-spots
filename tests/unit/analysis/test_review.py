import json
import re
from pathlib import Path
from xml.etree import ElementTree

import pytest

from regional_code_paper.reporting.editorial import editorial_checks, equal_values
from regional_code_paper.reporting.review import graph_svg

ROOT = Path(__file__).resolve().parents[3]


def test_review_graph_contains_each_node_and_edge():
    dot = '0[label = "inputs"];1[label = "model"];2[label = "audit"];0 -> 1;1 -> 2;'
    svg = graph_svg(dot)
    ElementTree.fromstring(svg)
    assert all(name in svg for name in ("inputs", "model", "audit"))
    assert len(re.findall("marker-end=", svg)) == 2
    with pytest.raises(ValueError, match="cycle"):
        graph_svg(dot + "2 -> 0;")


def test_json_source_comparison_checks_all_fields_and_values():
    assert equal_values({"score": [0.8, None]}, {"score": [0.8 + 1e-12, None]})
    assert not equal_values({"score": [0.8]}, {"score": [0.81]})
    assert not equal_values({"score": 0.8}, {"score": 0.8, "extra": 1})
    assert not equal_values(True, 1)


def test_editorial_panel_inventory_is_complete_and_assets_unchanged():
    manifest = json.loads((ROOT / "assets/review_figures/manifest.json").read_text())
    assert sum(len(f["panels"]) for f in manifest["figures"]) == 41
    checks = editorial_checks(ROOT, ROOT / "results")
    # These checks also run before a clean full build, when result sources do not exist.
    assets = [c for c in checks if c.get("method") == "asset SHA-256" or "check" in c]
    assert all(check["passed"] for check in assets), [c for c in assets if not c["passed"]]
