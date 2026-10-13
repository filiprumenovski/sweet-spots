"""Build an offline reviewer report from audit evidence, figures and their sources."""

from __future__ import annotations

import argparse
import base64
import html
import json
import os
import re
import subprocess
from pathlib import Path

import pandas as pd

from ..core.config import load_config
from ..core.io import atomic_write_text, sha256_file, write_json


def graph_svg(dot: str) -> str:
    nodes = dict(re.findall(r'(\d+)\[label = "([^"]+)"', dot))
    edges = re.findall(r"(\d+) -> (\d+)", dot)
    levels = {}
    while len(levels) < len(nodes):
        ready = [
            n for n in nodes if n not in levels and all(a in levels for a, b in edges if b == n)
        ]
        if not ready:
            raise ValueError("Workflow graph has a cycle or an unknown dependency")
        for node in ready:
            levels[node] = max((levels[a] + 1 for a, b in edges if b == node), default=0)
    positions, counts = {}, {}
    for node in sorted(nodes, key=lambda n: (levels[n], nodes[n])):
        level = levels[node]
        positions[node] = (20 + 275 * level, 20 + 70 * counts.get(level, 0))
        counts[level] = counts.get(level, 0) + 1
    width, height = (
        (max(levels.values(), default=0) + 1) * 275,
        max(counts.values(), default=1) * 70 + 20,
    )
    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {width} {height}" '
        'role="img" aria-label="Complete workflow dependency graph">'
        '<defs><marker id="arrow" markerWidth="7" markerHeight="7" refX="6" '
        'refY="3" orient="auto"><path d="M0,0 L0,6 L6,3z" fill="#99a9ae"/></marker></defs>'
    ]
    for a, b in edges:
        x, y = positions[a]
        xx, yy = positions[b]
        parts.append(
            f'<path d="M{x + 240},{y + 22} C{x + 265},{y + 22} {xx - 25},{yy + 22} '
            f'{xx},{yy + 22}" fill="none" stroke="#bbc9cf" marker-end="url(#arrow)"/>'
        )
    for node, (x, y) in positions.items():
        parts.append(
            f'<rect x="{x}" y="{y}" width="240" height="44" rx="5" '
            'fill="#edf4f5" stroke="#547a86"/>'
            f'<text x="{x + 8}" y="{y + 27}" font-size="10" '
            f'font-family="sans-serif">{html.escape(nodes[node])}</text>'
        )
    return "".join(parts) + "</svg>"


def source_files(project: Path) -> list[Path]:
    paths = sorted(
        p
        for folder in ("src", "workflow", "config")
        for p in (project / folder).rglob("*")
        if p.suffix in (".py", ".smk", ".yaml", ".json") or p.name == "Snakefile"
    )
    return [*paths, project / "pyproject.toml", project / "uv.lock", project / "Makefile"]


def render(project: Path, results: Path, output: Path) -> None:
    output.mkdir(parents=True, exist_ok=True)
    audit = json.loads((results / "audit/reproduction_report.json").read_text())
    manifest = json.loads((project / "assets/review_figures/manifest.json").read_text())
    captions = json.loads((project / "assets/review_figures/captions.json").read_text())
    dot = subprocess.run(
        [
            "uv",
            "run",
            "--frozen",
            "snakemake",
            "-s",
            "workflow/Snakefile",
            "all",
            "--rulegraph",
            "dot",
        ],
        cwd=project,
        text=True,
        capture_output=True,
        check=True,
    ).stdout
    atomic_write_text(output / "workflow.dot", dot)
    svg = graph_svg(dot)
    atomic_write_text(output / "workflow.svg", svg)
    write_json(
        output / "source_hashes.json",
        {str(p.relative_to(project)): sha256_file(p) for p in source_files(project)},
    )

    def link(path: Path, label: str | None = None) -> str:
        url = os.path.relpath(path, output)
        return f'<a href="{html.escape(url, quote=True)}">{html.escape(label or path.name)}</a>'

    def thumbnail(path: Path) -> str:
        if not path.is_file():
            return "<p>Figure unavailable; inspect the audit failure.</p>"
        data = base64.b64encode(path.read_bytes()).decode()
        return (
            f'<img loading="lazy" src="data:image/png;base64,{data}" '
            f'alt="{html.escape(path.stem)}">'
        )

    passed = audit["computational_pass"]
    status = "PASS" if passed else "FAILED"
    parts = [
        '<!doctype html><html lang="en"><meta charset="utf-8">'
        '<meta name="viewport" content="width=device-width, initial-scale=1">'
        "<title>Sweet spots — reproducibility review</title><style>"
        "body{font:16px/1.55 system-ui,sans-serif;color:#18313f;background:#f8fafb;"
        "max-width:1280px;margin:auto;padding:28px}a{color:#176b89}h1{font-size:36px;line-height:1.2}"
        "h2{margin-top:2.3em}table{border-collapse:collapse;width:100%;font-size:13px;background:white}"
        "th,td{padding:9px;text-align:left;border-bottom:1px solid #dce5e8;vertical-align:top}"
        "td{overflow-wrap:anywhere}th{background:#e7eff2}nav{display:flex;gap:20px;flex-wrap:wrap}"
        "article{background:white;border:1px solid "
        "#dae5e8;border-radius:8px;padding:20px;margin:24px 0}"
        ".pair{display:grid;grid-template-columns:1fr 1fr;gap:20px}img{width:100%;height:auto}"
        ".pass{color:#147460}.fail{color:#a33832}.scroll{overflow:auto}.dag{min-width:1800px}"
        "code,pre{font-size:12px}pre{white-space:pre-wrap;overflow-wrap:anywhere}"
        "details{margin:15px 0}summary{cursor:pointer;font-weight:600}"
        "@media(max-width:720px){body{padding:14px}.pair{grid-template-columns:1fr}}"
        "</style><header><p>Sweet spots / scientific reproducibility</p>"
        "<h1>Complete reviewer report "
        f'<span class="{"pass" if passed else "fail"}">{status}</span></h1>'
        "<p>Frozen inputs → protein-grouped analyses → 82 claims → nine generated figures. "
        "The revised editorial set is shown with panel-level provenance.</p>"
        '<nav><a href="#checks">Checks</a><a href="#claims">Claim comparison</a>'
        '<a href="#figures">Figures and sources</a><a '
        'href="#workflow">Workflow</a></nav></header>'
    ]
    parts.append(
        '<h2 id="checks">Reproduction checks</h2><p>Reference: '
        f"<code>{html.escape(audit.get('reference_commit', 'unknown'))}</code>. "
        "Numerical tolerance: rtol=1e-9, atol=1e-10. Both full prediction tables must "
        "also match byte for byte. These checks do not independently validate the "
        "scientific claims.</p>"
    )
    parts.append("<table><tr><th>Check group</th><th>Passed / total</th></tr>")
    for name, checks in audit.get("checks", {}).items():
        count = sum(bool(c["passed"]) for c in checks)
        parts.append(f"<tr><td>{html.escape(name)}</td><td>{count} / {len(checks)}</td></tr>")
    parts.append(
        "</table><p>"
        + link(results / "audit/reproduction_report.json", "Complete audit JSON")
        + " · "
        + link(output / "source_hashes.json", "Exact source-file hashes")
        + "</p>"
    )
    failures = [
        c for group in audit.get("checks", {}).values() for c in group if not c["passed"]
    ]
    if failures or audit.get("error"):
        parts.append(
            '<pre class="fail">'
            + html.escape(json.dumps(failures or audit["error"], indent=2))
            + "</pre>"
        )
    old = pd.read_csv(project / "reference/manuscript_result_ledger.csv").set_index("claim_id")
    new = pd.read_csv(results / "tables/manuscript_result_ledger.csv").set_index("claim_id")
    parts.append(
        '<h2 id="claims">Every archived claim, before and after</h2>'
        "<p>Panel references below are the archived figure numbering. Revised figure "
        "lineage follows in the next section; some supporting results are not displayed "
        "in the revised artwork. Confidence bounds and metadata are included in the audit.</p>"
        '<div class="scroll"><table><tr><th>Claim / archived panel</th><th>Archived '
        "value [CI]</th>"
        "<th>Reproduced value [CI]</th><th>Source</th></tr>"
    )

    def value(row):
        if row is None:
            return "MISSING"
        number = f"{row['value']:.12g}"
        if pd.notna(row.ci_low):
            number += f" [{row.ci_low:.12g}, {row.ci_high:.12g}]"
        return number

    for claim, row in old.iterrows():
        observed = new.loc[claim] if claim in new.index and new.index.is_unique else None
        source = results / row.source_file.removeprefix("results/")
        parts.append(
            f"<tr><td>{html.escape(claim)}<br>{html.escape(row.figure_panel)}</td>"
            f"<td>{value(row)}</td><td>{value(observed)}</td><td>{link(source)}</td></tr>"
        )
    parts.append("</table></div>")
    parts.append(
        '<h2 id="figures">Generated and revised figures</h2>'
        "<p>Left: analysis figures generated by the workflow. Right: selected editorial "
        "exports from 11 October, retained and verified against their recorded sources. "
        "They are not claimed to be newly rendered. Captions are review drafts; S3 "
        "still requires its original method-specific details during manuscript integration.</p>"
    )
    for figure in manifest["figures"]:
        name = figure["figure"]
        generated = {"Figure_3": "Figure_4", "Figure_4": "Figure_3"}.get(name, name)
        parts.append(
            f'<article id="{name}"><h3>{name.replace("_", " ")}</h3><div class="pair">'
            f"<div><h4>Generated {generated.replace('_', ' ')} "
            "(archived numbering)</h4>"
            + thumbnail(results / "figures" / f"{generated}.png")
            + link(results / "figures" / f"{generated}.pdf", "Generated PDF")
            + "</div><div><h4>Retained revised artwork</h4>"
            + thumbnail(project / "assets/review_figures" / f"{name}.png")
            + link(project / "assets/review_figures" / f"{name}.pdf", "Revised PDF")
            + " · "
            + link(project / "assets/review_figures" / f"{name}.svg", "Editable SVG")
            + "</div></div><p>"
            + html.escape(captions[name.removeprefix("Figure_")])
            + "</p>"
        )
        parts.append(
            "<details><summary>Panel sources and historical plotting code</summary><p>"
            + link(project / figure["recipe"], "Archived plotting recipe (reference only)")
            + "</p><table><tr><th>Panel</th><th>Origin</th><th>Exact source files</th></tr>"
        )
        for panel in figure["panels"]:
            links = [
                link(
                    results / s.removeprefix("results/")
                    if s.startswith("results/")
                    else project / s,
                    s,
                )
                for s in panel["sources"]
            ]
            parts.append(
                f"<tr><td>{panel['panel']}</td><td>{html.escape(panel['origin'])}</td>"
                "<td>" + "<br>".join(links) + "</td></tr>"
            )
        parts.append("</table></details></article>")
    parts.append(
        '<h2 id="workflow">Complete workflow graph</h2><p>Derived from Snakemake rule '
        "dependencies. Five outer folds and other shards are collapsed to their rule names.</p>"
        '<div class="scroll"><div class="dag">'
        + svg
        + "</div></div><p>"
        + link(output / "workflow.dot", "Graph source (DOT)")
        + " · "
        + link(output / "workflow.svg", "Standalone SVG")
        + "</p>"
        "<p>Run <code>make check</code> for synthetic CI; <code>make reproduce</code> "
        "for full scientific recomputation; <code>make review</code> to build this report "
        "and enforce the audit gate. Full reproduction requires the complete frozen "
        "input tree.</p>"
        "</html>"
    )
    atomic_write_text(output / "index.html", "".join(parts))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, required=True)
    args = parser.parse_args()
    config = load_config(args.config)
    render(config.project_root, config.results_root, config.results_root / "review")


if __name__ == "__main__":
    main()
