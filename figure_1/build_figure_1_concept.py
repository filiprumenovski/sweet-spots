"""Build Figure 1 (concept figure), version 2. Four panels: A enzyme and reaction,
B substrate map, C FOXK1, D cross-species transfer (letters run down each column).

Journal styling: white background, bold panel letters, no in-panel titles;
the legend carries the explanation. Starts from the vector re-assembly at
~/Documents/papers/Sweet_Spots_Figure_1_Recreated/Figure_1_Recreated.svg
(1500 x 1050 layout, vector text, embedded image assets) and changes:

  A  OGA and the cycling scheme are removed. Top to bottom: the full-length
     OGT dimer (PDB 7YEA, render_ogt_assembly.py), the reaction it catalyses
     as chemical structures (render_chemistry.py), and every validated
     substrate protein as a dot grouped by UniProt location
     (compute_substrate_map.py).
  B  FOXK1 throughout: the full AlphaFold model with every O-GlcNAc site and the three
     consensus regions (render_foxk1_model.py); inset of Thr683 with GlcNAc modeled from
     PDB 5UN8 (render_glcnac_site.py); residue track below (compute_track.py).
  C  Layout of the first iteration (species drawings, composition tiles, labelled
     transfer arrows) with transfer AUROCs and 95% intervals against the prespecified
     0.80 threshold; species named Homo sapiens and Oryza sativa.

Run from this folder: python3 build_figure_1_concept.py
Writes Figure_1_Concept_v2.svg, .pdf and .png (3000 x 2100 px, 10 x 7 in at 300 dpi).
"""
from __future__ import annotations

import base64
import io
import json
import math
import re
import subprocess
from pathlib import Path

from PIL import Image

HERE = Path(__file__).resolve().parent
ORGANISMS = HERE / "assets/organisms"   # human and rice drawings carried over from the earlier Figure 1 assembly
OUT = HERE / "Figure_1_Concept_v2"

INK = "#142b38"
ARROW = "#18313f"
MUTED = "#5d7380"
GLCNAC = "#22a6ed"
GLCNAC_EDGE = "#0875b8"

# Result ledger (results/tables/manuscript_result_ledger.csv), Figure 4A / Table S5
AUROC = {
    "h2r": (0.8065612089584693, 0.7677481780591642, 0.8426008963648411),  # transfer.human_to_target.rice
    "r2h": (0.8113448310305266, 0.7985542598101216, 0.8234664753668065),  # transfer.rice_to_human.human
}
THRESHOLD = 0.80
def esc(s: str) -> str:
    return s.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")


def png_uri(img: Image.Image) -> str:
    buf = io.BytesIO()
    img.save(buf, format="PNG", optimize=True)
    return "data:image/png;base64," + base64.b64encode(buf.getvalue()).decode()


def crop_alpha(path: Path) -> Image.Image:
    im = Image.open(path).convert("RGBA")
    return im.crop(im.split()[3].getbbox())


def recreated_assets() -> dict[str, str]:
    return {p.stem: "data:image/png;base64," + base64.b64encode(p.read_bytes()).decode()
            for p in sorted(ORGANISMS.glob("*.png"))}


def text(x, y, s, size=15, anchor="middle", weight=400, fill=INK, style="", cls="", raw=False):
    body = s if raw else esc(s)
    st = f' font-style="{style}"' if style else ""
    c = f' class="{cls}"' if cls else ""
    return (f'<text x="{x:.1f}" y="{y:.1f}" font-size="{size}" text-anchor="{anchor}" '
            f'font-weight="{weight}" fill="{fill}"{st}{c}>{body}</text>')


def text_w(s: str, size: float) -> float:
    """Rough Arial advance width."""
    narrow = set("ilIjtfr.,·:;'| ")
    wide = set("MWmw")
    w = 0.0
    for ch in s:
        w += 0.30 if ch in narrow else 0.86 if ch in wide else 0.70 if ch.isupper() else 0.56
    return w * size


def glcnac(x, y, s=11):
    return (f'<rect x="{x - s / 2:.1f}" y="{y - s / 2:.1f}" width="{s}" height="{s}" rx="2" '
            f'fill="{GLCNAC}" stroke="{GLCNAC_EDGE}" stroke-width="1"/>')


# ---------------------------------------------------------------- shared style
# 1500 x 1050 units print at 10 x 7 in, so 1 unit = 0.48 pt.
LETTER = 26      # ~12.5 pt bold panel letters
LABEL = 15       # ~7.2 pt
SMALL = 13       # ~6.2 pt

# Grid: column A x 28-640, column B/C x 700-1472; top 36, bottom 1026.
# Panel C and the substrate map in A share one baseline (ROW2).
COL_A = (28.0, 640.0)
COL_R = (700.0, 1472.0)
TOP = 40.0
ROW2 = 690.0


def letter(x, y, s):
    return text(x, y, s, LETTER, "start", weight=700)


# ---------------------------------------------------------------- panel A
COMPARTMENTS = [  # key in substrate_map.json, label, colour
    ("nucleus", "Nucleus", "#2f5d9e"),
    ("nucleocytoplasmic", "Nucleus and cytoplasm", "#5e8fc8"),
    ("cytoplasm", "Cytoplasm", "#2a8a8f"),
    ("membranes & secretory", "Membranes and secretory", "#7a62b8"),
    ("mitochondria", "Mitochondria", "#c47f17"),
    ("not annotated", "Not annotated", "#9fb0bb"),
]


def embed_svg(path: Path, x: float, y: float, k: float) -> tuple[str, float, float]:
    """Inline an RDKit SVG at (x, y) scaled by k; returns markup, width, height."""
    s = path.read_text()
    w = float(re.search(r"width='([\d.]+)px'", s).group(1))
    h = float(re.search(r"height='([\d.]+)px'", s).group(1))
    body = s[s.index(">", s.index("<svg")) + 1: s.rindex("</svg>")]
    body = re.sub(r"<rect[^>]*style='opacity:1.0;fill:#FFFFFF[^>]*/>", "", body)
    return (f'<svg x="{x:.1f}" y="{y:.1f}" width="{w * k:.1f}" height="{h * k:.1f}" '
            f'viewBox="0 0 {w} {h}" overflow="visible">{body}</svg>', w * k, h * k)


def panel_a(dimer: Image.Image) -> list[str]:
    L, R = COL_A
    cx = (L + R) / 2
    o: list[str] = [letter(L, TOP, "A")]

    # 1. Enzyme: full-length OGT dimer, filling the column
    dw = 548.0
    dh = dw * dimer.height / dimer.width
    dx, dy = cx - dw / 2, TOP + 6
    o.append(f'<image x="{dx:.1f}" y="{dy:.1f}" width="{dw:.1f}" height="{dh:.1f}" '
             f'preserveAspectRatio="xMidYMid meet" xlink:href="{png_uri(dimer)}"/>')
    o.append(text(cx + 150, dy + 26, "OGT homodimer (PDB 7YEA)", LABEL, "start"))

    # 2. Reaction at the active site, spanning the column
    ry0 = dy + dh + 26
    kc = 1.1
    acc, aw, ah = embed_svg(HERE / "assets/chem_acceptor.svg", L + 12, ry0 + 4, kc)
    px = R - 220.0
    prod, pw, ph = embed_svg(HERE / "assets/chem_product.svg", px, ry0, kc)
    o += [acc, prod]
    ay = ry0 + 4 + ah * 0.60
    a0, a1 = L + 12 + aw + 18, px - 14
    mid = (a0 + a1) / 2
    o.append(f'<path d="M{a0:.1f} {ay:.1f} L{a1:.1f} {ay:.1f}" stroke="{ARROW}" stroke-width="1.8" '
             f'fill="none" marker-end="url(#arrow)"/>')
    o.append(text(mid, ay + 20, "OGT", LABEL))
    o.append(f'<path d="M{a0 + 10:.1f} {ay - 34:.1f} Q{mid:.1f} {ay + 4:.1f} {a1 - 10:.1f} {ay - 34:.1f}" '
             f'stroke="{ARROW}" stroke-width="1.2" fill="none" marker-end="url(#arrowsm)"/>')
    o.append(text(a0 + 4, ay - 40, "UDP-GlcNAc", SMALL, "start"))
    o.append(text(a1 - 4, ay - 40, "UDP", SMALL, "end"))
    o.append(text(L + 12 + aw / 2, ry0 + 4 + ah + 16, "Ser", LABEL))
    o.append(text(px + pw / 2, ry0 + ph + 14, "O-GlcNAc-Ser", LABEL))
    gx, gy = px + pw + 14, ry0 + ph * 0.66
    o.append(glcnac(gx, gy, 14))
    o.append(text(gx + 12, gy + 5, "GlcNAc", SMALL, "start"))

    # 3. Substrates: every reported human O-GlcNAc protein (O-GlcNAc Database plus the
    # validated set), by location; solid = validated mapped sites, pale = none mapped
    sm = json.loads((HERE / "data/substrate_map.json").read_text())
    o.append(letter(L, ROW2, "B"))
    head_b = f"OGT, one enzyme with {sm['n_proteins']:,} human substrates"
    o.append(text(L + 30, ROW2, head_b, LABEL, "start"))
    k = 0.555
    g = sm["groups"]

    def tint(col, f=0.68):
        c = [int(col[i:i + 2], 16) for i in (1, 3, 5)]
        return "#" + "".join(f"{round(v + (255 - v) * f):02x}" for v in c)

    def cluster(key, label, col, ccx, ccy, ly):
        grp = g[key]
        out = []
        for fill, keep in ((tint(col), lambda c: c["sites"] == 0), (col, lambda c: c["sites"] > 0)):
            out.append(f'<g fill="{fill}">')
            for c in grp["circles"]:
                if keep(c):
                    out.append(f'<circle cx="{ccx + c["x"] * k:.2f}" cy="{ccy + c["y"] * k:.2f}" r="{c["r"] * k:.2f}"/>')
            out.append("</g>")
        out.append(text(ccx, ly, label, SMALL))
        out.append(text(ccx, ly + 15, f'{grp["n_proteins"]:,} ({grp["n_mapped"]:,} mapped)', SMALL, fill=MUTED))
        return out

    r = {key: g[key]["extent"] * k for key, _, _ in COMPARTMENTS}
    # each row ordered largest to smallest cluster
    row1 = sorted(COMPARTMENTS[:3], key=lambda c: -r[c[0]])
    row2 = sorted(COMPARTMENTS[3:], key=lambda c: -r[c[0]])
    gap = 22.0
    # row 1: the three large clusters, centred in the column
    m1 = max(r[key] for key, _, _ in row1)
    cy1 = ROW2 + 14 + m1
    ly1 = cy1 + m1 + 15
    w1 = sum(2 * r[key] for key, _, _ in row1) + gap * (len(row1) - 1)
    x = (L + R) / 2 - w1 / 2
    for key, label, col in row1:
        x += r[key]
        o += cluster(key, label, col, x, cy1, ly1)
        x += r[key] + gap
    # row 2: three small clusters plus the size key, packed and centred
    m2 = max(r[key] for key, _, _ in row2)
    cy2 = ly1 + 26 + m2
    ly2 = cy2 + m2 + 15
    widths = [max(2 * r[key], text_w(label, SMALL) * 0.92,
                  text_w(f'{g[key]["n_proteins"]:,} ({g[key]["n_mapped"]:,} mapped)', SMALL) * 0.92)
              for key, label, _ in row2]
    key_w = 104.0
    gap2 = 18.0
    w2 = sum(widths) + key_w + gap2 * len(row2)
    x = (L + R) / 2 - w2 / 2
    for (key, label, col), w in zip(row2, widths):
        o += cluster(key, label, col, x + w / 2, cy2, ly2)
        x += w + gap2
    kx = x + 14
    ky = cy2 + m2 - 2
    o.append(text(x + key_w / 2, ly2, "Mapped sites", SMALL))
    for n, dxk in ((0, -4), (1, 14), (10, 38), (100, 76)):
        rr = 1.45 * math.sqrt(max(n, 1)) * k
        if n == 0:
            o.append(f'<circle cx="{kx + dxk:.1f}" cy="{ky - rr:.1f}" r="{rr:.2f}" fill="#a9b5be"/>')
        else:
            o.append(f'<circle cx="{kx + dxk:.1f}" cy="{ky - rr:.1f}" r="{rr:.2f}" fill="none" '
                     f'stroke="{INK}" stroke-width="0.9"/>')
        o.append(text(kx + dxk, ly2 + 15, str(n), SMALL, fill=MUTED))
    print("panel B bottom", round(ly2 + 15), "row2 width", round(w2))
    return o


# ---------------------------------------------------------------- panel B
def tangents(c1, r1, c2, r2):
    """External tangent points between two circles: [(p1, p2), (p1, p2)]."""
    dx, dy = c2[0] - c1[0], c2[1] - c1[1]
    d = math.hypot(dx, dy)
    ux, uy = dx / d, dy / d
    k = (r1 - r2) / d
    h = math.sqrt(1 - k * k)
    res = []
    for s in (1, -1):
        nx, ny = k * ux - s * h * uy, k * uy + s * h * ux
        res.append(((c1[0] + r1 * nx, c1[1] + r1 * ny), (c2[0] + r2 * nx, c2[1] + r2 * ny)))
    return res


REGION_NAVY = "#1c5e8a"
CHAIN = "#3fa3a6"
GENE_B = "foxk1"
DOMAINS_B = [("FHA", 104, 206), ("Forkhead", 304, 396)]  # AlphaFold pLDDT > 70 segments
LENS_SITE = 683
LENS_R = 108.0
LENS_C = (COL_R[1] - LENS_R - 6, 168.0)


def protein_track(o: list[str], y0: float) -> None:
    """Residue track of the panel B protein (same protein as the drawing above)."""
    t = json.loads((HERE / f"data/{GENE_B}_track.json").read_text())
    L, R = COL_R
    tx0, tx1 = L + 92, R - 6
    n = t["length"]
    X = lambda p: tx0 + (p - 1) / (n - 1) * (tx1 - tx0)
    dy = 24.0
    top = y0 + 14
    o.append(text(L, y0, f'{t["gene"]} ({t["accession"]}), {n} aa', LABEL, "start"))
    for a, b in t["regions"]:
        o.append(f'<rect x="{X(a) - 2:.1f}" y="{top + dy - 4:.1f}" width="{X(b) - X(a) + 4:.1f}" '
                 f'height="{dy * 2 + 2:.1f}" fill="{REGION_NAVY}" fill-opacity="0.18"/>')
    # lane 0: domains
    yy = top + dy / 2
    o.append(f'<line x1="{tx0}" y1="{yy:.1f}" x2="{tx1}" y2="{yy:.1f}" stroke="{CHAIN}" stroke-width="2"/>')
    o.append(text(tx0 - 10, yy + 4.5, "Domains", SMALL, "end"))
    for name, a, b in DOMAINS_B:
        o.append(f'<rect x="{X(a):.1f}" y="{yy - 8:.1f}" width="{X(b) - X(a):.1f}" height="16" rx="3" fill="#2a8c8f"/>')
        o.append(f'<text x="{(X(a) + X(b)) / 2:.1f}" y="{yy + 4.5:.1f}" font-size="{SMALL - 1}" '
                 f'text-anchor="middle" fill="white">{name}</text>')
    for name, pos, col, i in (("Ser/Thr", t["st"], "#9aa7b0", 1), ("O-GlcNAc", t["oglcnac"], GLCNAC_EDGE, 2)):
        yy = top + i * dy + dy / 2
        o.append(f'<line x1="{tx0}" y1="{yy:.1f}" x2="{tx1}" y2="{yy:.1f}" stroke="#c9d3da" stroke-width="1"/>')
        o.append(text(tx0 - 10, yy + 4.5, name, SMALL, "end"))
        for p in pos:
            o.append(f'<line x1="{X(p):.1f}" y1="{yy - 7:.1f}" x2="{X(p):.1f}" y2="{yy + 7:.1f}" '
                     f'stroke="{col}" stroke-width="1.5"/>')
    ay = top + 3 * dy + 2
    o.append(f'<line x1="{tx0}" y1="{ay:.1f}" x2="{tx1}" y2="{ay:.1f}" stroke="{ARROW}" stroke-width="1"/>')
    for p in (1, 100, 200, 300, 400, 500, 600, n):
        o.append(f'<line x1="{X(p):.1f}" y1="{ay:.1f}" x2="{X(p):.1f}" y2="{ay + 4:.1f}" stroke="{ARROW}" stroke-width="1"/>')
        o.append(text(X(p), ay + 16, str(p), SMALL))


def panel_b(assets: dict[str, str], pep: Image.Image) -> list[str]:
    o: list[str] = [letter(COL_R[0], TOP, "C")]
    t = json.loads((HERE / f"data/{GENE_B}_track.json").read_text())
    model = crop_alpha(HERE / "assets/FOXK1_model.png")
    s2d = json.loads((HERE / "assets/FOXK1_model_sites2d.json").read_text())["sites"]
    bx, by, bw = COL_R[0] + 4, 52.0, 548.0
    bh = bw * model.height / model.width
    P = lambda site: (bx + s2d[str(site)][0] * bw, by + s2d[str(site)][1] * bh)

    c1, r1 = P(LENS_SITE), 11.0
    c2, r2 = LENS_C, LENS_R
    (a1, a2), (b1, b2) = tangents(c1, r1, c2, r2)
    o.append(f'<path d="M{a1[0]:.1f} {a1[1]:.1f} L{a2[0]:.1f} {a2[1]:.1f} '
             f'A{r2} {r2} 0 0 0 {b2[0]:.1f} {b2[1]:.1f} L{b1[0]:.1f} {b1[1]:.1f} Z" fill="#eaf3f9"/>')
    o.append(f'<image x="{bx}" y="{by}" width="{bw}" height="{bh:.1f}" preserveAspectRatio="xMidYMid meet" '
             f'xlink:href="{png_uri(model)}"/>')
    for p in t["oglcnac"]:
        x, y = P(p)
        o.append(f'<rect x="{x - 5:.1f}" y="{y - 5:.1f}" width="10" height="10" rx="1.5" '
                 f'fill="{GLCNAC}" stroke="white" stroke-width="1.1"/>')
    for p, q in ((a1, a2), (b1, b2)):
        o.append(f'<line x1="{p[0]:.1f}" y1="{p[1]:.1f}" x2="{q[0]:.1f}" y2="{q[1]:.1f}" '
                 f'stroke="{MUTED}" stroke-width="1" stroke-dasharray="4 3"/>')
    o.append(f'<circle cx="{c1[0]:.1f}" cy="{c1[1]:.1f}" r="{r1}" fill="none" stroke="{ARROW}" stroke-width="1.4"/>')
    site = crop_alpha(HERE / f"assets/FOXK1_Thr{LENS_SITE}_glcnac.png")
    o.append(f'<circle cx="{c2[0]}" cy="{c2[1]}" r="{r2}" fill="white"/>')
    pw = 2 * r2 - 22
    ph = pw * site.height / site.width
    o.append(f'<g clip-path="url(#lensclip)"><image x="{c2[0] - pw / 2:.1f}" y="{c2[1] - ph / 2:.1f}" '
             f'width="{pw:.1f}" height="{ph:.1f}" preserveAspectRatio="xMidYMid meet" '
             f'xlink:href="{png_uri(site)}"/></g>')
    o.append(f'<circle cx="{c2[0]}" cy="{c2[1]}" r="{r2}" fill="none" stroke="{ARROW}" stroke-width="1.4"/>')
    o.append(text(c2[0], c2[1] - r2 - 10, f'GlcNAc-Thr{LENS_SITE} <tspan fill="{MUTED}" font-size="{SMALL}">(modeled)</tspan>',
                  LABEL, raw=True))
    # key under the inset
    kx, ky = c2[0] - r2 + 20, c2[1] + r2 + 34
    o.append(f'<rect x="{kx - 6:.1f}" y="{ky - 6:.1f}" width="12" height="12" rx="1.5" fill="{GLCNAC}"/>')
    o.append(text(kx + 14, ky + 5, "O-GlcNAc site", LABEL, "start"))
    o.append(f'<line x1="{kx - 8:.1f}" y1="{ky + 28:.1f}" x2="{kx + 6:.1f}" y2="{ky + 28:.1f}" '
             f'stroke="{REGION_NAVY}" stroke-width="4" stroke-linecap="round"/>')
    o.append(text(kx + 14, ky + 33, "Consensus region", LABEL, "start"))

    protein_track(o, max(by + bh + 28, 470.0))
    return o


# ---------------------------------------------------------------- panel C
def panel_c(assets: dict[str, str]) -> list[str]:
    """Layout of the first iteration (icons, composition tiles, labelled transfer
    arrows, AUROC plot between them), restyled for the white-background figure."""
    L, R = COL_R
    o: list[str] = [letter(L, ROW2, "D")]
    blue, green = "#1d5f8a", "#4c8a4f"
    top = ROW2 + 22
    hh = 246.0
    hw = hh * 137 / 238
    rw = hh * 148 / 244
    o.append(f'<image x="{L:.1f}" y="{top:.1f}" width="{hw:.1f}" height="{hh:.1f}" '
             f'preserveAspectRatio="xMidYMid meet" xlink:href="{assets["Human"]}"/>')
    o.append(f'<image x="{R - rw:.1f}" y="{top:.1f}" width="{rw:.1f}" height="{hh:.1f}" '
             f'preserveAspectRatio="xMidYMid meet" xlink:href="{assets["Rice"]}"/>')
    o.append(text(L + hw / 2, top + hh + 22, "Homo sapiens", LABEL + 1, style="italic"))
    o.append(text(R - rw / 2, top + hh + 22, "Oryza sativa", LABEL + 1, style="italic"))

    m0, m1 = L + hw + 18, R - rw - 18          # middle column
    mc = (m0 + m1) / 2
    # composition tiles
    seq = [("W", "#f5dfad"), ("A", "#c8e4ee"), ("S", "#edc6c0"), ("A", "#f5dfad"), ("H", "#c8e4ee"),
           ("K", "#f5dfad"), ("A", "#edc6c0"), ("S", "#c8e4ee"), ("G", "#d9d0ed"), ("W", "#f5dfad"),
           ("A", "#c8e4ee"), ("R", "#dedede")]
    tw, step = 24.0, 27.0
    tx = mc - (step * len(seq) - (step - tw)) / 2
    for i, (aa, col) in enumerate(seq):
        x = tx + i * step
        o.append(f'<rect x="{x:.1f}" y="{top:.1f}" width="{tw}" height="27" rx="1" fill="{col}" stroke="#c3cdd5"/>')
        o.append(text(x + tw / 2, top + 19, aa, LABEL))

    # top arrow: human-trained model predicts rice
    ya = top + 74
    o.append(text(mc, ya - 12, 'Trained on <tspan font-style="italic">H. sapiens</tspan> → '
                  'predicts <tspan font-style="italic">O. sativa</tspan>', LABEL, raw=True))
    o.append(f'<path d="M{m0:.1f} {ya:.1f} L{m1:.1f} {ya:.1f}" stroke="{blue}" stroke-width="2.2" '
             f'fill="none" marker-end="url(#arrowblue)"/>')

    # AUROC plot
    x0, x1 = m0 + 14, m1 - 14

    def X(v):  # AUROC 0.5 .. 1.0
        return x0 + (v - 0.5) / 0.5 * (x1 - x0)

    pt = ya + 18
    ay = pt + 104
    o.append(f'<line x1="{X(0.5):.1f}" y1="{pt:.1f}" x2="{X(0.5):.1f}" y2="{ay:.1f}" stroke="#c3cdd5" stroke-width="1"/>')
    o.append(f'<line x1="{X(THRESHOLD):.1f}" y1="{pt:.1f}" x2="{X(THRESHOLD):.1f}" y2="{ay:.1f}" '
             f'stroke="{MUTED}" stroke-width="1.1" stroke-dasharray="4 3"/>')
    o.append(text(X(0.5) + 3, pt + 10, "chance", SMALL, "start", fill=MUTED))
    o.append(text(X(THRESHOLD) - 5, pt + 10, "prespecified 0.80", SMALL, "end", fill=MUTED))
    for key, col, y, lab in (("h2r", blue, pt + 40, "H. sapiens → O. sativa"),
                             ("r2h", green, pt + 76, "O. sativa → H. sapiens")):
        est, lo, hi = AUROC[key]
        o.append(text(X(0.5) + 8, y + 5, lab, SMALL + 1, "start", style="italic"))
        o.append(f'<line x1="{X(lo):.1f}" y1="{y:.1f}" x2="{X(hi):.1f}" y2="{y:.1f}" stroke="{col}" stroke-width="2.2"/>')
        for v in (lo, hi):
            o.append(f'<line x1="{X(v):.1f}" y1="{y - 6:.1f}" x2="{X(v):.1f}" y2="{y + 6:.1f}" stroke="{col}" stroke-width="1.8"/>')
        o.append(f'<circle cx="{X(est):.1f}" cy="{y:.1f}" r="6" fill="{col}" stroke="white" stroke-width="1.4"/>')
        o.append(text(X(hi) + 8, y + 5.5, f"{est:.2f}", LABEL + 1, "start", weight=700, fill=col))
        o.append(text(X(hi) + 46, y + 5, f"({lo:.2f}–{hi:.2f})", SMALL, "start", fill=MUTED))
    o.append(f'<line x1="{x0:.1f}" y1="{ay:.1f}" x2="{x1:.1f}" y2="{ay:.1f}" stroke="{ARROW}" stroke-width="1"/>')
    for v in (0.5, 0.6, 0.7, 0.8, 0.9, 1.0):
        o.append(f'<line x1="{X(v):.1f}" y1="{ay:.1f}" x2="{X(v):.1f}" y2="{ay + 4:.1f}" stroke="{ARROW}" stroke-width="1"/>')
        o.append(text(X(v), ay + 17, f"{v:.1f}", SMALL, fill=MUTED))
    o.append(text(x1, ay - 6, "AUROC", SMALL, "end", weight=700))

    # bottom arrow: rice-trained model predicts human
    yb = ay + 44
    o.append(f'<path d="M{m1:.1f} {yb:.1f} L{m0:.1f} {yb:.1f}" stroke="{green}" stroke-width="2.2" '
             f'fill="none" marker-end="url(#arrowgreen)"/>')
    o.append(text(mc, yb + 22, 'Trained on <tspan font-style="italic">O. sativa</tspan> → '
                  'predicts <tspan font-style="italic">H. sapiens</tspan>', LABEL, raw=True))
    o.append(text(mc, 1022, "Sequence composition only · alignment-free prediction", SMALL, fill=MUTED))
    return o


def defs() -> str:
    return f"""<defs>
<marker id="arrow" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="7" markerHeight="7" orient="auto-start-reverse"><path d="M0 0 L10 5 L0 10 Z" fill="{ARROW}"/></marker>
<marker id="arrowblue" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="6" markerHeight="6" orient="auto-start-reverse"><path d="M0 0 L10 5 L0 10 Z" fill="#1d5f8a"/></marker>
<marker id="arrowgreen" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="6" markerHeight="6" orient="auto-start-reverse"><path d="M0 0 L10 5 L0 10 Z" fill="#4c8a4f"/></marker>
<marker id="arrowsm" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="6" markerHeight="6" orient="auto-start-reverse"><path d="M0 0 L10 5 L0 10 Z" fill="{ARROW}"/></marker>
<clipPath id="lensclip"><circle cx="{LENS_C[0]}" cy="{LENS_C[1]}" r="{LENS_R - 1}"/></clipPath>
</defs>"""


def main() -> None:
    assets = recreated_assets()
    dimer = crop_alpha(HERE / "assets/OGT_7YEA_dimer.png")
    pep = crop_alpha(HERE / "structures/Glycopeptide_5UN8_chainE.png")

    parts = [
        '<svg xmlns="http://www.w3.org/2000/svg" xmlns:xlink="http://www.w3.org/1999/xlink" '
        'width="3000" height="2100" viewBox="0 0 1500 1050">',
        "<title>Figure 1</title>",
        "<desc>Figure 1 concept, version 2. OGT dimer rendered from PDB 7YEA with ligands placed by "
        "superposition of PDB 3PE4; glycopeptide rendered from PDB 5UN8 chain E; schematic protein, "
        "human and rice drawings carried over from the earlier assembly.</desc>",
        defs(),
        f"<style>text{{font-family:Arial,Helvetica,sans-serif;fill:{INK}}}</style>",
        '<rect width="1500" height="1050" fill="white"/>',
        *panel_a(dimer),
        *panel_b(assets, pep),
        *panel_c(assets),
        "</svg>",
    ]
    svg = "\n".join(parts)
    OUT.with_suffix(".svg").write_text(svg)
    subprocess.run(["rsvg-convert", "-f", "pdf", "--dpi-x", "300", "--dpi-y", "300", "-o", str(OUT.with_suffix(".pdf")),
                    str(OUT.with_suffix(".svg"))], check=True)
    subprocess.run(["rsvg-convert", "-b", "white", "-w", "3000", "-h", "2100", "-o",
                    str(OUT.with_suffix(".png")), str(OUT.with_suffix(".svg"))], check=True)
    (HERE / "Assembly_Record.json").write_text(json.dumps({
        "built_by": "build_figure_1_concept.py",
        "canvas_viewbox": [1500, 1050], "png_pixels": [3000, 2100], "print_inches": [10, 7],
        "panel_A_structure": "PDB 7YEA full-length OGT homodimer (render_ogt_assembly.py); UDP and "
                             "acceptor peptide from PDB 3PE4 superposed per protomer",
        "panel_A_map": "data/substrate_map.json from compute_substrate_map.py (3,532 proteins with "
                        "validated O-GlcNAc sites; UniProt reviewed human subcellular locations)",
        "panel_A_reaction": "assets/chem_*.svg from render_chemistry.py (RDKit; stereo checked by CIP)",
        "panel_B": "FOXK1 (P85037): AlphaFold AF-P85037-F1 v6 FHA (104-206) and forkhead (304-396) domains; "
                   "full model with site markers; GlcNAc on Thr683 fitted from PDB 5UN8",
        "panel_C_auroc": {k: {"estimate": v[0], "ci95": [v[1], v[2]]} for k, v in AUROC.items()},
        "panel_C_source": "manuscript_result_ledger.csv rows transfer.human_to_target.rice and "
                          "transfer.rice_to_human.human (Figure 4A, Table S5)",
        "carried_over_rasters": ["Human", "Rice"],
        "carried_over_from": "assets/organisms/",
    }, indent=2) + "\n")
    print("wrote", OUT.with_suffix(".svg"), OUT.with_suffix(".pdf"), OUT.with_suffix(".png"))


if __name__ == "__main__":
    main()
