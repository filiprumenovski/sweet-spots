"""Render the submitted Figures 2-6 and S1-S4 from assets/review_figures at Cell Press size with a standard type scale.
Print width 174 mm. Type scale at print size: panel letters 8 pt bold, titles 7 pt,
labels/ticks 6.5 pt, annotations 6 pt minimum. Arial throughout. Title band and
caption block are deleted (not just cropped). usage: render_cell_figs.py SRC_DIR OUT_DIR"""
import pathlib
import re
import subprocess
import sys

sys.path.insert(0, str(pathlib.Path(__file__).parent))
from despace import despace

SRC, OUT = map(pathlib.Path, sys.argv[1:3])
W_DESIGN = 540.0; W_MM = 174.0
k = (W_MM / 25.4 * 72) / W_DESIGN          # design pt -> print pt
def tier(sz):                                # design size -> target print pt
    if sz >= 8.75: return 8.0                # panel letters (9) / figure number (10)
    if sz >= 7.75: return 7.0                # panel titles (8, 8.5)
    if sz >= 6.9:  return 6.5                # axis labels, ticks, legends (7, 7.5)
    return 6.0                               # annotations: floor
# paper name: (review source, top cut, bottom keep) in design pt
PLAN = {"Figure_2": ("Figure_1", 29, 571), "Figure_3": ("Figure_2", 29, 552), "Figure_4": ("Figure_3", 29, 376.5),
        "Figure_5": ("Figure_4", 29, 496), "Figure_6": ("Figure_5", 29, 516.3),
        "Figure_S1": ("Figure_S1", 29, 490), "Figure_S2": ("Figure_S4", 29, 394),
        "Figure_S3": ("Figure_S3", 30, 502), "Figure_S4": ("Figure_S2", 29, 492.5)}
TEXT = re.compile(r'<text\b([^>]*)>(.*?)</text>', re.S)
report = {}
for name, (src, top, keep) in PLAN.items():
    s = (SRC / f"{src}.svg").read_text()
    s = s.replace('font-family="Liberation Sans"', 'font-family="Arial"')
    removed = resized = 0
    def fix(m):
        global removed, resized
        attrs, inner = m.group(1), m.group(2)
        tm = re.search(r'transform="matrix\(([^)]*)\)"', attrs)
        fs = re.search(r'font-size="([\d.]+)"', attrs)
        if not (tm and fs): return m.group(0)
        a, b, c, d, tx, ty = map(float, tm.group(1).replace(",", " ").split())
        sz = float(fs.group(1))
        if ty < top - 2 or ty > keep:        # title band or caption block: delete
            removed += 1; return ""
        new = tier(sz) / k
        r = new / sz
        # glyph x positions: scale about the label centre so centred labels stay centred
        xs_all = [[float(v) for v in xs.split()] for xs in re.findall(r'\bx="([^"]*)"', inner)]
        flat = [v for xs in xs_all for v in xs]
        width = (max(flat) + 0.55 * sz) if flat else 0
        cx = width / 2
        def scale_x(mx):
            vals = [float(v) for v in mx.group(1).split()]
            return 'x="' + " ".join(f"{cx + (v - cx) * r:.4f}" for v in vals) + '"'
        inner2 = re.sub(r'\bx="([^"]*)"', scale_x, inner)
        attrs2 = attrs.replace(f'font-size="{fs.group(1)}"', f'font-size="{new:.4f}"')
        if tier(sz) == 8.0 and len(re.sub(r'<[^>]+>', '', inner).strip()) == 1:
            attrs2 += ' font-weight="bold"'
        resized += 1
        return f'<text{attrs2}>{inner2}</text>'
    # exponents: a smaller text sharing the transform of a "10" base sits as its superscript
    texts = list(TEXT.finditer(s))
    plain = lambda m: re.sub(r'<[^>]+>', '', m.group(2)).strip()
    tf = lambda m: (re.search(r'transform="([^"]*)"', m.group(1)) or [None, None])[1] if re.search(r'transform="([^"]*)"', m.group(1)) else None
    fsz = lambda m: float(re.search(r'font-size="([\d.]+)"', m.group(1)).group(1))
    base_of = {}
    for i, m in enumerate(texts):
        if plain(m) == "10":
            for j in (i + 1, i - 1):
                if 0 <= j < len(texts) and tf(texts[j]) == tf(m) and fsz(texts[j]) < fsz(m):
                    base_of[j] = m
    out, last = [], 0
    for i, m in enumerate(texts):
        out.append(s[last:m.start()]); last = m.end()
        if i in base_of:
            b = base_of[i]
            bsz, esz = fsz(b), fsz(m)
            bnew, enew = tier(bsz) / k, tier(esz) / k
            bxs = [float(v) for v in re.search(r'\bx="([^"]*)"', b.group(2)).group(1).split()]
            adv = 0.556 * bsz
            bw = bxs[-1] + adv; cx = bw / 2; rb = bnew / bsz
            edge_new = cx + (bw - cx) * rb
            ex = [float(v) for v in re.search(r'\bx="([^"]*)"', m.group(2)).group(1).split()]
            ey = float(re.search(r'\by="([^"]*)"', m.group(2)).group(1))
            gap = ex[0] - bw
            nx = [edge_new + gap + (v - ex[0]) * (enew / esz) for v in ex]
            inner = re.sub(r'\bx="[^"]*"', 'x="' + " ".join(f"{v:.4f}" for v in nx) + '"', m.group(2), count=1)
            inner = re.sub(r'\by="[^"]*"', f'y="{ey * rb:.4f}"', inner, count=1)
            attrs = m.group(1).replace(re.search(r'font-size="[\d.]+"', m.group(1)).group(0), f'font-size="{enew:.4f}"')
            ty = float(re.search(r'matrix\(([^)]*)\)', attrs).group(1).replace(",", " ").split()[5])
            if ty < top - 2 or ty > keep: removed += 1; continue
            resized += 1; out.append(f'<text{attrs}>{inner}</text>')
        else:
            out.append(fix(m))
    out.append(s[last:])
    s2 = "".join(out)
    h = keep - top
    s2 = re.sub(r'<svg\b([^>]*?)\swidth="[^"]*"\s+height="[^"]*"\s+viewBox="[^"]*"',
                lambda m: f'<svg{m.group(1)} width="{W_MM}mm" height="{h * W_MM / W_DESIGN:.3f}mm" viewBox="0 {top} {W_DESIGN} {h}"',
                s2, count=1)
    assert f'viewBox="0 {top} ' in s2, name
    s2 = despace(s2)
    out_svg = OUT / f"{name}.svg"; out_svg.write_text(s2)
    subprocess.run(["rsvg-convert", "-f", "pdf", "-o", str(OUT / f"{name}.pdf"), str(out_svg)], check=True)
    subprocess.run(["rsvg-convert", "-f", "png", "--dpi-x", "1000", "--dpi-y", "1000", "-b", "white",
                    "-o", str(OUT / f"{name}_1000.png"), str(out_svg)], check=True)
    report[name] = (removed, resized, round(h * W_MM / W_DESIGN, 1))
for n, (rm, rs, hmm) in report.items(): print(f"{n}: deleted {rm} title/caption texts, restyled {rs}, size 174 x {hmm} mm")
