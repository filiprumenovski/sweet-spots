"""Remove space glyphs from SVG text so cairo/rsvg emits no Type 3 font.
Each space becomes a tspan shifted by one Arial space advance (0.2778 em)."""
import html
import re

SPACE_EM = 569 / 2048          # Arial U+0020 advance

def _fs(attrs, default):
    m = re.search(r'font-size="([\d.]+)"', attrs)
    return float(m.group(1)) if m else default

def _fix_text(m):
    attrs, inner = m.group(1), m.group(2)
    base = _fs(attrs, 16.0)
    # case 1: tspans carrying per-glyph x lists
    def per_glyph(t):
        tattrs, body = t.group(1), t.group(2)
        xm = re.search(r'\bx="([^"]*)"', tattrs)
        chars = list(html.unescape(body))
        xs = xm.group(1).split() if xm else []
        if not xm or len(xs) != len(chars) or "<" in body:
            return t.group(0)
        keep = [i for i, ch in enumerate(chars) if ch != " "]
        if len(keep) == len(chars): return t.group(0)
        new_body = html.escape("".join(chars[i] for i in keep), quote=False)
        new_x = " ".join(xs[i] for i in keep)
        return f'<tspan{tattrs.replace(xm.group(0), f"x={chr(34)}{new_x}{chr(34)}")}>{new_body}</tspan>'
    # case 2: plain text nodes; track the font size of the enclosing element
    out, stack = [], [base]
    carry = 0.0                                  # pending space width carried to the next text run
    for tok in re.split(r'(<[^>]+>)', inner):
        if tok.startswith("</tspan"):
            stack.pop(); out.append(tok)
        elif tok.startswith("<tspan"):
            stack.append(_fs(tok, stack[-1])); out.append(tok)
        elif tok.startswith("<") or tok == "":
            out.append(tok)
        else:
            w = SPACE_EM * stack[-1]
            parts = tok.split(" ")
            buf = ""
            for i, p in enumerate(parts):
                if i > 0: carry += w
                if p == "": continue
                buf += f'<tspan dx="{carry:.3f}">{p}</tspan>' if carry else p
                carry = 0.0
            out.append(buf)
    return f'<text{attrs}>{"".join(out)}</text>'

def despace(svg: str) -> str:
    return re.sub(r'<text\b([^>]*)>(.*?)</text>', _fix_text, svg, flags=re.S)
