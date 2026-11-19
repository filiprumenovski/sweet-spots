"""Compartment dot map of OGT substrate proteins (Figure 1B).

One circle per human protein in the O-GlcNAc Database (OVSlab, Wulff-Fuentes et al.;
$OVSLAB_CSV, Homo sapiens, isoform suffixes stripped).
Site counts are validated O-GlcNAc sites (12,491 sites on 3,532 proteins; same filter
as manuscript Methods) on those proteins. Circle area
is proportional to the validated site count; proteins without a validated site get a
1-site-sized circle flagged sites=0 and are packed last, forming the outer halo. Proteins are grouped by UniProt curated
subcellular location (data/uniprot_human_subcellular.tsv, reviewed human entries;
query in data/uniprot_human_subcellular.query.txt). Writes data/substrate_map.json.
"""
import csv, json, math, os, random, re, sys
from collections import Counter, defaultdict
from pathlib import Path
import pandas as pd

HERE = Path(__file__).resolve().parent
INPUTS = HERE.parent / "inputs"   # frozen repository inputs (fetched by ./reproduce)
fa = pd.read_parquet(INPUTS / "data/interim/fasta_human.parquet")
seqs = dict(zip(fa.accession, fa.sequence))
gene = dict(zip(fa.accession, fa.gene))

sites = set()
for r in csv.DictReader(open(INPUTS / "analysis/revalidation/data/atlas_unambiguous.csv", encoding="latin-1")):
    if r["species"].strip().lower() != "human" or r["site_residue"] not in ("S", "T"):
        continue
    s = seqs.get(r["accession"])
    try:
        p = int(float(r["position_in_protein"]))
    except ValueError:
        continue
    if s and 1 <= p <= len(s) and s[p - 1] == r["site_residue"]:
        sites.add((r["accession"], p))
n_sites = Counter(a for a, _ in sites)
print("sites", len(sites), "proteins", len(n_sites), file=sys.stderr)
_OVSLAB_CSV = os.environ.get("OVSLAB_CSV")
if not _OVSLAB_CSV:
    sys.exit("set OVSLAB_CSV to the path of the O-GlcNAc Database (OVSlab) all-species export (provider-restricted; see figure_1/README.md)")
OVS = Path(_OVSLAB_CSV)
ovs = pd.read_csv(OVS, encoding="latin-1", low_memory=False)
ovs_h = set(ovs.loc[ovs.organism == "Homo sapiens", "UniprotKB ID"].astype(str).str.strip().str.split("-").str[0])
universe = {a: n_sites.get(a, 0) for a in ovs_h}   # O-GlcNAc Database proteins only
print("OVSlab human", len(ovs_h), "universe", len(universe), file=sys.stderr)

loc = {}
with open(HERE / "data/uniprot_human_subcellular.tsv") as f:
    next(f)
    for line in f:
        acc, _, txt = line.rstrip("\n").partition("\t")
        loc[acc] = txt

def flags(txt: str) -> set[str]:
    t = re.sub(r"\{[^}]*\}", "", txt)          # evidence tags
    t = re.sub(r"Note=.*?(?=SUBCELLULAR LOCATION:|$)", "", t)  # free-text notes
    t = t.lower()
    f = set()
    if re.search(r"nucle(us|olus|oplasm|ar pore|ar speck|ar body|ar envelope)|chromosome|centromere|kinetochore|telomere|pml body", t):
        f.add("N")
    if re.search(r"cytoplasm|cytosol|cytoskeleton|centrosome|centriole|spindle|p-body|stress granule|microtubule|actin|myofibril|sarcomere|z line|cilium|cell projection|axon|dendrite|synapse|lamellipodium|filopodium|focal adhesion|cell cortex|perinuclear", t):
        f.add("C")
    if "mitochondri" in t:
        f.add("M")
    if re.search(r"endoplasmic reticulum|golgi|membrane|secreted|extracellular|endosome|lysosome|vesicle|vacuole|peroxisome|cell junction|cell surface", t):
        f.add("S")
    return f

def compartment(acc: str) -> str:
    f = flags(loc.get(acc, ""))
    if "N" in f and "C" in f:
        return "nucleocytoplasmic"
    if "N" in f:
        return "nucleus"
    if "M" in f:
        return "mitochondria"
    if "C" in f:
        return "cytoplasm"
    if "S" in f:
        return "membranes & secretory"
    return "not annotated"

groups = defaultdict(list)
for a, n in universe.items():
    groups[compartment(a)].append((a, n))
for g, v in sorted(groups.items(), key=lambda kv: -len(kv[1])):
    print(f"{g:24s} proteins {len(v):5d}  mapped {sum(n > 0 for _, n in v):5d}  sites {sum(n for _, n in v):6d}", file=sys.stderr)

# circle packing per group: radius = R0 * sqrt(sites); greedy placement on an
# Archimedean spiral, largest first, resuming the spiral near the last hit
R0, GAP = 1.45, 0.55
RMAX = R0 * math.sqrt(212) + GAP
CELL = 8.0

def pack(items):
    items = sorted(items, key=lambda t: -t[1])
    grid = defaultdict(list)
    def ok(x, y, r):
        gx, gy = int(x // CELL), int(y // CELL)
        k = int(math.ceil((r + RMAX) / CELL))
        for i in range(gx - k, gx + k + 1):
            for j in range(gy - k, gy + k + 1):
                for (px, py, pr) in grid.get((i, j), ()):
                    if (px - x) ** 2 + (py - y) ** 2 < (pr + r + GAP) ** 2:
                        return False
        return True
    out, t_last = [], 0.0
    for acc, n in items:
        r = R0 * math.sqrt(max(n, 1))
        t = max(0.0, t_last - 40.0)
        while True:
            rad = 0.55 * math.sqrt(t * 6)
            x, y = rad * math.cos(t), rad * math.sin(t)
            if ok(x, y, r):
                break
            t += 0.12
        t_last = max(t_last, t)
        grid[(int(x // CELL), int(y // CELL))].append((x, y, r))
        out.append({"acc": acc, "gene": gene.get(acc) or acc, "sites": n,
                    "x": round(x, 2), "y": round(y, 2), "r": round(r, 3)})
    ext = max(math.hypot(c["x"], c["y"]) + c["r"] for c in out)
    return out, ext

res = {"n_sites": sum(universe.values()), "n_proteins": len(universe),
       "n_mapped": sum(n > 0 for n in universe.values()),
       "n_ovslab_human": len(ovs_h), "radius_rule": f"r = {R0}*sqrt(sites)",
       "groups": {}}
for g, v in groups.items():
    circles, ext = pack(v)
    res["groups"][g] = {"n_proteins": len(v), "n_mapped": sum(n > 0 for _, n in v), "n_sites": sum(n for _, n in v),
                        "extent": round(ext, 1), "circles": circles}
    print(f"{g:24s} extent {ext:.1f}", file=sys.stderr)
(HERE / "data/substrate_map.json").write_text(json.dumps(res))
