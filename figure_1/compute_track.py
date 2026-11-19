"""Residue track for the panel B protein, Figure 1B (usage: python3 compute_track.py ACC GENE).

O-GlcNAc sites: validated atlas sites (same filter as compute_substrate_map.py).
Phosphosites: PhosphoSitePlus Phosphorylation_site_dataset (human, pS/pT),
validated against the canonical sequence. Regions: consensus O-GlcNAc regions
(oglcnac_consensus_regions.parquet). Writes data/<gene>_track.json.
"""
import csv, gzip, io, json, os, sys
from pathlib import Path
import pandas as pd

ACC, GENE = (sys.argv[1], sys.argv[2]) if len(sys.argv) > 2 else ("P85037", "FOXK1")
HERE = Path(__file__).resolve().parent
INPUTS = HERE.parent / "inputs"   # frozen repository inputs (fetched by ./reproduce)
if not os.environ.get("PSP_SITES"):
    sys.exit("set PSP_SITES to PhosphoSitePlus Phosphorylation_site_dataset.gz (provider-restricted; see figure_1/README.md)")
PSP = Path(os.environ["PSP_SITES"])

fa = pd.read_parquet(INPUTS / "data/interim/fasta_human.parquet")
seq = fa.loc[fa.accession == ACC, "sequence"].iloc[0]

og = set()
for r in csv.DictReader(open(INPUTS / "analysis/revalidation/data/atlas_unambiguous.csv", encoding="latin-1")):
    if r["accession"] != ACC or r["species"].strip().lower() != "human" or r["site_residue"] not in ("S", "T"):
        continue
    p = int(float(r["position_in_protein"]))
    if seq[p - 1] == r["site_residue"]:
        og.add(p)

raw = gzip.open(PSP, "rt", encoding="latin-1").read().splitlines()
start = next(i for i, l in enumerate(raw) if l.startswith("GENE") or l.startswith("PROTEIN"))
ph = set()
for r in csv.DictReader(io.StringIO("\n".join(raw[start:])), delimiter="\t"):
    if r.get("ACC_ID") != ACC or r.get("ORGANISM") != "human":
        continue
    m = r["MOD_RSD"].split("-")[0]
    aa, pos = m[0], int(m[1:])
    if aa in "ST" and pos <= len(seq) and seq[pos - 1] == aa:
        ph.add(pos)

reg = pd.read_parquet(INPUTS / "data/processed/regions/oglcnac_consensus_regions.parquet")
reg = reg[reg.accession == ACC].sort_values("start")
regions = [[int(a), int(b)] for a, b in zip(reg.start, reg.end)]
st = [i + 1 for i, c in enumerate(seq) if c in "ST"]
out = {"accession": ACC, "gene": GENE, "length": len(seq), "st": st,
       "oglcnac": sorted(og), "phospho": sorted(ph), "regions": regions}
(HERE / f"data/{GENE.lower()}_track.json").write_text(json.dumps(out))
print(f"length {len(seq)}  S/T {len(st)}  O-GlcNAc {len(og)}  phospho {len(ph)}  regions {regions}", file=sys.stderr)
