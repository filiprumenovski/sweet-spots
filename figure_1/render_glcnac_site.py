"""Model one O-GlcNAc site of FOXK1 at atomic scale for the Figure 1B inset.

beta-GlcNAc is taken from PDB 5UN8 chain E (NAG 200 on Ser149) and placed on the
FOXK1 acceptor by fitting the Ser149 N/CA/CB/OG atoms onto the acceptor's
N/CA/CB/OG (Ser) or N/CA/CB/OG1 (Thr). Residues site-4..site+4 of the AlphaFold
model (AF-P85037 v6) are drawn as sticks. The glycan position is a model.

Run from this folder: pymol -cq render_glcnac_site.py [-- SITE]
Writes assets/FOXK1_<Ser|Thr><SITE>_glcnac.png and a .txt with the fit RMSD.
"""
import sys
from pathlib import Path

from pymol import cmd

root = Path.cwd()
SITE = int(sys.argv[1]) if len(sys.argv) > 1 else 683
MODEL = root / "structures/AF-P85037-F1-model_v6.pdb"

cmd.reinitialize()
cmd.bg_color("white")
for k, v in {"ray_opaque_background": 0, "orthoscopic": 1, "antialias": 2,
             "ray_trace_mode": 1, "ray_trace_gain": 0.02, "ray_shadows": 0,
             "ambient": 0.45, "specular": 0.15, "stick_radius": 0.24}.items():
    cmd.set(k, v)
cmd.set_color("peptide_coral", [0.88, 0.36, 0.29])
cmd.set_color("glycan_blue", [0.12, 0.34, 0.70])

cmd.load(str(MODEL), "prot")
cmd.load(str(root / "structures/5UN8.pdb"), "g")
cmd.remove("g and not ((chain E and resi 149) or (chain E and resn NAG))")
resn = cmd.get_model(f"prot and resi {SITE} and name CA").atom[0].resn
oxygen = {"SER": "OG", "THR": "OG1"}[resn]
rms = cmd.pair_fit("g and resi 149 and name N", f"prot and resi {SITE} and name N",
                   "g and resi 149 and name CA", f"prot and resi {SITE} and name CA",
                   "g and resi 149 and name CB", f"prot and resi {SITE} and name CB",
                   "g and resi 149 and name OG", f"prot and resi {SITE} and name {oxygen}")
seg = f"prot and resi {SITE - 4}-{SITE + 4} and not hydro"
cmd.create("site", f"({seg}) or (g and resn NAG and not hydro)")  # one object so the bond can be drawn
cmd.bond(f"site and resi {SITE} and name {oxygen}", "site and resn NAG and name C1")
cmd.delete("prot")
cmd.delete("g")
cmd.hide("everything")
cmd.show("sticks", "site")
cmd.color("peptide_coral", "site and not resn NAG")
cmd.color("glycan_blue", "site and resn NAG")
cmd.orient("site")
cmd.zoom("site", 1.5)
cmd.ray(1100, 900)
tag = f"{resn.capitalize()}{SITE}"
cmd.png(str(root / f"assets/FOXK1_{tag}_glcnac.png"), dpi=300)
(root / f"assets/FOXK1_{tag}_graft.txt").write_text(
    f"GlcNAc of PDB 5UN8 chain E (NAG 200 on Ser149) fitted onto FOXK1 {tag} (AF-P85037 v6) "
    f"by N/CA/CB/{oxygen}; RMSD {rms:.2f} A\n")
print("pair_fit", tag, "RMSD", round(rms, 3))
