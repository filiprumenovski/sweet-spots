"""Render the full FOXK1 AlphaFold model (AF-P85037-F1 v6) for Figure 1B.

Cartoon of all 733 residues; consensus O-GlcNAc regions in navy. A second, flat-colour
pass draws one sphere per O-GlcNAc site with the same camera so the figure builder can
place site markers at their projected positions.

Outputs (assets/): FOXK1_model.png, FOXK1_model_sites2d.json
Run from this folder: pymol -cq render_foxk1_model.py
"""
import json
from pathlib import Path

import numpy as np
from PIL import Image
from pymol import cmd

root = Path.cwd()
out = root / "assets"
track = json.loads((root / "data/foxk1_track.json").read_text())
W = H = 1600

cmd.reinitialize()
cmd.bg_color("white")
for k, v in {"ray_opaque_background": 0, "orthoscopic": 1, "antialias": 2,
             "ray_trace_mode": 1, "ray_trace_gain": 0.02, "ray_shadows": 0,
             "ambient": 0.45, "specular": 0.15, "shininess": 18,
             "cartoon_fancy_helices": 1, "cartoon_flat_sheets": 1,
             "cartoon_loop_radius": 0.45, "cartoon_smooth_loops": 1,
             "cartoon_sampling": 14, "fog": 0}.items():
    cmd.set(k, v)
cmd.set_color("prot_teal", [0.33, 0.71, 0.71])
cmd.set_color("prot_helix", [0.16, 0.55, 0.58])
cmd.set_color("region_navy", [0.11, 0.37, 0.54])

cmd.load(str(root / "structures/AF-P85037-F1-model_v6.pdb"), "prot")
cmd.hide("everything")
cmd.show("cartoon", "prot")
cmd.color("prot_teal", "prot")
cmd.color("prot_helix", "prot and ss H")
reg_sel = " or ".join(f"resi {a}-{b}" for a, b in track["regions"])
cmd.color("region_navy", f"prot and ({reg_sel})")
cmd.orient("prot and b > 70")      # face the folded domains
cmd.zoom("prot", 0)
view = cmd.get_view()
cmd.ray(W, H)
cmd.png(str(out / "FOXK1_model.png"), dpi=300)

# marker pass: flat colours, same camera
cmd.hide("everything")
for k, v in {"ray_trace_mode": 0, "antialias": 0, "ambient": 1.0, "direct": 0.0,
             "reflect": 0.0, "specular": 0.0, "depth_cue": 0, "sphere_scale": 1.6}.items():
    cmd.set(k, v)
codes = {}
for i, p in enumerate(track["oglcnac"]):
    code = 8 * (i + 1)
    cmd.set_color(f"site{p}", [code / 255.0, 0.0, 1.0])
    cmd.show("spheres", f"prot and resi {p} and name CA")
    cmd.color(f"site{p}", f"prot and resi {p} and name CA")
    codes[p] = code
cmd.set_view(view)
cmd.ray(W, H)
tmp = out / "_markers.png"
cmd.png(str(tmp), dpi=300)

main = Image.open(out / "FOXK1_model.png")
bb = main.split()[3].getbbox()
m = np.array(Image.open(tmp).convert("RGBA")).astype(int)
r, g, b, a = m[..., 0], m[..., 1], m[..., 2], m[..., 3]
sites = {}
for p, code in codes.items():
    ys, xs = np.nonzero((a > 200) & (b > 150) & (g < 60) & (r == code))
    sites[str(p)] = [round(float((xs.mean() - bb[0]) / (bb[2] - bb[0])), 4),
                     round(float((ys.mean() - bb[1]) / (bb[3] - bb[1])), 4)]
tmp.unlink()
(out / "FOXK1_model_sites2d.json").write_text(json.dumps({"crop_bbox": list(bb), "sites": sites}))
print("sites placed", len(sites))
