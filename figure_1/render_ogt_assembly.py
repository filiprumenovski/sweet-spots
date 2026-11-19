"""Render the full-length human OGT dimer (PDB 7YEA) for Figure 1A.

7YEA (cryo-EM, 3.82 A) covers residues 23-1039 of both protomers but has no
ligands. To mark the active site, PDB 3PE4 (OGT 315-1028 with UDP and a
CKII-derived acceptor peptide) is superposed on the catalytic region of each
protomer, and only its UDP and peptide are drawn. Ligand placement is by
superposition, not observed in 7YEA.

Run from this folder: pymol -cq render_ogt_assembly.py
"""
from pathlib import Path
from pymol import cmd

root = Path.cwd()  # PyMOL rebinds __file__, so run from this folder
st = root / "structures"
out = root / "assets"
out.mkdir(exist_ok=True)

TPR_END = 465  # TPR superhelix ends near residue 465; catalytic domains follow

cmd.reinitialize()
cmd.bg_color("white")
for k, v in {
    "ray_opaque_background": 0, "orthoscopic": 1, "antialias": 2,
    "ray_trace_mode": 1, "ray_trace_gain": 0.02, "ray_shadows": 0,
    "ambient": 0.45, "specular": 0.15, "shininess": 18,
    "cartoon_fancy_helices": 1, "cartoon_flat_sheets": 1,
    "cartoon_loop_radius": 0.28, "cartoon_sampling": 14,
    "cartoon_transparency": 0, "stick_radius": 0.45,
    "sphere_scale": 0.9,
}.items():
    cmd.set(k, v)

cmd.set_color("catA", [0.13, 0.48, 0.53])
cmd.set_color("tprA", [0.47, 0.72, 0.72])
cmd.set_color("catB", [0.31, 0.45, 0.58])
cmd.set_color("tprB", [0.63, 0.74, 0.82])
cmd.set_color("peptide_coral", [0.88, 0.36, 0.29])
cmd.set_color("ligand_gold", [0.90, 0.64, 0.16])

cmd.load(str(st / "7YEA.cif"), "ogt")
cmd.remove("ogt and (solvent or hydro)")
cmd.hide("everything")
cmd.show("cartoon", "ogt")
cmd.color("catA", f"ogt and chain A and resi {TPR_END + 1}-1100")
cmd.color("tprA", f"ogt and chain A and resi 1-{TPR_END}")
cmd.color("catB", f"ogt and chain B and resi {TPR_END + 1}-1100")
cmd.color("tprB", f"ogt and chain B and resi 1-{TPR_END}")

rms = {}
for ch in ("A", "B"):
    name = f"lig{ch}"
    cmd.load(str(st / "3PE4.cif"), name)
    cmd.remove(f"{name} and not (chain A+B)")  # enzyme A, peptide B, UDP on A
    cmd.remove(f"{name} and (solvent or resn SO4 or hydro)")
    r = cmd.super(f"{name} and chain A and name CA and resi {TPR_END + 1}-1100",
                  f"ogt and chain {ch} and name CA and resi {TPR_END + 1}-1100")
    rms[ch] = (round(r[0], 2), r[1])
    cmd.hide("everything", name)
    cmd.show("sticks", f"{name} and (chain B or resn UDP)")
    cmd.show("spheres", f"{name} and resn UDP")
    cmd.color("peptide_coral", f"{name} and chain B")
    cmd.color("ligand_gold", f"{name} and resn UDP")
    cmd.util.cnc(f"{name} and chain B")
    cmd.color("peptide_coral", f"{name} and chain B and elem C")

print("superposition RMSD (A), atoms:", rms)

# Camera raised 40 degrees above the front (scissor) view, no sideways turn
# (option 3 of view_options.png). Pinned as an explicit matrix because
# cmd.orient() lands on a different handedness depending on the prior camera.
cmd.set_view((0.9998, -0.0122, 0.0147, -0.0191, -0.6425, 0.766, 0.0001, -0.7662, -0.6427,
              0.0, -0.0, -543.9324, 172.195, 172.1965, 181.4398, 428.8405, 659.0244, 20.0))
cmd.set("fog_start", 0.55)
cmd.ray(2800, 1800)
cmd.png(str(out / "OGT_7YEA_dimer.png"), dpi=300)
cmd.save(str(out / "OGT_7YEA_dimer.pse"))
(root / "assets" / "OGT_7YEA_superposition.txt").write_text(
    "3PE4 chain A superposed on 7YEA catalytic region (CA, resi >465)\n"
    + "\n".join(f"protomer {k}: RMSD {v[0]} A over {v[1]} atoms" for k, v in rms.items())
    + "\n")
