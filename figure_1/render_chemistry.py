"""2D chemical structures for the OGT reaction in Figure 1A (RDKit).

Acceptor: an L-serine residue within a peptide (wavy bonds = chain continues).
Product:  the same residue carrying beta-D-GlcNAc on its side-chain oxygen.
Stereocentres are checked against their CIP labels before drawing:
  Ser C-alpha (L) = S;  GlcNAc ring C1..C5 = R, R, R, S, R  (beta-D-gluco).
Writes assets/chem_acceptor.svg and assets/chem_product.svg.

Run with an RDKit-enabled Python, e.g.
  ~/Documents/GitHub/ogt-crackdown/.venv/bin/python render_chemistry.py
"""
from pathlib import Path

from rdkit import Chem
from rdkit.Chem import rdCIPLabeler, rdDepictor
from rdkit.Chem.Draw import rdMolDraw2D

HERE = Path(__file__).resolve().parent
ARIAL = "/System/Library/Fonts/Supplemental/Arial.ttf"

ACCEPTOR = "*N[C@@H](CO)C(*)=O"
# beta-D-GlcNAc (PubChem isomeric SMILES) with the anomeric OH replaced by Ser O-gamma
PRODUCT = "*N[C@@H](CO[C@@H]1O[C@H](CO)[C@@H](O)[C@H](O)[C@H]1NC(C)=O)C(*)=O"


def cip(mol):
    rdCIPLabeler.AssignCIPLabels(mol)
    return {a.GetIdx(): a.GetProp("_CIPCode") for a in mol.GetAtoms() if a.HasProp("_CIPCode")}


def check():
    acc = Chem.MolFromSmiles(ACCEPTOR)
    assert list(cip(acc).values()) == ["S"], cip(acc)
    prod = Chem.MolFromSmiles(PRODUCT)
    labels = cip(prod)
    ser_ca = 2
    assert labels[ser_ca] == "S", labels
    ring = prod.GetRingInfo().AtomRings()[0]
    ring_o = next(i for i in ring if prod.GetAtomWithIdx(i).GetSymbol() == "O")
    # walk ring from anomeric carbon (bonded to ring O and to Ser O-gamma)
    anomeric = next(i for i in ring if i != ring_o
                    and prod.GetBondBetweenAtoms(i, ring_o) is not None
                    and any(n.GetIdx() == 4 for n in prod.GetAtomWithIdx(i).GetNeighbors()))
    order = [anomeric]
    prev = ring_o
    while len(order) < 5:
        cur = order[-1]
        nxt = next(n.GetIdx() for n in prod.GetAtomWithIdx(cur).GetNeighbors()
                   if n.GetIdx() in ring and n.GetIdx() not in (prev, ring_o) and n.GetIdx() not in order)
        prev = cur
        order.append(nxt)
    got = [labels[i] for i in order]
    assert got == ["R", "R", "R", "S", "R"], got
    return prod, set(order) | {ring_o}


def orient(mol, n=1, ca=2, cb=3, c=None):
    """Rotate 2D coords: backbone N->C along +x, side chain pointing down."""
    import math
    from rdkit.Geometry import Point3D
    conf = mol.GetConformer()
    P = lambda i: conf.GetAtomPosition(i)
    vx, vy = P(c).x - P(n).x, P(c).y - P(n).y
    ang = -math.atan2(vy, vx)
    ca_, sa_ = math.cos(ang), math.sin(ang)
    for i in range(mol.GetNumAtoms()):
        p = P(i)
        conf.SetAtomPosition(i, Point3D(p.x * ca_ - p.y * sa_, p.x * sa_ + p.y * ca_, 0))
    if P(cb).y > P(ca).y:  # RDKit y points up; we want the side chain below
        for i in range(mol.GetNumAtoms()):
            p = P(i)
            conf.SetAtomPosition(i, Point3D(p.x, -p.y, 0))


def draw(mol, path):
    d = rdMolDraw2D.MolDraw2DSVG(-1, -1)
    o = d.drawOptions()
    o.fontFile = ARIAL
    o.bondLineWidth = 2
    o.fixedBondLength = 34
    o.padding = 0.04
    o.clearBackground = False
    o.dummiesAreAttachments = True
    o.additionalAtomLabelPadding = 0.08
    o.minFontSize = 15
    o.maxFontSize = 18
    o.setAtomPalette({-1: (0.08, 0.17, 0.22), 0: (0.45, 0.52, 0.56), 6: (0.08, 0.17, 0.22),
                      7: (0.18, 0.36, 0.62), 8: (0.75, 0.31, 0.30)})
    Chem.WedgeMolBonds(mol, mol.GetConformer())
    d.DrawMolecule(mol)
    d.FinishDrawing()
    path.write_text(d.GetDrawingText())


def main():
    prod, ring_atoms = check()
    rdDepictor.SetPreferCoordGen(True)
    rdDepictor.Compute2DCoords(prod)
    orient(prod, n=1, ca=2, cb=3, c=19)
    # acceptor reuses the product's serine coordinates so both drawings match
    acc = Chem.MolFromSmiles(ACCEPTOR)  # 0 * 1 N 2 CA 3 CB 4 OG 5 C 6 * 7 O
    rdDepictor.Compute2DCoords(acc)
    from rdkit.Geometry import Point3D
    pc, ac = prod.GetConformer(), acc.GetConformer()
    for ai, pi in zip(range(8), (0, 1, 2, 3, 4, 19, 20, 21)):
        p = pc.GetAtomPosition(pi)
        ac.SetAtomPosition(ai, Point3D(p.x, p.y, 0))
    out = HERE / "assets"
    out.mkdir(exist_ok=True)
    draw(acc, out / "chem_acceptor.svg")
    draw(prod, out / "chem_product.svg")
    print("stereo checks passed; wrote", out / "chem_acceptor.svg", out / "chem_product.svg")


if __name__ == "__main__":
    main()
