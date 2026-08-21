#!/usr/bin/env python3
"""Build capped 3D model compounds for the custom residues, and write mol2.

mol2 rather than SDF or PDB because antechamber must receive *both* the bond orders
and the atom names.  A PDB has names but no bond orders, so antechamber would have to
guess them from geometry -- and the two rings it would have to guess about are a
1,2,4,5-tetrazine and a 4,5-dihydropyridazine, which is exactly where a wrong guess
would go unnoticed.  An SDF has bond orders but drops the names, and the names are
what prepgen and every later atom selection key off.
"""
from __future__ import annotations

import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from residue_defs import RESIDUES, element_of  # noqa: E402

MOL2_ORDER = {1: "1", 2: "2", 3: "3", "ar": "ar"}

# Cap specifications: (name, element, bond order to the atom it hangs off).
# ACE/NME names follow aminont12/aminoct12 so the model reads like a real tripeptide.
ACE_CAP = [("C", "CAY", "C", 1), ("CAY", "OAY", "O", 2), ("CAY", "CAZ", "C", 1),
           ("CAZ", "HZ1", "H", 1), ("CAZ", "HZ2", "H", 1), ("CAZ", "HZ3", "H", 1)]
NME_CAP = [("C", "NNT", "N", 1), ("NNT", "HNT", "H", 1), ("NNT", "CNT", "C", 1),
           ("CNT", "HT1", "H", 1), ("CNT", "HT2", "H", 1), ("CNT", "HT3", "H", 1)]
NMETHYL_CAP = [("C", "NCP", "N", 1), ("NCP", "HCP", "H", 1), ("NCP", "CCP", "C", 1),
               ("CCP", "HC1", "H", 1), ("CCP", "HC2", "H", 1), ("CCP", "HC3", "H", 1)]
ACETYL_CAP = ACE_CAP
PHOSPHATE_CAP = [("C", "P", "P", 1), ("P", "OP1", "O", 2), ("P", "OP2", "O", 1),
                 ("P", "O5C", "O", 1), ("O5C", "C5C", "C", 1),
                 ("C5C", "H5A", "H", 1), ("C5C", "H5B", "H", 1), ("C5C", "H5C", "H", 1)]

# Which caps each residue gets, and to which of its own atoms they attach.
CAP_PLAN = {
    "TET": [("head", ACE_CAP), ("tail", NME_CAP)],
    "TDP": [("head", ACE_CAP), ("tail", NME_CAP), ("extra", NMETHYL_CAP)],
    "DNL": [("head", ACETYL_CAP), ("tail", PHOSPHATE_CAP)],
}
FORMAL_CHARGE = {"TET": 0, "TDP": 0, "DNL": -1}   # DNL's cap carries a phosphodiester


def build(name: str, outdir: Path, seed: int = 0xC0FFEE):
    from rdkit import Chem
    from rdkit.Chem import AllChem

    res = RESIDUES[name]
    bt = {1: Chem.BondType.SINGLE, 2: Chem.BondType.DOUBLE, 3: Chem.BondType.TRIPLE,
          "ar": Chem.BondType.AROMATIC}
    rw, idx, elem = Chem.RWMol(), {}, {}

    def add(nm, el):
        a = Chem.Atom(el)
        a.SetNoImplicit(True)
        idx[nm] = rw.AddAtom(a)
        elem[nm] = el
        return idx[nm]

    for a in res["heavy"]:
        add(a, element_of(a))
    for parent, hs in res["hydrogens"].items():
        for h in hs:
            add(h, "H")
    for a, b, o in res["bonds"]:
        rw.AddBond(idx[a], idx[b], bt[o])
    for parent, hs in res["hydrogens"].items():
        for h in hs:
            rw.AddBond(idx[parent], idx[h], bt[1])

    anchors = {"head": res["head"], "tail": res.get("tail"),
               "extra": res.get("extra_bond", {}).get("this_atom")}
    caps = []
    for where, spec in CAP_PLAN[name]:
        anchor = anchors[where]
        for src, nm, el, order in spec:
            add(nm, el)
            parent = anchor if src == "C" else src
            rw.AddBond(idx[parent], idx[nm], bt[order])
            caps.append(nm)
        # the first cap atom bonds to the residue's anchor, later ones to each other
        # (handled above by treating "C" as "the anchor")

    if name == "DNL":
        rw.GetAtomWithIdx(idx["OP2"]).SetFormalCharge(-1)

    mol = rw.GetMol()
    Chem.SanitizeMol(mol)
    fc = Chem.GetFormalCharge(mol)
    if fc != FORMAL_CHARGE[name]:
        sys.exit(f"{name}: model formal charge {fc}, expected {FORMAL_CHARGE[name]}")

    if AllChem.EmbedMolecule(mol, randomSeed=seed, useRandomCoords=False,
                             maxAttempts=500) != 0:
        if AllChem.EmbedMolecule(mol, randomSeed=seed, useRandomCoords=True,
                                 maxAttempts=1000) != 0:
            sys.exit(f"{name}: 3D embedding failed")
    props = AllChem.MMFFGetMoleculeProperties(mol)
    ff = (AllChem.MMFFGetMoleculeForceField(mol, props) if props
          else AllChem.UFFGetMoleculeForceField(mol))
    ff.Minimize(maxIts=5000)

    names = list(idx)
    inv = {v: k for k, v in idx.items()}
    conf = mol.GetConformer()
    lines = ["@<TRIPOS>MOLECULE", name,
             f"{mol.GetNumAtoms()} {mol.GetNumBonds()} 1 0 0", "SMALL", "USER_CHARGES",
             "", "@<TRIPOS>ATOM"]
    for i in range(mol.GetNumAtoms()):
        p = conf.GetAtomPosition(i)
        nm = inv[i]
        q = mol.GetAtomWithIdx(i).GetFormalCharge()
        lines.append(f"{i + 1:>7} {nm:<8}{p.x:10.4f}{p.y:10.4f}{p.z:10.4f} "
                     f"{elem[nm]:<6}{1:>5} {name:<8}{float(q):>10.4f}")
    lines.append("@<TRIPOS>BOND")
    for k, b in enumerate(mol.GetBonds(), start=1):
        o = {Chem.BondType.SINGLE: "1", Chem.BondType.DOUBLE: "2",
             Chem.BondType.TRIPLE: "3", Chem.BondType.AROMATIC: "ar"}[b.GetBondType()]
        lines.append(f"{k:>6}{b.GetBeginAtomIdx() + 1:>6}{b.GetEndAtomIdx() + 1:>6} {o}")
    out = outdir / f"{name}_model.mol2"
    out.write_text("\n".join(lines) + "\n")

    smi = Chem.MolToSmiles(Chem.RemoveHs(Chem.Mol(mol)))
    print(f"  {name}: {mol.GetNumAtoms()} atoms ({len(caps)} in caps), "
          f"formal charge {fc}")
    print(f"    SMILES (from the tables, heavy atoms): {smi}")
    return out, caps, fc, smi


if __name__ == "__main__":
    outdir = Path(sys.argv[1]) if len(sys.argv) > 1 else HERE
    for n in ("TET", "TDP", "DNL"):
        d = outdir / n
        d.mkdir(parents=True, exist_ok=True)
        build(n, d)
