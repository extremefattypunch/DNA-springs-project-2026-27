#!/usr/bin/env python3
"""The covalent tether from the protein backbone to the DNA phosphate, as one molecule.

TDP + DNL are built together, embedded and MMFF-optimised as a single molecule, then
posed by rigid-body superposition onto the target residue's backbone followed by
rotations about the tether's genuinely rotatable bonds.

Why not let tleap build it
--------------------------
An Amber prep file stores a *tree* of internal coordinates, so tleap builds a
residue's missing atoms by walking that tree -- and ring-closing bonds simply end up
wherever the walk leaves them.  Asked to build TDP from backbone + CB, tleap returned
a cyclopropane with sides of 5.28, 4.70 and 1.52 A instead of three of 1.51: the
bicyclo[6.1.0]nonane cage was torn open.  The cage is the mechanically load-bearing
part of this tether, so its geometry has to come from a real 3D optimisation.

Rotatable bonds
---------------
Found from the bond graph rather than listed by hand: a bond is rotatable if it is
single, not in a ring, and both ends carry at least one other heavy atom.  For this
tether that yields the three side-chain torsions, the carbamate's three, the amide,
and the six of the hexyl arm -- the fused ring system correctly contributes none.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "02_params"))
from residue_defs import RESIDUES, element_of  # noqa: E402

BACKBONE = ("N", "CA", "C", "O", "H", "HA")


def kabsch(P, Q):
    """Rotation+translation taking P onto Q (both n x 3)."""
    pc, qc = P.mean(axis=0), Q.mean(axis=0)
    H = (P - pc).T @ (Q - qc)
    U, S, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    R = Vt.T @ np.diag([1.0, 1.0, d]) @ U.T
    return R, qc - R @ pc


def rodrigues(axis, angle):
    a = axis / np.linalg.norm(axis)
    K = np.array([[0, -a[2], a[1]], [a[2], 0, -a[0]], [-a[1], a[0], 0]])
    return np.eye(3) + np.sin(angle) * K + (1 - np.cos(angle)) * (K @ K)


def dihedral(p0, p1, p2, p3):
    b1 = p2 - p1
    b1 = b1 / np.linalg.norm(b1)
    v = (p1 - p0) - ((p1 - p0) @ b1) * b1
    w = (p3 - p2) - ((p3 - p2) @ b1) * b1
    return np.arctan2(np.cross(b1, v) @ w, v @ w)


class Arm:
    """TDP side chain + DNL, with coordinates and a torsion-rotation interface."""

    def __init__(self, names, coords, bonds, elements):
        self.names = list(names)
        self.index = {n: i for i, n in enumerate(self.names)}
        self.xyz = np.asarray(coords, float)
        self.bonds = [(self.index[a], self.index[b]) for a, b in bonds]
        self.elements = list(elements)
        self._adj = {i: set() for i in range(len(self.names))}
        for a, b in self.bonds:
            self._adj[a].add(b)
            self._adj[b].add(a)

    def moving_set(self, i, j):
        """Atoms on j's side when bond i-j is cut, or None if the bond is in a ring.

        The traversal must be allowed to *reach* i by another route -- that is the
        whole ring test.  Skipping i explicitly, as an earlier version did, makes the
        test unfalsifiable: every aromatic and cyclopropane bond then reports as
        rotatable with the entire molecule downstream of it.
        """
        seen, stack = {j}, [j]
        while stack:
            x = stack.pop()
            for y in self._adj[x]:
                if (x, y) in ((i, j), (j, i)):
                    continue                  # the cut bond itself
                if y in seen:
                    continue
                if y == i:
                    return None               # reached i another way: it is a ring
                seen.add(y)
                stack.append(y)
        return sorted(seen)

    def rotatable(self, exclude_backbone=True):
        """[(i, j, moving atoms)] for the tether's rotatable single bonds.

        Backbone bonds are excluded (except CA-CB, which is chi1): the backbone
        belongs to the protein and is not ours to turn.  A bond whose distal side
        carries no heavy atom is also excluded -- rotating a terminal hydroxyl or
        methyl changes nothing about where the tether reaches.
        """
        bb = {self.index[n] for n in BACKBONE if n in self.index}
        out = []
        for a, b in self.bonds:
            if self.elements[a] == "H" or self.elements[b] == "H":
                continue
            if exclude_backbone and a in bb and b in bb:
                continue
            heavy_a = [x for x in self._adj[a] if self.elements[x] != "H" and x != b]
            heavy_b = [x for x in self._adj[b] if self.elements[x] != "H" and x != a]
            if not heavy_a or not heavy_b:
                continue
            mov = self.moving_set(a, b)
            if mov is None:
                continue
            if not any(self.elements[m] != "H" for m in mov):
                continue
            # always rotate the side that does *not* contain the backbone
            if any(m in bb for m in mov):
                a, b = b, a
                mov = self.moving_set(a, b)
                if mov is None or any(m in bb for m in mov):
                    continue
            out.append((a, b, mov))
        return out

    def rotate(self, i, j, mov, delta):
        R = rodrigues(self.xyz[j] - self.xyz[i], delta)
        pivot = self.xyz[j]
        self.xyz[mov] = (R @ (self.xyz[mov] - pivot).T).T + pivot

    def place_on_backbone(self, target: dict):
        """Superimpose (N, CA, CB) onto the target residue's own N, CA, CB."""
        src = np.array([self.xyz[self.index[n]] for n in ("N", "CA", "CB")])
        dst = np.array([target[n] for n in ("N", "CA", "CB")])
        R, t = kabsch(src, dst)
        self.xyz = (R @ self.xyz.T).T + t
        rms = np.sqrt(((np.array([self.xyz[self.index[n]] for n in ("N", "CA", "CB")])
                        - dst) ** 2).sum(axis=1).mean())
        return rms


def build_arm(dna_residue: str = "DNL") -> Arm:
    """TDP side chain + the DNA-side arm as one MMFF-optimised molecule.

    ``dna_residue`` is DNL for a chimera (the arm continues into the duplex) or DNH
    for the clicked, unloaded control (the arm ends in a free alcohol).
    """
    from rdkit import Chem
    from rdkit.Chem import AllChem

    bt = {1: Chem.BondType.SINGLE, 2: Chem.BondType.DOUBLE, 3: Chem.BondType.TRIPLE,
          "ar": Chem.BondType.AROMATIC}
    tdp, dnl = RESIDUES["TDP"], RESIDUES[dna_residue]
    rw, idx, names, elems = Chem.RWMol(), {}, [], []

    def add(nm, el, tag):
        a = Chem.Atom(el)
        a.SetNoImplicit(True)
        idx[(tag, nm)] = rw.AddAtom(a)
        names.append(nm if tag == "TDP" else f"{nm}*")
        elems.append(el)
        return idx[(tag, nm)]

    for res, tag in ((tdp, "TDP"), (dnl, "DNL")):
        for a in res["heavy"]:
            add(a, element_of(a), tag)
        for parent, hs in res["hydrogens"].items():
            for h in hs:
                add(h, "H", tag)
    bonds_named = []
    for res, tag in ((tdp, "TDP"), (dnl, "DNL")):
        for a, b, o in res["bonds"]:
            rw.AddBond(idx[(tag, a)], idx[(tag, b)], bt[o])
            bonds_named.append(((tag, a), (tag, b)))
        for parent, hs in res["hydrogens"].items():
            for h in hs:
                rw.AddBond(idx[(tag, parent)], idx[(tag, h)], bt[1])
                bonds_named.append(((tag, parent), (tag, h)))
    # the carbamate C-N bond joining the two residues
    rw.AddBond(idx[("TDP", "CN")], idx[("DNL", "N")], bt[1])
    bonds_named.append((("TDP", "CN"), ("DNL", "N")))

    # cap the free valences so RDKit can sanitise: an H on the backbone N and C=O
    # already exist; the backbone C needs one more substituent, and DNL's OL one.
    cap_h = []
    tail_caps = [("TDP", "C")]
    if dnl.get("tail"):
        tail_caps.append(("DNL", dnl["tail"]))
    for tag, nm in tail_caps:
        h = add(f"HCAP{len(cap_h)}", "H", tag)
        rw.AddBond(idx[(tag, nm)], h, bt[1])
        bonds_named.append(((tag, nm), (tag, f"HCAP{len(cap_h)}")))
        cap_h.append(h)

    mol = rw.GetMol()
    Chem.SanitizeMol(mol)
    if AllChem.EmbedMolecule(mol, randomSeed=0xBEEF, maxAttempts=1000) != 0:
        if AllChem.EmbedMolecule(mol, randomSeed=0xBEEF, useRandomCoords=True,
                                 maxAttempts=2000) != 0:
            sys.exit("arm: 3D embedding failed")
    props = AllChem.MMFFGetMoleculeProperties(mol)
    ff = (AllChem.MMFFGetMoleculeForceField(mol, props) if props
          else AllChem.UFFGetMoleculeForceField(mol))
    ff.Minimize(maxIts=20000)
    conf = mol.GetConformer()
    xyz = np.array([[conf.GetAtomPosition(i).x, conf.GetAtomPosition(i).y,
                     conf.GetAtomPosition(i).z] for i in range(mol.GetNumAtoms())])

    # drop the two capping hydrogens
    keep = [i for i in range(mol.GetNumAtoms()) if i not in cap_h]
    remap = {old: new for new, old in enumerate(keep)}
    names2 = [names[i] for i in keep]
    elems2 = [elems[i] for i in keep]
    bonds2 = [(names[a], names[b]) for (ta, na), (tb, nb) in bonds_named
              for a, b in [(idx[(ta, na)], idx[(tb, nb)])]
              if a in remap and b in remap]
    return Arm(names2, xyz[keep], bonds2, elems2)


def build_residue_mol(name: str) -> Arm:
    """A single custom residue as an MMFF-optimised Arm (capped at its tail)."""
    from rdkit import Chem
    from rdkit.Chem import AllChem

    bt = {1: Chem.BondType.SINGLE, 2: Chem.BondType.DOUBLE, 3: Chem.BondType.TRIPLE,
          "ar": Chem.BondType.AROMATIC}
    res = RESIDUES[name]
    rw, idx, names, elems = Chem.RWMol(), {}, [], []

    def add(nm, el):
        a = Chem.Atom(el)
        a.SetNoImplicit(True)
        idx[nm] = rw.AddAtom(a)
        names.append(nm)
        elems.append(el)
        return idx[nm]

    for a in res["heavy"]:
        add(a, element_of(a))
    for parent, hs in res["hydrogens"].items():
        for h in hs:
            add(h, "H")
    bonds = []
    for a, b, o in res["bonds"]:
        rw.AddBond(idx[a], idx[b], bt[o])
        bonds.append((a, b))
    for parent, hs in res["hydrogens"].items():
        for h in hs:
            rw.AddBond(idx[parent], idx[h], bt[1])
            bonds.append((parent, h))
    # cap the backbone carbonyl carbon, and the tail if it is not the backbone C
    caps = []
    for anchor in {"C", res.get("tail")} - {None}:
        nm = f"HCAP{len(caps)}"
        add(nm, "H")
        rw.AddBond(idx[anchor], idx[nm], bt[1])
        bonds.append((anchor, nm))
        caps.append(nm)
    mol = rw.GetMol()
    Chem.SanitizeMol(mol)
    if AllChem.EmbedMolecule(mol, randomSeed=0xF00D, maxAttempts=1000) != 0:
        AllChem.EmbedMolecule(mol, randomSeed=0xF00D, useRandomCoords=True,
                              maxAttempts=2000)
    props = AllChem.MMFFGetMoleculeProperties(mol)
    ff = (AllChem.MMFFGetMoleculeForceField(mol, props) if props
          else AllChem.UFFGetMoleculeForceField(mol))
    ff.Minimize(maxIts=20000)
    conf = mol.GetConformer()
    xyz = np.array([[conf.GetAtomPosition(i).x, conf.GetAtomPosition(i).y,
                     conf.GetAtomPosition(i).z] for i in range(mol.GetNumAtoms())])
    return Arm(names, xyz, bonds, elems)


def validate_cage(arm: Arm) -> dict:
    g = lambda n: arm.xyz[arm.index[n]]  # noqa: E731
    sides = [np.linalg.norm(g("CO1") - g("CO8")), np.linalg.norm(g("CO8") - g("CO9")),
             np.linalg.norm(g("CO9") - g("CO1"))]
    rep = {"cyclopropane_sides_A": [round(float(s), 3) for s in sides],
           "azadiene_torsion_deg": round(float(np.degrees(
               dihedral(g("CP3"), g("NP2"), g("NP1"), g("CP6")))), 2),
           "CB_to_OL_A": round(float(np.linalg.norm(g("OL*") - g("CB"))), 2),
           "n_rotatable_bonds": len(arm.rotatable())}
    if max(sides) > 1.7:
        sys.exit(f"cyclopropane is not closed: sides {rep['cyclopropane_sides_A']}")
    return rep


if __name__ == "__main__":
    a = build_arm()
    print(f"arm: {len(a.names)} atoms "
          f"({sum(1 for e in a.elements if e != 'H')} heavy)")
    rep = validate_cage(a)
    for k, v in rep.items():
        print(f"  {k}: {v}")
    print("\n  rotatable bonds (distal side moves):")
    for i, j, mov in a.rotatable():
        print(f"    {a.names[i]:<5}-{a.names[j]:<5}  {len(mov):>2} atoms move")
