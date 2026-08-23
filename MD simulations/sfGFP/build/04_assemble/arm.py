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

    def __init__(self, names, coords, bonds, elements,
                 dna_atom: str = "OL*", linker: str = "DNL", chain=("TDP", "DNL")):
        # `dna_atom` is the arm atom that bonds the duplex.  It is not always the same
        # atom: the 5' linker ends in the oxygen that becomes nucleotide 1's bridging
        # O5' (DNL.OL), while the 3' linker ends in its own phosphorus, which bonds
        # nucleotide N's existing O3' (DNP.P).  Callers must not hard-code "OL*".
        self.dna_atom = dna_atom
        self.linker = linker
        self.chain = tuple(chain)
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


def build_arm(chain=("TDP", "DNL")) -> Arm:
    """A whole tether -- protein-side residue through to the DNA -- as one molecule.

    ``chain`` is the residue names in order from the protein outwards.  Two residues is
    the sfGFP tether (``("TDP", "DNL")``, or DNP at the 3' end, or DNH for the clicked
    control that stops at a free alcohol); three is the ubiquitin one
    (``("CMS", "CDB", "DNL")``), where the crosslinker is split at an amide so that
    neither half is an unreasonable AM1-BCC problem.  A bare string is accepted as the
    old two-residue form for convenience: ``build_arm("DNP")``.

    The whole thing is built and MMFF-optimised as a single molecule rather than
    residue by residue, because the geometry that matters is across the junctions --
    and because tleap cannot build the polycyclic pieces at all: handed the cyclopropane
    of the sTCO cage it returned sides of 5.28, 4.70 and 1.52 A, and the DBCO cage is
    worse.  Consecutive residues are joined using each one's own declared extra_bond,
    so the connectivity comes from the same tables the parameters were fitted against.
    """
    if isinstance(chain, str):
        chain = ("TDP", chain)
    from rdkit import Chem
    from rdkit.Chem import AllChem

    bt = {1: Chem.BondType.SINGLE, 2: Chem.BondType.DOUBLE, 3: Chem.BondType.TRIPLE,
          "ar": Chem.BondType.AROMATIC}
    residues = [RESIDUES[n] for n in chain]
    first, dnl = residues[0], residues[-1]
    rw, idx, names, elems = Chem.RWMol(), {}, [], []

    # Atoms of the protein-side residue keep their bare names; each residue further out
    # gets one more trailing star.  The stars are how the assembler knows which atoms
    # belong to which residue when it writes the PDB, and they are stripped there.
    def suffix(tag):
        return "*" * tag

    def add(nm, el, tag):
        if (tag, nm) in idx:
            raise ValueError(f"duplicate atom {chain[tag]}.{nm} in the arm")
        a = Chem.Atom(el)
        a.SetNoImplicit(True)
        idx[(tag, nm)] = rw.AddAtom(a)
        names.append(f"{nm}{suffix(tag)}")
        elems.append(el)
        return idx[(tag, nm)]

    for tag, res in enumerate(residues):
        for a in res["heavy"]:
            add(a, element_of(a), tag)
        for parent, hs in res["hydrogens"].items():
            for h in hs:
                add(h, "H", tag)
    bonds_named = []
    for tag, res in enumerate(residues):
        for a, b, o in res["bonds"]:
            rw.AddBond(idx[(tag, a)], idx[(tag, b)], bt[o])
            bonds_named.append(((tag, a), (tag, b)))
        for parent, hs in res["hydrogens"].items():
            for h in hs:
                rw.AddBond(idx[(tag, parent)], idx[(tag, h)], bt[1])
                bonds_named.append(((tag, parent), (tag, h)))
    # Join consecutive residues using the upstream one's declared extra_bond, and check
    # it points where we are actually attaching: a mismatch here would build a molecule
    # whose connectivity differs from the tables the charges were fitted against.
    for tag in range(len(residues) - 1):
        eb = residues[tag].get("extra_bond")
        if eb is None:
            sys.exit(f"{chain[tag]} declares no extra_bond, so it cannot be chained "
                     f"to {chain[tag + 1]}")
        want, got = eb["partner_residue"], chain[tag + 1]
        # DNL, DNP and DNH are interchangeable at this junction -- all three present an
        # amide nitrogen named N -- so only the amino-arm family is treated as a match.
        family = {"DNL", "DNP", "DNH"}
        if want != got and not ({want, got} <= family):
            sys.exit(f"{chain[tag]}.extra_bond points at {want}, not {got}")
        rw.AddBond(idx[(tag, eb["this_atom"])],
                   idx[(tag + 1, eb["partner_atom"])], bt[1])
        bonds_named.append(((tag, eb["this_atom"]),
                            (tag + 1, eb["partner_atom"])))

    # Cap the free valences so RDKit can sanitise.  The backbone N and C=O already
    # carry their hydrogens; the backbone C needs one more substituent, and so does
    # whichever linker terminal the DNA will take over.  Which terminal that is
    # differs between the two arms, and so does the cap:
    #   DNL hands its TAIL oxygen to nucleotide 1's phosphorus  -> one H (a free
    #       alcohol, which is exactly what the DNH control keeps);
    #   DNP hands its HEAD phosphorus to nucleotide N's O3'      -> OH, which makes the
    #       fragment a phosphate monoester.  A bare H on phosphorus would satisfy the
    #       valence but MMFF has no business typing a P-H here.
    last = len(residues) - 1
    cap_ids, cap_chain = [], []
    # An amino-acid residue is recognised by its head/tail being the backbone N and C
    # rather than by a flag, because that flag lives in build_residues' own table and
    # not in the residue definitions; the backbone carbonyl then needs one substituent.
    if first["head"] == "N" and first.get("tail") == "C":
        cap_chain.append((0, "C", ["H"]))
    if dnl.get("tail"):
        cap_chain.append((last, dnl["tail"], ["H"]))
        dna_atom = f"{dnl['tail']}{suffix(last)}"
    elif dnl["head"] != "N":
        cap_chain.append((last, dnl["head"], ["O", "H"]))
        dna_atom = f"{dnl['head']}{suffix(last)}"
    else:
        # DNH, the clicked control: its head N is already satisfied by the carbamate
        # and its far end is a real free alcohol with its hydrogen in the tables, so
        # there is no free valence to cap and nothing is handed to a duplex.
        dna_atom = f"OL{suffix(last)}"
    for tag, nm, elements in cap_chain:
        prev = (tag, nm)
        for el in elements:
            cname = f"{el}CAP{len(cap_ids)}"
            i = add(cname, el, tag)
            rw.AddBond(idx[prev], i, bt[1])
            bonds_named.append((prev, (tag, cname)))
            cap_ids.append(i)
            prev = (tag, cname)

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

    # drop the capping atoms
    keep = [i for i in range(mol.GetNumAtoms()) if i not in cap_ids]
    remap = {old: new for new, old in enumerate(keep)}
    names2 = [names[i] for i in keep]
    elems2 = [elems[i] for i in keep]
    bonds2 = [(names[a], names[b]) for (ta, na), (tb, nb) in bonds_named
              for a, b in [(idx[(ta, na)], idx[(tb, nb)])]
              if a in remap and b in remap]
    return Arm(names2, xyz[keep], bonds2, elems2,
               dna_atom=dna_atom, linker=chain[-1], chain=tuple(chain))


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


# Rings whose geometry is worth asserting, by the atom names that make them up.  Each
# arm has only some of these, so each is checked only if its atoms are present -- the
# sfGFP tether has the strained sTCO cyclopropane and the dihydropyridazine, the
# ubiquitin one has the DBCO cage's two benzo rings and a triazole.  Names carry the
# star suffix of whichever residue they belong to, so the lookup tries every suffix.
# Each entry is (atom names, maximum rms deviation from the ring's own best-fit plane).
# The tolerances differ because the rings differ: an aromatic ring should be flat to
# within a few hundredths of an angstrom, but the clicked product is a 4,5-DIHYDRO
# pyridazine, whose two sp3 carbons make it a half-chair -- it is not supposed to be
# planar, and holding it to an aromatic tolerance rejects a perfectly good structure.
# What matters there is the azadiene, checked separately as a torsion.
RINGS = {
    "sTCO cyclopropane": (("CO1", "CO8", "CO9"), None),
    "dihydropyridazine": (("CP3", "NP2", "NP1", "CP6", "CQ5", "CQ4"), 0.45),
    "tetrazine": (("CT3", "NT2", "NT1", "CT6", "NT5", "NT4"), 0.10),
    "DBCO benzo 1": (("CP1", "CP2", "CP3", "CP4", "CP5", "CP6"), 0.10),
    "DBCO benzo 2": (("CQ1", "CQ2", "CQ3", "CQ4", "CQ5", "CQ6"), 0.10),
    "SPAAC triazole": (("NJ1", "NJ2", "NJ3", "CJ1", "CJ2"), 0.10),
}
# The conjugated C=N-N=C of the clicked adduct, which should be close to flat even
# though the ring around it is not.
AZADIENE = ("CP3", "NP2", "NP1", "CP6")
AZADIENE_MAX_DEV_DEG = 35.0


def _find(arm, name):
    """Resolve an atom name that may carry any number of trailing stars."""
    for s in ("", "*", "**", "***"):
        if name + s in arm.index:
            return arm.xyz[arm.index[name + s]]
    return None


def validate_cage(arm: Arm) -> dict:
    """Check the rings the arm actually has, and report its reach.

    Written to be adaptive rather than to assume one tether's atom names, because the
    two chemistries share no rings: asking the ubiquitin arm for the sTCO cyclopropane
    simply raised a KeyError.  Ring geometry is worth asserting at all because RDKit's
    embedding plus MMFF is the only thing standing between us and a broken cage --
    handed the same job, tleap returned cyclopropane sides of 5.28, 4.70 and 1.52 A.
    """
    rep = {"chain": " + ".join(arm.chain), "linker": arm.linker,
           "dna_atom": arm.dna_atom, "n_rotatable_bonds": len(arm.rotatable())}
    cb = _find(arm, "CB")
    dna = arm.xyz[arm.index[arm.dna_atom]]
    if cb is not None:
        rep["CB_to_dna_atom_A"] = round(float(np.linalg.norm(dna - cb)), 2)

    problems = []
    for label, (names, max_rms) in RINGS.items():
        pts = [_find(arm, n) for n in names]
        if any(q is None for q in pts):
            continue
        pts = np.array(pts)
        if len(pts) == 3:                       # a three-ring has no planarity to test
            sides = [float(np.linalg.norm(pts[i] - pts[(i + 1) % 3])) for i in range(3)]
            rep[f"{label} sides_A"] = [round(s, 3) for s in sides]
            if max(sides) > 1.7:
                problems.append(f"{label} is not closed: {sides}")
            continue
        c = pts.mean(axis=0)
        normal = np.linalg.svd(pts - c)[2][-1]
        rms = float(np.sqrt((((pts - c) @ normal) ** 2).mean()))
        rep[f"{label} planarity_rms_A"] = round(rms, 3)
        if max_rms is not None and rms > max_rms:
            problems.append(f"{label} rms {rms:.3f} A exceeds {max_rms} A")

    az = [_find(arm, n) for n in AZADIENE]
    if all(q is not None for q in az):
        d = float(np.degrees(dihedral(*az)))
        rep["azadiene_torsion_deg"] = round(d, 2)
        dev = min(abs(abs(d) - 180.0), abs(d))
        if dev > AZADIENE_MAX_DEV_DEG:
            problems.append(f"azadiene twisted {dev:.1f} deg out of plane")
    if problems:
        sys.exit("arm geometry: " + "; ".join(problems))
    rep["rings_checked"] = sum(1 for k in rep if "sides_A" in k or "planarity" in k)
    return rep


if __name__ == "__main__":
    for chain in (("TDP", "DNL"), ("TDP", "DNP"), ("TDP", "DNH"),
                  ("CMS", "CDB", "DNL"), ("CMS", "CDB", "DNP")):
        a = build_arm(chain)
        print(f"\n=== {' + '.join(chain)} ===")
        print(f"  {len(a.names)} atoms "
              f"({sum(1 for e in a.elements if e != 'H')} heavy)")
        for k, v in validate_cage(a).items():
            print(f"  {k}: {v}")
        print(f"  rotatable bonds: {len(a.rotatable())}")
