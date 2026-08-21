#!/usr/bin/env python3
"""Mutate a residue to one of the custom residues, using real 3D geometry.

Coordinates come from an MMFF-optimised model of the residue, superimposed on the
target's own N/CA/CB, then its torsions are scanned to relieve clashes with the rest
of the protein.  tleap is not asked to build the side chain: an Amber prep file is a
tree of internal coordinates, so tleap places ring-closing bonds wherever the tree
walk happens to leave them -- for TDP that tore the bicyclononane cage open
(cyclopropane sides 5.28/4.70/1.52 A instead of three of 1.51).
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
from scipy.spatial import cKDTree

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from arm import Arm, build_residue_mol  # noqa: E402
from assemble_chimera import WATER, clash_score, fmt, read_pdb  # noqa: E402


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--protein", required=True)
    ap.add_argument("--residue", required=True, help="TET or TDP")
    ap.add_argument("--sites", nargs="+", type=int, required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--seed", type=int, default=11)
    args = ap.parse_args()
    rng = np.random.default_rng(args.seed)

    prot = read_pdb(args.protein)
    solute = [a for a in prot if a["resn"] not in WATER]
    waters = [a for a in prot if a["resn"] in WATER]
    pres = {}
    for a in solute:
        pres.setdefault(a["resi"], {})[a["name"]] = a

    others = np.array([a["xyz"] for a in solute
                       if a["resi"] not in args.sites])
    tree = cKDTree(others)
    placed = {}
    for site in args.sites:
        mol = build_residue_mol(args.residue)
        mol.place_on_backbone({k: pres[site][k]["xyz"] for k in ("N", "CA", "CB")})
        rot = mol.rotatable()
        heavy = [i for i, e in enumerate(mol.elements) if e != "H"]
        heavy_tree = tree
        rot = [r for r in rot]

        def score(x):
            return clash_score(x[heavy], heavy_tree)

        best_xyz, best = mol.xyz.copy(), score(mol.xyz)
        for start_i in range(24):
            mol.xyz = best_xyz.copy() if start_i == 0 else mol.xyz
            if start_i:
                mol.xyz = best_xyz.copy()
                for (i, j, mv) in rot:
                    mol.rotate(i, j, mv, rng.uniform(0, 2 * np.pi))
            cur = score(mol.xyz)
            for coarse in (60, 20, 6):
                for _ in range(4):
                    improved = False
                    for (i, j, mv) in rot:
                        for sgn in (+1, -1):
                            mol.rotate(i, j, mv, np.radians(sgn * coarse))
                            new = score(mol.xyz)
                            if new < cur - 1e-9:
                                cur, improved = new, True
                            else:
                                mol.rotate(i, j, mv, np.radians(-sgn * coarse))
                    if not improved:
                        break
            if cur < best:
                best, best_xyz = cur, mol.xyz.copy()
        mol.xyz = best_xyz
        print(f"  {args.residue} at {site}: clash score {best:.2f}")
        placed[site] = mol

    lines, serial = [], 1
    bb = {"N", "CA", "C", "O", "H", "HA"}
    for resi in sorted(pres):
        if resi in placed:
            m = placed[resi]
            for nm in m.names:
                if nm.startswith("HCAP"):
                    continue
                xyz = (pres[resi][nm]["xyz"] if nm in bb and nm in pres[resi]
                       else m.xyz[m.index[nm]])
                lines.append(fmt(serial, nm, args.residue, "A", resi, xyz,
                                 m.elements[m.index[nm]]))
                serial += 1
        else:
            for nm, at in pres[resi].items():
                lines.append(fmt(serial, nm, at["resn"], "A", resi, at["xyz"],
                                 at["elem"]))
                serial += 1
    lines.append(f"TER   {serial:5d}")
    serial += 1
    newxyz = np.array([m.xyz[i] for m in placed.values() for i in range(len(m.names))])
    ntree = cKDTree(newxyz)
    wres = {}
    for w in waters:
        wres.setdefault(w["resi"], []).append(w)
    kept = 0
    for resi, ws in wres.items():
        if min(ntree.query(np.array([w["xyz"] for w in ws]), k=1)[0]) < 2.6:
            continue
        kept += 1
        for w in ws:
            lines.append(fmt(serial, w["name"], "WAT", "W", resi, w["xyz"], w["elem"]))
            serial += 1
        lines.append(f"TER   {serial:5d}")
        serial += 1
    Path(args.out).write_text("\n".join(lines) + "\nEND\n")
    print(f"wrote {args.out} ({serial - 1} atoms, {kept}/{len(wres)} waters kept)")


if __name__ == "__main__":
    main()
