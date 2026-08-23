#!/usr/bin/env python3
"""Assemble the clicked-but-unloaded control: sfGFP + two full tethers, no DNA.

This is the control that makes the chimera comparison mean something.  The
force-clamp ladder runs on TET -- the unclicked tetrazine -- so it cannot serve as
the zero-force reference for a chimera built from TDP: the two differ by the entire
sTCO adduct plus a six-carbon arm, and that difference moves the attachment sites by
about as much as the spring does.  Here the chemistry is identical to a chimera's and
only the duplex is missing, so a chimera minus this control isolates the load.

The tethers are posed to point away from the barrel and away from each other, which
is the unloaded state a free arm would adopt on average -- not curled against the
surface, which would bias the sites inwards and manufacture the very effect we are
trying to measure.
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
from arm import build_arm  # noqa: E402
sys.path.insert(0, str(HERE.parent / "02_params"))
from residue_defs import RESIDUES  # noqa: E402
from assemble_chimera import WATER, clash_score, fmt, read_pdb  # noqa: E402


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--protein", required=True)
    ap.add_argument("--sites", nargs=2, type=int, default=[133, 149])
    ap.add_argument("--out", required=True)
    ap.add_argument("--seed", type=int, default=23)
    ap.add_argument("--starts", type=int, default=32)
    ap.add_argument("--arm", nargs="+", default=["TDP", "DNH"], metavar="RESIDUE",
                    help="the tether, protein outwards, terminating in something inert. "
                         "sfGFP: TDP DNH.  Ubiquitin: CMS CDB DNH.")
    a = ap.parse_args()
    rng = np.random.default_rng(a.seed)

    prot = read_pdb(a.protein)
    solute = [x for x in prot if x["resn"] not in WATER]
    waters = [x for x in prot if x["resn"] in WATER]
    pres = {}
    for x in solute:
        pres.setdefault(x["resi"], {})[x["name"]] = x
    others = np.array([x["xyz"] for x in solute if x["resi"] not in a.sites])
    tree = cKDTree(others)
    centroid = np.array([x["xyz"] for x in solute]).mean(axis=0)

    arms = {}
    for site in a.sites:
        mol = build_arm(tuple(a.arm))
        mol.place_on_backbone({k: pres[site][k]["xyz"] for k in ("N", "CA", "CB")})
        rot = mol.rotatable()
        heavy = [i for i, e in enumerate(mol.elements) if e != "H"]
        ol = mol.index[mol.dna_atom]
        cb = mol.xyz[mol.index["CB"]]
        outward = cb - centroid
        outward /= np.linalg.norm(outward)

        def score(x):
            # reward the arm reaching outward along the local surface normal, and
            # penalise contacts; no target point, because there is nothing to reach
            return clash_score(x[heavy], tree) - 2.0 * float((x[ol] - cb) @ outward)

        best_xyz, best = mol.xyz.copy(), score(mol.xyz)
        for s in range(a.starts):
            mol.xyz = best_xyz.copy()
            if s:
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
        arms[site] = mol
        reach = float(np.linalg.norm(mol.xyz[ol] - cb))
        print(f"  tether at {site}: extends {reach:.1f} A from CB, "
              f"clash {clash_score(mol.xyz[heavy], tree):.2f}")

    lines, serial, ordinal = [], 1, 0
    tdp_ord, link_ord = {}, {}
    bb = {"N", "CA", "C", "O", "H", "HA"}
    for resi in sorted(pres):
        ordinal += 1
        if resi in arms:
            tdp_ord[resi] = ordinal
            m = arms[resi]
            for nm in m.names:
                if nm.endswith("*") or nm.startswith("HCAP"):
                    continue
                xyz = (pres[resi][nm]["xyz"] if nm in bb and nm in pres[resi]
                       else m.xyz[m.index[nm]])
                lines.append(fmt(serial, nm, a.arm[0], "A", resi, xyz,
                                 m.elements[m.index[nm]]))
                serial += 1
        else:
            for nm, at in pres[resi].items():
                lines.append(fmt(serial, nm, at["resn"], "A", resi, at["xyz"],
                                 at["elem"]))
                serial += 1
    lines.append(f"TER   {serial:5d}")
    serial += 1
    # One residue per star level, so a three-residue tether writes its two non-protein
    # residues separately.  Each arm gets its own chain letter and its residues are
    # numbered 1..k within it; nothing chains head-to-tail here, since every link is
    # declared explicitly below.
    for k, site in enumerate(a.sites):
        m = arms[site]
        for level in range(1, len(a.arm)):
            star = "*" * level
            ordinal += 1
            link_ord[(site, level)] = ordinal
            n = 0
            for nm in m.names:
                if not (nm.endswith(star) and not nm.endswith(star + "*")):
                    continue
                lines.append(fmt(serial, nm[:-level], a.arm[level], "EF"[k], level,
                                 m.xyz[m.index[nm]], m.elements[m.index[nm]]))
                serial += 1
                n += 1
            if n == 0:
                sys.exit(f"no atoms at star level {level} for the arm at {site}")
            # See the note in assemble_chimera.py: tleap chains a prep residue's last
            # main-chain atom to the next residue's head on its own, so two linker
            # residues written back to back get an extra bond nobody asked for.
            if level < len(a.arm) - 1:
                lines.append(f"TER   {serial:5d}")
                serial += 1
        lines.append(f"TER   {serial:5d}")
        serial += 1

    newxyz = np.array([m.xyz[i] for m in arms.values() for i in range(len(m.names))])
    ntree = cKDTree(newxyz)
    wres = {}
    for w in waters:
        wres.setdefault(w["resi"], []).append(w)
    kept = 0
    for resi, ws in wres.items():
        if min(ntree.query(np.array([w["xyz"] for w in ws]), k=1)[0]) < 2.6:
            continue
        kept += 1
        ordinal += 1
        for w in ws:
            lines.append(fmt(serial, w["name"], "WAT", "W", resi, w["xyz"], w["elem"]))
            serial += 1
        lines.append(f"TER   {serial:5d}")
        serial += 1
    Path(a.out).write_text("\n".join(lines) + "\nEND\n")
    # Read from each residue's own extra_bond, so the declared topology and the built
    # geometry come from the same tables.
    bonds = []
    for s in a.sites:
        for level in range(len(a.arm) - 1):
            eb = RESIDUES[a.arm[level]]["extra_bond"]
            up = tdp_ord[s] if level == 0 else link_ord[(s, level)]
            bonds.append(f"{up}@{eb['this_atom']}:"
                         f"{link_ord[(s, level + 1)]}@{eb['partner_atom']}")
    Path(a.out).with_suffix(".json").write_text(json.dumps(
        {"sites": a.sites, "tleap_bonds": bonds, "waters_kept": kept,
         "waters_dropped": len(wres) - kept}, indent=2))
    print(f"wrote {a.out} ({serial - 1} atoms, {kept}/{len(wres)} waters kept)")
    print(f"tleap bonds: {' '.join(bonds)}")


if __name__ == "__main__":
    main()
