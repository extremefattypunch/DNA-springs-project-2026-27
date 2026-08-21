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
from assemble_chimera import WATER, clash_score, fmt, read_pdb  # noqa: E402


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--protein", required=True)
    ap.add_argument("--sites", nargs=2, type=int, default=[133, 149])
    ap.add_argument("--out", required=True)
    ap.add_argument("--seed", type=int, default=23)
    ap.add_argument("--starts", type=int, default=32)
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
        mol = build_arm("DNH")
        mol.place_on_backbone({k: pres[site][k]["xyz"] for k in ("N", "CA", "CB")})
        rot = mol.rotatable()
        heavy = [i for i, e in enumerate(mol.elements) if e != "H"]
        ol = mol.index["OL*"]
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
    tdp_ord, dnh_ord = {}, {}
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
                lines.append(fmt(serial, nm, "TDP", "A", resi, xyz,
                                 m.elements[m.index[nm]]))
                serial += 1
        else:
            for nm, at in pres[resi].items():
                lines.append(fmt(serial, nm, at["resn"], "A", resi, at["xyz"],
                                 at["elem"]))
                serial += 1
    lines.append(f"TER   {serial:5d}")
    serial += 1
    for k, site in enumerate(a.sites):
        m = arms[site]
        ordinal += 1
        dnh_ord[site] = ordinal
        for nm in m.names:
            if not nm.endswith("*"):
                continue
            lines.append(fmt(serial, nm[:-1], "DNH", "EF"[k], 1,
                             m.xyz[m.index[nm]], m.elements[m.index[nm]]))
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
    bonds = [f"{tdp_ord[s]}@CN:{dnh_ord[s]}@N" for s in a.sites]
    Path(a.out).with_suffix(".json").write_text(json.dumps(
        {"sites": a.sites, "tleap_bonds": bonds, "waters_kept": kept,
         "waters_dropped": len(wres) - kept}, indent=2))
    print(f"wrote {a.out} ({serial - 1} atoms, {kept}/{len(wres)} waters kept)")
    print(f"tleap bonds: {' '.join(bonds)}")


if __name__ == "__main__":
    main()
