#!/usr/bin/env python3
"""Assemble the sfGFP-DNA spring chimera: protein + two tethers + a bent duplex.

Geometry, in order
------------------
1. The two attachment sites A = CB(Asp133) and B = CB(Asn149) are 31.3 A apart.  The
   duplex's own 5'P-to-5'P distance D is fixed by how far it was bent.  Target
   phosphate positions are placed symmetrically:

       P_A = A - a*u + b*n        P_B = B + a*u + b*n
       u = unit(B - A),  a = (D - |AB|)/2,  n = outward surface normal

   The u components differ by |AB| + 2a = D by construction, and the n components
   cancel, so the target separation matches the duplex exactly for any b.  b then
   lifts both attachment points off the protein surface -- without it the tethers
   would have to run flat along the barrel -- and is capped so that each tether's
   span stays inside its reach.

2. The duplex is placed by aligning its P-to-P vector to the target vector.  Two
   points leave one degree of freedom, a spin about that axis, which is scanned and
   chosen to push the arc away from the protein and minimise contacts.

3. Each tether is then posed by rigid-body superposition onto its residue's backbone
   followed by torsion optimisation to bring its bridging oxygen to the phosphate.
   Coordinate descent from several random starts, scoring distance-to-target plus a
   clash penalty.

The DNL linker leads each DNA strand in the output rather than being a separate
molecule.  That is not cosmetic: tleap applies the 5'-terminal template to whichever
residue starts a nucleic-acid chain, which strips the phosphate the tether needs to
bond to.  With DNL first, the nucleotide is no longer chain-initial, keeps its
phosphate, and tleap makes the DNL-to-phosphorus bond itself through the residues'
head/tail connectivity -- leaving only the two carbamate bonds to declare explicitly.
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
from arm import Arm, build_arm, validate_cage  # noqa: E402

WATER = {"WAT", "HOH"}


def read_pdb(path):
    out = []
    for l in Path(path).read_text().splitlines():
        if l[:6] not in ("ATOM  ", "HETATM"):
            continue
        el = (l[76:78].strip() or l[12:16].strip()[0]).upper()
        out.append({"name": l[12:16].strip(), "resn": l[17:20].strip(),
                    "chain": l[21], "resi": int(l[22:26]), "elem": el,
                    "xyz": np.array([float(l[30:38]), float(l[38:46]),
                                     float(l[46:54])])})
    return out


def fmt(serial, name, resn, chain, resi, xyz, elem):
    fld = f" {name:<3}" if len(name) < 4 else name
    return (f"ATOM  {serial:5d} {fld}{resn:>4}{chain:>2}{resi:4d}    "
            f"{xyz[0]:8.3f}{xyz[1]:8.3f}{xyz[2]:8.3f}"
            f"{1.00:6.2f}{0.00:6.2f}          {elem:>2}")


def align_two_points(src_a, src_b, dst_a, dst_b, spin_deg, extra=None):
    """Rigid transform taking src_a,src_b onto dst_a,dst_b, plus a spin about the axis."""
    from arm import rodrigues
    v1, v2 = src_b - src_a, dst_b - dst_a
    n1, n2 = np.linalg.norm(v1), np.linalg.norm(v2)
    axis = np.cross(v1 / n1, v2 / n2)
    s = np.linalg.norm(axis)
    if s < 1e-9:
        R = np.eye(3) if (v1 / n1) @ (v2 / n2) > 0 else -np.eye(3)
    else:
        ang = np.arctan2(s, (v1 / n1) @ (v2 / n2))
        R = rodrigues(axis, ang)
    R = rodrigues(v2, np.radians(spin_deg)) @ R
    mid_s, mid_d = 0.5 * (src_a + src_b), 0.5 * (dst_a + dst_b)
    return R, mid_d - R @ mid_s


def clash_score(pts, tree, cutoff=3.0):
    if len(pts) == 0:
        return 0.0
    d, _ = tree.query(pts, k=1)
    viol = np.clip(cutoff - d, 0, None)
    return float((viol ** 2).sum())


def pose_arm(arm: Arm, site: dict, ol_target, protein_tree, rng, n_starts=8,
             sweeps=4, step_deg=(60, 20, 6)):
    """Coordinate descent over the tether's torsions to reach ol_target."""
    arm.place_on_backbone(site)
    rot = arm.rotatable()
    ol = arm.index["OL*"]
    heavy = [i for i, e in enumerate(arm.elements) if e != "H"]
    base = arm.xyz.copy()

    def score():
        d = np.linalg.norm(arm.xyz[ol] - ol_target)
        return d * d + 0.5 * clash_score(arm.xyz[heavy], protein_tree)

    best_xyz, best = None, np.inf
    for start in range(n_starts):
        arm.xyz = base.copy()
        if start:
            for (i, j, mov) in rot:
                arm.rotate(i, j, mov, rng.uniform(0, 2 * np.pi))
        cur = score()
        for coarse in step_deg:
            for _ in range(sweeps):
                improved = False
                for (i, j, mov) in rot:
                    for sgn in (+1, -1):
                        arm.rotate(i, j, mov, np.radians(sgn * coarse))
                        new = score()
                        if new < cur - 1e-9:
                            cur, improved = new, True
                        else:
                            arm.rotate(i, j, mov, np.radians(-sgn * coarse))
                if not improved:
                    break
        if cur < best:
            best, best_xyz = cur, arm.xyz.copy()
    arm.xyz = best_xyz
    return {"score": round(float(best), 3),
            "OL_to_target_A": round(float(np.linalg.norm(arm.xyz[ol] - ol_target)), 3),
            "clash": round(clash_score(arm.xyz[heavy], protein_tree), 3)}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--protein", required=True)
    ap.add_argument("--dna", required=True)
    ap.add_argument("--sites", nargs=2, type=int, default=[133, 149])
    ap.add_argument("--out", required=True)
    ap.add_argument("--extension-frac", type=float, default=0.88,
                    help="fraction of the tether's full reach to use for the span")
    ap.add_argument("--water-clash", type=float, default=2.6)
    ap.add_argument("--seed", type=int, default=7)
    args = ap.parse_args()
    rng = np.random.default_rng(args.seed)

    prot = read_pdb(args.protein)
    dna = read_pdb(args.dna)
    solute = [a for a in prot if a["resn"] not in WATER]
    waters = [a for a in prot if a["resn"] in WATER]
    pres = {}
    for a in solute:
        pres.setdefault(a["resi"], {})[a["name"]] = a
    for s in args.sites:
        if s not in pres:
            sys.exit(f"site {s} not present in {args.protein}")

    arm = build_arm()
    cage = validate_cage(arm)
    # full reach: fully extend every torsion away from CB, then measure
    reach_arm = build_arm()
    reach_arm.place_on_backbone({k: pres[args.sites[0]][k]["xyz"]
                                for k in ("N", "CA", "CB")})
    cb = reach_arm.xyz[reach_arm.index["CB"]]
    far = cb + 1000.0 * (cb - np.mean([a["xyz"] for a in solute], axis=0))
    pose_arm(reach_arm, {k: pres[args.sites[0]][k]["xyz"] for k in ("N", "CA", "CB")},
             far, cKDTree(np.array([[1e6, 1e6, 1e6]])), rng, n_starts=4)
    reach = float(np.linalg.norm(reach_arm.xyz[reach_arm.index["OL*"]] - cb))
    print(f"tether: {cage['n_rotatable_bonds']} torsions, cyclopropane "
          f"{cage['cyclopropane_sides_A']}, extended CB->OL reach {reach:.2f} A")

    # --- target phosphate positions ---
    A = pres[args.sites[0]]["CB"]["xyz"]
    B = pres[args.sites[1]]["CB"]["xyz"]
    u = (B - A) / np.linalg.norm(B - A)
    centroid = np.mean([a["xyz"] for a in solute], axis=0)
    n = 0.5 * (A + B) - centroid
    n = n - (n @ u) * u
    n /= np.linalg.norm(n)

    dres = {}
    for a in dna:
        dres.setdefault((a["chain"], a["resi"]), {})[a["name"]] = a
    dchains = sorted({c for c, _ in dres})
    PA = dres[(dchains[0], 1)]["P"]["xyz"]
    PB = dres[(dchains[1], 1)]["P"]["xyz"]
    D = float(np.linalg.norm(PA - PB))
    ab = float(np.linalg.norm(B - A))
    a_off = 0.5 * (D - ab)
    span_target = args.extension_frac * reach
    b_off = float(np.sqrt(max(0.0, span_target ** 2 - a_off ** 2)))
    tA = A - a_off * u + b_off * n
    tB = B + a_off * u + b_off * n
    print(f"sites {ab:.2f} A apart; duplex 5'-5' {D:.2f} A -> a = {a_off:.2f}, "
          f"b = {b_off:.2f} A; tether span {np.linalg.norm(tA - A):.2f} A "
          f"({100 * np.linalg.norm(tA - A) / reach:.0f}% of reach)")
    if np.linalg.norm(tA - A) > reach:
        sys.exit(f"required tether span {np.linalg.norm(tA - A):.2f} A exceeds the "
                 f"reach {reach:.2f} A: bend the duplex further")

    # --- place the duplex: align P-P, then scan the spin about that axis ---
    ptree = cKDTree(np.array([a["xyz"] for a in solute]))
    dxyz = np.array([a["xyz"] for a in dna])
    best = None
    for spin in range(0, 360, 5):
        R, t = align_two_points(PA, PB, tA, tB, spin)
        moved = (R @ dxyz.T).T + t
        cl = clash_score(moved, ptree, cutoff=4.0)
        # reward the arc bulging away from the protein
        bulge = float((moved.mean(axis=0) - 0.5 * (tA + tB)) @ n)
        sc = cl - 2.0 * bulge
        if best is None or sc < best[0]:
            best = (sc, spin, R, t, cl, bulge)
    sc, spin, R, t, cl, bulge = best
    for a in dna:
        a["xyz"] = R @ a["xyz"] + t
    dres = {}
    for a in dna:
        dres.setdefault((a["chain"], a["resi"]), {})[a["name"]] = a
    print(f"duplex placed: spin {spin} deg, clash score {cl:.1f}, "
          f"arc bulge {bulge:+.1f} A along the outward normal")

    # --- pose the two tethers ---
    arms, reports = {}, {}
    for site, tgt, ch in ((args.sites[0], tA, dchains[0]),
                          (args.sites[1], tB, dchains[1])):
        P = dres[(ch, 1)]["P"]["xyz"]
        ol_t = P + 1.61 * (A + B) / 2.0 * 0  # placeholder, replaced below
        ol_t = P - 1.61 * (P - pres[site]["CB"]["xyz"]) / np.linalg.norm(
            P - pres[site]["CB"]["xyz"])
        this = build_arm()
        rep = pose_arm(this, {k: pres[site][k]["xyz"] for k in ("N", "CA", "CB")},
                       ol_t, ptree, rng)
        ol = this.xyz[this.index["OL*"]]
        rep["OL_to_P_A"] = round(float(np.linalg.norm(ol - P)), 3)
        arms[site] = this
        reports[str(site)] = rep
        print(f"tether at {site}: OL-to-target {rep['OL_to_target_A']} A, "
              f"OL-P bond {rep['OL_to_P_A']} A (ideal 1.61), clash {rep['clash']}")

    # --- write the assembled PDB ---
    lines, serial = [], 1
    ordinal, tdp_ord, dnl_ord = 0, {}, {}
    resi_order = sorted(pres)
    for resi in resi_order:
        ordinal += 1
        if resi in args.sites:
            tdp_ord[resi] = ordinal
            a = arms[resi]
            bb = {"N", "CA", "C", "O", "H", "HA"}
            for nm in a.names:
                if nm.endswith("*"):
                    continue
                if nm in bb and nm not in pres[resi]:
                    continue
                xyz = (pres[resi][nm]["xyz"] if nm in pres[resi] and nm in bb
                       else a.xyz[a.index[nm]])
                lines.append(fmt(serial, nm, "TDP", "A", resi, xyz,
                                 a.elements[a.index[nm]]))
                serial += 1
        else:
            for nm, at in pres[resi].items():
                lines.append(fmt(serial, nm, at["resn"], "A", resi, at["xyz"],
                                 at["elem"]))
                serial += 1
    lines.append(f"TER   {serial:5d}")
    serial += 1

    dna_chain_ids = {dchains[0]: "C", dchains[1]: "D"}
    n_bp = max(r for _, r in dres)
    for src, site in ((dchains[0], args.sites[0]), (dchains[1], args.sites[1])):
        a = arms[site]
        ordinal += 1
        dnl_ord[site] = ordinal
        for nm in a.names:
            if not nm.endswith("*"):
                continue
            lines.append(fmt(serial, nm[:-1], "DNL", dna_chain_ids[src], 0,
                             a.xyz[a.index[nm]], a.elements[a.index[nm]]))
            serial += 1
        for i in range(1, n_bp + 1):
            ordinal += 1
            rn = next(iter(dres[(src, i)].values()))["resn"]
            for nm, at in dres[(src, i)].items():
                lines.append(fmt(serial, nm, rn, dna_chain_ids[src], i, at["xyz"],
                                 at["elem"]))
                serial += 1
        lines.append(f"TER   {serial:5d}")
        serial += 1

    # crystal waters, minus any that the DNA or tethers now occupy
    newxyz = np.array([a.xyz[i] for a in arms.values()
                       for i in range(len(a.names))] + [x["xyz"] for x in dna])
    ntree = cKDTree(newxyz)
    kept = 0
    wres = {}
    for w in waters:
        wres.setdefault(w["resi"], []).append(w)
    for resi, ws in wres.items():
        if min(ntree.query(np.array([w["xyz"] for w in ws]), k=1)[0]) < args.water_clash:
            continue
        ordinal += 1
        kept += 1
        for w in ws:
            lines.append(fmt(serial, w["name"], "WAT", "W", resi, w["xyz"], w["elem"]))
            serial += 1
        lines.append(f"TER   {serial:5d}")
        serial += 1
    Path(args.out).write_text("\n".join(lines) + "\nEND\n")

    bonds = [f"{tdp_ord[s]}@CN:{dnl_ord[s]}@N" for s in args.sites]
    meta = {"protein": args.protein, "dna": args.dna, "sites": args.sites,
            "tether": cage, "reach_A": round(reach, 2),
            "sites_CB_CB_A": round(ab, 2), "duplex_5p_5p_A": round(D, 2),
            "a_offset_A": round(a_off, 2), "b_offset_A": round(b_off, 2),
            "duplex_spin_deg": spin, "arm_reports": reports,
            "waters_kept": kept, "waters_dropped": len(wres) - kept,
            "tleap_bonds": bonds}
    Path(args.out).with_suffix(".json").write_text(json.dumps(meta, indent=2))
    print(f"\nwrote {args.out}  ({serial - 1} atoms, {kept} waters kept, "
          f"{len(wres) - kept} dropped)")
    print(f"tleap bonds to declare: {' '.join(bonds)}")


if __name__ == "__main__":
    main()
