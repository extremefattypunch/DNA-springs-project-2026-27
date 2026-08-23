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
import math
import sys
from pathlib import Path

import numpy as np
from scipy.spatial import cKDTree

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from arm import Arm, build_arm, validate_cage  # noqa: E402

WATER = {"WAT", "HOH"}
PRIME3 = "O3'"


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
             sweeps=4, step_deg=(60, 20, 6), bond_length=None):
    """Coordinate descent over the tether's torsions.

    With ``bond_length`` set, ``ol_target`` is the *anchor atom* and the objective is
    the length of the bond about to be formed to it.  Aiming at a fixed point instead
    -- one bond short of the anchor along some chosen direction -- overconstrains the
    problem: it fixes the bond's direction as well as its length, and missing that
    point by half an angstrom sideways can leave the bond itself half an angstrom too
    long.  The longer 3' arm, with 17 torsions to settle, is where that showed.
    """
    arm.place_on_backbone(site)
    rot = arm.rotatable()
    ol = arm.index[arm.dna_atom]
    heavy = [i for i, e in enumerate(arm.elements) if e != "H"]
    base = arm.xyz.copy()

    def score():
        d = np.linalg.norm(arm.xyz[ol] - ol_target)
        if bond_length is not None:
            d = d - bond_length
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
    return {"score": round(float(best), 3), "dna_atom": arm.dna_atom,
            "to_target_A": round(float(np.linalg.norm(arm.xyz[ol] - ol_target)), 3),
            "bond_objective": bond_length,
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

    # The two arms are chemically different, so they are built and measured separately.
    # Site 0 carries the 5' linker DNL, whose terminal oxygen becomes nucleotide 1's
    # bridging O5'.  Site 1 carries the 3' linker DNP, which brings its own phosphate
    # and bonds nucleotide N's existing O3'.  DNP is the longer arm by roughly the
    # width of that phosphate, so the placement below cannot assume symmetry.
    LINKERS = {args.sites[0]: "DNL", args.sites[1]: "DNP"}
    cage, reach = {}, {}
    for site, linker in LINKERS.items():
        a0 = build_arm(linker)
        cage[linker] = validate_cage(a0)
        # full reach: extend every torsion away from CB and away from the protein
        ra = build_arm(linker)
        bb = {k: pres[site][k]["xyz"] for k in ("N", "CA", "CB")}
        ra.place_on_backbone(bb)
        cb = ra.xyz[ra.index["CB"]]
        far = cb + 1000.0 * (cb - np.mean([a["xyz"] for a in solute], axis=0))
        pose_arm(ra, bb, far, cKDTree(np.array([[1e6, 1e6, 1e6]])), rng, n_starts=4)
        reach[site] = float(np.linalg.norm(ra.xyz[ra.index[ra.dna_atom]] - cb))
        print(f"tether at {site}: {linker}, {cage[linker]['n_rotatable_bonds']} "
              f"torsions, cyclopropane {cage[linker]['cyclopropane_sides_A']}, "
              f"extended CB->{ra.dna_atom[:-1]} reach {reach[site]:.2f} A")

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
    n_bp = max(r for _, r in dres)
    mod, comp = dchains[0], dchains[1]          # modified strand, complement
    # BOTH anchors are on the modified strand: its 5' phosphate and its 3' O3'.  This
    # is the whole point of the correction -- with a linker on each of two different
    # strands the tether bonds are bridges in the bond graph and no covalent cycle
    # exists, so the "loop" is held shut by base pairing alone and the construct can
    # swap arms with another chimera to form a dimer.
    PA = dres[(mod, 1)]["P"]["xyz"]                    # DNL.OL will bond this P
    PB = dres[(mod, n_bp)][PRIME3]["xyz"]              # DNP.P will bond this O3'
    D = float(np.linalg.norm(PA - PB))
    ab = float(np.linalg.norm(B - A))

    # Place the two anchor targets outboard of the two C-beta, both pushed the same
    # distance b along the outward normal so the arc bulges away from the barrel:
    #   tA = A - aA*u + b*n,  tB = B + aB*u + b*n
    # Requiring |tA-tB| = D gives aA + aB = D - ab, and requiring each arm to use the
    # same fraction of its own reach gives aA = sqrt(sA^2 - b^2), aB = sqrt(sB^2 - b^2).
    # The old code fixed aA = aB = (D-ab)/2, which is only right when the two arms are
    # the same length.  DNP is ~5 A longer than DNL, so solve for b instead.
    sA = args.extension_frac * reach[args.sites[0]]
    sB = args.extension_frac * reach[args.sites[1]]
    need = D - ab
    if sA + sB < need:
        sys.exit(f"the two tethers reach {sA:.2f} + {sB:.2f} = {sA + sB:.2f} A but the "
                 f"duplex anchors are {D:.2f} A apart and the sites {ab:.2f} A, so "
                 f"{need:.2f} A must be spanned: bend the duplex further")
    lo, hi = 0.0, min(sA, sB)
    for _ in range(100):                          # g(b) decreases monotonically
        b = 0.5 * (lo + hi)
        g = math.sqrt(max(0.0, sA ** 2 - b * b)) + math.sqrt(max(0.0, sB ** 2 - b * b))
        if g > need:
            lo = b
        else:
            hi = b
    b_off = 0.5 * (lo + hi)
    aA = math.sqrt(max(0.0, sA ** 2 - b_off ** 2))
    aB = math.sqrt(max(0.0, sB ** 2 - b_off ** 2))
    tA = A - aA * u + b_off * n
    tB = B + aB * u + b_off * n
    print(f"sites {ab:.2f} A apart; duplex anchors (5'P nt1 -> 3'O3' nt{n_bp}) "
          f"{D:.2f} A apart")
    print(f"  targets: a = {aA:.2f}/{aB:.2f}, b = {b_off:.2f} A; tether spans "
          f"{np.linalg.norm(tA - A):.2f} A ({100 * np.linalg.norm(tA - A) / reach[args.sites[0]]:.0f}% "
          f"of reach) and {np.linalg.norm(tB - B):.2f} A "
          f"({100 * np.linalg.norm(tB - B) / reach[args.sites[1]]:.0f}%)")

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
    # ideal bond lengths at the two junctions: DNL.OL-P is a phosphodiester ester bond,
    # DNP.P-O3' the same bond seen from the other side.  Both are 1.61 A (parm10 OS-P).
    for site, key, atom in ((args.sites[0], (mod, 1), "P"),
                            (args.sites[1], (mod, n_bp), PRIME3)):
        anchor = dres[key][atom]["xyz"]
        this = build_arm(LINKERS[site])
        rep = pose_arm(this, {k: pres[site][k]["xyz"] for k in ("N", "CA", "CB")},
                       anchor, ptree, rng, n_starts=16, bond_length=1.61)
        rep["bond_length_A"] = round(float(np.linalg.norm(
            this.xyz[this.index[this.dna_atom]] - anchor)), 3)
        rep["linker"] = LINKERS[site]
        arms[site] = this
        reports[str(site)] = rep
        print(f"tether at {site} ({LINKERS[site]}): to-target {rep['to_target_A']} A, "
              f"{this.dna_atom[:-1]}-{atom} bond {rep['bond_length_A']} A "
              f"(ideal 1.61), clash {rep['clash']}")

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

    dna_chain_ids = {mod: "C", comp: "D"}

    def write_linker(site, resi):
        """One linker residue, its atoms being the star-tagged half of the arm."""
        nonlocal serial, ordinal
        a = arms[site]
        ordinal += 1
        dnl_ord[site] = ordinal
        for nm in a.names:
            if not nm.endswith("*"):
                continue
            lines.append(fmt(serial, nm[:-1], LINKERS[site], "C", resi,
                             a.xyz[a.index[nm]], a.elements[a.index[nm]]))
            serial += 1

    # The modified strand, in tleap order.  DNL comes FIRST so that its tail oxygen
    # links to nucleotide 1's phosphorus by ordinary head-to-tail chaining, and DNP
    # comes LAST so nucleotide N's O3' links to its phosphorus the same way.  Writing
    # them inside the chain like this is also what keeps nucleotides 1 and N on the
    # *internal* templates -- they are no longer chain-initial or chain-final -- which
    # is exactly the charge bookkeeping DNL's -0.3079 and DNP's -0.6921 assume.
    write_linker(args.sites[0], 0)
    for i in range(1, n_bp + 1):
        ordinal += 1
        rn = next(iter(dres[(mod, i)].values()))["resn"]
        for nm, at in dres[(mod, i)].items():
            lines.append(fmt(serial, nm, rn, "C", i, at["xyz"], at["elem"]))
            serial += 1
    write_linker(args.sites[1], n_bp + 1)
    lines.append(f"TER   {serial:5d}")
    serial += 1

    # The complement carries no linker at all.  Its own 5'/3' termini were already
    # written as DA5/DA3-style names by build_dna.py, phosphate stripped.
    for i in range(1, n_bp + 1):
        ordinal += 1
        rn = next(iter(dres[(comp, i)].values()))["resn"]
        for nm, at in dres[(comp, i)].items():
            lines.append(fmt(serial, nm, rn, "D", i, at["xyz"], at["elem"]))
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

    # Only the two carbamates need declaring: they cross from the protein chain to the
    # DNA chain, so tleap cannot infer them.  The two DNA-side bonds (DNL.OL to
    # nucleotide 1's P, and nucleotide N's O3' to DNP.P) are ordinary head-to-tail
    # links within chain C and tleap makes them from the residue order alone.
    bonds = [f"{tdp_ord[s]}@CN:{dnl_ord[s]}@N" for s in args.sites]
    meta = {"protein": args.protein, "dna": args.dna, "sites": args.sites,
            "linkers": {str(k): v for k, v in LINKERS.items()},
            "tether": cage, "reach_A": {str(k): round(v, 2) for k, v in reach.items()},
            "sites_CB_CB_A": round(ab, 2), "duplex_anchor_span_A": round(D, 2),
            "a_offsets_A": [round(aA, 2), round(aB, 2)], "b_offset_A": round(b_off, 2),
            "duplex_spin_deg": spin, "arm_reports": reports,
            "waters_kept": kept, "waters_dropped": len(wres) - kept,
            "tleap_bonds": bonds}
    Path(args.out).with_suffix(".json").write_text(json.dumps(meta, indent=2))
    print(f"\nwrote {args.out}  ({serial - 1} atoms, {kept} waters kept, "
          f"{len(wres) - kept} dropped)")
    print(f"tleap bonds to declare: {' '.join(bonds)}")


if __name__ == "__main__":
    main()
