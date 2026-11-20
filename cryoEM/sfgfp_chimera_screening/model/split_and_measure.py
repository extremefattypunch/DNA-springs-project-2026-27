#!/usr/bin/env python3
"""Split the assembled L6_loop50 chimera into solute / protein / DNA PDBs and measure it.

Input : L6_loop50/chimera_loop50.pdb (+ .json) from build_L6_loop50.sh, and optionally
        L6_loop50/tleap/build_report.json + system.prmtop (covalent-loop assertion, full-H MW).
Output: L6_loop50_solute.pdb, L6_protein.pdb, L6_dna.pdb, geometry.json in L6_loop50/.

Protein file = chain A (sfGFP incl. the two TDP click-adduct residues) + the DNL/DNP linker
residues (chain C resi 0 and N+1).  DNA file = nucleotides only (chain C 1..N, chain D).
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

D = Path(__file__).resolve().parent / "L6_loop50"
SRC = D / "chimera_loop50.pdb"
SITES = (133, 149)
LINKERS = {"DNL", "DNP"}
WATER = {"WAT", "HOH"}
MASS = {"H": 1.008, "C": 12.011, "N": 14.007, "O": 15.999, "P": 30.974, "S": 32.06}
BASE = {"N1", "C2", "N3", "C4", "C5", "C6", "N7", "C8", "N9", "O2", "O4", "O6", "N2",
        "N4", "N6", "C7"}


def parse(path):
    recs = []
    for l in path.read_text().splitlines():
        if l[:6] in ("ATOM  ", "HETATM"):
            el = l[76:78].strip() or l[12:16].strip()[0]
            recs.append(dict(line=l, name=l[12:16].strip(), resn=l[17:20].strip(),
                             chain=l[21], resi=int(l[22:26]), el=el,
                             xyz=np.array([float(l[30:38]), float(l[38:46]),
                                           float(l[46:54])])))
        elif l.startswith("TER"):
            recs.append(dict(line="TER"))
    return recs


def write(path, recs, keep, remark):
    out, prev_ter = [f"REMARK   {remark}"], True
    for r in recs:
        if r["line"] == "TER":
            if not prev_ter:
                out.append("TER")
                prev_ter = True
        elif keep(r):
            out.append(r["line"])
            prev_ter = False
    if not prev_ter:
        out.append("TER")
    path.write_text("\n".join(out) + "\nEND\n")


def rg(xyz, m):
    c = (m[:, None] * xyz).sum(0) / m.sum()
    return float(np.sqrt((m * ((xyz - c) ** 2).sum(1)).sum() / m.sum()))


def max_dim(xyz):
    from scipy.spatial import ConvexHull
    h = xyz[ConvexHull(xyz).vertices]
    d = np.linalg.norm(h[:, None] - h[None], axis=-1)
    return float(d.max())


def circle_fit(pts):
    c0 = pts.mean(0)
    _, s, vt = np.linalg.svd(pts - c0)
    e1, e2, nrm = vt
    p2 = np.c_[(pts - c0) @ e1, (pts - c0) @ e2]
    A = np.c_[2 * p2, np.ones(len(p2))]
    sol, *_ = np.linalg.lstsq(A, (p2 ** 2).sum(1), rcond=None)
    cx, cy, k = sol
    R = float(np.sqrt(k + cx ** 2 + cy ** 2))
    centre = c0 + cx * e1 + cy * e2
    resid = np.linalg.norm(p2 - [cx, cy], axis=1) - R
    planarity = float(s[2] / np.sqrt(len(pts)))      # rms out-of-plane
    ang = np.unwrap(np.arctan2(p2[:, 1] - cy, p2[:, 0] - cx))
    return R, centre, float(np.sqrt((resid ** 2).mean())), planarity, \
        float(np.degrees(abs(ang[-1] - ang[0])))


def main():
    recs = parse(SRC)
    atoms = [r for r in recs if r["line"] != "TER"]
    solute = [a for a in atoms if a["resn"] not in WATER]
    is_prot = lambda a: a["chain"] == "A" or a["resn"] in LINKERS
    is_dna = lambda a: a["chain"] in "CD" and a["resn"] not in LINKERS
    remark = "L6_loop50 sfGFP(133/149)-50bp closed-loop DNA chimera; ideal starting model " \
             "for fitting/validation ONLY - never a picking/refinement reference"
    write(D / "L6_loop50_solute.pdb", recs, lambda a: a["resn"] not in WATER, remark)
    write(D / "L6_protein.pdb", recs, is_prot, remark + " | protein+TDP adducts+DNL/DNP linkers")
    write(D / "L6_dna.pdb", recs, is_dna, remark + " | DNA duplex only")

    def xyz(sel):
        return np.array([a["xyz"] for a in sel])

    def mass(sel):
        return np.array([MASS[a["el"].capitalize()] for a in sel])

    idx = {(a["chain"], a["resi"], a["name"]): a["xyz"] for a in atoms}
    prot = [a for a in solute if is_prot(a)]
    dna = [a for a in solute if is_dna(a)]
    n_bp = max(a["resi"] for a in dna if a["chain"] == "C")

    cbcb = float(np.linalg.norm(idx[("A", SITES[0], "CB")] - idx[("A", SITES[1], "CB")]))
    anchor = float(np.linalg.norm(idx[("C", 1, "P")] - idx[("C", n_bp, "O3'")]))
    # base-pair centres (C i pairs D n+1-i), then a one-turn running mean = helix axis
    bpc = []
    for i in range(1, n_bp + 1):
        sel = [a["xyz"] for a in dna if a["name"] in BASE and
               ((a["chain"] == "C" and a["resi"] == i) or
                (a["chain"] == "D" and a["resi"] == n_bp + 1 - i))]
        bpc.append(np.mean(sel, 0))
    bpc = np.array(bpc)
    w = 10
    axis = np.array([bpc[i:i + w].mean(0) for i in range(n_bp - w + 1)])
    contour_axis = float(np.linalg.norm(np.diff(axis, axis=0), axis=1).sum())
    R, centre, rms, planar, arc_axis_deg = circle_fit(axis)
    R_bp, _, rms_bp, _, _ = circle_fit(bpc)
    contour_full = (n_bp - 1) * 3.375                   # ideal B rise along the arc
    bend_deg = float(np.degrees(contour_full / R))
    t0 = axis[3] - axis[0]
    t1 = axis[-1] - axis[-4]
    tangent_deg = float(np.degrees(np.arccos(t0 @ t1 / np.linalg.norm(t0) /
                                             np.linalg.norm(t1))))
    e2e_bp = float(np.linalg.norm(bpc[-1] - bpc[0]))
    chord_from_R = float(2 * R * np.sin(contour_full / (2 * R)))

    # protein-DNA contacts (heavy atoms), excluding the bonded ends (nt 1 and nt N, linkers)
    from scipy.spatial import cKDTree
    pa = [a for a in prot if a["el"] != "H" and a["chain"] == "A"]
    da = [a for a in dna if a["el"] != "H"]
    dt = cKDTree(xyz(da))
    dmin, _ = dt.query(xyz(pa), k=1)
    clash25 = int((dmin < 2.5).sum())
    dmin_all = float(dmin.min())
    clash30 = int((dmin < 3.0).sum())
    allh = [a for a in solute if a["el"] != "H"]
    st = cKDTree(xyz(allh))
    pairs = st.query_pairs(2.2)
    # non-bonded close pairs = within 2.2 A but in residues more than one apart in sequence
    def rk(a):
        return (a["chain"], a["resi"])
    # TDP149 vs Tyr151/His199 (1.9-2.2 A) is inherited unchanged from L4_loop40 (same
    # seed, same tether poses); minimisation relieves it.
    # expected covalent links between non-sequence-adjacent residues (bond <= 1.85 A; the CRO peptide links are built long, ~1.8 A):
    # CRO-64/68, TDP-DNL/DNP carbamates.  Anything else under 2.2 A is a clash.
    close = []
    for i, j in pairs:
        a, b = allh[i], allh[j]
        if rk(a) == rk(b) or (a["chain"] == b["chain"] and abs(a["resi"] - b["resi"]) == 1):
            continue
        d = float(np.linalg.norm(a["xyz"] - b["xyz"]))
        close.append((d, f"{a['chain']}{a['resi']}{a['resn']}.{a['name']}",
                      f"{b['chain']}{b['resi']}{b['resn']}.{b['name']}"))
    bonded = [c for c in close if c[0] <= 1.85]
    nb_close = [c for c in close if c[0] > 1.85]

    junction = {
        "DNL.OL-P(nt1)": float(np.linalg.norm(idx[("C", 0, "OL")] - idx[("C", 1, "P")])),
        f"DNP.P-O3'(nt{n_bp})": float(np.linalg.norm(idx[("C", n_bp + 1, "P")]
                                                      - idx[("C", n_bp, "O3'")])),
    }

    meta = json.loads((D / "chimera_loop50.json").read_text())
    geo = {
        "system": "L6_loop50",
        "description": "sfGFP (2B3P numbering, D133/N149 -> Tet2-Et) clicked to a 50 bp "
                       "duplex modified on ONE strand at both ends (5' DNL, 3' DNP); "
                       "ideal B-DNA bent onto a circular arc; as-built, unminimised",
        "use": "rigid-body fitting / validation only; NOT a cryo-EM picking or refinement "
               "reference",
        "n_bp": n_bp, "sequence_modified_strand":
            json.loads((D / "dna" / "spring_50bp.json").read_text())["sequence"],
        "target_anchor_span_A": 66.0,
        "anchor_CB133_CB149_A": round(cbcb, 2),
        "dna_end_to_end": {
            "anchor_span_P1_to_O3prime50_A": round(anchor, 2),
            "basepair_centre_bp1_to_bp50_A": round(e2e_bp, 2),
            "arc_chord_from_fit_A": round(chord_from_R, 2),
            "note": "anchor span = the two atoms the linkers bond (same strand); arc chord "
                    "= 2R sin(L/2R) for the fitted axis circle and the 49-step B-DNA contour",
        },
        "arc": {
            "radius_fit_axis_A": round(R, 2),
            "radius_fit_bp_centres_A": round(R_bp, 2),
            "radius_as_built_A": json.loads((D / "dna" / "spring_50bp.json").read_text())
                                 ["bend"]["radius_A"],
            "circle_fit_rms_A": round(rms, 2),
            "out_of_plane_rms_A": round(planar, 2),
            "contour_length_A": round(contour_full, 1),
            "bend_angle_contour_over_R_deg": round(bend_deg, 1),
            "end_tangent_angle_deg": round(tangent_deg, 1),
            "deg_per_bp_step": round(bend_deg / (n_bp - 1), 2),
        },
        "max_dimension_A": {"solute": round(max_dim(xyz(solute)), 1),
                            "protein": round(max_dim(xyz(prot)), 1),
                            "dna": round(max_dim(xyz(dna)), 1)},
        "radius_of_gyration_A": {"solute": round(rg(xyz(solute), mass(solute)), 2),
                                 "protein": round(rg(xyz(prot), mass(prot)), 2),
                                 "dna": round(rg(xyz(dna), mass(dna)), 2),
                                 "weighting": "atomic mass, atoms present in the PDB "
                                              "(heavy atoms + linker H)"},
        "atoms_in_pdb": {"solute": len(solute), "protein": len(prot), "dna": len(dna),
                         "hydrogens": sum(a["el"] == "H" for a in solute)},
        "covalent_junctions_A": {k: round(v, 3) for k, v in junction.items()},
        "tether_reports": meta["arm_reports"],
        "duplex_spin_deg": meta["duplex_spin_deg"],
        "protein_dna_contacts": {"min_sfGFP_to_DNA_heavy_atom_A": round(dmin_all, 2),
                                 "protein_heavy_atoms_within_2.5A_of_DNA": clash25,
                                 "within_3.0A": clash30,
                                 "nonbonded_heavy_pairs_1.85_to_2.2A": len(nb_close),
                                 "nonbonded_pairs": [[round(d, 2), x, y]
                                                     for d, x, y in sorted(nb_close)],
                                 "cross_residue_bonds_under_1.85A":
                                     [[round(d, 2), x, y] for d, x, y in sorted(bonded)]},
    }

    # molecular weight: prefer full-H topology from tleap, fall back to heavy-atom + H estimate
    rep_p = D / "tleap" / "build_report.json"
    prm = D / "tleap" / "system.prmtop"
    mw = {}
    if prm.exists():
        import parmed
        p = parmed.load_file(str(prm))
        sol = [r for r in p.residues if r.name not in ("WAT", "Na+", "Cl-", "MG")]
        dna_names = {"DA", "DT", "DG", "DC", "DA5", "DT5", "DG5", "DC5", "DA3", "DT3",
                     "DG3", "DC3"}
        m_d = sum(a.mass for r in sol if r.name in dna_names for a in r.atoms)
        m_p = sum(a.mass for r in sol if r.name not in dna_names for a in r.atoms)
        mw = {"source": "tleap prmtop (all H, Amber protonation)",
              "protein_plus_adducts_and_linkers_kDa": round(m_p / 1000, 2),
              "dna_kDa": round(m_d / 1000, 2), "total_kDa": round((m_p + m_d) / 1000, 2)}
    else:
        mw = {"source": "PDB heavy atoms only (no H) -- underestimate",
              "protein_plus_adducts_and_linkers_kDa": round(mass(prot).sum() / 1000, 2),
              "dna_kDa": round(mass(dna).sum() / 1000, 2)}
    geo["molecular_weight"] = mw
    if rep_p.exists():
        r = json.loads(rep_p.read_text())
        loop = r.get("covalent_loop", {})
        geo["covalent_loop"] = {"closed": bool(loop.get("closed")),
                                "source": "build_system.py --expect-loop TDP (bond-graph "
                                          "cycle test on the tleap prmtop)",
                                "detail": loop, "tleap_problems": r.get("problems"),
                                "bonds_verified": r.get("bonds_verified"),
                                "net_charge_solvated": r.get("net_charge"),
                                "atoms_solvated": r.get("atoms")}
    else:
        geo["covalent_loop"] = {"closed": None, "source": "tleap not run",
                                "junction_distances_A": junction}
    (D / "geometry.json").write_text(json.dumps(geo, indent=2))
    print(json.dumps({k: geo[k] for k in ("anchor_CB133_CB149_A", "dna_end_to_end", "arc",
                                          "max_dimension_A", "radius_of_gyration_A",
                                          "molecular_weight", "protein_dna_contacts",
                                          "covalent_junctions_A")}, indent=1))
    print("covalent loop closed:", geo["covalent_loop"]["closed"])


if __name__ == "__main__":
    sys.exit(main())
