#!/usr/bin/env python
"""Plan section 5F - Fig 3 geometry from the deposited models 9MIP (Aalpha) and 9MF5 (Abeta).

Reproduces, with the paper's numbers alongside:
  Fig 3A  HEAT1-HEAT15 centre-of-mass distance (paper 62.7 A Aalpha, 59.6 A Abeta)
  Fig 3B  per-HEAT-repeat CoM shift after superposing on HEAT 12-15 (paper: largest at HEAT 3-6, ~2-3 A)
  Fig 3C  HEAT3-HEAT5 salt-bridge backbone distance R113-E191 (Aalpha) vs R125-E203 (Abeta) (paper 13.2 / 12.4 A)
  Fig 3D/E interface residues of the A subunit within 4 A of B56g and of Calpha (paper lists 11 conserved pairs)
Also records what is needed for the Fig 2B test: is 9MIP chain B 404-449 (absent from 9MF5) helical, and
how far it sits from the A N-terminus.

HEAT repeat boundaries: two definitions are computed and both reported -
  (i) DSSP-derived: helices >= 7 residues on chain A, paired consecutively into 15 (A,B) helix pairs;
  (ii) nominal 39-residue grid  HEAT k = 8+39(k-1) ... 8+39k-1  (Aalpha numbering; Abeta = +12),
      anchored by R113 in HEAT3, E191 in HEAT5, W257 in HEAT7, K416/R418 in HEAT11, D572/D574 in HEAT15.
Outputs: $PP2A_DATA/em/geometry/geometry.json and $PP2A_REPO/results/geometry_metrics.csv
(rows in the metrics.csv schema: track,system,quantity,selection_or_method,value,unit,paper_value,notes).
"""
import csv
import json
import os
import sys

import mdtraj as md
import numpy as np

REPO = os.environ.get("PP2A_REPO", "/n/home05/ianpoon/DNA-springs-project-2026-27/cryoEM/pp2a_day2025_reproduction")
DATA = os.environ.get("PP2A_DATA", "/n/netscratch/hekstra_lab/Lab/ian_poon/pp2a_day2025")
MODELS = os.path.join(DATA, "em", "models")
OUT = os.path.join(DATA, "em", "geometry")
RES_DIR = os.path.join(REPO, "results")

SYSTEMS = {"Aa": {"pdb": "9MIP", "offset": 0}, "Ab": {"pdb": "9MF5", "offset": 12}}
PAPER = {
    "heat1_heat15_com": {"Aa": 62.7, "Ab": 59.6},
    "salt_bridge_ca": {"Aa": 13.2, "Ab": 12.4},
    "salt_bridge_res": {"Aa": (113, 191), "Ab": (125, 203)},
    # Fig 3D/E residue lists in Aalpha numbering (Abeta = +12)
    "iface_B56": [100, 101, 183, 257],
    "iface_C": [416, 418, 456, 495, 498, 572, 574],
}
NOMINAL_START = 8
HEAT_LEN = 39


def kabsch(P, Q):
    """Rotation R, translation t minimising |R P + t - Q| (rows = points)."""
    Pc, Qc = P.mean(0), Q.mean(0)
    H = (P - Pc).T @ (Q - Qc)
    U, S, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    D = np.diag([1, 1, d])
    R = Vt.T @ D @ U.T
    t = Qc - R @ Pc
    return R, t


def chain_by_id(top, cid):
    for ch in top.chains:
        if ch.chain_id == cid:
            return ch
    raise KeyError(cid)


def ca_xyz(traj, chain, resseqs):
    """Cα coordinates (Å) for the given author residue numbers; missing residues skipped."""
    out, found = [], []
    lut = {r.resSeq: r for r in chain.residues}
    for rs in resseqs:
        r = lut.get(rs)
        if r is None:
            continue
        ca = [a.index for a in r.atoms if a.name == "CA"]
        if ca:
            out.append(traj.xyz[0, ca[0]] * 10.0)
            found.append(rs)
    return np.array(out), found


def dssp_helices(traj, chain, min_len=7):
    sel = traj.topology.select(f"chainid {chain.index}")
    sub = traj.atom_slice(sel)
    ss = md.compute_dssp(sub, simplified=True)[0]
    resseq = [r.resSeq for r in sub.topology.residues]
    segs, start = [], None
    for i, s in enumerate(ss + ["X"]):
        if s == "H" and start is None:
            start = i
        elif s != "H" and start is not None:
            if i - start >= min_len:
                segs.append((resseq[start], resseq[i - 1]))
            start = None
    return segs


def heat_ranges_nominal(offset):
    return {k: (NOMINAL_START + HEAT_LEN * (k - 1) + offset, NOMINAL_START + HEAT_LEN * k - 1 + offset)
            for k in range(1, 16)}


def heat_ranges_dssp(segs, offset):
    """Pair consecutive helices into HEAT repeats; returns None if the count is not 30."""
    # keep helices inside the solenoid (drop anything N-terminal of the first repeat start)
    segs = [s for s in segs if s[1] >= NOMINAL_START + offset]
    if len(segs) < 30:
        return None, len(segs)
    segs = segs[:30]
    return {k: (segs[2 * (k - 1)][0], segs[2 * (k - 1) + 1][1]) for k in range(1, 16)}, len(segs)


def com_per_heat(traj, chain, ranges):
    coms = {}
    for k, (a, b) in ranges.items():
        xyz, found = ca_xyz(traj, chain, range(a, b + 1))
        coms[k] = (xyz.mean(0), len(found))
    return coms


def interface_residues(traj, chainA, chainX, cutoff=4.0):
    top = traj.topology
    A = top.select(f"chainid {chainA.index} and not element H")
    X = top.select(f"chainid {chainX.index} and not element H")
    xa, xx = traj.xyz[0, A] * 10, traj.xyz[0, X] * 10
    # brute force is fine: ~5k x ~3k
    d2 = ((xa[:, None, :] - xx[None, :, :]) ** 2).sum(-1)
    close = (d2.min(1) <= cutoff ** 2)
    res = sorted({top.atom(int(a)).residue.resSeq for a in A[close]})
    return res


def main():
    os.makedirs(OUT, exist_ok=True)
    os.makedirs(RES_DIR, exist_ok=True)
    result, rows = {}, []

    def row(system, quantity, method, value, unit, paper="", notes=""):
        rows.append(["em", system, quantity, method, f"{value:.3f}" if isinstance(value, float) else value,
                     unit, paper, notes])

    trajs, chains = {}, {}
    for sysname, meta in SYSTEMS.items():
        t = md.load(os.path.join(MODELS, meta["pdb"] + ".pdb"))
        trajs[sysname] = t
        chains[sysname] = {cid: chain_by_id(t.topology, cid) for cid in "ABC"}
        A = chains[sysname]["A"]
        rs = [r.resSeq for r in A.residues]
        result[sysname] = {"pdb": meta["pdb"], "chainA_range": (min(rs), max(rs)), "n_res_A": len(rs)}

    # ---- HEAT boundaries and centres of mass
    heat = {}
    for sysname, meta in SYSTEMS.items():
        t, A = trajs[sysname], chains[sysname]["A"]
        off = meta["offset"]
        segs = dssp_helices(t, A)
        dssp_r, nseg = heat_ranges_dssp(segs, off)
        nom_r = heat_ranges_nominal(off)
        result[sysname]["dssp_helices"] = segs
        result[sysname]["dssp_helix_count"] = nseg
        result[sysname]["heat_nominal"] = nom_r
        result[sysname]["heat_dssp"] = dssp_r
        heat[sysname] = {"nominal": com_per_heat(t, A, nom_r)}
        if dssp_r:
            heat[sysname]["dssp"] = com_per_heat(t, A, dssp_r)
        # sanity anchors (Aalpha numbering)
        anchors = {"R113": 113, "E191": 191, "W257": 257, "K416": 416, "R418": 418, "D572": 572, "D574": 574}
        lut = {r.resSeq: r.name for r in A.residues}
        result[sysname]["anchors"] = {k: lut.get(v + off) for k, v in anchors.items()}

    # ---- Fig 3A: HEAT1-HEAT15 CoM distance (intramolecular)
    for sysname in SYSTEMS:
        for method, coms in heat[sysname].items():
            d = float(np.linalg.norm(coms[1][0] - coms[15][0]))
            result[sysname][f"heat1_heat15_com_{method}"] = d
            row(sysname, "HEAT1-HEAT15 CoM distance", f"Calpha CoM, {method} HEAT bounds", d, "A",
                PAPER["heat1_heat15_com"][sysname], "Fig 3A")

    # ---- Fig 3B: per-HEAT CoM shift after superposing Abeta on Aalpha using HEAT 12-15 Calpha
    for method in ("nominal", "dssp"):
        if method not in heat["Aa"] or method not in heat["Ab"]:
            continue
        ra, rb = result["Aa"][f"heat_{method}"], result["Ab"][f"heat_{method}"]
        # matched Calpha: same Aalpha-numbered residue present in both
        resA = [r for k in range(12, 16) for r in range(ra[k][0], ra[k][1] + 1)]
        xa, fa = ca_xyz(trajs["Aa"], chains["Aa"]["A"], resA)
        xb, fb = ca_xyz(trajs["Ab"], chains["Ab"]["A"], [r + 12 for r in resA])
        common = sorted(set(fa) & set(r - 12 for r in fb))
        xa = np.array([xa[fa.index(r)] for r in common]); xb = np.array([xb[fb.index(r + 12)] for r in common])
        R, tvec = kabsch(xb, xa)
        fit_rmsd = float(np.sqrt((((xb @ R.T + tvec) - xa) ** 2).sum(1).mean()))
        shifts = {}
        for k in range(1, 16):
            cb = heat["Ab"][method][k][0] @ R.T + tvec
            shifts[k] = float(np.linalg.norm(cb - heat["Aa"][method][k][0]))
            row("Ab-vs-Aa", f"HEAT{k} CoM shift", f"after HEAT12-15 Calpha fit ({method}), n_fit={len(common)}",
                shifts[k], "A", "2-3 (HEAT 3-6 largest)" if 3 <= k <= 6 else "", "Fig 3B")
        result[f"com_shift_{method}"] = {"fit_rmsd": fit_rmsd, "n_fit_atoms": len(common), "shift": shifts,
                                          "largest": sorted(shifts, key=shifts.get, reverse=True)[:4]}

    # ---- Fig 3C: salt bridge R113-E191 / R125-E203, Calpha-Calpha and closest sidechain N-O
    for sysname in SYSTEMS:
        t, A = trajs[sysname], chains[sysname]["A"]
        r1, r2 = PAPER["salt_bridge_res"][sysname]
        lut = {r.resSeq: r for r in A.residues}
        R_, E_ = lut[r1], lut[r2]
        assert R_.name == "ARG" and E_.name == "GLU", (R_.name, E_.name)
        ca1 = [a.index for a in R_.atoms if a.name == "CA"][0]; ca2 = [a.index for a in E_.atoms if a.name == "CA"][0]
        dca = float(np.linalg.norm(t.xyz[0, ca1] - t.xyz[0, ca2]) * 10)
        n_at = [a.index for a in R_.atoms if a.name in ("NE", "NH1", "NH2")]
        o_at = [a.index for a in E_.atoms if a.name in ("OE1", "OE2")]
        dno = float(min(np.linalg.norm(t.xyz[0, i] - t.xyz[0, j]) * 10 for i in n_at for j in o_at)) if n_at and o_at else float("nan")
        result[sysname]["salt_bridge"] = {"pair": f"R{r1}-E{r2}", "ca_ca": dca, "min_N_O": dno}
        row(sysname, f"salt bridge R{r1}-E{r2} backbone distance", "Calpha-Calpha", dca, "A",
            PAPER["salt_bridge_ca"][sysname], "Fig 3C")
        row(sysname, f"salt bridge R{r1}-E{r2} side-chain distance", "min Arg N - Glu O", dno, "A", "", "Fig 3C")

    # ---- Fig 3D/E: interface residues
    for sysname, meta in SYSTEMS.items():
        t, ch = trajs[sysname], chains[sysname]
        off = meta["offset"]
        for partner, key in (("B", "iface_B56"), ("C", "iface_C")):
            res = interface_residues(t, ch["A"], ch[partner])
            expected = [r + off for r in PAPER[key]]
            hits = [r for r in expected if r in res]
            lut = {r.resSeq: r.name for r in ch["A"].residues}
            result[sysname][f"interface_{partner}"] = {"n_residues": len(res), "residues": res,
                                                       "paper_residues": {f"{lut.get(r, '?')}{r}": (r in res) for r in expected}}
            row(sysname, f"A-{'B56g' if partner == 'B' else 'Calpha'} interface residues within 4 A", "heavy atoms",
                len(res), "count", "", f"paper's {len(expected)} residues found: {len(hits)}/{len(expected)}  (Fig 3{'D' if partner == 'B' else 'E'})")

    # ---- Fig 2B support: 9MIP chain B 404-449 (not in 9MF5) - helical? where relative to A N-terminus?
    tA, chA = trajs["Aa"], chains["Aa"]
    segsB = dssp_helices(tA, chA["B"], min_len=5)
    tail = [s for s in segsB if s[1] >= 404]
    xb, fb = ca_xyz(tA, chA["B"], range(404, 450))
    xn, fn = ca_xyz(tA, chA["A"], range(7, 21))
    dmin = float(np.min(np.linalg.norm(xb[:, None] - xn[None], axis=-1))) if len(xb) and len(xn) else float("nan")
    result["fig2b_helix_candidate"] = {"chainB_range_9MIP": [min(r.resSeq for r in chA["B"].residues), max(r.resSeq for r in chA["B"].residues)],
                                       "chainB_range_9MF5": [min(r.resSeq for r in chains["Ab"]["B"].residues), max(r.resSeq for r in chains["Ab"]["B"].residues)],
                                       "helices_in_404_449": tail, "n_CA_404_449": len(fb),
                                       "min_CA_distance_to_A_res7_20": dmin}
    row("Aa", "B56g 404-449 min Calpha distance to A N-terminus (7-20)", "9MIP", dmin, "A", "", "Fig 2B region check")

    json.dump(result, open(os.path.join(OUT, "geometry.json"), "w"), indent=2, default=lambda o: list(o) if hasattr(o, "__iter__") else str(o))
    with open(os.path.join(RES_DIR, "geometry_metrics.csv"), "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["track", "system", "quantity", "selection_or_method", "value", "unit", "paper_value", "notes"])
        w.writerows(rows)

    # console summary
    for s in SYSTEMS:
        r = result[s]
        print(f"{s} ({r['pdb']}): chain A {r['chainA_range']}, DSSP helices>=7: {r['dssp_helix_count']}, anchors {r['anchors']}")
        for m in ("nominal", "dssp"):
            k = f"heat1_heat15_com_{m}"
            if k in r:
                print(f"   HEAT1-15 CoM [{m}]: {r[k]:.1f} A   (paper {PAPER['heat1_heat15_com'][s]})")
        print(f"   salt bridge {r['salt_bridge']['pair']}: Ca-Ca {r['salt_bridge']['ca_ca']:.1f} A (paper {PAPER['salt_bridge_ca'][s]}), min N-O {r['salt_bridge']['min_N_O']:.1f} A")
        for p in "BC":
            print(f"   interface A-{p}: {r[f'interface_{p}']['n_residues']} residues; paper's: {r[f'interface_{p}']['paper_residues']}")
    for m in ("nominal", "dssp"):
        k = f"com_shift_{m}"
        if k in result:
            sh = result[k]["shift"]
            print(f"CoM shift [{m}] fit RMSD {result[k]['fit_rmsd']:.2f} A on {result[k]['n_fit_atoms']} CA; per HEAT: " +
                  " ".join(f"{i}:{sh[i]:.1f}" for i in range(1, 16)) + f"  largest {result[k]['largest']}")
    print("Fig 2B candidate:", result["fig2b_helix_candidate"])


if __name__ == "__main__":
    sys.exit(main())
