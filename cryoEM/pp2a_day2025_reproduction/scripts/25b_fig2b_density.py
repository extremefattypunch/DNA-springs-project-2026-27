#!/usr/bin/env python
"""Plan section 5F - quantify the Fig 2B claim: at 6 A the B56g C-terminal helix is visible in the
Aalpha map and absent from the Abeta map.

Region: 9MIP chain B residues 402-436 (helical in 9MIP, sitting 6.9 A from the A N-terminus; not
modelled in 9MF5, which stops at B 403).  For the Abeta map the helix atoms are carried into the 9MF5
frame by a Kabsch fit of the B56g core (9MIP B 38-403 onto 9MF5 B 38-403).
Readout in each 6 A low-passed raw map (RELION relion_image_handler --lowpass 6 of the half-map average;
also the low-passed deposited sharpened map for comparison):
    r = <rho>(helix atoms, 3 A dilation) / <rho>(B56g core atoms 38-403, 3 A dilation)
plus the same ratio for two control regions (A-subunit HEAT 7 and the C-subunit core) and the fraction
of helix voxels above the map's B56g-core median.  Pre-registered expectation: r_Aa >= 0.6, r_Ab <= 0.25.
Outputs $PP2A_DATA/em/geometry/fig2b_density.json and rows in results/fig2b_metrics.csv.
"""
import csv
import json
import os
import sys

import mdtraj as md
import mrcfile
import numpy as np
from scipy.ndimage import binary_dilation, generate_binary_structure

REPO = os.environ.get("PP2A_REPO", "/n/home05/ianpoon/DNA-springs-project-2026-27/cryoEM/pp2a_day2025_reproduction")
DATA = os.environ.get("PP2A_DATA", "/n/netscratch/hekstra_lab/Lab/ian_poon/pp2a_day2025")
M = os.path.join(DATA, "em", "models")
R = os.path.join(DATA, "em", "relion")
OUT = os.path.join(DATA, "em", "geometry")
HELIX = (402, 436)
CORE = (38, 403)
ISO = {"Aa": dict(pdb="9MIP", off=0), "Ab": dict(pdb="9MF5", off=12)}


def kabsch(P, Q):
    Pc, Qc = P.mean(0), Q.mean(0)
    U, S, Vt = np.linalg.svd((P - Pc).T @ (Q - Qc))
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    Rm = Vt.T @ np.diag([1, 1, d]) @ U.T
    return Rm, Qc - Rm @ Pc


def atoms(traj, chain_id, lo, hi, ca_only=False):
    top = traj.topology
    sel = [a.index for a in top.atoms if a.residue.chain.chain_id == chain_id and lo <= a.residue.resSeq <= hi
           and a.element.symbol != "H" and (not ca_only or a.name == "CA")]
    return traj.xyz[0, sel] * 10.0, [top.atom(i).residue.resSeq for i in sel]


def open_map(path):
    m = mrcfile.open(path, permissive=True)
    data = np.asarray(m.data, dtype=np.float32)          # axes: z, y, x (mrcfile default for mapc=1,mapr=2,maps=3)
    vs = float(m.voxel_size.x)
    origin = np.array([float(m.header.origin.x), float(m.header.origin.y), float(m.header.origin.z)])
    nstart = np.array([int(m.header.nxstart), int(m.header.nystart), int(m.header.nzstart)])
    if np.allclose(origin, 0) and np.any(nstart != 0):
        origin = nstart * vs
    return data, vs, origin


def mask_from_atoms(xyz, shape, vs, origin, dilate_A=3.0):
    ijk = np.round((xyz - origin) / vs).astype(int)          # x, y, z voxel indices
    mask = np.zeros(shape, dtype=bool)
    ok = np.all((ijk >= 0) & (ijk < np.array(shape)[::-1]), axis=1)
    ijk = ijk[ok]
    mask[ijk[:, 2], ijk[:, 1], ijk[:, 0]] = True
    n = int(round(dilate_A / vs))
    if n > 0:
        mask = binary_dilation(mask, structure=generate_binary_structure(3, 1), iterations=n)
    return mask, ok.mean()


def main():
    os.makedirs(OUT, exist_ok=True)
    os.makedirs(os.path.join(REPO, "results"), exist_ok=True)
    t_a = md.load(os.path.join(M, "9MIP.pdb"))
    t_b = md.load(os.path.join(M, "9MF5.pdb"))
    helix_xyz, _ = atoms(t_a, "B", *HELIX)
    result, rows = {}, []
    for iso, meta in ISO.items():
        traj = t_a if iso == "Aa" else t_b
        # helix coordinates in this model's frame
        if iso == "Aa":
            hx = helix_xyz
            fit_rmsd = 0.0
        else:
            ca_a, ra = atoms(t_a, "B", *CORE, ca_only=True)
            ca_b, rb = atoms(t_b, "B", *CORE, ca_only=True)
            common = sorted(set(ra) & set(rb))
            A = np.array([ca_a[ra.index(r)] for r in common])
            B = np.array([ca_b[rb.index(r)] for r in common])
            Rm, tv = kabsch(A, B)
            fit_rmsd = float(np.sqrt((((A @ Rm.T + tv) - B) ** 2).sum(1).mean()))
            hx = helix_xyz @ Rm.T + tv
        core_xyz, _ = atoms(traj, "B", *CORE)
        ctrl_heat7, _ = atoms(traj, "A", 242 + meta["off"], 280 + meta["off"])
        ctrl_c, _ = atoms(traj, "C", 30, 280)
        result[iso] = {"b56_core_fit_rmsd_A": fit_rmsd, "maps": {}}
        for tag in ("raw_lp6", "sharp_lp6", "raw"):
            path = os.path.join(R, iso, f"{tag}.mrc")
            if not os.path.exists(path):
                continue
            data, vs, origin = open_map(path)
            m_hel, cov = mask_from_atoms(hx, data.shape, vs, origin)
            m_core, _ = mask_from_atoms(core_xyz, data.shape, vs, origin)
            m_h7, _ = mask_from_atoms(ctrl_heat7, data.shape, vs, origin)
            m_c, _ = mask_from_atoms(ctrl_c, data.shape, vs, origin)
            m_core_only = m_core & ~m_hel
            solvent = data[~binary_dilation(m_core | m_hel | m_h7 | m_c, iterations=6)]
            mu0, sd0 = float(solvent.mean()), float(solvent.std())
            core_med = float(np.median(data[m_core_only]))
            r = float((data[m_hel].mean() - mu0) / (data[m_core_only].mean() - mu0))
            res = {"ratio_helix_over_B56core": r,
                   "ratio_HEAT7_over_B56core": float((data[m_h7].mean() - mu0) / (data[m_core_only].mean() - mu0)),
                   "ratio_Ccore_over_B56core": float((data[m_c].mean() - mu0) / (data[m_core_only].mean() - mu0)),
                   "helix_frac_above_core_median": float((data[m_hel] > core_med).mean()),
                   "helix_mean_sigma": float((data[m_hel].mean() - mu0) / sd0),
                   "B56core_mean_sigma": float((data[m_core_only].mean() - mu0) / sd0),
                   "n_helix_voxels": int(m_hel.sum()), "atom_coverage": float(cov), "voxel_A": vs}
            result[iso]["maps"][tag] = res
            rows.append(["em", iso, "Fig2B helix/B56-core density ratio", f"{tag} map, 9MIP B402-436 mask 3 A", f"{r:.3f}", "",
                         ">=0.6 (Aa) / <=0.25 (Ab) pre-registered", f"HEAT7 control {res['ratio_HEAT7_over_B56core']:.2f}, C-core control {res['ratio_Ccore_over_B56core']:.2f}"])
            print(f"{iso} {tag:9s} helix/B56core = {r:.3f}   helix {res['helix_mean_sigma']:.2f} sigma, core {res['B56core_mean_sigma']:.2f} sigma, "
                  f"frac>core-median {res['helix_frac_above_core_median']:.2f}  (controls HEAT7 {res['ratio_HEAT7_over_B56core']:.2f}, C {res['ratio_Ccore_over_B56core']:.2f}; fit RMSD {fit_rmsd:.2f} A)")
    json.dump(result, open(os.path.join(OUT, "fig2b_density.json"), "w"), indent=2)
    with open(os.path.join(REPO, "results", "fig2b_metrics.csv"), "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["track", "system", "quantity", "selection_or_method", "value", "unit", "paper_value", "notes"])
        w.writerows(rows)


if __name__ == "__main__":
    sys.exit(main())
