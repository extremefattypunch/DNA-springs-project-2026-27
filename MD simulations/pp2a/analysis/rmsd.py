#!/usr/bin/env python3
"""Backbone rigidity: RMSD against the first frame, per-residue RMSF, secondary structure.

Follows Day et al. 2025 so the apo baseline can be set against their Fig 1: C-alpha RMSD
after superposition on the core (residues 14-589), reported for the whole chain, for the
core alone, and for the N-terminal tail 1-13 after the core fit (their "N-terminus"
panel, 10-15 A for Aalpha).  The reference is frame 0 of the analysed window, as theirs
was ("frame 1").  RMSF is about the mean structure after the same core fit, with the six
residues either side of each attachment site excluded from the fit, as the sfGFP work
did, so the tethers cannot steer the superposition.  DSSP is simplified H/E/C per
residue, averaged over frames -- the test that the modification changes no secondary
structure.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import mdtraj as md
import numpy as np
import pandas as pd

from common import CORE, N_TAIL, SITES, ca_by_resnum, heat_of, load, load_resmap, summ


def analyse(traj, top_path, resmap_path, out_prefix, stride=1, skip=0, sites=SITES):
    t = load(traj, top_path, stride, skip)
    top = t.topology
    resmap = load_resmap(top, resmap_path)
    ca = ca_by_resnum(top, resmap)
    nums = sorted(ca)
    all_idx = np.array([ca[n] for n in nums])
    core_nums = [n for n in nums if CORE[0] <= n <= CORE[1]
                 and all(abs(n - s) > 6 for s in sites)]
    core = np.array([ca[n] for n in core_nums])
    tail = np.array([ca[n] for n in nums if N_TAIL[0] <= n <= N_TAIL[1]])
    rows = {"time_ps": t.time if t.time is not None else np.arange(t.n_frames)}

    rows["rmsd_all_A"] = md.rmsd(t, t, 0, atom_indices=all_idx) * 10.0
    rows["rmsd_core_A"] = md.rmsd(t, t, 0, atom_indices=core) * 10.0
    t.superpose(t, 0, atom_indices=core)
    d = t.xyz[:, tail] - t.xyz[0, tail]
    rows["rmsd_tail_corefit_A"] = np.sqrt((d ** 2).sum(-1).mean(-1)) * 10.0
    d = t.xyz[:, all_idx] - t.xyz[0, all_idx]
    rows["rmsd_all_corefit_A"] = np.sqrt((d ** 2).sum(-1).mean(-1)) * 10.0

    rmsf = md.rmsf(t, t, 0, atom_indices=all_idx) * 10.0
    per_res = pd.DataFrame({"resnum": nums, "heat": [heat_of(n) for n in nums],
                            "resname": [top.atom(ca[n]).residue.name for n in nums],
                            "rmsf_A": np.round(rmsf, 4)})

    # DSSP on the protein only (DNA/water confuse mdtraj's implementation)
    prot = t.atom_slice(top.select("protein"))
    codes = md.compute_dssp(prot, simplified=True)
    ss = {}
    for col, r in enumerate(prot.topology.residues):
        c = codes[:, col]
        ss[r.index] = (float((c == "H").mean()), float((c == "E").mean()), float((c == "C").mean()))
    # map back: protein residues keep their order, so pair by sequence position
    prot_res = [r for r in prot.topology.residues]
    dssp_by_num = {}
    for num, (pr) in zip(sorted(n for n in nums), prot_res):
        dssp_by_num[num] = ss[pr.index]
    per_res["helix_frac"] = [round(dssp_by_num.get(n, (np.nan,) * 3)[0], 3) for n in nums]
    per_res.to_csv(f"{out_prefix}_rmsf.csv", index=False)

    core_mask = per_res.resnum.between(*CORE)
    summary = {
        "n_frames": int(t.n_frames),
        "reference": "frame 0 of the analysed window; core fit = CA 14-589 minus sites +/-6",
        "rmsd_all_A": summ(rows["rmsd_all_A"]),
        "rmsd_core_A": summ(rows["rmsd_core_A"]),
        "rmsd_tail_corefit_A": summ(rows["rmsd_tail_corefit_A"]),
        "rmsd_all_corefit_A": summ(rows["rmsd_all_corefit_A"]),
        "rmsf_A": {"mean_all": round(float(rmsf.mean()), 3),
                   "mean_core": round(float(per_res.rmsf_A[core_mask].mean()), 3),
                   "mean_tail": round(float(per_res.rmsf_A[~core_mask].mean()), 3),
                   "max": round(float(rmsf.max()), 3),
                   "argmax_resnum": int(per_res.resnum[per_res.rmsf_A.idxmax()]),
                   "at_sites": {str(s): round(float(per_res.rmsf_A[per_res.resnum == s].iloc[0]), 3)
                                for s in sites if (per_res.resnum == s).any()},
                   "per_heat_mean": {int(k): round(float(g.rmsf_A.mean()), 3)
                                     for k, g in per_res.groupby("heat")}},
        "helix_fraction_core": round(float(per_res.helix_frac[core_mask].mean()), 4),
        "helix_fraction_at_sites": {str(s): per_res.helix_frac[per_res.resnum == s].iloc[0]
                                    for s in sites if (per_res.resnum == s).any()},
        "day2025_reference": "Aalpha: N-term 1-13 RMSD ~10-15 A, rest ~2-4 A over 100 ns; RMSF 1-3 A beyond residue 20",
    }
    print(f"  RMSD core {rows['rmsd_core_A'].mean():.2f} A (max {rows['rmsd_core_A'].max():.2f}), "
          f"tail (core fit) {rows['rmsd_tail_corefit_A'].mean():.2f} A; "
          f"RMSF core {summary['rmsf_A']['mean_core']:.2f} A, at sites {summary['rmsf_A']['at_sites']}; "
          f"helix fraction {summary['helix_fraction_core']:.3f}")
    pd.DataFrame(rows).to_csv(f"{out_prefix}_rmsd.csv", index=False)
    Path(f"{out_prefix}_rmsd.json").write_text(json.dumps(summary, indent=2))
    return summary


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--traj", required=True)
    ap.add_argument("--top", required=True)
    ap.add_argument("--resmap", default=None)
    ap.add_argument("--out-prefix", required=True)
    ap.add_argument("--stride", type=int, default=1)
    ap.add_argument("--skip", type=int, default=0)
    a = ap.parse_args()
    analyse(a.traj, a.top, a.resmap, a.out_prefix, a.stride, a.skip)
