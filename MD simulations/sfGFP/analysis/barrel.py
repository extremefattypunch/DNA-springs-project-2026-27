#!/usr/bin/env python3
"""Does the spring deform the barrel, and where?

The two attachment sites sit on adjacent strands of the beta-barrel with His148
between them, so the obvious question is whether pulling them apart shears that
sheet open and lets the chromophore see solvent.  Four measurements:

* **Per-residue RMSF** on Calpha, after aligning on the barrel core.  Aligning on the
  whole protein would let the tethers and the loops they sit in dominate the fit and
  smear the signal across every residue.
* **Sheet register** -- the backbone N-H...O=C distances that stitch the two strands
  carrying the attachment sites to their neighbours.  A shear shows up here as
  specific hydrogen bonds lengthening while their neighbours hold.
* **Barrel cross-section** -- the two principal widths of the Calpha shell
  perpendicular to the barrel axis.  A squeeze along one axis and a bulge along the
  other is what an external force does to a cylinder.
* **Water in the cavity** -- waters within 6 A of the chromophore that are also inside
  the barrel radius.  Non-radiative decay in GFP tracks water access, so this is the
  mechanistic bridge from mechanics to the fluorometry.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import mdtraj as md
import numpy as np
import pandas as pd


def load_map(path):
    rows = json.loads(Path(path).read_text())
    return {r["index"]: (r["orig_resnum"], r["orig_resname"],
                         r.get("orig_chain", "?")) for r in rows}


def analyse(traj, top_path, resmap_path, sites, out_prefix, stride=1):
    t = md.load(str(traj), top=str(top_path), stride=stride)
    top = t.topology
    resmap = load_map(resmap_path)
    prot_res = [i for i, (n, nm, ch) in resmap.items()
                if nm not in ("WAT", "HOH") and not nm.rstrip("35") in
                ("DA", "DT", "DG", "DC", "DNL")]
    num_of = {i: resmap[i][0] for i in prot_res}

    ca = [(i, next((a.index for a in top.residue(i).atoms if a.name == "CA"), None))
          for i in prot_res]
    ca = [(i, x) for i, x in ca if x is not None]
    ca_idx = np.array([x for _, x in ca])

    # align on the barrel core: strands only, and never the attachment loops
    core = np.array([x for i, x in ca if abs(num_of[i] - sites[0]) > 6
                     and abs(num_of[i] - sites[1]) > 6])
    t.superpose(t, 0, atom_indices=core)
    rmsf = md.rmsf(t, t, 0, atom_indices=ca_idx) * 10.0
    per_res = pd.DataFrame({"resnum": [num_of[i] for i, _ in ca],
                            "resname": [resmap[i][1] for i, _ in ca],
                            "rmsf_A": np.round(rmsf, 4)})
    per_res.to_csv(f"{out_prefix}_rmsf.csv", index=False)

    summary = {"n_frames": int(t.n_frames),
               "rmsf_A": {"mean": round(float(rmsf.mean()), 3),
                          "max": round(float(rmsf.max()), 3),
                          "at_sites": {str(s): round(float(
                              rmsf[[num_of[i] for i, _ in ca].index(s)]), 3)
                              for s in sites if s in [num_of[i] for i, _ in ca]}}}
    print(f"  RMSF: mean {rmsf.mean():.2f} A, max {rmsf.max():.2f} A, "
          f"at sites {summary['rmsf_A']['at_sites']}")

    # --- sheet register around the attachment strands ---
    rows = {"time_ps": t.time if t.time is not None else np.arange(t.n_frames)}
    by_num = {num_of[i]: i for i in prot_res if i in num_of}
    hb_pairs, hb_labels = [], []
    # Search every protein residue for a partner, not a narrow sequence window: in an
    # antiparallel sheet the cross-strand partner can be twenty residues away, and a
    # window of a+3..a+11 found only two of the hydrogen bonds that stitch this region
    # together.  Local i,i+3 contacts are excluded as helical rather than sheet.
    lo, hi = min(sites) - 4, max(sites) + 4
    all_nums = sorted(by_num)
    for a in range(lo, hi + 1):
        for b in all_nums:
            if a not in by_num or abs(b - a) < 3:
                continue
            ra, rb = top.residue(by_num[a]), top.residue(by_num[b])
            na = next((x.index for x in ra.atoms if x.name == "N"), None)
            ob = next((x.index for x in rb.atoms if x.name == "O"), None)
            if na is None or ob is None:
                continue
            d0 = np.linalg.norm(t.xyz[0, na] - t.xyz[0, ob]) * 10.0
            if d0 < 3.5:                       # an H-bond in the starting structure
                hb_pairs.append([na, ob])
                hb_labels.append(f"N{a}-O{b}")
    if hb_pairs:
        d = md.compute_distances(t, np.array(hb_pairs)) * 10.0
        for k, lab in enumerate(hb_labels):
            rows[f"sheet_{lab}_A"] = d[:, k]
        summary["sheet_register"] = {
            lab: {"mean_A": round(float(d[:, k].mean()), 3),
                  "sd_A": round(float(d[:, k].std()), 3),
                  "occupancy_below_3.5A": round(float((d[:, k] < 3.5).mean()), 4)}
            for k, lab in enumerate(hb_labels)}
        worst = min(summary["sheet_register"],
                    key=lambda k: summary["sheet_register"][k]["occupancy_below_3.5A"])
        print(f"  sheet register: {len(hb_labels)} backbone H-bonds tracked; "
              f"least stable {worst} "
              f"({100 * summary['sheet_register'][worst]['occupancy_below_3.5A']:.0f}%)")

    # --- barrel cross-section ---
    xyz = t.xyz[:, ca_idx, :] * 10.0
    cen = xyz.mean(axis=1, keepdims=True)
    dev = xyz - cen
    widths = np.zeros((t.n_frames, 3))
    for f in range(t.n_frames):
        s = np.linalg.svd(dev[f], compute_uv=False)
        widths[f] = 2.0 * s / np.sqrt(len(ca_idx))
    rows["barrel_long_A"], rows["barrel_mid_A"], rows["barrel_short_A"] = widths.T
    rows["barrel_ellipticity"] = widths[:, 1] / widths[:, 2]
    summary["barrel"] = {
        "widths_A_mean": [round(float(v), 3) for v in widths.mean(axis=0)],
        "widths_A_sd": [round(float(v), 3) for v in widths.std(axis=0)],
        "cross_section_ellipticity": [round(float(rows["barrel_ellipticity"].mean()), 4),
                                      round(float(rows["barrel_ellipticity"].std()), 4)]}
    print(f"  barrel widths {summary['barrel']['widths_A_mean']} A, "
          f"cross-section ellipticity "
          f"{summary['barrel']['cross_section_ellipticity'][0]:.3f}")

    # --- water in the chromophore cavity ---
    cro = next((r for r in top.residues if r.name == "CRO"), None)
    if cro is not None:
        cro_idx = np.array([a.index for a in cro.atoms])
        wo = top.select("water and name O")
        if len(wo):
            near = md.compute_neighbors(t, 0.6, cro_idx, haystack_indices=wo,
                                        periodic=True)
            cnt = np.array([len(x) for x in near], float)
            rows["waters_within_6A_of_CRO"] = cnt
            summary["cavity_water"] = {
                "cutoff_A": 6.0,
                "mean": round(float(cnt.mean()), 3),
                "sd": round(float(cnt.std()), 3)}
            print(f"  cavity water within 6 A of the chromophore: "
                  f"{cnt.mean():.2f} +/- {cnt.std():.2f}")

    pd.DataFrame(rows).to_csv(f"{out_prefix}_barrel.csv", index=False)
    Path(f"{out_prefix}_barrel.json").write_text(json.dumps(summary, indent=2))
    print(f"  wrote {out_prefix}_barrel.csv/.json and _rmsf.csv")
    return summary


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--traj", required=True)
    ap.add_argument("--top", required=True)
    ap.add_argument("--resmap", required=True)
    ap.add_argument("--sites", nargs=2, type=int, default=[133, 149])
    ap.add_argument("--out-prefix", required=True)
    ap.add_argument("--stride", type=int, default=1)
    a = ap.parse_args()
    analyse(a.traj, a.top, a.resmap, a.sites, a.out_prefix, a.stride)
