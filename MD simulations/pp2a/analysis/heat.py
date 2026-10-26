#!/usr/bin/env python3
"""Does the horseshoe kink?  Bend, twist and aperture of the HEAT solenoid.

The literature reports where PR65 hinges (HEAT 10-11 from crystal comparison, 11-12 in
Day et al. 2025, 12-13 in the Yilmaz/Bahar MD) but gives no bend-angle definition.  This
uses the same construction the sfGFP work used for the DNA helical axis, applied to the
solenoid: the axis is the chain of HEAT-repeat C-alpha centroids c_1 ... c_15, the
tangent at step k is t_k = c_{k+1} - c_k, and

    bend_k   = angle(t_{k-1}, t_k)           the kink at repeat k, degrees
    twist_k  = dihedral(c_{k-1}, c_k, c_{k+1}, c_{k+2})   the superhelical twist
    total    = angle(t_1, t_14)              how far the solenoid turns end to end

plus the quantities other groups quote: the HEAT1-HEAT15 centroid distance (Day:
62.7 A in the B56 holoenzyme), the Calpha Asn30-Phe578 distance (Yilmaz "N29-F577":
66.5 +/- 9.7 A apo, 47.7 A compact 6NTS, 70.3 A extended 1B3U), the Arg113-Glu191
inter-repeat salt bridge (Day Fig 3C), and the internal Calpha RMSD of the N / central / C
blocks after superposing each block on itself (Yilmaz Fig 8b: 2.1 / 2.7 / 1.4 A apo).

A hinge is a repeat whose bend angle moves -- in mean against the clicked control, or
in spread -- not a repeat with a large static bend; a solenoid is bent everywhere.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import mdtraj as md
import numpy as np
import pandas as pd

from common import (LITERATURE_PAIR, N_HEAT, REGIONS, SALT_BRIDGE, SITES, angle_deg,
                    atom_by_resnum, ca_by_resnum, dihedral_deg, heat_ranges, load,
                    load_resmap, summ)


def centroids(t, ca, ranges):
    """(n_frames, n_heat, 3) C-alpha centroid of each repeat, in A."""
    out = np.zeros((t.n_frames, len(ranges), 3))
    for j, (k, (lo, hi)) in enumerate(sorted(ranges.items())):
        idx = [ca[n] for n in range(lo, hi + 1) if n in ca]
        out[:, j] = t.xyz[:, idx].mean(axis=1) * 10.0
    return out


def analyse(traj, top_path, resmap_path, out_prefix, stride=1, skip=0):
    t = load(traj, top_path, stride, skip)
    top = t.topology
    resmap = load_resmap(top, resmap_path)
    ca = ca_by_resnum(top, resmap)
    ranges = heat_ranges()
    c = centroids(t, ca, ranges)                     # frames x 15 x 3
    tang = c[:, 1:] - c[:, :-1]                      # frames x 14 x 3, t_1..t_14
    rows = {"time_ps": t.time if t.time is not None else np.arange(t.n_frames)}
    summary = {"n_frames": int(t.n_frames), "axis": "HEAT C-alpha centroids, nominal 39-residue grid"}

    bend = angle_deg(tang[:, :-1], tang[:, 1:])      # frames x 13: bend at repeats 2..14
    twist = dihedral_deg(c[:, :-3], c[:, 1:-2], c[:, 2:-1], c[:, 3:])   # frames x 12: at steps 2..13
    for j in range(bend.shape[1]):
        rows[f"bend_H{j + 2}_deg"] = bend[:, j]
    for j in range(twist.shape[1]):
        rows[f"twist_H{j + 2}_deg"] = twist[:, j]
    rows["bend_total_deg"] = angle_deg(tang[:, 0], tang[:, -1])
    rows["bend_sum_deg"] = bend.sum(axis=1)
    rows["twist_sum_deg"] = twist.sum(axis=1)
    rows["aperture_H1_H15_A"] = np.linalg.norm(c[:, 0] - c[:, -1], axis=1)
    # radius of curvature of the solenoid from the end-to-end chord and the arc length
    arc = np.linalg.norm(tang, axis=2).sum(axis=1)
    chord = rows["aperture_H1_H15_A"]
    rows["axis_arc_A"] = arc
    rows["axis_straightness"] = chord / arc

    profile = pd.DataFrame({"heat": np.arange(2, 2 + bend.shape[1]),
                            "bend_deg_mean": bend.mean(0).round(3),
                            "bend_deg_sd": bend.std(0).round(3)})
    tw = pd.DataFrame({"heat": np.arange(2, 2 + twist.shape[1]),
                       "twist_deg_mean": twist.mean(0).round(3),
                       "twist_deg_sd": twist.std(0).round(3)})
    profile = profile.merge(tw, on="heat", how="left")
    profile.to_csv(f"{out_prefix}_heat_profile.csv", index=False)
    summary["bend"] = {
        "per_repeat_mean_deg": {int(h): round(float(m), 2) for h, m in zip(profile.heat, profile.bend_deg_mean)},
        "per_repeat_sd_deg": {int(h): round(float(s), 2) for h, s in zip(profile.heat, profile.bend_deg_sd)},
        "most_variable_repeat": int(profile.heat[profile.bend_deg_sd.idxmax()]),
        "total_deg": summ(rows["bend_total_deg"]),
        "sum_deg": summ(rows["bend_sum_deg"]),
    }
    summary["twist"] = {
        "per_step_mean_deg": {int(h): round(float(m), 2) for h, m in zip(tw.heat, tw.twist_deg_mean)},
        "sum_deg": summ(rows["twist_sum_deg"]),
    }
    summary["aperture_H1_H15_A"] = summ(rows["aperture_H1_H15_A"])
    summary["axis_straightness"] = summ(rows["axis_straightness"])
    print(f"  solenoid: total bend {rows['bend_total_deg'].mean():.1f} deg, "
          f"H1-H15 {rows['aperture_H1_H15_A'].mean():.1f} +/- {rows['aperture_H1_H15_A'].std():.1f} A, "
          f"most variable repeat H{summary['bend']['most_variable_repeat']}")

    # --- distances other groups quote ---
    la, lb = LITERATURE_PAIR
    if la in ca and lb in ca:
        d = md.compute_distances(t, [[ca[la], ca[lb]]], periodic=False)[:, 0] * 10.0
        rows[f"d_CA{la}_CA{lb}_A"] = d
        summary["literature_pair"] = {"atoms": f"CA{la}-CA{lb} (Yilmaz 'N29-F577')",
                                      "A": summ(d),
                                      "reference_apo_MD_A": "66.5 +/- 9.7 (CHARMM36, 3 x 704 ns)",
                                      "reference_6NTS_A": 47.7, "reference_1B3U_A": 70.3}
        print(f"  CA{la}-CA{lb}: {d.mean():.1f} +/- {d.std():.1f} A (lit. apo 66.5 +/- 9.7)")
    (ra, rn), (ea, en) = SALT_BRIDGE
    i = atom_by_resnum(top, resmap, ra, "CA")
    j = atom_by_resnum(top, resmap, ea, "CA")
    if i is not None and j is not None:
        d = md.compute_distances(t, [[i, j]], periodic=False)[:, 0] * 10.0
        rows[f"d_CA{ra}_CA{ea}_A"] = d
        summary["salt_bridge_backbone"] = {"atoms": f"CA{ra}-CA{ea} ({rn}{ra}-{en}{ea}, Day Fig 3C)",
                                           "A": summ(d), "reference_9MIP_A": 13.2}
        i2 = atom_by_resnum(top, resmap, ra, "CZ")
        j2 = atom_by_resnum(top, resmap, ea, "CD")
        if i2 is not None and j2 is not None:
            d2 = md.compute_distances(t, [[i2, j2]], periodic=False)[:, 0] * 10.0
            rows[f"d_R{ra}CZ_E{ea}CD_A"] = d2
            summary["salt_bridge_sidechain"] = {"atoms": f"CZ{ra}-CD{ea}", "A": summ(d2),
                                                "fraction_formed_below_5A": round(float((d2 < 5.0).mean()), 3)}
    sa, sb = SITES
    if sa in ca and sb in ca:
        d = md.compute_distances(t, [[ca[sa], ca[sb]]], periodic=False)[:, 0] * 10.0
        rows[f"d_CA{sa}_CA{sb}_A"] = d
        summary["sites_CA_CA_A"] = summ(d)

    # --- region internal RMSD (each block superposed on itself) ---
    region = {}
    for name, (k0, k1) in REGIONS.items():
        lo, hi = ranges[k0][0], ranges[k1][1]
        idx = np.array([ca[n] for n in range(lo, hi + 1) if n in ca])
        sub = t.atom_slice(idx)
        r = md.rmsd(sub, sub, 0) * 10.0
        key = f"rmsd_internal_{name.split()[0]}_A"
        rows[key] = r
        region[name] = {"residues": [lo, hi], "A": summ(r)}
    summary["region_internal_rmsd"] = region
    summary["region_reference"] = "Yilmaz 2025 apo: N 2.1, central 2.7, C 1.4 A"
    print("  region internal RMSD: " + ", ".join(
        f"{k.split()[0]} {v['A']['mean']:.2f} A" for k, v in region.items()))

    # --- per-repeat RMSF of centroids (which repeats move as blocks) ---
    core = np.array([ca[n] for n in range(14, 590) if n in ca])
    t.superpose(t, 0, atom_indices=core)
    c2 = centroids(t, ca, ranges)
    disp = np.linalg.norm(c2 - c2.mean(0), axis=2)
    summary["repeat_centroid_rmsf_A"] = {int(k): round(float(np.sqrt((disp[:, k - 1] ** 2).mean())), 2)
                                         for k in range(1, N_HEAT + 1)}

    pd.DataFrame(rows).to_csv(f"{out_prefix}_heat.csv", index=False)
    Path(f"{out_prefix}_heat.json").write_text(json.dumps(summary, indent=2))
    print(f"  wrote {out_prefix}_heat.csv/.json and _heat_profile.csv")
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
