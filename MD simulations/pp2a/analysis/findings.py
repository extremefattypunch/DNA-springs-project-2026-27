#!/usr/bin/env python3
"""Reduce summary.csv and the per-frame CSVs to the numbers the wiki page quotes.

Single replicates, so every comparison between systems is made on block means: each
100 ns series is cut into ten 9.5 ns blocks after the 5 ns discard, and Welch's t is
run on the two sets of ten block means.  That is an honest error for one trajectory
and a crude one -- it treats blocks as independent, which they are for anything that
decorrelates faster than ~10 ns and are not for the slow horseshoe motions.  The page
says so where it matters.

The zero reference for every chimera number is P2_clicked: the protein with the same
tether chemistry and no DNA.  Comparing against the unmodified protein is the mistake
the sfGFP review caught.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

from common import CONTROL, SITES, SYSTEMS, analysis_prefix, heat_ranges

HERE = Path(__file__).resolve().parent
N_BLOCKS = 10


def blocks(series, n=N_BLOCKS):
    x = np.asarray(series, float)
    x = x[~np.isnan(x)]
    if len(x) < 2 * n:
        return x
    return np.array([b.mean() for b in np.array_split(x, n)])


def welch(a, b):
    a, b = np.asarray(a, float), np.asarray(b, float)
    if len(a) < 2 or len(b) < 2:
        return None
    va, vb = a.var(ddof=1) / len(a), b.var(ddof=1) / len(b)
    if va + vb == 0:
        return None
    t = (a.mean() - b.mean()) / np.sqrt(va + vb)
    dof = (va + vb) ** 2 / (va ** 2 / (len(a) - 1) + vb ** 2 / (len(b) - 1))
    p = float(2.0 * stats.t.sf(abs(t), dof))
    return {"delta": float(a.mean() - b.mean()), "sem_delta": float(np.sqrt(va + vb)),
            "t": float(t), "dof": float(dof), "p": p, "n_blocks": [len(a), len(b)],
            "resolved": bool(p < 0.05)}


def per_frame(system, kind):
    f = Path(str(analysis_prefix(system)) + f"_{kind}.csv")
    return pd.read_csv(f) if f.exists() else None


def main(summary=None, out=None):
    df = pd.read_csv(summary or HERE / "summary.csv")
    present = [s for s in SYSTEMS if s in set(df.system)]
    f = {"total_ns": float(df.ns.sum()), "n_runs": int(len(df)), "systems_present": present,
         "labels": {s: SYSTEMS[s]["label"] for s in present}, "control": CONTROL,
         "sites": list(SITES), "block_method": f"{N_BLOCKS} blocks per run after the discard; Welch t on block means"}

    def col(s, c):
        if c not in df.columns:
            return None
        v = df.loc[df.system == s, c]
        return None if v.empty or pd.isna(v.iloc[0]) else v.iloc[0]

    # ---- end-to-end and its displacement against the clicked zero ----
    mech = {s: per_frame(s, "mechanics") for s in present}
    f["extension"] = {}
    for s in present:
        m = mech.get(s)
        if m is None or "anchor_CB_CB_A" not in m:
            continue
        x = m["anchor_CB_CB_A"].to_numpy()
        b = blocks(x)
        f["extension"][s] = {"mean_A": float(x.mean()), "sd_A": float(x.std()),
                             "sem_block_A": float(b.std(ddof=1) / np.sqrt(len(b))) if len(b) > 1 else None,
                             "min_A": float(x.min()), "max_A": float(x.max()),
                             "blocks_A": [round(float(v), 3) for v in b]}
    f["extension_vs_control"] = {}
    if CONTROL in f["extension"]:
        base = f["extension"][CONTROL]["blocks_A"]
        for s in present:
            if s == CONTROL or s not in f["extension"]:
                continue
            w = welch(f["extension"][s]["blocks_A"], base)
            if w:
                w["label"] = SYSTEMS[s]["label"]
                f["extension_vs_control"][s] = w
    f["clicked_vs_tet"] = (welch(f["extension"][CONTROL]["blocks_A"], f["extension"]["P1_tet"]["blocks_A"])
                           if CONTROL in f["extension"] and "P1_tet" in f["extension"] else None)

    # ---- the strand ----
    if "P3_ss50" in present and col("P3_ss50", "ssdna.n_nt") is not None:
        s = "P3_ss50"
        f["ss_spring"] = {
            "n_nt": int(col(s, "ssdna.n_nt")),
            "contour_nm": col(s, "ssdna.marko_siggia.contour_nm"),
            "natural_rms_nm": col(s, "ssdna.marko_siggia.natural_rms_end_to_end_nm"),
            "state": col(s, "ssdna.marko_siggia.state_at_mean_x"),
            "x_A": {k: col(s, f"ssdna.end_to_end.A.{k}") for k in ("mean", "sd", "min", "max")},
            "tension_pN_clipped": col(s, "ssdna.marko_siggia.tension_pN_clipped.mean"),
            "tension_pN_unclipped": col(s, "ssdna.marko_siggia.tension_pN_unclipped.mean"),
            "ideal_chain_pN": col(s, "ssdna.marko_siggia.ideal_chain_pN.mean"),
            "fraction_of_frames_stretched": col(s, "ssdna.marko_siggia.fraction_of_frames_stretched"),
            "rg_A": col(s, "ssdna.rg_A.mean"),
            "sharp_turns_per_frame": col(s, "ssdna.trace_bend.sharp_turns_per_frame.mean"),
            "hottest_nt": col(s, "ssdna.trace_bend.hottest_nt"),
            "stacked_fraction": col(s, "ssdna.stacking.fraction.mean"),
            "nt_in_contact": col(s, "ssdna.protein_contact.nt_in_contact_per_frame.mean"),
            "fraction_of_strand_on_protein": col(s, "ssdna.protein_contact.fraction_of_strand"),
            "end_to_end_pbc": col(s, "ssdna.end_to_end.pbc.convention"),
        }
        if col(s, "contact.residues_over_half") is not None:
            f["ss_spring"]["protein_residues_touched_over_half"] = str(col(s, "contact.residues_over_half"))
            f["ss_spring"]["protein_residues_ever_touched"] = col(s, "contact.residues_ever_touched")

    # ---- measured force windows ----
    fxp = HERE / "force_extension.csv"
    if fxp.exists():
        fx = pd.read_csv(fxp)
        f["windows"] = [{k: (None if pd.isna(r[k]) else (float(r[k]) if isinstance(r[k], (int, float, np.floating)) else r[k]))
                         for k in ("tag", "x0_A", "mean_r_A", "sd_r_A", "force_pN", "force_sem_pN", "box_ok",
                                   "box_width_A", "ms_tension_pN", "ms_tension_unclipped_pN", "ideal_chain_pN",
                                   "ms_state", "ns", "independent_samples")}
                        for _, r in fx.iterrows()]

    # ---- rigidity ----
    f["rigidity"] = {}
    for s in present:
        f["rigidity"][s] = {k: col(s, k) for k in (
            "rmsd.rmsd_core_A.mean", "rmsd.rmsd_core_A.max", "rmsd.rmsd_tail_corefit_A.mean",
            "rmsd.rmsd_all_A.mean", "rmsd.rmsd_core_to_mean_A.mean", "rmsd.rmsf_A.mean_core", "rmsd.rmsf_A.mean_tail", "rmsd.rmsf_A.max",
            "rmsd.rmsf_A.argmax_resnum", f"rmsd.rmsf_A.at_sites.{SITES[0]}", f"rmsd.rmsf_A.at_sites.{SITES[1]}",
            "rmsd.helix_fraction_core", f"rmsd.helix_fraction_at_sites.{SITES[0]}",
            f"rmsd.helix_fraction_at_sites.{SITES[1]}")}
        f["rigidity"][s]["rmsf_per_heat"] = {k: col(s, f"rmsd.rmsf_A.per_heat_mean.{k}") for k in range(0, 16)}
    # block comparison of the core RMSF is not possible from a single value; compare the core RMSD series
    f["rigidity_vs_control"] = {}
    rm = {s: per_frame(s, "rmsd") for s in present}
    if rm.get(CONTROL) is not None:
        for s in present:
            if s == CONTROL or rm.get(s) is None:
                continue
            # per-block core RMSD about each block's own mean: a reference-free amplitude,
            # unlike RMSD-to-frame-0, which compares drift from two different starting frames
            if "rmsd_core_to_mean_A" not in rm[s] or "rmsd_core_to_mean_A" not in rm[CONTROL]:
                continue
            w = welch(blocks(rm[s]["rmsd_core_to_mean_A"]), blocks(rm[CONTROL]["rmsd_core_to_mean_A"]))
            if w:
                w["quantity"] = "core RMSD about the run mean, block means"
                f["rigidity_vs_control"][s] = w

    # ---- the solenoid ----
    f["heat"] = {}
    ht = {s: per_frame(s, "heat") for s in present}
    for s in present:
        d = {k: col(s, k) for k in (
            "heat.bend.total_deg.mean", "heat.bend.total_deg.sd", "heat.bend.sum_deg.mean",
            "heat.twist.sum_deg.mean", "heat.aperture_H1_H15_A.mean", "heat.aperture_H1_H15_A.sd",
            "heat.axis_straightness.mean", "heat.literature_pair.A.mean", "heat.literature_pair.A.sd",
            "heat.salt_bridge_backbone.A.mean", "heat.salt_bridge_sidechain.fraction_formed_below_5A",
            "heat.sites_CA_CA_A.mean", "heat.bend.most_variable_repeat")}
        d["region_rmsd_A"] = {r: col(s, f"heat.region_internal_rmsd.{r}.A.mean")
                              for r in ("N (HEAT 1-3)", "central (HEAT 4-12)", "C (HEAT 13-15)")}
        d["bend_per_repeat_mean"] = {k: col(s, f"heat.bend.per_repeat_mean_deg.{k}") for k in range(2, 15)}
        d["bend_per_repeat_sd"] = {k: col(s, f"heat.bend.per_repeat_sd_deg.{k}") for k in range(2, 15)}
        d["repeat_centroid_rmsf"] = {k: col(s, f"heat.repeat_centroid_rmsf_A.{k}") for k in range(1, 16)}
        f["heat"][s] = d
    f["hinge_vs_control"] = {}
    if ht.get(CONTROL) is not None:
        for s in present:
            if s == CONTROL or ht.get(s) is None:
                continue
            shifts = {}
            for k in range(2, 15):
                c = f"bend_H{k}_deg"
                if c in ht[s] and c in ht[CONTROL]:
                    w = welch(blocks(ht[s][c]), blocks(ht[CONTROL][c]))
                    if w:
                        shifts[k] = {"delta_deg": w["delta"], "p": w["p"],
                                     "sd_ratio": float(ht[s][c].std() / max(ht[CONTROL][c].std(), 1e-9))}
            if shifts:
                top = sorted(shifts, key=lambda k: shifts[k]["p"])[:3]
                f["hinge_vs_control"][s] = {"per_repeat": shifts, "largest_shifts": top,
                                            "aperture": welch(blocks(ht[s]["aperture_H1_H15_A"]),
                                                              blocks(ht[CONTROL]["aperture_H1_H15_A"]))}

    # ---- the sites ----
    f["sites"] = {}
    try:
        att = json.loads((HERE.parent / "build" / "01_protein" / "attachment_sites.json").read_text())
        pair = next(iter(att["pairs"].values()))
        f["sites"]["native_AF"] = {str(SITES[0]): {"phi": pair["site_a"].get("phi_deg"), "psi": pair["site_a"].get("psi_deg"),
                                                   "resname": pair["site_a"]["resname"], "plddt": pair["site_a"].get("plddt")},
                                   str(SITES[1]): {"phi": pair["site_b"].get("phi_deg"), "psi": pair["site_b"].get("psi_deg"),
                                                   "resname": pair["site_b"]["resname"], "plddt": pair["site_b"].get("plddt")}}
        f["sites"]["CB_CB_AF_A"] = pair["d_anchor_anchor_A"]
        f["sites"]["numbering_note"] = att["numbering_note"]
    except Exception:                                             # noqa: BLE001
        pass
    for s in present:
        f["sites"][s] = {"label": SYSTEMS[s]["label"]}
        for site in SITES:
            f["sites"][s][str(site)] = {"residue": col(s, f"site.resname_{site}"),
                                        "phi": col(s, f"site.phi_{site}.mean"), "phi_sd": col(s, f"site.phi_{site}.sd"),
                                        "psi": col(s, f"site.psi_{site}.mean"), "psi_sd": col(s, f"site.psi_{site}.sd"),
                                        "ring_neighbours_per_frame": col(s, f"site.ring_neighbours_{site}.contacts_per_frame")}

    out = Path(out or HERE / "findings.json")
    out.write_text(json.dumps(f, indent=2, default=float))
    print(f"wrote {out}")
    for s, w in f["extension_vs_control"].items():
        print(f"  {s:<11} Cb-Cb vs {CONTROL}: {w['delta']:+.2f} +/- {w['sem_delta']:.2f} A, p = {w['p']:.3f}")
    return f


if __name__ == "__main__":
    main()
