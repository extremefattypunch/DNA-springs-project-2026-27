#!/usr/bin/env python3
"""Reduce analysis/summary.csv to the numbers the report quotes.

Kept separate from the report so the report renders numbers it did not compute, and
so every figure in the report can be checked against the same JSON.  Anything the
statistics do not support is marked here, not softened in the prose: ``resolved``
flags whether a difference survives its own spread, and ``bistable`` flags an
observable whose replicates sit at two separated values, where a mean is meaningless.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
SPRINGS = [("S3_spring27", "27 bp"), ("S4_spring40", "40 bp"),
           ("S5_spring40nick", "40 bp, nicked")]
PROTEINS = [("S0_wt", "WT sfGFP"), ("S1_tet", "2× Tet2-Et, unloaded"),
            ("S2_clicked", "clicked, no DNA"), ("S3_spring27", "27 bp spring"),
            ("S4_spring40", "40 bp spring"), ("S5_spring40nick", "40 bp nicked")]
T95_DF4 = 2.776          # two-sided 95% t for 4 degrees of freedom (n=3 vs n=3)


def agg(sub, col):
    v = sub[col].dropna().values.astype(float)
    if not len(v):
        return None
    return {"mean": float(v.mean()), "sd": float(v.std(ddof=1)) if len(v) > 1 else 0.0,
            "sem": float(v.std(ddof=1) / np.sqrt(len(v))) if len(v) > 1 else 0.0,
            "min": float(v.min()), "max": float(v.max()), "n": int(len(v)),
            "values": [round(float(x), 4) for x in v]}


def welch(a, b):
    """Two-sample t and whether it clears 95% with these tiny n."""
    a, b = np.asarray(a, float), np.asarray(b, float)
    if len(a) < 2 or len(b) < 2:
        return None
    sp = np.sqrt(((len(a) - 1) * a.var(ddof=1) + (len(b) - 1) * b.var(ddof=1))
                 / (len(a) + len(b) - 2))
    if sp == 0:
        return None
    t = (a.mean() - b.mean()) / (sp * np.sqrt(1 / len(a) + 1 / len(b)))
    return {"delta": float(a.mean() - b.mean()), "t": float(t),
            "resolved": bool(abs(t) > T95_DF4)}


def main(summary=None, out=None):
    summary = Path(summary or HERE / "summary.csv")
    df = pd.read_csv(summary)
    f = {"total_ns": float(df["ns"].sum()), "n_replicates": int(len(df)),
         "systems_present": sorted(set(df.system))}

    # ---- springs ----
    f["springs"] = {}
    for k, nm in SPRINGS:
        sub = df[df.system == k]
        if not len(sub):
            continue
        f["springs"][k] = {
            "label": nm, "n_bp": int(sub["mech.spring.n_bp"].iloc[0]),
            "nicked": bool(sub["mech.spring.nicked"].iloc[0]),
            "tau_c": float(sub["mech.spring.tau_c_pN_nm"].iloc[0]),
            "regime": str(sub["mech.spring.regime"].iloc[0]),
            "contour_nm": float(sub["mech.spring.contour_nm"].iloc[0]),
            "x_A": agg(sub, "mech.spring.x_mean_A"),
            "x_wander_A": agg(sub, "mech.spring.x_sd_A"),
            "force_pN": agg(sub, "mech.spring.force_pN_mean"),
            "energy_kT": agg(sub, "mech.spring.energy_kT_mean"),
            "bend_total_deg": agg(sub, "mech.bend.total_bend_deg_mean"),
            "bend_per_step_deg": agg(sub, "mech.bend.per_step_mean_deg"),
            "hottest_bp": agg(sub, "mech.bend.hottest_bp"),
            "bp_open": agg(sub, "mech.base_pair_opening.mean_open_count"),
            "extension_A": agg(sub, "mech.anchor_CB_CB_A.mean"),
        }

    # ---- the clamp ladder, per force and per topology ----
    f["clamp"] = {}
    for base in ("S6_clamp", "S2_clicked"):
        sub = df[(df.system == base) & df.clamp_pN.notna()]
        if not len(sub):
            continue
        ladder = []
        for fx, g in sub.groupby("clamp_pN"):
            a = agg(g, "mech.anchor_CB_CB_A.mean")
            ladder.append({"pN": float(fx), **a})
        x = sub["clamp_pN"].values.astype(float)
        y = sub["mech.anchor_CB_CB_A.mean"].values.astype(float)
        entry = {"ladder": sorted(ladder, key=lambda r: r["pN"]), "n": int(len(x))}
        if len(x) > 3 and len(set(x)) > 2:
            k, b = np.polyfit(x, y, 1)
            se = (np.sqrt(((y - (k * x + b)) ** 2).sum() / (len(x) - 2))
                  / np.sqrt(((x - x.mean()) ** 2).sum()))
            tcrit = 2.04
            entry.update({
                "slope_A_per_pN": float(k), "slope_se": float(se),
                "ci95": [float(k - tcrit * se), float(k + tcrit * se)],
                "slope_resolved": bool(k - tcrit * se > 0),
                "kappa_lower_bound_kT_nm2": float(10.0 / ((k + tcrit * se) * 4.2))
                if k + tcrit * se > 0 else None})
        f["clamp"][base] = entry

    # ---- protein extension, every system, against the matched zero-force control ----
    zero = df[(df.system == "S6_clamp") & (df.clamp_pN == 0)]
    f["extension"] = {"zero_force_ref": agg(zero, "mech.anchor_CB_CB_A.mean")}
    for k, nm in PROTEINS:
        sub = df[df.system == k]
        if len(sub):
            f["extension"][k] = {"label": nm,
                                 **agg(sub, "mech.anchor_CB_CB_A.mean")}
    # does the ncAA substitution alone move the sites?
    w = welch(df[df.system == "S1_tet"]["mech.anchor_CB_CB_A.mean"].values,
              df[df.system == "S0_wt"]["mech.anchor_CB_CB_A.mean"].values)
    f["extension"]["tet_vs_wt"] = w
    # strongest against weakest spring, matched chemistry
    w2 = welch(df[df.system == "S3_spring27"]["mech.anchor_CB_CB_A.mean"].values,
               df[df.system == "S5_spring40nick"]["mech.anchor_CB_CB_A.mean"].values)
    f["extension"]["strongest_vs_weakest_spring"] = w2
    if f["springs"]:
        fs = [f["springs"][k]["force_pN"]["mean"] for k, _ in SPRINGS
              if k in f["springs"]]
        es = [f["springs"][k]["extension_A"]["mean"] for k, _ in SPRINGS
              if k in f["springs"]]
        if len(fs) > 2:
            f["extension"]["chimera_slope_A_per_pN"] = float(np.polyfit(fs, es, 1)[0])

    # ---- chromophore ----
    hb_cols = [c for c in df.columns
               if c.startswith("chromo.hbonds.") and c.endswith(".occupancy")]
    f["hbonds"] = {}
    for c in hb_cols:
        name = c[len("chromo.hbonds."):-len(".occupancy")]
        per = {}
        bistable = False
        for k, nm in PROTEINS:
            sub = df[df.system == k]
            a = agg(sub, c)
            if a is None:
                continue
            per[k] = a
            v = np.array(a["values"])
            # two separated clusters with nothing in between: a mean is meaningless
            if len(v) > 2 and (v.max() - v.min()) > 0.5 and \
                    not ((v > 0.25) & (v < 0.75)).any():
                bistable = True
        f["hbonds"][name] = {"per_system": per, "bistable": bistable}
    f["twist"] = {}
    for d in ("tau_phenol_bridge", "phi_bridge_imidazolinone"):
        f["twist"][d] = {k: agg(df[df.system == k],
                                f"chromo.dihedrals.{d}.twist_mean_deg")
                         for k, _ in PROTEINS if len(df[df.system == k])}
    f["cavity_water"] = {k: agg(df[df.system == k], "barrel.cavity_water.mean")
                         for k, _ in PROTEINS if len(df[df.system == k])}
    f["rmsf"] = {k: agg(df[df.system == k], "barrel.rmsf_A.mean")
                 for k, _ in PROTEINS if len(df[df.system == k])}
    f["labels"] = {k: nm for k, nm in PROTEINS}

    out = Path(out or HERE / "findings.json")
    out.write_text(json.dumps(f, indent=2))
    print(f"wrote {out}")
    print(f"  {f['n_replicates']} replicates, {f['total_ns']:.0f} ns")
    print(f"  springs: {list(f['springs'])}")
    print(f"  bistable H-bonds: "
          f"{[n for n, v in f['hbonds'].items() if v['bistable']]}")
    return f


if __name__ == "__main__":
    main()
