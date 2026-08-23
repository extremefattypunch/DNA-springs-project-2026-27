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
from scipy import stats
import pandas as pd

HERE = Path(__file__).resolve().parent
# Two families, each with its OWN zero-force control, because a chimera has to be
# compared against a protein carrying the same tether chemistry and nothing else.  For
# sfGFP that is S2_clicked; for ubiquitin it is U1_clicked.  Comparing a ubiquitin
# chimera against the sfGFP control would be meaningless, and comparing either against
# an unclicked protein is the mistake that produced the retracted claim.
FAMILIES = {
    "sfGFP": {
        "control": "S2_clicked",
        "springs": [
            # the superseded two-arm builds, kept and reported as what they are
            ("S3_spring27", "27 bp, two-arm"), ("S4_spring40", "40 bp, two-arm"),
            ("S5_spring40nick", "40 bp nicked, two-arm"),
            # the corrected closed-loop builds
            ("L3_loop27", "27 bp, closed loop"), ("L4_loop40", "40 bp, closed loop"),
            ("L5_loop40nick", "40 bp nicked, closed loop"),
            ("C1_ss50", "50 nt ssDNA, closed loop"), ("C2_ssarms", "two 25 nt arms"),
        ],
        "proteins": [("S0_wt", "WT sfGFP"), ("S1_tet", "2× Tet2-Et, unloaded"),
                     ("S2_clicked", "clicked, no DNA")],
    },
    "ubiquitin": {
        "control": "U1_clicked",
        "springs": [("U2_ss50", "50 nt ssDNA, closed loop"),
                    ("U3_ds50", "50 bp dsDNA, closed loop"),
                    ("U4_ssarms", "two 25 nt arms")],
        "proteins": [("U0_apo", "ubiquitin, apo"),
                     ("U1_clicked", "ubiquitin, clicked, no DNA")],
    },
}
SPRINGS = [s for f in FAMILIES.values() for s in f["springs"]]
PROTEINS = [s for f in FAMILIES.values() for s in f["proteins"] + f["springs"]]


def agg(sub, col):
    v = sub[col].dropna().values.astype(float)
    if not len(v):
        return None
    return {"mean": float(v.mean()), "sd": float(v.std(ddof=1)) if len(v) > 1 else 0.0,
            "sem": float(v.std(ddof=1) / np.sqrt(len(v))) if len(v) > 1 else 0.0,
            "min": float(v.min()), "max": float(v.max()), "n": int(len(v)),
            "values": [round(float(x), 4) for x in v]}


def welch(a, b):
    """Welch's two-sample t, with a p-value and the Welch-Satterthwaite d.o.f.

    Actually Welch's now: the previous version pooled the two variances, which is
    Student's test and assumes equal variance -- not safe when one group is three
    replicates of one construct and the other is six of another.  It also returned no
    p-value, so comparisons could only be called resolved or not against a fixed
    critical t for 4 degrees of freedom, which is wrong whenever the group sizes differ.
    """
    a, b = np.asarray(a, float), np.asarray(b, float)
    if len(a) < 2 or len(b) < 2:
        return None
    va, vb = a.var(ddof=1) / len(a), b.var(ddof=1) / len(b)
    if va + vb == 0:
        return None
    t = (a.mean() - b.mean()) / np.sqrt(va + vb)
    dof = (va + vb) ** 2 / (va ** 2 / (len(a) - 1) + vb ** 2 / (len(b) - 1))
    p = float(2.0 * stats.t.sf(abs(t), dof))
    return {"delta": float(a.mean() - b.mean()), "t": float(t),
            "dof": float(dof), "p": p, "n": [len(a), len(b)],
            "resolved": bool(p < 0.05)}


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

    # ---- protein extension, against the MATCHED control ----
    # The zero has to be the clicked, unloaded protein.  Using the unclicked one, which
    # is what the clamp ladder runs on, is what produced the retracted claim that the
    # springs extend the attachment sites: installing the two Tet2-Et residues moves the
    # sites 0.62 A on its own and clicking gives it back, so every chimera looks ~0.65 A
    # "extended" from the chemistry alone -- 15 to 30 times any force effect.  Both
    # references are recorded, and which is which is stated in the key.
    unclicked = df[(df.system == "S6_clamp") & (df.clamp_pN == 0)]
    clicked = df[(df.system == "S2_clicked") & (df.clamp_pN.fillna(0) == 0)]
    f["extension"] = {
        "zero_force_ref": agg(unclicked, "mech.anchor_CB_CB_A.mean"),
        "zero_force_ref_note": "UNCLICKED protein at 0 pN -- not a valid baseline for a "
                               "chimera; kept only because the clamp ladder runs on it",
        "matched_control": agg(clicked, "mech.anchor_CB_CB_A.mean"),
        "matched_control_note": "clicked, unloaded (S2_clicked at 0 pN) -- the right "
                                "zero for a chimera",
    }
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
    # Each chimera against the matched control, with a Bonferroni correction for the
    # three comparisons -- without it a single p = 0.02 out of three looks decisive.
    f["extension"]["vs_matched_control"] = {}
    for fam, spec in FAMILIES.items():
        ctrl = df[(df.system == spec["control"]) & (df.clamp_pN.fillna(0) == 0)]
        base = ctrl["mech.anchor_CB_CB_A.mean"].dropna().values
        if not len(base):
            continue
        present = [(k, nm) for k, nm in spec["springs"]
                   if len(df[df.system == k]["mech.anchor_CB_CB_A.mean"].dropna())]
        # Bonferroni over the comparisons actually made in THIS family, not over a
        # global count -- correcting an sfGFP comparison for ubiquitin's tests would be
        # both wrong and conveniently conservative.
        for k, nm in present:
            arr = df[df.system == k]["mech.anchor_CB_CB_A.mean"].dropna().values
            w3 = welch(arr, base)
            if w3 is None:
                continue
            w3["p_bonferroni"] = min(1.0, w3["p"] * len(present))
            w3["label"] = nm
            w3["family"] = fam
            w3["control"] = spec["control"]
            w3["n_comparisons_in_family"] = len(present)
            f["extension"]["vs_matched_control"][k] = w3

    if f["springs"]:
        fs = [f["springs"][k]["force_pN"]["mean"] for k, _ in SPRINGS
              if k in f["springs"]]
        es = [f["springs"][k]["extension_A"]["mean"] for k, _ in SPRINGS
              if k in f["springs"]]
        if len(fs) > 2:
            f["extension"]["chimera_slope_A_per_pN"] = float(np.polyfit(fs, es, 1)[0])

    # ---- backbone conformation at the two modified sites ----
    # Circular means across replicates, because these are angles: site 133 sits near
    # +60 deg, and a residue that occasionally samples the other side of the branch cut
    # would otherwise drag the mean somewhere it never visits.
    if "site.phi_133.mean" in df.columns:
        f["sites"] = {}
        for k, nm in PROTEINS:
            s = df[df.system == k]
            if not len(s):
                continue
            entry = {"label": nm}
            for site in (133, 149):
                rn = s.get(f"site.resname_{site}")
                d = {"residue": (rn.iloc[0] if rn is not None and len(rn) else None)}
                for ang in ("phi", "psi"):
                    c = f"site.{ang}_{site}.mean"
                    if c not in s or s[c].isna().all():
                        continue
                    v = np.radians(s[c].dropna().values)
                    d[ang] = round(float(np.degrees(np.arctan2(
                        np.sin(v).mean(), np.cos(v).mean()))), 2)
                st = f"site.stack_Tyr200_{site}.contacts_per_frame"
                if st in s and not s[st].isna().all():
                    d["Tyr200_contacts_per_frame"] = round(
                        float(s[st].dropna().mean()), 2)
                entry[str(site)] = d
            f["sites"][k] = entry

    # ---- how much of the protein the DNA is actually touching ----
    # This is the measurement that decides between "the spring applies a force" and
    # "the DNA leans on the barrel", and the second is what the data supports.
    if "contact.atoms_in_contact_mean" in df.columns:
        f["contacts"] = {}
        for k, nm in SPRINGS + [("C1_ss50", "50 nt ssDNA spring"),
                                ("C2_ssarms", "two 25 nt arms"),
                                ("L3_loop27", "closed loop, 27 bp"),
                                ("L4_loop40", "closed loop, 40 bp"),
                                ("L5_loop40nick", "closed loop, 40 bp nicked")]:
            s = df[df.system == k]
            s = s[s["contact.atoms_in_contact_mean"].notna()]
            if not len(s):
                continue
            f["contacts"][k] = {
                "label": nm,
                "atoms_within_4.5A": agg(s, "contact.atoms_in_contact_mean"),
                "frames_with_contact": float(s["contact.frames_with_contact"].mean()),
                "residues_touched": float(s["contact.residues_ever_touched"].mean())}

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
