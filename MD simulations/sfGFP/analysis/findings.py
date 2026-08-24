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
    # Two kinds, and they must not be mixed.  A duplex is an elastic rod that pushes its
    # ends apart and gets Zocchi's f(x); a single strand is an entropic chain that pulls
    # them together and gets Marko-Siggia.  mechanics.py writes them under different keys
    # for exactly that reason, so a system with no mech.spring.* block is not missing data
    # -- it is a single strand, and it belongs in ss_springs.
    f["springs"], f["ss_springs"] = {}, {}
    for k, nm in SPRINGS:
        sub = df[df.system == k]
        if not len(sub):
            continue
        if "mech.ss_spring.n_nt" in sub and sub["mech.ss_spring.n_nt"].notna().any():
            f["ss_springs"][k] = {
                "label": nm,
                "n_nt": int(sub["mech.ss_spring.n_nt"].dropna().iloc[0]),
                "contour_nm": float(sub["mech.ss_spring.contour_nm"].dropna().iloc[0]),
                "persistence_nm": float(
                    sub["mech.ss_spring.persistence_nm"].dropna().iloc[0]),
                "natural_rms_nm": float(
                    sub["mech.ss_spring.natural_rms_end_to_end_nm"].dropna().iloc[0]),
                "state": str(sub["mech.ss_spring.state"].dropna().iloc[0]),
                "x_A": agg(sub, "mech.ss_spring.x_mean_A"),
                "tension_pN": agg(sub, "mech.ss_spring.tension_pN_mean"),
                "extension_A": agg(sub, "mech.anchor_CB_CB_A.mean"),
            }
            continue
        if "mech.spring.n_bp" not in sub or sub["mech.spring.n_bp"].isna().all():
            continue
        f["springs"][k] = {
            "label": nm, "n_bp": int(sub["mech.spring.n_bp"].dropna().iloc[0]),
            "nicked": bool(sub["mech.spring.nicked"].dropna().iloc[0]),
            "tau_c": float(sub["mech.spring.tau_c_pN_nm"].dropna().iloc[0]),
            "regime": str(sub["mech.spring.regime"].dropna().iloc[0]),
            "contour_nm": float(sub["mech.spring.contour_nm"].dropna().iloc[0]),
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
        # agg returns None when the column is entirely absent or NaN, which happens for
        # any system with no protein -- the bare duplexes.  Spreading None into a dict
        # raises, so ask first.
        a_ = agg(sub, "mech.anchor_CB_CB_A.mean") if len(sub) else None
        if a_:
            f["extension"][k] = {"label": nm, **a_}
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
    # The step the point-3 retraction turns on: what CLICKING alone does, unclicked to
    # clicked, no DNA either side.  Distinct from tet_vs_wt (installing the residues) and
    # opposite in sign -- installing them shortens the separation, clicking lengthens it
    # again -- so the two must not be quoted for each other.
    for nm, a, b in (("clicked_vs_tet", "S1_tet", "S2_clicked"),
                     ("clicked_vs_wt", "S0_wt", "S2_clicked")):
        x = df[(df.system == a) & (df.clamp_pN.fillna(0) == 0)][
            "mech.anchor_CB_CB_A.mean"].dropna().values
        y = df[(df.system == b) & (df.clamp_pN.fillna(0) == 0)][
            "mech.anchor_CB_CB_A.mean"].dropna().values
        w = welch(y, x)
        if w:
            f["extension"][nm] = w

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

    # ---- does closing the covalent loop matter? ----
    # The pairs below differ in ONE thing: whether a covalent path runs from one
    # attachment site out along the tether, through the DNA, and back to the other.
    # Same protein, same chemistry, same spring, same length -- so this isolates
    # topology from everything else, which is what the whole point-2 rebuild was for.
    TOPOLOGY_PAIRS = [("L3_loop27", "S3_spring27", "27 bp"),
                      ("L4_loop40", "S4_spring40", "40 bp"),
                      ("L5_loop40nick", "S5_spring40nick", "40 bp nicked"),
                      ("C1_ss50", "C2_ssarms", "50 nt ssDNA, sfGFP"),
                      ("U2_ss50", "U4_ssarms", "50 nt ssDNA, ubiquitin")]
    col = "mech.anchor_CB_CB_A.mean"
    pairs, diffs = {}, []
    for closed, open_, nm in TOPOLOGY_PAIRS:
        a = df[df.system == closed][col].dropna().values
        b = df[df.system == open_][col].dropna().values
        w = welch(a, b)
        if w is None:
            continue
        w["label"] = nm
        w["closed"], w["open"] = closed, open_
        pairs[nm] = w
        diffs.append(w["delta"])
    if pairs:
        entry = {"pairs": pairs, "n_pairs": len(pairs),
                 "closed_minus_open_A": [round(float(x), 3) for x in diffs],
                 "mean_closed_minus_open_A": round(float(np.mean(diffs)), 3),
                 "all_same_sign": bool(len(set(np.sign(diffs))) == 1)}
        def _paired(d):
            # Each pair is its own control, so the quantity with a sampling distribution
            # is the set of differences, tested against zero.
            d = np.array(d, float)
            tt = stats.ttest_1samp(d, 0.0)
            return {"n_pairs": len(d), "mean_A": round(float(d.mean()), 3),
                    "t": round(float(tt[0]), 3), "dof": len(d) - 1,
                    "p": round(float(tt[1]), 4), "resolved": bool(tt[1] < 0.05)}
        if len(diffs) > 2:
            entry["paired_across_pairs"] = _paired(diffs)
        # The three duplex pairs are the comparison the point-2 rebuild was designed to
        # make, so they are reported as their own planned test rather than as a subset
        # chosen after the fact.  The ssDNA pairs came later and are kept separate for
        # the same reason -- and they do NOT reproduce the duplex ordering, which is
        # stated rather than buried.
        ds = [w["delta"] for nm, w in pairs.items() if "bp" in nm]
        ss = [w["delta"] for nm, w in pairs.items() if "ssDNA" in nm]
        if len(ds) > 2:
            entry["duplex_pairs_planned"] = _paired(ds)
        if ss:
            entry["ssDNA_pairs"] = {"deltas_A": [round(float(x), 3) for x in ss],
                                    "reproduces_duplex_ordering": False}
        f["topology"] = entry

    # ---- is the compression a spring effect at all? ----
    # The whole project rests on the DNA acting as a spring, so the question deserves a
    # number rather than an impression: does the displacement of the attachment sites
    # track the force the spring is supposed to deliver, or does it not?
    vm = f["extension"].get("vs_matched_control") or {}
    if len(vm) > 3:
        force, delta, contact = [], [], []
        for k, v in vm.items():
            force.append((f["springs"].get(k, {}).get("force_pN", {}) or {})
                         .get("mean", 0.0) or 0.0)
            delta.append(v["delta"])
            c = ((f.get("contacts", {}).get(k, {}) or {})
                 .get("atoms_within_4.5A", {}) or {}).get("mean")
            contact.append(c)
        force, delta = np.array(force), np.array(delta)
        has_c = [i for i, c in enumerate(contact) if c is not None]
        rf = stats.pearsonr(force, delta)
        entry = {"n_constructs": len(vm),
                 "vs_modelled_force": {"r": round(float(rf[0]), 3),
                                       "p": round(float(rf[1]), 3)},
                 "mean_delta_no_spring": round(float(delta[force == 0].mean()), 3)
                 if (force == 0).any() else None,
                 "mean_delta_with_spring": round(float(delta[force > 0].mean()), 3)
                 if (force > 0).any() else None}
        if len(has_c) > 3:
            rc = stats.pearsonr([contact[i] for i in has_c], [delta[i] for i in has_c])
            entry["vs_contact"] = {"r": round(float(rc[0]), 3),
                                   "p": round(float(rc[1]), 3)}
        if (force == 0).any() and (force > 0).any():
            w = welch(delta[force == 0], delta[force > 0])
            if w:
                entry["no_spring_vs_spring"] = w
        f["is_it_a_spring"] = entry

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
