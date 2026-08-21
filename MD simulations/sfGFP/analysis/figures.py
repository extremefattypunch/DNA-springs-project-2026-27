#!/usr/bin/env python3
"""Report figures, drawn from analysis/summary.csv and the per-run CSVs.

Form choices follow from what each panel's data has to do:

* force response -- a relationship between an applied force and a deformation, so a
  line with markers and error bars across replicates.  The three chimeras are
  overlaid, but distinguished by *marker shape and direct label* rather than by three
  more hues: this is a scatter form, where only the first three palette slots are
  validated on the all-pairs list, and shape is a secondary encoding that keeps
  identity off colour entirely.
* RMSF -- change along a sequence, so a line per system with direct labels.
* H-bond occupancy -- identity across a handful of categories, so a dot plot: a
  grouped bar with five systems x seven bonds would be unreadable.
* bend profile -- change along the duplex, so a line per spring.

Every panel ships the CSV it was drawn from, which is also what discharges the
relief rule for the palette slots that fall below 3:1 on a light surface.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt          # noqa: E402
import numpy as np                        # noqa: E402
import pandas as pd                       # noqa: E402

import sys
sys.path.insert(0, str(Path(__file__).resolve().parent))
import vizstyle as VS                     # noqa: E402
from spring_model import Spring, TAU_C_NICKED, TAU_C_CONTINUOUS  # noqa: E402

VS.apply(matplotlib)
MARKERS = {"S3_spring27": "o", "S4_spring40": "s", "S5_spring40nick": "D"}


def label_right(ax, x, y, text, color):
    ax.annotate(text, (x, y), xytext=(4, 0), textcoords="offset points",
                va="center", ha="left", fontsize=8.5, color=VS.INK_2)


def fig_force_response(df, out):
    fig, ax = plt.subplots(figsize=(6.2, 4.0))
    clamp = df[df.system == "S6_clamp"].dropna(subset=["clamp_pN"])
    tbl = []
    if len(clamp):
        g = clamp.groupby("clamp_pN")["mech.anchor_CB_CB_A.mean"]
        f = np.array(sorted(g.groups))
        m = g.mean().reindex(f).values
        n = g.count().reindex(f).values
        sd = g.std().reindex(f).fillna(0.0).values
        sem = np.where(n > 1, sd / np.sqrt(n), np.nan)
        ax.errorbar(f, m, yerr=sem, color=VS.SERIES[0], marker="o", markersize=6,
                    capsize=3, lw=2.0, zorder=3, label="force clamp on 2x Tet2-Et")
        for fi, mi, ni in zip(f, m, n):
            tbl.append({"applied_pN": fi, "extension_A": round(mi, 3), "n_rep": int(ni)})
        if len(f) > 2 and np.isfinite(m).all():
            k, b = np.polyfit(f, m, 1)
            ax.plot(f, k * f + b, color=VS.SERIES[0], lw=1.0, ls=":", zorder=2)
            kappa = 1.0 / k if k else np.nan          # A per pN
            ax.annotate(f"compliance {k:.3f} " + r"$\mathrm{\AA\,pN^{-1}}$"
                        + f"\nstiffness {kappa:.1f} " + r"$\mathrm{pN\,\AA^{-1}}$",
                        (0.03, 0.95), xycoords="axes fraction", va="top",
                        fontsize=8.5, color=VS.INK_2)
    for sysname, sub in df[df.system.isin(MARKERS)].groupby("system"):
        fcol = "mech.spring.force_pN_mean"
        if fcol not in sub or sub[fcol].isna().all():
            continue
        fx = sub[fcol].mean()
        ext = sub["mech.anchor_CB_CB_A.mean"].mean()
        ax.plot([fx], [ext], marker=MARKERS[sysname], ms=9, color=VS.SERIES[1],
                mew=1.4, mfc="none", zorder=4)
        label_right(ax, fx, ext, VS.SYSTEM_LABEL[sysname], VS.SERIES[1])
        tbl.append({"applied_pN": round(float(fx), 3),
                    "extension_A": round(float(ext), 3),
                    "n_rep": int(len(sub)), "label": VS.SYSTEM_LABEL[sysname]})
    for name, y in (("WT sfGFP", df[df.system == "S0_wt"]["mech.anchor_CB_CB_A.mean"]),
                    ("2x Tet2-Et, no force",
                     df[df.system == "S1_tet"]["mech.anchor_CB_CB_A.mean"])):
        if len(y.dropna()):
            ax.axhline(y.mean(), color=VS.INK_MUTED, lw=0.9, ls="--", zorder=1)
            ax.annotate(name, (ax.get_xlim()[1], y.mean()), xytext=(-2, 3),
                        textcoords="offset points", ha="right", fontsize=8,
                        color=VS.INK_MUTED)
    ax.set_xlabel("force applied between the attachment sites (pN)")
    ax.set_ylabel(r"C$\beta$–C$\beta$ separation, Asp133–Asn149 ($\mathrm{\AA}$)")
    ax.set_title("Force response of the attachment sites")
    ax.grid(axis="y", alpha=0.6)
    ax.legend(loc="lower right")
    fig.savefig(out / "fig1_force_response.png")
    pd.DataFrame(tbl).to_csv(out / "fig1_force_response.csv", index=False)
    plt.close(fig)
    return "fig1_force_response"


def fig_rmsf(df, runs, out):
    fig, ax = plt.subplots(figsize=(7.6, 3.8))
    frames = []
    # Drop the ACE/NME caps and the two residues either side of them.  Their RMSF is
    # 4-11 A -- they are free ends, not structure -- and leaving them in squeezes every
    # real feature into the bottom tenth of the axis.
    lo, hi = 8, 229
    for sysname in ["S0_wt", "S1_tet", "S3_spring27", "S4_spring40", "S5_spring40nick"]:
        files = sorted(Path(runs).glob(f"{sysname}/rep*/03_production/analysis_rmsf.csv"))
        if not files:
            continue
        d = pd.concat([pd.read_csv(f) for f in files])
        d = d[(d.resnum >= lo) & (d.resnum <= hi)]
        g = d.groupby("resnum")["rmsf_A"].mean()
        ax.plot(g.index, g.values, color=VS.SYSTEM_COLOR[sysname], lw=1.6,
                label=f"{VS.SYSTEM_LABEL[sysname]} (n={len(files)})")
        label_right(ax, g.index[-1], g.values[-1], VS.SYSTEM_LABEL[sysname],
                    VS.SYSTEM_COLOR[sysname])
        frames.append(g.rename(sysname))
    # Stagger the marker labels: 148 and 149 are one residue apart and their labels
    # printed on top of each other as "H1498".
    top = ax.get_ylim()[1]
    for x, txt, dy in ((66, "CRO 66", -2), (133, "Asp133", -2),
                       (148, "His148", -13), (149, "Asn149", -24)):
        ax.axvline(x, color=VS.INK_MUTED, lw=0.8, ls=":", zorder=0)
        ax.annotate(txt, (x, top), xytext=(3, dy), textcoords="offset points",
                    va="top", fontsize=7.5, color=VS.INK_MUTED)
    ax.set_xlim(lo - 6, hi + 22)
    ax.set_xlabel("residue (2B3P numbering; ACE/NME caps and residues 2–7, 230–232 omitted)")
    ax.set_ylabel(r"C$\alpha$ RMSF ($\mathrm{\AA}$)")
    ax.set_title("Backbone mobility, aligned on the barrel core")
    ax.grid(axis="y", alpha=0.6)
    ax.legend(loc="upper right", ncols=1, bbox_to_anchor=(1.0, 1.02))
    fig.savefig(out / "fig2_rmsf.png")
    if frames:
        pd.concat(frames, axis=1).to_csv(out / "fig2_rmsf.csv")
    plt.close(fig)
    return "fig2_rmsf"


def fig_hbonds(df, out):
    cols = [c for c in df.columns
            if c.startswith("chromo.hbonds.") and c.endswith(".occupancy")]
    if not cols:
        return None
    labels = [c[len("chromo.hbonds."):-len(".occupancy")] for c in cols]
    systems = [s for s in ["S0_wt", "S1_tet", "S3_spring27", "S4_spring40",
                           "S5_spring40nick"] if s in set(df.system)]
    fig, ax = plt.subplots(figsize=(7.0, 0.55 * len(labels) + 2.0))
    y = np.arange(len(labels))
    tbl = {}
    for k, s in enumerate(systems):
        sub = df[df.system == s]
        v = np.array([sub[c].mean() * 100 for c in cols], float)
        ax.plot(v, y + (k - (len(systems) - 1) / 2) * 0.13, MARKERS.get(s, "o"),
                ms=7, color=VS.SYSTEM_COLOR[s], mfc=VS.SYSTEM_COLOR[s],
                mec=VS.SURFACE, mew=1.2, label=VS.SYSTEM_LABEL[s], zorder=3)
        tbl[VS.SYSTEM_LABEL[s]] = np.round(v, 2)
    for yy in y:
        ax.axhline(yy, color=VS.GRID, lw=0.8, zorder=0)
    ax.set_yticks(y, [l.replace(" -> ", " → ") for l in labels], fontsize=8)
    ax.set_xlim(-3, 103)
    ax.set_xlabel("occupancy (%)  —  H···A < 2.5 Å and D–H···A > 120°")
    ax.set_title("Chromophore hydrogen-bond network")
    ax.grid(axis="x", alpha=0.6)
    ax.legend(loc="lower left", ncols=2)
    fig.savefig(out / "fig3_chromophore_hbonds.png")
    pd.DataFrame(tbl, index=labels).to_csv(out / "fig3_chromophore_hbonds.csv")
    plt.close(fig)
    return "fig3_chromophore_hbonds"


def fig_bend(runs, out):
    fig, ax = plt.subplots(figsize=(6.6, 3.6))
    any_data, frames = False, []
    for sysname in ["S3_spring27", "S4_spring40", "S5_spring40nick"]:
        files = sorted(Path(runs).glob(
            f"{sysname}/rep*/03_production/analysis_bend_profile.csv"))
        if not files:
            continue
        d = pd.concat([pd.read_csv(f) for f in files]).groupby("bp").mean()
        ax.plot(d.index, d["bend_deg_mean"], color=VS.SYSTEM_COLOR[sysname], lw=1.8,
                label=f"{VS.SYSTEM_LABEL[sysname]} (n={len(files)})")
        ax.fill_between(d.index, d["bend_deg_mean"] - d["bend_deg_sd"],
                        d["bend_deg_mean"] + d["bend_deg_sd"],
                        color=VS.SYSTEM_COLOR[sysname], alpha=0.14, lw=0)
        label_right(ax, d.index[-1], d["bend_deg_mean"].iloc[-1],
                    VS.SYSTEM_LABEL[sysname], VS.SYSTEM_COLOR[sysname])
        frames.append(d["bend_deg_mean"].rename(sysname))
        any_data = True
    if not any_data:
        plt.close(fig)
        return None
    ax.set_xlabel("base pair along the duplex")
    ax.set_ylabel("axis bend per step (deg)")
    ax.set_title("Where the spring bends — a kink would be a local spike")
    ax.grid(axis="y", alpha=0.6)
    ax.legend(loc="upper left")
    fig.savefig(out / "fig4_dna_bend.png")
    pd.concat(frames, axis=1).to_csv(out / "fig4_dna_bend.csv")
    plt.close(fig)
    return "fig4_dna_bend"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--summary", default="analysis/summary.csv")
    ap.add_argument("--runs", default="data/runs")
    ap.add_argument("--out", default="figures")
    a = ap.parse_args()
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    df = pd.read_csv(a.summary)
    made = [fig_force_response(df, out), fig_rmsf(df, a.runs, out),
            fig_hbonds(df, out), fig_bend(a.runs, out)]
    made = [m for m in made if m]
    print("wrote: " + ", ".join(f"{m}.png" for m in made))
    (out / "figures.json").write_text(json.dumps(made, indent=2))


if __name__ == "__main__":
    main()
