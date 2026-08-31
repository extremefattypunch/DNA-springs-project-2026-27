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
MEAN_MARKERS = ["o", "s", "D"]


def label_right(ax, x, y, text, color, dx=9, dy=0):
    ax.annotate(text, (x, y), xytext=(dx, dy), textcoords="offset points",
                va="center", ha="left", fontsize=8.5, color=VS.INK_2)


def labels_right(ax, items, dx=9, min_gap_pts=11.5):
    """Right-edge direct labels, pushed apart so near-equal series stay readable.

    Series that end within a few tenths of each other print on top of one another
    otherwise -- the two 40 bp springs did exactly that.  Positions are nudged in
    display points, so the spacing holds whatever the data range is.
    """
    if not items:
        return
    items = sorted(items, key=lambda it: it[1])
    fig = ax.figure
    fig.canvas.draw()
    trans = ax.transData
    pts = [trans.transform((x, y))[1] for x, y, _, _ in items]
    for i in range(1, len(pts)):
        if pts[i] - pts[i - 1] < min_gap_pts:
            pts[i] = pts[i - 1] + min_gap_pts
    for (x, y, text, color), py in zip(items, pts):
        y_disp = trans.transform((x, y))[1]
        ax.annotate(text, (x, y), xytext=(dx, py - y_disp),
                    textcoords="offset points", va="center", ha="left",
                    fontsize=8.5, color=VS.INK_2)


def fig_force_response(df, out):
    fig, ax = plt.subplots(figsize=(6.2, 4.0))
    clamp = df[df.system.isin(("S6_clamp", "S2_clicked"))].dropna(subset=["clamp_pN"])
    tbl = []
    for k_sys, (sysname, style) in enumerate(
            (("S6_clamp", dict(color=VS.SERIES[0], marker="o",
                               label="axial clamp on 2× Tet2-Et")),
             ("S2_clicked", dict(color=VS.SERIES[2], marker="^",
                                 label="axial clamp on the clicked tether")))):
        sub = clamp[clamp.system == sysname]
        if not len(sub):
            continue
        # every replicate, faint -- the spread is the story here, not the mean
        ax.plot(sub["clamp_pN"], sub["mech.anchor_CB_CB_A.mean"], style["marker"],
                ms=4, color=style["color"], alpha=0.35, mew=0, zorder=2)
        g = sub.groupby("clamp_pN")["mech.anchor_CB_CB_A.mean"]
        f = np.array(sorted(g.groups))
        m = g.mean().reindex(f).values
        n = g.count().reindex(f).values
        sd = g.std().reindex(f).fillna(0.0).values
        sem = np.where(n > 1, sd / np.sqrt(n), np.nan)
        ax.errorbar(f, m, yerr=sem, color=style["color"], marker=style["marker"],
                    markersize=7, capsize=3, lw=2.0, zorder=3, label=style["label"])
        for fi, mi, ni, si in zip(f, m, n, sem):
            tbl.append({"series": sysname, "applied_pN": fi,
                        "extension_A": round(mi, 3), "sem_A": round(float(si), 4),
                        "n_rep": int(ni)})
        if sysname == "S6_clamp" and len(sub) > 3:
            x, y = sub["clamp_pN"].values.astype(float), \
                sub["mech.anchor_CB_CB_A.mean"].values
            kk, bb = np.polyfit(x, y, 1)
            se = (np.sqrt(((y - (kk * x + bb)) ** 2).sum() / (len(x) - 2))
                  / np.sqrt(((x - x.mean()) ** 2).sum()))
            lo, hi = kk - 2.04 * se, kk + 2.04 * se
            ax.plot(f, kk * f + bb, color=style["color"], lw=1.0, ls=":", zorder=2)
            kap = 10.0 / (hi * 4.2) if hi > 0 else np.inf
            ax.annotate(
                f"compliance {kk:.4f} ± {se:.4f} " + r"$\mathrm{\AA\,pN^{-1}}$"
                + f"\n95% CI [{lo:.4f}, {hi:.4f}] — includes zero"
                + f"\nso only a bound: " + r"$\kappa$" + f" > {kap:.0f} "
                + r"$k_\mathrm{B}T\,\mathrm{nm^{-2}}$",
                (0.38, 0.04), xycoords="axes fraction", va="bottom",
                fontsize=8.2, color=VS.INK_2)
    for sysname, sub in df[df.system.isin(MARKERS)].groupby("system"):
        fcol = "mech.spring.force_pN_mean"
        if fcol not in sub or sub[fcol].isna().all():
            continue
        fx = sub[fcol].mean()
        e = sub["mech.anchor_CB_CB_A.mean"]
        sem = e.std() / np.sqrt(len(e)) if len(e) > 1 else np.nan
        ax.plot(sub[fcol], e, MARKERS[sysname], ms=4, color=VS.SERIES[1],
                alpha=0.35, mew=0, zorder=3)
        ax.errorbar([fx], [e.mean()], yerr=[sem], marker=MARKERS[sysname], ms=10,
                    color=VS.SERIES[1], mew=1.6, mfc="none", capsize=3, lw=0,
                    elinewidth=1.4, zorder=5)
        label_right(ax, fx, e.mean(), VS.SYSTEM_LABEL[sysname], VS.SERIES[1], dx=11)
        tbl.append({"series": sysname, "applied_pN": round(float(fx), 3),
                    "extension_A": round(float(e.mean()), 3),
                    "sem_A": round(float(sem), 4), "n_rep": int(len(sub))})
    for name, y in (("WT sfGFP", df[df.system == "S0_wt"]["mech.anchor_CB_CB_A.mean"]),
                    ("2× Tet2-Et, unloaded",
                     df[df.system == "S1_tet"]["mech.anchor_CB_CB_A.mean"])):
        if len(y.dropna()):
            ax.axhline(y.mean(), color=VS.INK_MUTED, lw=0.9, ls="--", zorder=1)
            ax.annotate(name, (ax.get_xlim()[1], y.mean()), xytext=(-3, 4),
                        textcoords="offset points", ha="right", fontsize=8,
                        color=VS.INK_MUTED, zorder=6)
    ax.set_xlabel("force applied between the attachment sites (pN)")
    ax.set_ylabel(r"C$\beta$–C$\beta$ separation, Asp133–Asn149 ($\mathrm{\AA}$)")
    ax.set_title("Force response of the attachment sites")
    ax.set_xlim(left=-1.2)
    ax.grid(axis="y", alpha=0.6)
    ax.legend(loc="upper left")
    fig.savefig(out / "fig1_force_response.png")
    pd.DataFrame(tbl).to_csv(out / "fig1_force_response.csv", index=False)
    plt.close(fig)
    return "fig1_force_response"


def fig_rmsf(df, runs, out):
    fig, ax = plt.subplots(figsize=(7.6, 3.8))
    frames, pending = [], []
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
        pending.append((g.index[-1], float(g.values[-1]), VS.SYSTEM_LABEL[sysname],
                        VS.SYSTEM_COLOR[sysname]))
        frames.append(g.rename(sysname))
    labels_right(ax, pending)
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
    systems = [s for s in ["S0_wt", "S1_tet", "S2_clicked", "S3_spring27",
                           "S4_spring40", "S5_spring40nick"] if s in set(df.system)]
    fig, ax = plt.subplots(figsize=(8.4, 0.62 * len(labels) + 2.2))
    y = np.arange(len(labels))
    tbl = {}
    for k, s in enumerate(systems):
        sub = df[df.system == s]
        off = (k - (len(systems) - 1) / 2) * 0.15
        v = np.array([sub[c].mean() * 100 for c in cols], float)
        # every replicate as a faint tick: two of these contacts are bistable, either
        # ~0% or ~90% in a given replicate, and a mean of three coin flips is not a
        # measurement.  Showing the replicates makes that visible instead of hiding it.
        for r_i, (_, row) in enumerate(sub.iterrows()):
            rv = np.array([row[c] * 100 for c in cols], float)
            ax.plot(rv, y + off, "|", ms=7, color=VS.SYSTEM_COLOR[s], alpha=0.45,
                    mew=1.2, zorder=2)
        ax.plot(v, y + off, MARKERS.get(s, "o"),
                ms=6.5, color=VS.SYSTEM_COLOR[s], mfc=VS.SYSTEM_COLOR[s],
                mec=VS.SURFACE, mew=1.2, label=VS.SYSTEM_LABEL[s], zorder=3)
        tbl[VS.SYSTEM_LABEL[s]] = np.round(v, 2)
        tbl[VS.SYSTEM_LABEL[s] + " range"] = [
            f"{sub[c].min()*100:.0f}–{sub[c].max()*100:.0f}" for c in cols]
    for yy in y:
        ax.axhline(yy, color=VS.GRID, lw=0.8, zorder=0)
    ax.set_yticks(y, [l.replace(" -> ", " → ").replace(" <- ", " ← ")
                      for l in labels], fontsize=8.2)
    ax.set_xlim(-3, 103)
    ax.set_xlabel("occupancy (%)   ·   ticks = replicates, filled = mean")
    ax.set_title("Chromophore hydrogen-bond network")
    ax.grid(axis="x", alpha=0.6)
    ax.legend(loc="lower center", ncols=3, bbox_to_anchor=(0.5, 1.03),
              borderaxespad=0)
    fig.savefig(out / "fig3_chromophore_hbonds.png", bbox_inches="tight")
    pd.DataFrame(tbl, index=labels).to_csv(out / "fig3_chromophore_hbonds.csv")
    plt.close(fig)
    return "fig3_chromophore_hbonds"


def fig_bend(runs, out):
    fig, ax = plt.subplots(figsize=(7.0, 3.8))
    any_data, frames, pending = False, [], []
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
        pending.append((d.index[-1], float(d["bend_deg_mean"].iloc[-1]),
                        VS.SYSTEM_LABEL[sysname], VS.SYSTEM_COLOR[sysname]))
        frames.append(d["bend_deg_mean"].rename(sysname))
        any_data = True
    if not any_data:
        plt.close(fig)
        return None
    labels_right(ax, pending)
    ax.set_xlabel("base pair along the duplex")
    ax.set_ylabel("axis bend per step (deg)")
    ax.set_title("Where the spring bends — a kink would be a local spike")
    ax.grid(axis="y", alpha=0.6)
    ax.legend(loc="upper left")
    fig.savefig(out / "fig4_dna_bend.png")
    pd.concat(frames, axis=1).to_csv(out / "fig4_dna_bend.csv")
    plt.close(fig)
    return "fig4_dna_bend"



def fig_kink_threshold(runs, outdir):
    """How often a bent duplex breaks an internal base pair, against curvature.

    Form.  The per-run measurement is bimodal, not continuous: a run either keeps every
    interior pair shut for its whole trajectory or it opens one and holds it open 30-100%
    of frames.  Averaging that gives a number no run ever exhibits, and joining those
    averages gave a sawtooth that read as violent structure where the truth is a smooth
    rise in a probability.  So the measure plotted is the one that actually varies with
    curvature -- the FRACTION OF RUNS in which any interior pair opens -- with the
    individual runs shown underneath as the binary data they are.  One measure, one axis.

    The band is a Wilson 68% interval, which is the right interval for a proportion from
    four runs: a normal approximation would put the lower edge below zero and the upper
    edge above one, and at n = 4 that is not a rounding detail.

    Curvature is the nominal circular-arc equivalent of the separation the duplex
    actually adopted.  The windows load it at its two ends, where the chimera's tethers
    distribute the bend, so this coordinate is comparable rather than identical -- the
    caption says so.

    One series, because the claim is about the force field rather than about any one
    duplex; which construct each run came from is in the shipped CSV.
    """
    f = Path(__file__).resolve().parent / "force_extension.csv"
    if not f.exists():
        print("  kink threshold: no force_extension.csv yet")
        return None
    d = pd.read_csv(f)
    if "wc.max_internal_occupancy" not in d.columns:
        print("  kink threshold: no internal-opening column yet")
        return None
    if "box_ok" in d.columns:
        d = d[d.box_ok.astype(str).str.lower().isin(("true", "1"))]
    if "box_ok_on_measured_r" not in d.columns and "box_width_A" in d.columns:
        d = d[(d.mean_r_A / d.box_width_A) <= 0.85]
    d = d[d["wc.max_internal_occupancy"].notna()]
    if len(d) < 4:
        print("  kink threshold: too few valid windows")
        return None

    RISE = 3.376
    # n_bp comes from wc.pairs, which force_extension.py MEASURES from the
    # topology, rather than from a map keyed on system name.  The map had to
    # be edited in five files whenever a system was added, and a system
    # missing from it silently fell back to 27.
    NBP = (d.groupby("system")["wc.pairs"].max().dropna()
           .astype(int).to_dict() if "wc.pairs" in d else {})

    def per_bp(chord, n_bp):
        L = (n_bp - 1) * RISE
        lo, hi = L / (2 * np.pi), 1e7
        for _ in range(300):
            mid = 0.5 * (lo + hi)
            lo, hi = ((mid, hi) if 2 * mid * np.sin(L / (2 * mid)) < chord
                      else (lo, mid))
        return np.degrees(RISE / (0.5 * (lo + hi)))

    xcol = "mean_r_A" if "mean_r_A" in d.columns else "x0_A"
    d = d.assign(deg_bp=[per_bp(float(x), NBP.get(s, 27))
                         for x, s in zip(d[xcol], d.system)],
                 opened=d["wc.max_internal_occupancy"] > 0)

    def wilson(k, n, z=1.0):
        """68% Wilson interval -- correct at n = 4, where the normal one is not."""
        if not n:
            return (0.0, 1.0)
        ph, z2 = k / n, z * z
        c = (ph + z2 / (2 * n)) / (1 + z2 / n)
        h = (z / (1 + z2 / n)) * np.sqrt(ph * (1 - ph) / n + z2 / (4 * n * n))
        return max(0.0, c - h), min(1.0, c + h)

    lo_e, hi_e = np.floor(d.deg_bp.min()), np.ceil(d.deg_bp.max())
    edges = np.arange(lo_e, hi_e + 1e-9, 1.0)
    xs, ys, los, his, ns = [], [], [], [], []
    for a, b in zip(edges[:-1], edges[1:]):
        g = d[(d.deg_bp >= a) & (d.deg_bp < b)]
        if not len(g):
            continue
        k, n = int(g.opened.sum()), len(g)
        l, h = wilson(k, n)
        xs.append(0.5 * (a + b)); ys.append(k / n); los.append(l); his.append(h)
        ns.append(n)

    fig, ax = plt.subplots(figsize=(7.6, 4.5))
    ax.axvspan(4.8, 5.2, color=VS.SERIES[0], alpha=0.13, lw=0, zorder=0)
    ax.axvspan(8.6, 9.3, color=VS.GRID, alpha=0.9, lw=0, zorder=0)
    ax.annotate("our chimera\n4.8–5.2°/bp\nnone open in 3 × 50 ns", (5.0, 0.55),
                ha="center", va="top", fontsize=8.5, color=VS.SERIES[0], zorder=5)
    ax.annotate("Curuksu et al.\ntype II kink\n8.6–9.3", (8.95, 0.42), ha="center",
                va="top", fontsize=8.5, color=VS.INK_MUTED, zorder=5)

    # the individual runs, as the binary observations they are
    ax.plot(d.deg_bp, d.opened.astype(float), "|", ms=11, mew=1.6,
            color=VS.SERIES[0], alpha=0.4, zorder=2)

    ax.fill_between(xs, los, his, color=VS.SERIES[0], alpha=0.18, lw=0, zorder=1)
    ax.plot(xs, ys, "-o", lw=2.2, ms=8, color=VS.SERIES[0], mec="white", mew=1.6,
            zorder=3)
    for x, y, n in zip(xs, ys, ns):
        ax.annotate(f"{n}", (x, y), xytext=(0, 11), textcoords="offset points",
                    ha="center", fontsize=7.5, color=VS.INK_MUTED, zorder=5)

    ax.set_xlabel("curvature of the duplex (degrees per base-pair step, nominal "
                  "circular arc)")
    ax.set_ylabel("fraction of runs in which an\ninterior base pair opens")
    ax.set_ylim(-0.09, 1.16)
    ax.set_xlim(4.2, 10.4)
    ax.set_yticks([0, 0.25, 0.5, 0.75, 1.0])
    ax.set_title("Base-pair opening becomes certain only at the published kink curvature")
    ax.grid(axis="y", alpha=0.6)
    fig.tight_layout()
    fig.savefig(Path(outdir) / "fig5_kink_threshold.png", dpi=200, bbox_inches="tight")
    pd.DataFrame({"deg_bp_bin_centre": xs, "fraction_opening": ys,
                  "wilson68_lo": los, "wilson68_hi": his, "n_runs": ns}).to_csv(
        Path(outdir) / "fig5_kink_threshold.csv", index=False)
    keep = [c for c in ("system", "x0_A", "replicate", "k_kJ_mol_nm2", "mean_r_A",
                        "deg_bp", "wc.max_internal_occupancy",
                        "wc.most_open_internal_bp") if c in d.columns]
    d[keep].to_csv(Path(outdir) / "fig5_kink_threshold_runs.csv", index=False)
    plt.close(fig)
    return "fig5_kink_threshold"


def fig_force_extension(runs, outdir):
    """Measured force against measured extension, for the three bare duplexes.

    This is the panel that answers "is the force a model or a measurement".  Every other
    force in this project is Zocchi's f(x) evaluated at a simulated distance; these come
    out of the force field itself, as the mean restraint force <k(r - x0)>.

    Form: force is a magnitude against a continuum, so points with a line, one measure on
    one axis.  Three constructs is three series, so a legend AND direct labels -- identity
    is never carried by colour alone.  The two reference curves are annotation, not
    series: they wear grey and are labelled in place.

    The two restraint stiffnesses are the internal control rather than a series.  A soft
    restraint lets the duplex settle several angstrom past its target, so if the force
    were an artefact of how hard we pulled, the k = 33 and k = 200 points would sit on
    different curves.  Plotted against the distance the duplex actually adopted, they
    have to fall on one curve, and whether they do is visible here.  Marker shape carries
    which is which; the reader is not asked to compare them by colour.
    """
    f = Path(__file__).resolve().parent / "force_extension.csv"
    if not f.exists():
        print("  force-extension: no force_extension.csv yet")
        return None
    d = pd.read_csv(f)
    need = {"mean_r_A", "force_pN", "box_ok", "k_kJ_mol_nm2"}
    if not need <= set(d.columns):
        print("  force-extension: csv predates the box check")
        return None
    d = d[d.box_ok.astype(str).str.lower().isin(("true", "1"))]
    # box_ok in an older csv was judged on the TARGET separation only.  The measured
    # ratio is mean_r/box_width exactly, so derive it rather than skipping the check --
    # a missing column must not silently relax a guard.
    # box_ok already carries the measured-separation test when the csv is current;
    # only older files need it derived, and deriving it from rounded columns is what
    # made one boundary run appear in some tables and not others.
    if "box_ok_on_measured_r" not in d.columns:
        d = d[(d.mean_r_A / d.box_width_A) <= 0.85]
    if not len(d):
        print("  force-extension: no valid windows")
        return None

    # n_bp comes from wc.pairs, which force_extension.py MEASURES from the
    # topology, rather than from a map keyed on system name.  The map had to
    # be edited in five files whenever a system was added, and a system
    # missing from it silently fell back to 27.
    NBP = (d.groupby("system")["wc.pairs"].max().dropna()
           .astype(int).to_dict() if "wc.pairs" in d else {})
    B, RISE_NM = 200.0, 0.33
    # Grouped by the PHYSICAL duplex, not the system name: D_27bp and B_27bp are the
    # same 27 bp molecule in a 78 A and an 85 A box, so their force curves are one
    # curve, and the windows they share are a cross-check between boxes rather than two
    # unrelated series.  `duplex` is written by force_extension.py from the measured
    # base-pair count, so a new box needs no edit here.
    if "duplex" not in d.columns:
        d = d.assign(duplex=d.system.map(
            lambda s: f"{NBP.get(s, 27)} bp" + (", nicked" if "nick" in s.lower() else "")))
    order = sorted(set(d.duplex), key=lambda s: (int(s.split()[0]), "nick" in s))
    NICE = {s: s for s in order}

    fig, ax = plt.subplots(figsize=(7.6, 4.6))

    # --- annotation: the two things the measurement can be compared against ---
    for n_bp in sorted({int(s.split()[0]) for s in order}):
        eul = np.pi ** 2 * B / (((n_bp - 1) * RISE_NM) ** 2)
        ax.axhline(eul, color=VS.INK_MUTED, lw=1.1, ls=(0, (4, 3)), alpha=0.8, zorder=1)
        ax.annotate(f"Euler load, {n_bp} bp — {eul:.0f} pN", (0.012, eul),
                    xycoords=("axes fraction", "data"), xytext=(0, 4),
                    textcoords="offset points", ha="left", va="bottom", fontsize=8.5,
                    color=VS.INK_MUTED, zorder=5)
    pass

    # Zocchi's kinked-branch range per duplex length, keyed on base pairs so it follows
    # the duplex rather than the system name.
    zoc = {27: (2.4, 6.6), 40: (1.8, 3.2)}
    lens = {int(s.split()[0]) for s in order}
    lo = min(zoc[L][0] for L in lens if L in zoc)
    hi = max(zoc[L][1] for L in lens if L in zoc)
    ax.axhspan(lo, hi, color=VS.GRID, alpha=0.85, lw=0, zorder=0)
    ax.annotate("what Zocchi’s kinked branch predicts for these same duplexes",
                (0.012, lo + 0.35), xycoords=("axes fraction", "data"), ha="left",
                va="bottom", fontsize=8.5, color=VS.INK_MUTED, zorder=5)

    # --- the data ---
    # Drawn as a PLATEAU, not as a curve through the points.  Regressing force on
    # separation gives p = 0.29-0.54 for all three constructs: there is no resolvable
    # slope, which is the signature of post-buckling -- a buckled rod carries a roughly
    # constant load however far you compress it.  A line joining consecutive points
    # would invent a trend out of which restraint stiffness happened to land where, so
    # each construct gets the flat mean it actually measures, banded by its standard
    # error, spanning only the separations sampled.
    for i, s in enumerate(order):
        g = d[d.duplex == s]
        col = VS.SERIES[i]
        for k, mk in ((33.0, "s"), (200.0, "o")):
            gg = g[g.k_kJ_mol_nm2 == k]
            if len(gg):
                ax.plot(gg.mean_r_A, gg.force_pN, mk, ms=7, color=col, alpha=0.5,
                        lw=0, mec="white", mew=0.9, zorder=2,
                        label=NICE[s] if k == 200.0 else None)
        m, e = float(g.force_pN.mean()), float(g.force_pN.sem())
        r0, r1 = float(g.mean_r_A.min()), float(g.mean_r_A.max())
        ax.fill_between([r0, r1], m - e, m + e, color=col, alpha=0.20, lw=0, zorder=1)
        ax.plot([r0, r1], [m, m], "-", lw=2.4, color=col, zorder=3,
                solid_capstyle="round")
        ax.annotate(f"{NICE[s]}   {m:.1f} pN", (r1, m), xytext=(8, 0),
                    textcoords="offset points", ha="left", va="center",
                    fontsize=9, color=col, zorder=5)

    ax.plot([], [], "o", color=VS.INK_MUTED, ms=6, lw=0, label="stiff restraint (k = 200)")
    ax.plot([], [], "s", color=VS.INK_MUTED, ms=6, lw=0, label="soft restraint (k = 33)")
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.17), frameon=False,
              fontsize=8.5, ncol=5, handletextpad=0.4, columnspacing=1.4)
    ax.set_xlabel("end-to-end separation the duplex actually adopted (Å)")
    ax.set_ylabel("force the duplex exerts (pN)\nmeasured as ⟨k(r − x₀)⟩")
    ax.set_title("Force measured from the force field, not from a model")
    ax.set_xlim(float(d.mean_r_A.min()) - 3, float(d.mean_r_A.max()) + 17)
    ax.set_ylim(0, max(40.0, float(d.force_pN.max()) + 4))
    ax.grid(axis="y", alpha=0.6)
    fig.tight_layout()
    fig.savefig(Path(outdir) / "fig6_force_extension.png", dpi=200,
                bbox_inches="tight")
    keep = [c for c in ("duplex", "system", "box_A", "x0_A", "replicate",
                        "k_kJ_mol_nm2", "mean_r_A", "force_pN", "force_sem_pN")
            if c in d.columns]
    d[keep].to_csv(Path(outdir) / "fig6_force_extension.csv", index=False)
    plt.close(fig)
    return "fig6_force_extension"


def fig_headline(runs, outdir):
    """Site displacement against the force that is supposed to cause it.

    This is the project's central claim and until now it existed only as a table, while
    the ubiquitin half of the evidence -- five systems, fifteen replicates -- appeared in
    no figure at all.  Figure 1 shows three constructs, and they are the superseded
    two-arm builds.

    Form: the question is whether two continuous quantities are related, so a scatter,
    one measure per axis.  The point of the panel is a NON-relationship, which is a hard
    thing to show honestly -- a fitted line through noise invites the eye to see a trend.
    So the fit is drawn only as a faint band with its correlation stated in words, and
    the zero line is the reference that matters: everything below it is compression.

    Two proteins is two series, so a legend, and each is also direct-labelled -- identity
    is never carried by colour alone.  The constructs carrying no spring at all sit at
    x = 0 and are ringed, because they are the comparison that decides the question: if
    displacement were a force effect they should show none.
    """
    fj = Path(__file__).resolve().parent / "findings.json"
    if not fj.exists():
        print("  headline: no findings.json")
        return None
    f = json.loads(fj.read_text())
    vm = f.get("extension", {}).get("vs_matched_control") or {}
    sp = f.get("springs") or {}
    lab = f.get("labels") or {}
    if len(vm) < 4:
        print("  headline: too few constructs")
        return None

    pts = []
    for k, v in vm.items():
        force = (sp.get(k, {}).get("force_pN", {}) or {}).get("mean", 0.0) or 0.0
        # sem of the difference, recovered from Welch's t: delta / t
        sem = abs(v["delta"] / v["t"]) if v.get("t") else 0.0
        pts.append({"k": k, "label": lab.get(k, k), "force": force,
                    "delta": v["delta"], "sem": sem,
                    "fam": v.get("family", "sfGFP"), "p": v.get("p", 1.0)})

    fams = sorted({p_["fam"] for p_ in pts})
    colour = {fam: VS.SERIES[i] for i, fam in enumerate(fams)}
    fig, ax = plt.subplots(figsize=(7.8, 4.8))

    ax.axhline(0.0, color=VS.INK_MUTED, lw=1.1, zorder=1)
    # kept left: the right-hand end of this line is where the one positive construct
    # sits, and the label was landing on its error bar
    ax.annotate("no change from the construct's own clicked control",
                (0.30, 0.0), xycoords=("axes fraction", "data"), xytext=(0, 6),
                textcoords="offset points", ha="center", va="bottom", fontsize=8.5,
                color=VS.INK_MUTED, zorder=5)

    sq = f.get("is_it_a_spring") or {}
    r = (sq.get("vs_modelled_force") or {}).get("r")
    pv = (sq.get("vs_modelled_force") or {}).get("p")
    if r is not None and len(pts) > 2:
        xs = np.array([p_["force"] for p_ in pts])
        ys = np.array([p_["delta"] for p_ in pts])
        b, a0 = np.polyfit(xs, ys, 1)
        gx = np.linspace(-0.4, xs.max() * 1.08, 50)
        # Faint, and unlabelled as a "fit": the honest reading is that there is none.
        ax.plot(gx, a0 + b * gx, "-", lw=1.4, color=VS.GRID, zorder=1)

    for fam in fams:
        g = [p_ for p_ in pts if p_["fam"] == fam]
        c = colour[fam]
        ax.errorbar([p_["force"] for p_ in g], [p_["delta"] for p_ in g],
                    yerr=[p_["sem"] for p_ in g], fmt="o", ms=9, color=c, ecolor=c,
                    elinewidth=1.3, capsize=3, mec="white", mew=1.5, zorder=3,
                    label=fam, lw=0)
        # ring the springless constructs: they are the ones that decide the question
        noload = [p_ for p_ in g if p_["force"] == 0]
        if noload:
            ax.plot([p_["force"] for p_ in noload], [p_["delta"] for p_ in noload],
                    "o", ms=17, mfc="none", mec=c, mew=1.6, alpha=0.7, zorder=2)

    # label the three constructs the argument actually names: the strongest spring, the
    # most compressed, and the single one that is not compressed at all
    hi = max(pts, key=lambda p_: p_["force"])
    lo = min(pts, key=lambda p_: p_["delta"])
    top = max(pts, key=lambda p_: p_["delta"])
    labelled = [(hi, -10, 14, "right"), (lo, 12, -4, "left")]
    if top is not hi and top is not lo:
        labelled.append((top, -12, 6, "right"))
    for p_, dx, dy, ha in labelled:
        ax.annotate(p_["label"], (p_["force"], p_["delta"]), xytext=(dx, dy),
                    textcoords="offset points", ha=ha, fontsize=8.5,
                    color=colour[p_["fam"]], zorder=5)
    ax.annotate("ringed: no spring at all\n(slack strand, or arms that cannot close)",
                (0.02, 0.04), xycoords="axes fraction", fontsize=8.5,
                color=VS.INK_2, zorder=5)
    if r is not None:
        ax.annotate(f"r = {r:+.2f},  p = {pv:.2f}\nacross {len(pts)} constructs "
                    f"in {len(fams)} proteins",
                    (0.985, 0.045), xycoords="axes fraction", ha="right",
                    fontsize=9.5, color=VS.INK, zorder=5)

    ax.legend(loc="upper right", frameon=False, fontsize=9)
    ax.set_xlabel("force the spring model says the DNA delivers (pN)")
    ax.set_ylabel("displacement of the attachment sites (Å)\nagainst each construct's "
                  "own clicked control")
    ax.set_title("The sites move together, and not by more when the force is larger")
    ax.set_xlim(-0.6, max(p_["force"] for p_ in pts) * 1.12 + 0.4)
    ax.grid(axis="y", alpha=0.6)
    fig.tight_layout()
    fig.savefig(Path(outdir) / "fig7_headline.png", dpi=200, bbox_inches="tight")
    pd.DataFrame(pts).to_csv(Path(outdir) / "fig7_headline.csv", index=False)
    plt.close(fig)
    return "fig7_headline"


def fig_ramachandran(runs, outdir):
    """phi/psi at the two modified sites, WT against Tet2-Et against clicked.

    A scatter is the right form here because the claim is about *where in the
    Ramachandran plane* a residue sits, and the interesting fact -- that site 133 lives
    at positive phi, in the left-handed alpha region -- is a statement about position,
    not about a mean.  The generously drawn regions are the conventional ones; they are
    guides, not contours from a database.
    """
    sysdirs = {"S0_wt": "wild type", "S1_tet": "2x Tet2-Et (unclicked)",
               "S2_clicked": "clicked adduct"}
    data = {}
    for s in sysdirs:
        for f in sorted(Path(runs).glob(f"{s}/rep*/03_production/analysis_sites.csv")):
            d = pd.read_csv(f)
            for site in (133, 149):
                a, b = f"phi_{site}_deg", f"psi_{site}_deg"
                if a in d and b in d:
                    data.setdefault((s, site), []).append(d[[a, b]].to_numpy())
    if not data:
        print("  ramachandran: no site CSVs yet"); return None
    for k in data:
        data[k] = np.vstack(data[k])

    fig, axes = plt.subplots(1, 2, figsize=(9.6, 4.9), sharex=True, sharey=True)
    regions = [(-180, -30, -75, 30, r"$\beta$ / PPII"),
               (-160, -40, -70, -5, None),
               (-100, -30, -80, -5, r"right-handed $\alpha$"),
               (25, 100, -20, 90, r"left-handed $\alpha$")]
    for ax, site in zip(axes, (133, 149)):
        for x0, x1, y0, y1, lab in regions:
            ax.add_patch(plt.Rectangle((x0, y0), x1 - x0, y1 - y0, fc="none",
                                       ec=VS.GRID, lw=1.2, zorder=0))
            if lab:
                ax.text(x0 + 3, y1 - 6, lab, fontsize=7.5, color=VS.INK_MUTED,
                        va="top", zorder=1)
        for i, (s, label) in enumerate(sysdirs.items()):
            v = data.get((s, site))
            if v is None:
                continue
            ax.scatter(v[:, 0], v[:, 1], s=4, alpha=0.20, lw=0,
                       color=VS.SERIES[i], zorder=2)
            # circular mean, so a residue near the wrap does not average to nonsense
            mx = np.degrees(np.arctan2(np.sin(np.radians(v[:, 0])).mean(),
                                       np.cos(np.radians(v[:, 0])).mean()))
            my = np.degrees(np.arctan2(np.sin(np.radians(v[:, 1])).mean(),
                                       np.cos(np.radians(v[:, 1])).mean()))
            # Marker SHAPE as well as colour, because the three means sit within a few
            # degrees of one another at site 133 and overlapping discs of different
            # hues are not distinguishable.
            ax.plot(mx, my, marker=MEAN_MARKERS[i], ms=10, mfc=VS.SERIES[i],
                    mec="white", mew=1.8, ls="none", zorder=4)
            # Labels go in a fixed stack in an empty corner rather than beside each
            # mean: at site 133 all three land on top of each other otherwise, which is
            # what the first draft did.
            ax.plot(-165, -108 - 22 * i, marker=MEAN_MARKERS[i], ms=8,
                    mfc=VS.SERIES[i], mec="white", mew=1.4, ls="none",
                    clip_on=False, zorder=5)
            ax.text(-150, -108 - 22 * i, f"{label}   {mx:+.0f}, {my:+.0f}",
                    fontsize=8.5, color=VS.INK_2, va="center", zorder=5)
        ax.axhline(0, color=VS.INK_MUTED, lw=0.6, zorder=1)
        ax.axvline(0, color=VS.INK_MUTED, lw=0.6, zorder=1)
        ax.set_xlim(-180, 180); ax.set_ylim(-180, 180)
        ax.set_xticks([-180, -90, 0, 90, 180]); ax.set_yticks([-180, -90, 0, 90, 180])
        ax.set_xlabel(r"$\phi$ (degrees)")
        ax.set_title(f"site {site}", fontsize=10)
    axes[0].set_ylabel(r"$\psi$ (degrees)")
    fig.suptitle("Backbone conformation at the two modified sites", fontsize=11.5)
    fig.tight_layout()
    out = Path(outdir) / "ramachandran.png"
    fig.savefig(out, dpi=200, bbox_inches="tight")
    plt.close(fig)
    rows = [{"system": s, "site": site,
             "phi_deg": float(np.degrees(np.arctan2(
                 np.sin(np.radians(v[:, 0])).mean(),
                 np.cos(np.radians(v[:, 0])).mean()))),
             "psi_deg": float(np.degrees(np.arctan2(
                 np.sin(np.radians(v[:, 1])).mean(),
                 np.cos(np.radians(v[:, 1])).mean()))),
             "frames": len(v)}
            for (s, site), v in sorted(data.items())]
    pd.DataFrame(rows).to_csv(out.with_suffix(".csv"), index=False)
    return "ramachandran"


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
            fig_hbonds(df, out), fig_bend(a.runs, out),
             fig_ramachandran(a.runs, out),
             fig_kink_threshold(a.runs, out),
             fig_force_extension(a.runs, out),
             fig_headline(a.runs, out)]
    made = [m for m in made if m]
    print("wrote: " + ", ".join(f"{m}.png" for m in made))
    (out / "figures.json").write_text(json.dumps(made, indent=2))


if __name__ == "__main__":
    main()
