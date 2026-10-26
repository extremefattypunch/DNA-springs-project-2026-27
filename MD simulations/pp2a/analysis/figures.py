#!/usr/bin/env python3
"""Wiki figures for the PP2A spring, drawn from summary.csv and the per-run CSVs.

Every panel ships the CSV it was drawn from (vizstyle's relief rule).  Figures are
prefixed ``pp2a_`` because finalize.sh copies them into the shared docs/figures/.
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

from common import CONTROL, SCRATCH, SITES, SYSTEMS, analysis_prefix, heat_ranges  # noqa: E402
SCRATCH_RUNS = SCRATCH / "runs"
import vizstyle as VS                     # noqa: E402  (sfGFP)
from figures import labels_right          # noqa: E402  (sfGFP helper)
from spring_model import SsSpring, KT, XI_SS  # noqa: E402

VS.apply(matplotlib)
HERE = Path(__file__).resolve().parent


def color(s):
    return VS.SERIES[SYSTEMS[s]["slot"]]


def label(s):
    return SYSTEMS[s]["label"]


def per_frame(s, kind):
    f = Path(str(analysis_prefix(s)) + f"_{kind}.csv")
    return pd.read_csv(f) if f.exists() else None


def present(df):
    return [s for s in SYSTEMS if s in set(df.system)]


_DF = {"df": None}


def time_ns(d, s):
    """Frame times from the collector's bookkeeping: DCD files carry no clock, so mdtraj's
    time_ps is a frame index.  skip_ns + i x stride x report_ps."""
    df = _DF["df"]
    r = df[df.system == s].iloc[0]
    skip, stride, rep = float(r.get("skip_ns", 5.0)), float(r.get("stride", 2)), float(r.get("report_ps", 20.0))
    return skip + np.arange(len(d)) * stride * rep / 1000.0


def heat_boxes(ax, y, lo=0):
    for k, (a, b) in heat_ranges().items():
        if k % 2 == 0:
            ax.axvspan(a - 0.5, b + 0.5, color=VS.GRID, alpha=0.5, lw=0, zorder=0)
        ax.text((a + b) / 2, y, str(k), ha="center", va="top", fontsize=6.5, color=VS.INK_MUTED)


def fig_rmsd(df, out):
    fig, axes = plt.subplots(1, 3, figsize=(11, 3.4), sharex=True)
    frames = []
    for s in present(df):
        d = per_frame(s, "rmsd")
        if d is None:
            continue
        t = time_ns(d, s)
        for ax, c in zip(axes, ("rmsd_core_A", "rmsd_tail_corefit_A", "rmsd_all_A")):
            ax.plot(t, d[c], color=color(s), lw=1.1, alpha=0.9)
        frames.append(pd.DataFrame({"time_ns": t, **{f"{s}.{c}": d[c] for c in
                                                     ("rmsd_core_A", "rmsd_tail_corefit_A", "rmsd_all_A")}}))
    for ax, ttl in zip(axes, ("core 14–589, self-fit", "N-tail 1–13 after the core fit", "all Cα, self-fit")):
        ax.set_title(ttl, fontsize=9.5)
        ax.set_xlabel("time (ns)")
        ax.grid(axis="y", alpha=0.6)
    axes[0].set_ylabel(r"Cα RMSD to frame 0 ($\mathrm{\AA}$)")
    h = [plt.Line2D([], [], color=color(s), lw=2, label=label(s)) for s in present(df)]
    axes[2].legend(handles=h, loc="upper left", fontsize=7.5)
    fig.suptitle("Backbone rigidity: Day et al.'s split, N-tail against the rest", fontsize=11)
    fig.savefig(out / "pp2a_fig1_rmsd.png")
    if frames:
        pd.concat(frames, axis=1).to_csv(out / "pp2a_fig1_rmsd.csv", index=False)
    plt.close(fig)
    return "pp2a_fig1_rmsd"


def fig_rmsf(df, out):
    fig, ax = plt.subplots(figsize=(9, 3.8))
    frames, pending = [], []
    for s in present(df):
        d = per_frame(s, "rmsf")
        if d is None:
            continue
        ax.plot(d.resnum, d.rmsf_A, color=color(s), lw=1.3, label=label(s))
        pending.append((int(d.resnum.iloc[-1]), float(d.rmsf_A.iloc[-1]), label(s), color(s)))
        frames.append(d.set_index("resnum")["rmsf_A"].rename(s))
    labels_right(ax, pending)
    ax.set_ylim(0, min(ax.get_ylim()[1], 16))
    heat_boxes(ax, ax.get_ylim()[1] * 0.98)
    for x, txt in ((SITES[0], "Asp7 (D6X)"), (SITES[1], "Asn553 (N552X)")):
        ax.axvline(x, color=VS.INK_MUTED, lw=0.8, ls=":", zorder=1)
        ax.annotate(txt, (x, ax.get_ylim()[1] * 0.85), xytext=(3, 0), textcoords="offset points",
                    fontsize=7.5, color=VS.INK_MUTED)
    ax.set_xlim(0, 640)
    ax.set_xlabel("residue (UniProt P30153 numbering; HEAT repeats shaded, numbered at the top)")
    ax.set_ylabel(r"Cα RMSF ($\mathrm{\AA}$)")
    ax.set_title("Backbone mobility, aligned on the core")
    ax.grid(axis="y", alpha=0.6)
    fig.savefig(out / "pp2a_fig2_rmsf.png")
    if frames:
        pd.concat(frames, axis=1).to_csv(out / "pp2a_fig2_rmsf.csv")
    plt.close(fig)
    return "pp2a_fig2_rmsf"


def fig_end_to_end(df, out):
    fig, axes = plt.subplots(1, 3, figsize=(12, 3.6), gridspec_kw={"width_ratios": [2.2, 1, 1.4]})
    frames = []
    for s in present(df):
        m = per_frame(s, "mechanics")
        if m is None or "anchor_CB_CB_A" not in m:
            continue
        t = time_ns(m, s)
        axes[0].plot(t, m.anchor_CB_CB_A, color=color(s), lw=1.0, alpha=0.9)
        axes[1].hist(m.anchor_CB_CB_A, bins=40, orientation="horizontal", histtype="step",
                     color=color(s), lw=1.4, density=True)
        frames.append(pd.DataFrame({"time_ns": t, f"{s}.anchor_CB_CB_A": m.anchor_CB_CB_A}))
        if "spring_end_to_end_A" in m:
            axes[0].plot(t, m.spring_end_to_end_A, color=color(s), lw=0.9, ls="--", alpha=0.7)
            frames[-1][f"{s}.spring_end_to_end_A"] = m.spring_end_to_end_A
    axes[0].set_xlabel("time (ns)")
    axes[0].set_ylabel(r"Cβ7–Cβ553 ($\mathrm{\AA}$); dashed: DNA end-to-end")
    axes[0].set_title("Distance between the attachment points")
    axes[1].set_xlabel("density"); axes[1].set_yticklabels([])
    axes[1].sharey(axes[0])
    # the literature coordinate
    for s in present(df):
        h = per_frame(s, "heat")
        if h is None:
            continue
        c = [c for c in h.columns if c.startswith("d_CA30_CA578")]
        if c:
            axes[2].hist(h[c[0]], bins=40, histtype="step", color=color(s), lw=1.4, density=True, label=label(s))
    for x, txt in ((47.7, "6NTS compact"), (66.5, "apo MD, Yilmaz"), (70.3, "1B3U extended")):
        axes[2].axvline(x, color=VS.INK_MUTED, lw=0.8, ls=":")
        axes[2].annotate(txt, (x, axes[2].get_ylim()[1]), xytext=(2, -2), textcoords="offset points",
                         fontsize=6.5, color=VS.INK_MUTED, rotation=90, va="top")
    axes[2].set_xlabel(r"Cα Asn30–Phe578 ($\mathrm{\AA}$), the literature end-to-end")
    axes[2].set_ylabel("density")
    axes[2].legend(fontsize=7, loc="upper right")
    for ax in axes:
        ax.grid(axis="y", alpha=0.5)
    fig.savefig(out / "pp2a_fig3_end_to_end.png")
    if frames:
        pd.concat(frames, axis=1).to_csv(out / "pp2a_fig3_end_to_end.csv", index=False)
    plt.close(fig)
    return "pp2a_fig3_end_to_end"


def fig_force_extension(df, out):
    fxp = HERE / "force_extension.csv"
    fig, ax = plt.subplots(figsize=(7.2, 4.2))
    sp = SsSpring(50)
    x = np.linspace(20, 180, 400)
    s = np.clip(x / 10.0 / sp.L0, 0, 0.999)
    unclipped = (KT / XI_SS) * (1.0 / (4.0 * (1.0 - s) ** 2) - 0.25 + s)
    clipped = np.array([sp.tension(v / 10.0) for v in x])
    ideal = 3.0 * KT * (x / 10.0) / (2.0 * XI_SS * sp.L0)
    ax.plot(x, -unclipped, color=VS.INK_MUTED, lw=1.2, ls="--", label="Marko–Siggia (Zocchi eq. 2.207)")
    ax.plot(x, -clipped, color=VS.INK, lw=1.6, label="Marko–Siggia, sfGFP convention (zero below the rms size)")
    ax.plot(x, -ideal, color=VS.INK_MUTED, lw=1.0, ls=":", label="ideal chain 3kT x / (2ξL)")
    ax.axvline(sp.rms_end_to_end * 10, color=VS.GRID, lw=1.2)
    ax.text(sp.rms_end_to_end * 10 + 1, ax.get_ylim()[0] + 0.1, "rms size 94 Å", fontsize=7, color=VS.INK_MUTED)
    rows = None
    if fxp.exists():
        rows = pd.read_csv(fxp)
        ok = rows[rows.box_ok.astype(bool)]
        bad = rows[~rows.box_ok.astype(bool)]
        ax.errorbar(ok.mean_r_A, ok.force_pN, yerr=ok.force_sem_pN, fmt="o", color=VS.SERIES[2], ms=6,
                    capsize=3, label="measured: bare 50 nt strand, k = 33 kJ mol⁻¹ nm⁻²")
        if len(bad):
            ax.errorbar(bad.mean_r_A, bad.force_pN, yerr=bad.force_sem_pN, fmt="o", mfc="none",
                        color=VS.SERIES[2], ms=6, capsize=3, label="measured, box too small")
        for _, r in rows.iterrows():
            ax.annotate(f"x₀ {r.x0_A:.0f}", (r.mean_r_A, r.force_pN), xytext=(5, 5),
                        textcoords="offset points", fontsize=7, color=VS.INK_2)
    # where the chimera actually sits
    m = per_frame("P3_ss50", "mechanics")
    if m is not None and "spring_end_to_end_A" in m:
        v = m.spring_end_to_end_A
        ax.axvspan(v.quantile(0.05), v.quantile(0.95), color=VS.SERIES[2], alpha=0.12, lw=0)
        ax.text(v.mean(), ax.get_ylim()[1] * 0.9, "chimera's own\nend-to-end (5–95 %)", ha="center",
                fontsize=7.5, color=VS.INK_2)
    ax.axhline(0, color=VS.INK_MUTED, lw=0.6)
    ax.set_xlabel(r"end-to-end distance of the strand ($\mathrm{\AA}$)")
    ax.set_ylabel("force on the anchors (pN; + pushes apart, − pulls together)")
    ax.set_title("What the single strand delivers: model against force field")
    ax.legend(fontsize=7.5, loc="lower left")
    ax.grid(alpha=0.5)
    fig.savefig(out / "pp2a_fig4_force_extension.png")
    pd.DataFrame({"x_A": x, "ms_unclipped_pull_pN": unclipped, "ms_clipped_pull_pN": clipped,
                  "ideal_chain_pull_pN": ideal}).to_csv(out / "pp2a_fig4_force_extension_model.csv", index=False)
    if rows is not None:
        rows.to_csv(out / "pp2a_fig4_force_extension.csv", index=False)
    plt.close(fig)
    return "pp2a_fig4_force_extension"


def fig_heat_bend(df, out):
    fig, axes = plt.subplots(1, 3, figsize=(12, 3.6), gridspec_kw={"width_ratios": [1.4, 1.4, 1.6]})
    frames = []
    for s in present(df):
        p = Path(str(analysis_prefix(s)) + "_heat_profile.csv")
        h = per_frame(s, "heat")
        if not p.exists() or h is None:
            continue
        d = pd.read_csv(p)
        axes[0].errorbar(d.heat, d.bend_deg_mean, yerr=d.bend_deg_sd, color=color(s), lw=1.3,
                         marker=SYSTEMS[s]["marker"], ms=4, capsize=2, label=label(s))
        axes[1].errorbar(d.heat, d.twist_deg_mean, yerr=d.twist_deg_sd, color=color(s), lw=1.3,
                         marker=SYSTEMS[s]["marker"], ms=4, capsize=2)
        axes[2].plot(time_ns(h, s), h.aperture_H1_H15_A, color=color(s), lw=1.0, alpha=0.9)
        frames.append(d.set_index("heat").add_prefix(f"{s}."))
    axes[0].set_xlabel("HEAT repeat"); axes[0].set_ylabel("bend at the repeat (°)")
    axes[0].set_title("Where the solenoid bends (mean ± sd)")
    axes[0].legend(fontsize=7, loc="upper left")
    axes[1].set_xlabel("HEAT repeat"); axes[1].set_ylabel("superhelical twist at the step (°)")
    axes[1].set_title("…and how it twists")
    axes[2].set_xlabel("time (ns)"); axes[2].set_ylabel(r"HEAT1–HEAT15 centroid distance ($\mathrm{\AA}$)")
    axes[2].set_title("Horseshoe aperture (Day: 62.7 Å in the holoenzyme)")
    axes[2].axhline(62.7, color=VS.INK_MUTED, lw=0.8, ls=":")
    for ax in axes:
        ax.grid(alpha=0.5)
    for k in (5, 7, 10, 11, 12, 13):
        axes[0].axvline(k, color=VS.GRID, lw=0.8)
    fig.savefig(out / "pp2a_fig5_heat_bend.png")
    if frames:
        pd.concat(frames, axis=1).to_csv(out / "pp2a_fig5_heat_bend.csv")
    plt.close(fig)
    return "pp2a_fig5_heat_bend"


def fig_ssdna(df, out):
    """Does the DNA kink?  Where the strand bends, how sharply, and whether stretching it
    (the bare-strand windows) straightens it -- the single-strand analogue of the sfGFP
    base-pair-opening-against-curvature figure."""
    s = "P3_ss50"
    prof = Path(str(analysis_prefix(s)) + "_ssdna_profile.csv")
    d = per_frame(s, "ssdna")
    if d is None or not prof.exists():
        print("  ssdna: no P3_ss50 analysis yet"); return None
    pr = pd.read_csv(prof)
    fig, axes = plt.subplots(2, 2, figsize=(11, 6.6))
    c = color(s)
    ax = axes[0, 0]
    ax.plot(pr.nt, pr.bend_deg_mean, color=c, lw=1.4, label="C4′-trace bend, chimera")
    ax.fill_between(pr.nt, pr.bend_deg_mean - pr.bend_deg_sd, pr.bend_deg_mean + pr.bend_deg_sd, color=c, alpha=0.15, lw=0)
    ax.plot(pr.nt, pr.axis_curvature_deg_per_nt, color=VS.INK_2, lw=1.2, ls="--", label="smoothed-axis curvature (°/nt)")
    # the bare-strand windows, if analysed: same profile at each held separation
    fxp = HERE / "force_extension.csv"
    win_rows = pd.read_csv(fxp) if fxp.exists() else None
    if win_rows is not None:
        for i, (_, r) in enumerate(win_rows.sort_values("x0_A").iterrows()):
            wp = Path(SCRATCH_RUNS / r.tag / "rep1" / "03_production" / "analysis_ssdna_profile.csv")
            if wp.exists():
                w = pd.read_csv(wp)
                ax.plot(w.nt, w.bend_deg_mean, color=VS.SERIES[3 + i % 4], lw=0.9, alpha=0.8,
                        label=f"bare strand held at ⟨r⟩ = {r.mean_r_A:.0f} Å")
    ax2 = ax.twinx()
    ax2.bar(pr.nt, pr.sharp_turn_fraction, color=VS.SERIES[7], alpha=0.35, width=0.8)
    ax2.set_ylabel("fraction of frames with a > 90° turn", fontsize=8, color=VS.SERIES[7])
    ax2.set_ylim(0, max(0.05, float(pr.sharp_turn_fraction.max() or 0) * 1.3))
    ax.set_xlabel("nucleotide (5′ → 3′)"); ax.set_ylabel("bend (°)")
    ax.set_title("Where the strand bends, and where it turns sharply")
    ax.legend(fontsize=6.5, loc="upper left")
    ax = axes[0, 1]
    ax.plot(pr.nt, pr.bii_fraction, color=c, lw=1.3, label="BII fraction (ε − ζ > 0)")
    ax.plot(pr.nt, pr.unstacked_fraction_step_in, color=VS.INK_2, lw=1.1, ls="--", label="incoming step unstacked")
    ax.set_ylim(0, 1); ax.set_xlabel("nucleotide (5′ → 3′)"); ax.set_ylabel("fraction of frames")
    ax.set_title("Backbone torsion state and base stacking along the strand")
    ax.legend(fontsize=7.5)
    ax = axes[1, 0]
    t = time_ns(d, s)
    ax.plot(t, d.sharp_turns, color=c, lw=0.8, label="sharp turns per frame")
    ax.plot(t, d.axis_curvature_max_deg_per_nt / 10.0, color=VS.INK_2, lw=0.8, ls="--", label="max axis curvature / 10 (°/nt)")
    ax.set_xlabel("time (ns)"); ax.set_ylabel("count  |  °/nt ÷ 10"); ax.legend(fontsize=7.5)
    ax.set_title("Kinks over time")
    ax = axes[1, 1]
    ax.scatter(d.end_to_end_A, d.axis_curvature_mean_deg_per_nt, s=5, alpha=0.3, color=c, lw=0, label="chimera, per frame")
    if win_rows is not None and "kink.axis_curvature_deg_per_nt" in win_rows:
        ok = win_rows.dropna(subset=["kink.axis_curvature_deg_per_nt"])
        ax.scatter(ok.mean_r_A, ok["kink.axis_curvature_deg_per_nt"], s=60, marker="s", color=VS.SERIES[3],
                   edgecolor="white", zorder=4, label="bare strand, restraint windows")
    for y, lab in ((2.55, "Zocchi τc/B, nicked 2.6°/bp"), (3.21, "continuous 3.2°/bp"), (5.5, "sfGFP duplex opening onset 5.5°/bp")):
        ax.axhline(y, color=VS.INK_MUTED, lw=0.7, ls=":")
        ax.text(ax.get_xlim()[0] if False else d.end_to_end_A.min(), y + 0.3, lab, fontsize=6.5, color=VS.INK_MUTED)
    ax.set_xlabel(r"strand end-to-end ($\mathrm{\AA}$)"); ax.set_ylabel("mean axis curvature (°/nt)")
    ax.set_title("Curvature against extension — the duplex kink thresholds for scale")
    ax.legend(fontsize=7.5, loc="upper right")
    for ax in axes.flat:
        ax.grid(alpha=0.5)
    fig.suptitle("Does the DNA kink?", fontsize=11.5)
    fig.savefig(out / "pp2a_fig6_dna_kink.png")
    d.assign(time_ns=t).to_csv(out / "pp2a_fig6_dna_kink.csv", index=False)
    pr.to_csv(out / "pp2a_fig6_dna_kink_profile.csv", index=False)
    plt.close(fig)
    return "pp2a_fig6_dna_kink"


def fig_strand(df, out):
    s = "P3_ss50"
    d = per_frame(s, "ssdna")
    if d is None:
        return None
    fig, axes = plt.subplots(1, 2, figsize=(10, 3.4))
    t = time_ns(d, s)
    c = color(s)
    axes[0].plot(t, d.end_to_end_A, color=c, lw=1.0, label="end-to-end")
    axes[0].plot(t, d.rg_A, color=VS.INK_2, lw=1.0, label="radius of gyration")
    axes[0].axhline(94.1, color=VS.INK_MUTED, lw=0.8, ls=":")
    axes[0].text(t[0], 95, "rms size of a free 50-mer (model)", fontsize=7, color=VS.INK_MUTED)
    axes[0].set_ylabel(r"$\mathrm{\AA}$"); axes[0].legend(fontsize=7.5); axes[0].set_title("The strand's size")
    axes[1].plot(t, d.ms_tension_unclipped_pN, color=VS.INK_MUTED, lw=1.0, ls="--", label="Marko–Siggia (Zocchi eq. 2.207)")
    axes[1].plot(t, d.ms_tension_pN, color=c, lw=1.2, label="Marko–Siggia, sfGFP zero-below-rms convention")
    axes[1].plot(t, d.stacked_fraction, color=VS.SERIES[7], lw=0.8, alpha=0.7, label="stacked fraction (0–1)")
    axes[1].set_ylabel("modelled pull (pN)  |  fraction"); axes[1].legend(fontsize=7.5)
    axes[1].set_title("Modelled tension at that end-to-end distance")
    for ax in axes:
        ax.set_xlabel("time (ns)"); ax.grid(alpha=0.5)
    fig.savefig(out / "pp2a_fig8_strand.png")
    d.assign(time_ns=t).to_csv(out / "pp2a_fig8_strand.csv", index=False)
    plt.close(fig)
    return "pp2a_fig8_strand"


def fig_zocchi(df, out):
    """Zocchi's f(x) for the duplexes that would push at this span, with the single
    strand's Marko-Siggia pull and the chimera's own span distribution."""
    zp = HERE / "zocchi.json"
    cp = HERE / "zocchi_curves.csv"
    if not (zp.exists() and cp.exists()):
        return None
    z = json.loads(zp.read_text())
    cv = pd.read_csv(cp)
    fig, ax = plt.subplots(figsize=(7.4, 4.4))
    for i, n in enumerate((27, 40, 50)):
        col = VS.SERIES[i]
        ax.plot(cv.x_A, cv[f"f_{n}bp_nicked_pN"], color=col, lw=1.6, label=f"{n} bp duplex, nicked τc (full model: kinked below x_c, bent above)")
        ax.plot(cv.x_A, cv[f"f_{n}bp_continuous_pN"], color=col, lw=1.0, ls="--", alpha=0.8)
        ax.plot(cv.x_A, cv[f"f_{n}bp_nicked_unkinked_pN"], color=col, lw=0.8, ls=":", alpha=0.8)
    ax.plot(cv.x_A, -cv["ss50_pull_pN"], color=VS.INK, lw=1.8, label="50 nt single strand, Marko–Siggia (pulls)")
    sp = z.get("span", {})
    if sp.get("sd_A"):
        ax.axvspan(sp.get("p5_A", sp["A"] - sp["sd_A"]), sp.get("p95_A", sp["A"] + sp["sd_A"]),
                   color=VS.SERIES[2], alpha=0.12, lw=0)
    ax.axvline(sp.get("A", 90), color=VS.SERIES[2], lw=1.0)
    ax.text(sp.get("A", 90) + 1, ax.get_ylim()[1] * 0.92 if False else 20, f"span {sp.get('A', 90):.0f} Å\n({'simulated' if 'simulated' in sp.get('source', '') else 'build target'})",
            fontsize=7.5, color=VS.INK_2)
    for n, r in z["kink"]["softening_points"].items():
        ax.axvline(r["x_c_A"], color=VS.INK_MUTED, lw=0.6, ls=":")
        ax.text(r["x_c_A"] + 1, -4.5, f"x_c {n}", fontsize=6.5, color=VS.INK_MUTED, rotation=90, va="bottom")
    ax.axhline(0, color=VS.INK_MUTED, lw=0.6)
    ax.set_ylim(-6, 28); ax.set_xlim(20, 170)
    ax.set_xlabel(r"end-to-end distance the spring must span ($\mathrm{\AA}$)")
    ax.set_ylabel("force on the anchors (pN; + pushes apart, − pulls together)")
    ax.set_title("Zocchi's spring at this span — solid: f(x) with kinking (nicked τc), dashed: continuous τc, dotted: bent branch without a kink")
    ax.legend(fontsize=7, loc="upper right")
    ax.grid(alpha=0.5)
    fig.savefig(out / "pp2a_fig9_zocchi_model.png")
    cv.to_csv(out / "pp2a_fig9_zocchi_model.csv", index=False)
    plt.close(fig)
    return "pp2a_fig9_zocchi_model"


def fig_ramachandran(df, out):
    native = None
    try:
        att = json.loads((HERE.parent / "build/01_protein/attachment_sites.json").read_text())
        pair = next(iter(att["pairs"].values()))
        native = {SITES[0]: (pair["site_a"]["phi_deg"], pair["site_a"]["psi_deg"]),
                  SITES[1]: (pair["site_b"]["phi_deg"], pair["site_b"]["psi_deg"])}
    except Exception:                                             # noqa: BLE001
        pass
    data = {}
    for s in present(df):
        d = per_frame(s, "sites")
        if d is None:
            continue
        for site in SITES:
            a, b = f"phi_{site}_deg", f"psi_{site}_deg"
            if a in d and b in d:
                data[(s, site)] = d[[a, b]].to_numpy()
    if not data:
        print("  ramachandran: no site CSVs yet"); return None
    fig, axes = plt.subplots(1, 2, figsize=(9.6, 4.9), sharex=True, sharey=True)
    regions = [(-180, -30, -75, 30, r"$\beta$ / PPII"), (-100, -30, -80, -5, r"right-handed $\alpha$"),
               (25, 100, -20, 90, r"left-handed $\alpha$")]
    rows = []
    for ax, site in zip(axes, SITES):
        for x0, x1, y0, y1, lab in regions:
            ax.add_patch(plt.Rectangle((x0, y0), x1 - x0, y1 - y0, fc="none", ec=VS.GRID, lw=1.2, zorder=0))
            ax.text(x0 + 3, y1 - 6, lab, fontsize=7.5, color=VS.INK_MUTED, va="top", zorder=1)
        for i, s in enumerate(present(df)):
            v = data.get((s, site))
            if v is None:
                continue
            ax.scatter(v[:, 0], v[:, 1], s=4, alpha=0.18, lw=0, color=color(s), zorder=2)
            mx = np.degrees(np.arctan2(np.sin(np.radians(v[:, 0])).mean(), np.cos(np.radians(v[:, 0])).mean()))
            my = np.degrees(np.arctan2(np.sin(np.radians(v[:, 1])).mean(), np.cos(np.radians(v[:, 1])).mean()))
            ax.plot(mx, my, marker=SYSTEMS[s]["marker"], ms=10, mfc=color(s), mec="white", mew=1.8, ls="none", zorder=4)
            ax.plot(-165, -100 - 20 * i, marker=SYSTEMS[s]["marker"], ms=8, mfc=color(s), mec="white", mew=1.4,
                    ls="none", clip_on=False, zorder=5)
            ax.text(-150, -100 - 20 * i, f"{label(s)}   {mx:+.0f}, {my:+.0f}", fontsize=8, color=VS.INK_2, va="center")
            rows.append({"system": s, "site": site, "phi_deg": mx, "psi_deg": my, "frames": len(v)})
        if native and site in native and native[site][0] is not None:
            ax.plot(*native[site], marker="*", ms=13, mfc="white", mec=VS.INK, mew=1.2, ls="none", zorder=6)
            ax.text(native[site][0] + 6, native[site][1], "AlphaFold model", fontsize=7.5, color=VS.INK_2, va="center")
            rows.append({"system": "AF_model", "site": site, "phi_deg": native[site][0], "psi_deg": native[site][1], "frames": 1})
        ax.axhline(0, color=VS.INK_MUTED, lw=0.6, zorder=1); ax.axvline(0, color=VS.INK_MUTED, lw=0.6, zorder=1)
        ax.set_xlim(-180, 180); ax.set_ylim(-180, 180)
        ax.set_xticks([-180, -90, 0, 90, 180]); ax.set_yticks([-180, -90, 0, 90, 180])
        ax.set_xlabel(r"$\phi$ (degrees)")
        ax.set_title({SITES[0]: "Asp7 (construct D6X)", SITES[1]: "Asn553 (construct N552X)"}[site], fontsize=10)
    axes[0].set_ylabel(r"$\psi$ (degrees)")
    fig.suptitle("Backbone conformation at the two modified sites", fontsize=11.5)
    fig.savefig(out / "pp2a_ramachandran.png", dpi=200, bbox_inches="tight")
    pd.DataFrame(rows).to_csv(out / "pp2a_ramachandran.csv", index=False)
    plt.close(fig)
    return "pp2a_ramachandran"


def fig_headline(df, out):
    fp = HERE / "findings.json"
    if not fp.exists():
        return None
    f = json.loads(fp.read_text())
    ext = f.get("extension", {})
    if not ext:
        return None
    fig, ax = plt.subplots(figsize=(6.4, 3.6))
    rows = []
    for i, s in enumerate(present(df)):
        if s not in ext:
            continue
        e = ext[s]
        ax.errorbar(i, e["mean_A"], yerr=e["sem_block_A"] or 0, fmt=SYSTEMS[s]["marker"], color=color(s), ms=9,
                    capsize=4, mec="white", mew=1.2)
        ax.plot([i] * 2, [e["min_A"], e["max_A"]], color=color(s), lw=0.8, alpha=0.5)
        txt = f"{e['mean_A']:.1f} Å"
        w = f.get("extension_vs_control", {}).get(s)
        if w:
            txt += f"\nΔ {w['delta']:+.1f} Å, p = {w['p']:.2f}"
        ax.annotate(txt, (i, e["mean_A"]), xytext=(12, 0), textcoords="offset points", fontsize=7.5, va="center")
        rows.append({"system": s, **e, **({"delta_vs_control": w["delta"], "p": w["p"]} if w else {})})
    ax.set_xticks(range(len(present(df))))
    ax.set_xticklabels([label(s) for s in present(df)], fontsize=8)
    ax.set_ylabel(r"Cβ7–Cβ553 ($\mathrm{\AA}$), mean ± block SEM (range as a thin line)")
    if CONTROL in ext:
        ax.set_title(f"The attachment-point separation, against the clicked zero ({label(CONTROL)})")
    else:
        ax.set_title("The attachment-point separation (Cβ7–Cβ553); clicked control pending")
    ax.grid(axis="y", alpha=0.5)
    fig.savefig(out / "pp2a_fig7_headline.png")
    pd.DataFrame(rows).to_csv(out / "pp2a_fig7_headline.csv", index=False)
    plt.close(fig)
    return "pp2a_fig7_headline"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--summary", default=str(HERE / "summary.csv"))
    ap.add_argument("--out", default=str(HERE.parent / "figures"))
    a = ap.parse_args()
    df = pd.read_csv(a.summary)
    _DF["df"] = df
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    made = []
    for fn in (fig_rmsd, fig_rmsf, fig_end_to_end, fig_force_extension, fig_heat_bend, fig_ssdna,
               fig_strand, fig_zocchi, fig_ramachandran, fig_headline):
        try:
            r = fn(df, out)
            if r:
                made.append(r)
                print(f"  {r}")
        except Exception as e:                                    # noqa: BLE001
            print(f"  {fn.__name__} failed: {type(e).__name__}: {e}")
    (out / "figures.json").write_text(json.dumps(made, indent=2))


if __name__ == "__main__":
    main()
