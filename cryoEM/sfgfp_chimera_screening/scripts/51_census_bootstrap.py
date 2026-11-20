#!/usr/bin/env python
"""Population census of picked objects with micrograph-bootstrap uncertainty, pooled and per EPU grid square.

Inputs (all read here; nothing typed in):
  exports/W2_2D_rods__particles.npz      J27 round-1 2D of the rod stream (all 8,036 curated micrographs)
  exports/W2_2D_anchors__particles.npz   J28 round-1 2D of the anchor stream (random 25 % subset: 2,009 micrographs)
  exports/W2_curated__exposures.npz      J17 curated exposures (micrograph uid -> EPU movie name)
  results/class2d_labels_W2_2D_{rods,anchors}.csv   blind-review consensus label per 2D class (classes >= 270 ptcls)
  raw/Images-Disc1/GridSquare_*/Data/*_Fractions.tiff   EPU movie -> grid square
Method: per micrograph, count particles per consensus label (rod-stream particles in the 96 unreviewed small classes
are excluded from the denominator, as in report §8, and counted separately).  Pooled fraction = sum over micrographs
of label counts / sum of labelled particles (particle-weighted).  95 % CI: percentile bootstrap over micrographs
(2,000 resamples with replacement, seed 20261115), the ratio recomputed per resample.  Per grid square: the same,
resampling micrographs within the square.  Heterogeneity across squares per label: Q = sum_s (p_s - p)^2 / se_s^2
(se_s = bootstrap SD), compared with chi-square on S-1 df; Q/df >> 1 means the squares differ by more than
micrograph-to-micrograph sampling within squares.
Crowding covariates per square.  'labelled_particles_per_micrograph' counts the stream's own labelled particles, so it
CONTAINS the category being correlated (e.g. anchor dot pairs are 62 % of anchors): its rho is partly circular.  Two
covariates independent of the 2D labels are therefore added: rod-stream particles per micrograph (all J27 particles,
labelled or not; independent of the anchor stream) and G0-detected dots per micrograph (rod_end_coloc/dots from
scripts/17_rod_end_colocalisation.py, divided by the valid-area fraction 1 - frac_bad; independent of both pickers'
2D labels), each the median over the square's micrographs.  Median ice and crowding co-vary across squares, so the
two cannot be separated with 20 squares; the rho between them is recorded.
Outputs: results/final/census.csv, results/final/census_by_square.csv, figures/final/census.png (+ census.json).
"""
import csv
import json
import os
import re
import sys
from collections import Counter

import numpy as np
from scipy import stats

REPO = "/n/home05/ianpoon/DNA-springs-project-2026-27/cryoEM/sfgfp_chimera_screening"
DATA = "/n/netscratch/hekstra_lab/Lab/ian_poon/sfgfp_chimera_screening_20261112"
EXP = f"{DATA}/exports"
RAW = f"{DATA}/raw/Images-Disc1"
OUT_CSV = f"{REPO}/results/final/census.csv"
OUT_SQ = f"{REPO}/results/final/census_by_square.csv"
FIG = f"{REPO}/figures/final/census.png"
NBOOT, SEED = 2000, 20261115
DOTS = f"{DATA}/rod_end_coloc/dots"

STREAMS = {
    "rods": dict(npz=f"{EXP}/W2_2D_rods__particles.npz", labels=f"{REPO}/results/class2d_labels_W2_2D_rods.csv",
                 job="J27", mics="all curated micrographs (J17)",
                 picker="J18 ellipse blob 20-180 A, top 50 % local power per defocus bin (J22), 200 A 2D window"),
    "anchors": dict(npz=f"{EXP}/W2_2D_anchors__particles.npz", labels=f"{REPO}/results/class2d_labels_W2_2D_anchors.csv",
                    job="J28", mics="random 25 % of curated micrographs (J25, seed 20261113)",
                    picker="J19 circle blob 30-50 A, top 50 % local power per defocus bin (J24), 220 A 2D window"),
}
NICE = {"dot_long_rod": "dot at end of ≥100 Å rod", "ice_contaminant": "ice / contamination",
        "rod_no_dot": "rod, no dot (free DNA)", "dot_pair_thread": "two dots + faint thread",
        "dot_pair": "two dots, unconnected", "noise_junk": "noise / junk", "unclear": "unclear",
        "dot_short_rod": "dot + short (<100 Å) rod", "single_dot": "single dot (free sfGFP)",
        "closed_loop": "closed loop", "dot_cluster": "dot cluster"}
CATEGORIES = ["closed_loop", "dot_long_rod", "dot_short_rod", "rod_no_dot", "dot_pair_thread", "dot_pair",
              "single_dot", "dot_cluster", "ice_contaminant", "noise_junk", "unclear"]


def grid_square_map():
    """EPU movie stem (FoilHole_..._Data_..._<date>_<time>) -> GridSquare_<id>, from the raw tree."""
    m = {}
    for gs in sorted(os.listdir(RAW)):
        d = os.path.join(RAW, gs, "Data")
        if not gs.startswith("GridSquare_") or not os.path.isdir(d):
            continue
        for f in os.listdir(d):
            if f.endswith("_Fractions.tiff"):
                m[f[:-len("_Fractions.tiff")]] = gs
    return m


def square_order():
    """Acquisition order of grid squares from results/raw_inventory.csv (first movie mtime)."""
    rows = list(csv.DictReader(open(f"{REPO}/results/raw_inventory.csv")))
    rows.sort(key=lambda r: r["first_movie_mtime"])
    return [r["grid_square"] for r in rows]


def boot_frac(C, rng, nboot=NBOOT):
    """C: (n_mic, n_lab) counts.  Returns pooled fraction, 2.5/97.5 percentiles and SD over micrograph resamples."""
    tot = C.sum()
    p = C.sum(0) / tot if tot else np.full(C.shape[1], np.nan)
    n = C.shape[0]
    B = np.empty((nboot, C.shape[1]))
    for b in range(nboot):
        w = np.bincount(rng.integers(0, n, n), minlength=n).astype(float)
        num = w @ C
        B[b] = num / num.sum()
    return p, np.percentile(B, 2.5, axis=0), np.percentile(B, 97.5, axis=0), B.std(0, ddof=1)


def main():
    rng = np.random.default_rng(SEED)
    ex = np.load(f"{EXP}/W2_curated__exposures.npz")
    gsmap = grid_square_map()
    stem = [re.sub(r"_Fractions\.tiff$", "", re.sub(r"^\d+_", "", os.path.basename(p))) for p in ex["movie_blob__path"]]
    mic_sq = {int(u): gsmap.get(s, "unmapped") for u, s in zip(ex["uid"], stem)}
    n_unmapped = sum(v == "unmapped" for v in mic_sq.values())
    order = [s for s in square_order() if s in set(mic_sq.values())]
    sq_rank = {s: i + 1 for i, s in enumerate(order)}
    ex_sq = np.array([mic_sq[int(u)] for u in ex["uid"]])
    ice = ex["ctf_stats__ice_thickness_rel"]
    dfm = 0.5 * (ex["ctf__df1_A"] + ex["ctf__df2_A"]) / 1e4
    sq_cov = {s: dict(median_ice_rel=round(float(np.median(ice[ex_sq == s])), 4),
                      median_defocus_um=round(float(np.median(dfm[ex_sq == s])), 3),
                      n_curated_micrographs=int((ex_sq == s).sum())) for s in order}
    # independent crowding covariates (see docstring)
    zr = np.load(STREAMS["rods"]["npz"])
    rmic = Counter(int(u) for u in zr["location__micrograph_uid"])
    rods_pm = np.array([rmic.get(int(u), 0) for u in ex["uid"]], float)
    dot_sq, dot_pm = [], []
    for f in sorted(os.listdir(DOTS)):
        if not f.endswith(".npz"):
            continue
        dz = np.load(os.path.join(DOTS, f))
        for pth, n, fb in zip(dz["mic_path"], dz["n"], dz["frac_bad"]):
            st_ = re.sub(r"_Fractions_patch_aligned_doseweighted\.mrc$", "", re.sub(r"^\d+_", "", os.path.basename(str(pth))))
            dot_sq.append(gsmap.get(st_, "unmapped"))
            dot_pm.append(float(n) / max(1e-6, 1.0 - float(fb)))
    dot_sq, dot_pm = np.array(dot_sq), np.array(dot_pm)
    for s in order:
        sq_cov[s]["median_rod_stream_particles_per_micrograph"] = float(np.median(rods_pm[ex_sq == s]))
        sq_cov[s]["median_detected_dots_per_micrograph"] = round(float(np.median(dot_pm[dot_sq == s])), 1)
        sq_cov[s]["n_micrographs_with_dot_detection"] = int((dot_sq == s).sum())
    indep_src = dict(rod_stream_particles_per_micrograph=f"{STREAMS['rods']['npz']} location__micrograph_uid "
                                                         "(all J27 particles, labelled or not), median over the square's curated micrographs",
                     detected_dots_per_micrograph=f"{DOTS}/task_*.npz n / (1 - frac_bad) (G0 dot detector, "
                                                  "scripts/17_rod_end_colocalisation.py), median over the square's micrographs",
                     n_micrographs_with_dot_detection=int(len(dot_pm)),
                     n_dot_micrographs_unmapped=int((dot_sq == "unmapped").sum()))

    rows, sq_rows, meta = [], [], dict(script="scripts/51_census_bootstrap.py", n_boot=NBOOT, seed=SEED,
                                       exposures=f"{EXP}/W2_curated__exposures.npz", n_curated=len(ex["uid"]),
                                       grid_square_source=f"{RAW}/GridSquare_*/Data/*_Fractions.tiff",
                                       n_micrographs_unmapped=n_unmapped, square_order_source="results/raw_inventory.csv",
                                       streams={})
    plot = {}
    for sname, S in STREAMS.items():
        z = np.load(S["npz"])
        cls = z["alignments2D__class"].astype(int)
        mic = z["location__micrograph_uid"].astype(np.uint64)
        lab_of = {int(r["class"]): r["consensus_label"] for r in csv.DictReader(open(S["labels"]))}
        labs = [l for l in CATEGORIES if l in set(lab_of.values())]
        li = {l: i for i, l in enumerate(labs)}
        pl = np.array([li.get(lab_of.get(c, ""), -1) for c in cls])
        n_unrev = int((pl < 0).sum())
        if sname == "rods":
            mics = np.array(sorted(int(u) for u in ex["uid"]), dtype=np.uint64)   # every curated micrograph
        else:
            mics = np.unique(mic)                                                  # the 2,009-micrograph subset
        mi = {int(u): i for i, u in enumerate(mics)}
        ok = pl >= 0
        C = np.zeros((len(mics), len(labs)))
        np.add.at(C, (np.array([mi[int(u)] for u in mic[ok]]), pl[ok]), 1)
        U = np.zeros(len(mics))
        np.add.at(U, np.array([mi[int(u)] for u in mic[~ok]], dtype=int), 1)
        p, lo, hi, sd = boot_frac(C, rng)
        per_mic = C.sum(1)
        for j, l in enumerate(labs):
            rows.append(dict(stream=sname, label=l, description=NICE[l], particles=int(C[:, j].sum()),
                             labelled_particles=int(C.sum()), fraction=round(p[j], 5), ci95_lo=round(lo[j], 5),
                             ci95_hi=round(hi[j], 5), boot_sd=round(sd[j], 5), n_micrographs=len(mics),
                             particles_per_micrograph_mean=round(float(C[:, j].mean()), 3)))
        sqs = np.array([mic_sq[int(u)] for u in mics])
        het = {}
        sq_p = {}
        for s in order:
            sel = sqs == s
            if sel.sum() == 0:
                continue
            ps, los, his, sds = boot_frac(C[sel], rng)
            sq_p[s] = (ps, los, his, sds, int(sel.sum()), int(C[sel].sum()))
            ppm = float(C[sel].sum() / sel.sum())
            for j, l in enumerate(labs):
                sq_rows.append(dict(stream=sname, grid_square=s, acquisition_rank=sq_rank[s], label=l,
                                    n_micrographs=int(sel.sum()), particles=int(C[sel, j].sum()),
                                    labelled_particles=int(C[sel].sum()), fraction=round(ps[j], 5),
                                    ci95_lo=round(los[j], 5), ci95_hi=round(his[j], 5), boot_sd=round(sds[j], 5),
                                    labelled_particles_per_micrograph=round(ppm, 2),
                                    square_median_ice_rel=sq_cov[s]["median_ice_rel"],
                                    square_median_defocus_um=sq_cov[s]["median_defocus_um"],
                                    square_median_rod_stream_particles_per_micrograph=sq_cov[s][
                                        "median_rod_stream_particles_per_micrograph"],
                                    square_median_detected_dots_per_micrograph=sq_cov[s][
                                        "median_detected_dots_per_micrograph"]))
        for j, l in enumerate(labs):
            ps = np.array([v[0][j] for v in sq_p.values()])
            ss = np.array([v[3][j] for v in sq_p.values()])
            good = ss > 0
            Q = float((((ps - p[j]) / ss) ** 2)[good].sum())
            df = int(good.sum()) - 1
            het[l] = dict(square_min=round(float(ps.min()), 4), square_max=round(float(ps.max()), 4),
                          square_min_id=list(sq_p)[int(ps.argmin())], square_max_id=list(sq_p)[int(ps.argmax())],
                          between_square_sd=round(float(ps.std(ddof=1)), 4),
                          median_within_square_boot_sd=round(float(np.median(ss)), 4),
                          Q=round(Q, 1), df=df, Q_over_df=round(Q / df, 2), p_chi2=float(stats.chi2.sf(Q, df)))
            sl = list(sq_p)
            ppm = np.array([sq_p[x][5] / sq_p[x][4] for x in sl])
            for cov, vals in (("labelled_particles_per_micrograph", ppm),
                              ("rod_stream_particles_per_micrograph",
                               np.array([sq_cov[x]["median_rod_stream_particles_per_micrograph"] for x in sl])),
                              ("detected_dots_per_micrograph",
                               np.array([sq_cov[x]["median_detected_dots_per_micrograph"] for x in sl])),
                              ("median_ice_rel", np.array([sq_cov[x]["median_ice_rel"] for x in sl])),
                              ("median_defocus_um", np.array([sq_cov[x]["median_defocus_um"] for x in sl])),
                              ("acquisition_rank", np.array([sq_rank[x] for x in sl]))):
                r, pv = stats.spearmanr(ps, vals)
                het[l][f"spearman_vs_{cov}"] = dict(rho=round(float(r), 3), p=float(pv))
        meta["streams"][sname] = dict(job=S["job"], source=S["npz"], labels_source=S["labels"], micrographs=S["mics"],
                                      picker=S["picker"], n_particles=int(len(cls)), n_labelled=int(C.sum()),
                                      n_unreviewed_small_classes=n_unrev, n_micrographs=len(mics),
                                      n_micrographs_with_labelled_particles=int((per_mic > 0).sum()),
                                      n_squares=len(sq_p), heterogeneity=het)
        plot[sname] = dict(labs=labs, p=p, lo=lo, hi=hi, sq=sq_p)
        print(sname, len(cls), int(C.sum()), n_unrev, len(mics), {l: round(p[j], 4) for j, l in enumerate(labs)}, flush=True)

    for path, rr in ((OUT_CSV, rows), (OUT_SQ, sq_rows)):
        with open(path, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=list(rr[0]))
            w.writeheader()
            w.writerows(rr)
    meta["squares_in_order"] = order
    meta["square_covariates"] = sq_cov
    meta["independent_crowding_covariates"] = indep_src
    meta["covariate_notes"] = {
        "labelled_particles_per_micrograph": ("the stream's own labelled particles per micrograph; INCLUDES the category "
                                              "being correlated (e.g. anchor 'two dots' is ~62 % of anchors), so its "
                                              "rho is partly circular; kept for reference only"),
        "rod_stream_particles_per_micrograph": "independent of the anchor stream's labels (not of the rod stream's)",
        "detected_dots_per_micrograph": "independent of both streams' 2D labels",
        "ice_vs_crowding": "median ice and crowding co-vary across squares, so their effects cannot be separated"}
    iv = np.array([sq_cov[x]["median_ice_rel"] for x in order])
    conf = {}
    for k in ("median_detected_dots_per_micrograph", "median_rod_stream_particles_per_micrograph"):
        r, pv = stats.spearmanr(iv, [sq_cov[x][k] for x in order])
        conf[k] = dict(rho=round(float(r), 3), p=float(pv))
    for sname in plot:
        sl = list(plot[sname]["sq"])
        ppm = [plot[sname]["sq"][x][5] / plot[sname]["sq"][x][4] for x in sl]
        r, pv = stats.spearmanr([sq_cov[x]["median_ice_rel"] for x in sl], ppm)
        conf[f"{sname}_labelled_particles_per_micrograph"] = dict(rho=round(float(r), 3), p=float(pv))
    meta["ice_vs_crowding_spearman_over_squares"] = conf
    meta["square_median_ice_range"] = [float(iv.min()), float(iv.max())]
    # do the two streams' per-square compositions co-vary? (same squares, independent pickers)
    cross = {}
    for la in plot["anchors"]["labs"]:
        ja = plot["anchors"]["labs"].index(la)
        for lr in plot["rods"]["labs"]:
            jr = plot["rods"]["labs"].index(lr)
            sl = [x for x in order if x in plot["anchors"]["sq"] and x in plot["rods"]["sq"]]
            r, pv = stats.spearmanr([plot["anchors"]["sq"][x][0][ja] for x in sl], [plot["rods"]["sq"][x][0][jr] for x in sl])
            cross[f"anchors:{la} vs rods:{lr}"] = dict(rho=round(float(r), 3), p=float(pv))
    meta["cross_stream_square_spearman"] = cross
    figure(plot, meta)
    meta.update(figure=os.path.relpath(FIG, REPO), tables=[os.path.relpath(OUT_CSV, REPO), os.path.relpath(OUT_SQ, REPO)],
                plotted_numbers="bars, whiskers and per-square dots are the fraction / ci95_lo / ci95_hi columns of "
                                "results/final/census.csv and results/final/census_by_square.csv")
    json.dump(meta, open(FIG.replace(".png", ".json"), "w"), indent=1)


def figure(plot, meta):
    import textwrap
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    sys.path.insert(0, "/n/home05/ianpoon/DNA-springs-project-2026-27/cryoEM/pp2a_day2025_reproduction/scripts")
    import vizstyle
    vizstyle.apply(matplotlib)
    S = vizstyle.SERIES
    col = {"rods": S[0], "anchors": S[1]}
    # per-square panels: colour follows the label, not its rank
    lab_col = {"dot_long_rod": S[0], "ice_contaminant": S[2], "rod_no_dot": S[6],
               "dot_pair": S[1], "single_dot": S[5], "dot_short_rod": S[4]}
    title = {"rods": "Rod stream (J27)\n{n:,} labelled particles, all {m:,} curated micrographs",
             "anchors": "Anchor stream (J28)\n{n:,} labelled particles, {m:,}-micrograph random subset"}
    fig = plt.figure(figsize=(12.5, 12.0))
    sf_top, sf_bot, sf_cap = fig.subfigures(3, 1, height_ratios=[1.0, 0.92, 0.27])
    ax_top = sf_top.subplots(1, 2, gridspec_kw=dict(wspace=0.08))
    ax_bot = sf_bot.subplots(1, 2, gridspec_kw=dict(wspace=0.08))
    for k, sname in enumerate(("rods", "anchors")):
        P = plot[sname]
        st = meta["streams"][sname]
        het = st["heterogeneity"]
        order = np.argsort(P["p"])            # largest at top
        ax = ax_top[k]
        y = np.arange(len(order))
        ax.barh(y, 100 * P["p"][order], color=col[sname], height=0.62, zorder=2)
        ax.errorbar(100 * P["p"][order], y, xerr=[100 * (P["p"] - P["lo"])[order], 100 * (P["hi"] - P["p"])[order]],
                    fmt="none", ecolor=vizstyle.INK, elinewidth=1.2, capsize=3, zorder=4)
        xmax = 0
        for jj, j in enumerate(order):
            vals = np.array([v[0][j] for v in P["sq"].values()])
            ax.scatter(100 * vals, np.full(len(vals), jj) + np.linspace(-0.22, 0.22, len(vals)), s=9,
                       color=vizstyle.INK_MUTED, alpha=0.85, zorder=3, linewidths=0)
            xr = 100 * max(P["hi"][j], vals.max()) + 1.5
            ax.text(xr, jj, f"{100 * P['p'][j]:.1f} % [{100 * P['lo'][j]:.1f}–{100 * P['hi'][j]:.1f}]",
                    va="center", fontsize=8, color=vizstyle.INK)
            xmax = max(xmax, xr)
        ax.set_yticks(y, [NICE[P["labs"][j]] for j in order])
        ax.set_xlabel("% of labelled picked objects")
        ax.set_xlim(0, xmax + 22)
        ax.grid(axis="x", zorder=0)
        ax.set_title(title[sname].format(n=st["n_labelled"], m=st["n_micrographs"]), fontsize=10)
        ax2 = ax_bot[k]
        top3 = list(order[::-1][:3])
        sqs = list(P["sq"])
        x = np.arange(1, len(sqs) + 1)
        for t, j in enumerate(top3):
            l = P["labs"][j]
            ps = np.array([P["sq"][s_][0][j] for s_ in sqs])
            lo = np.array([P["sq"][s_][1][j] for s_ in sqs])
            hi = np.array([P["sq"][s_][2][j] for s_ in sqs])
            off = (t - 1) * 0.24
            ax2.errorbar(x + off, 100 * ps, yerr=[100 * (ps - lo), 100 * (hi - ps)], fmt="o", ms=3.5, color=lab_col[l],
                         elinewidth=1, capsize=0, label=f"{NICE[l]} (Q/df {het[l]['Q_over_df']:.0f})")
            ax2.axhline(100 * P["p"][j], color=lab_col[l], lw=0.9, ls="--", zorder=0)
        ax2.set_xticks(x, [str(i) for i in x], fontsize=8)
        ax2.set_xlim(0.4, len(sqs) + 0.6)
        ax2.set_xlabel("grid square, in acquisition order (IDs in census_by_square.csv)")
        ax2.set_ylabel("% of labelled picked objects")
        lo_, hi_ = ax2.get_ylim()
        ax2.set_ylim(lo_, hi_ + 0.28 * (hi_ - lo_))
        ax2.grid(axis="y")
        ax2.set_title("Per grid square: three largest categories, 95 % CI (dashed: pooled)", fontsize=9.5)
        ax2.legend(loc="upper left", ncol=1, fontsize=8, handletextpad=0.3, borderaxespad=0.2)
    fig.suptitle("Census of picked objects by blind-consensus 2D label (W2 round-1 2D)", x=0.01, ha="left", fontsize=12)
    qd = [h["Q_over_df"] for st_ in meta["streams"].values() for h in st_["heterogeneity"].values()]
    dfs = sorted({h["df"] for st_ in meta["streams"].values() for h in st_["heterogeneity"].values()})
    nsq = meta["streams"]["rods"]["n_squares"]
    ha, hr = meta["streams"]["anchors"]["heterogeneity"], meta["streams"]["rods"]["heterogeneity"]
    ir, cf = meta["square_median_ice_range"], meta["ice_vs_crowding_spearman_over_squares"]
    cap = ("These are fractions of PICKED objects, not molar ratios: each stream's picker and 2D window biases what it "
           "sees (rods: ellipse blob 20–180 Å + top-50 % local power; anchors: 30–50 Å circle blob + top-50 % local "
           "power, random 25 % of micrographs). Bars: pooled particle-weighted fraction, with value and 95 % percentile "
           f"bootstrap CI over micrographs (2,000 resamples); grey dots: the {nsq} grid squares. Bottom: per-square "
           "fractions with within-square bootstrap CIs. Q/df = between-square scatter relative to within-square "
           f"sampling error (chi-square, {'/'.join(map(str, dfs))} df); every category varies between squares far beyond "
           f"sampling (Q/df {min(qd):.0f}–{max(qd):.0f}), so the pooled CI understates square-to-square variation. "
           "The variation tracks ice and crowding across squares (Spearman over squares), and the two are confounded: "
           f"square median relative ice spans only {ir[0]:.2f}–{ir[1]:.2f} and falls as crowding rises (ρ ice vs detected "
           f"dots per micrograph {cf['median_detected_dots_per_micrograph']['rho']:+.2f}), so their effects cannot be "
           f"separated. Anchor 'two dots' rises with detected dots per micrograph (ρ "
           f"{ha['dot_pair']['spearman_vs_detected_dots_per_micrograph']['rho']:+.2f}; independent of the 2D labels) and "
           f"with rod-stream particles per micrograph (ρ "
           f"{ha['dot_pair']['spearman_vs_rod_stream_particles_per_micrograph']['rho']:+.2f}), and falls with median ice "
           f"(ρ {ha['dot_pair']['spearman_vs_median_ice_rel']['rho']:+.2f}); free rods rise with ice (ρ "
           f"{hr['rod_no_dot']['spearman_vs_median_ice_rel']['rho']:+.2f}; vs detected dots ρ "
           f"{hr['rod_no_dot']['spearman_vs_detected_dots_per_micrograph']['rho']:+.2f}). (Anchors per micrograph, ρ "
           f"{ha['dot_pair']['spearman_vs_labelled_particles_per_micrograph']['rho']:+.2f}, is not used: it includes "
           "the dot-pair category itself.) "
           "Rod-stream particles in unreviewed "
           f"(<270-particle) classes ({meta['streams']['rods']['n_unreviewed_small_classes']:,}) are excluded.")
    meta["caption"] = cap
    sf_cap.text(0.01, 0.97, "\n".join(textwrap.wrap(cap, 205)), fontsize=8.2, color=vizstyle.INK_2, ha="left", va="top")
    fig.savefig(FIG, dpi=200)
    plt.close(fig)


if __name__ == "__main__":
    sys.exit(main())
