"""A5 edge-on ring test (hypothesis 5 of results/gel_em_reconciliation.json). CPU only; real coordinates only.

Hypothesis 5: intact ~106 A DNA rings with sfGFP on the rim sit edge-on on the grid. Edge-on, a ring projects as a
<= 110 A bar; the far rim (where the arc runs along the beam) adds a second faint "dot". They would then have been
filed as "two dots + faint thread" = dot_pair_thread classes (rods 186, 101, 187, 188, 196; anchors 5, 19, 33, 59,
101). Chance neighbour pairs (expected NN ~75 A at ~4,500 dots/um^2) give the same look. Two tests:

(1) PAIR SEPARATION. Each class average is reduced to two dot positions: A = brightest local maximum, B = the
    brightest other local maximum >= 35 A away (Gaussian-smoothed average, robust z). They are mapped onto the
    micrograph per particle with the 2D pose/shift and the micrograph-orientation conventions determined
    empirically by scripts/17_rod_end_colocalisation.py (WORK/conventions.json; re-checked here per stream by the
    rate at which a detected dot lies within 20 A of the mapped A position, for all 16 candidate conventions).
    A_dot / B_dot = nearest detected dot (scripts/12 detector, dots from script 17's detect stage) within 20 A of
    the mapped A / B (B_dot != A_dot). Unique (A_dot, B_dot) pairs give the detected separation histogram
    (40-160 A, 5 A bins). Rotation null: B rotated about A by a uniform angle (20 replicates).
    Pair correlation g(r), 40-160 A, 5 A bins, for (i) ALL detected dots and (ii) the unique A dots of each class
    group: neighbour counts around each centre divided by the counts around the same centre shifted 150-300 A in a
    random direction on the same micrograph (shift null, 50 replicates; matches local density on the ~300 A scale).
    Only centres whose whole 160 A disc lies inside the detector's valid mask (bin-4 EDT) are used, observed and
    null alike (redrawn up to 50 times), so no edge correction is needed. Uncertainty of g: bootstrap over
    micrographs (1,000 resamples). Edge-on rings predict a narrow excess at 85-110 A in g_all beyond chance; the
    global excess number of pairs bounds how many class pairs can be real (excess_pairs_all = excess neighbours
    per eligible dot x total dots / 2, compared with the number of unique class pairs at 85-110 A).
    Comparators: rod-stream dot_pair (199, 105) and anchor-stream dot_pair (37 classes; separated dots, no thread)
    and single_dot / dot_long_rod (A only).
(2) THREAD DENSITY on the class averages. Integrated density per A along the A->B segment, middle part only
    (22 A from each dot centre excluded): strip integral |t| <= 15 A (also 25 A) minus the per-s median of the
    side strips 25 A <= |t| <= 40 A; and a polar version around A (sector +/-25 deg minus the same-radius median of
    the perpendicular sectors), which is robust to rotational blur about the aligned dot. Single-duplex reference,
    same stream, same procedure: rods = rod_no_dot averages (strip along the rod axis, ends inset 20 A) and
    dot_long_rod averages (polar, from 25 A beyond the dot to the half-max end minus 15 A, matched radii);
    anchors (no rod_no_dot class exists) = dot_short_rod averages (same polar/strip from the dot). Cross-check:
    density per A divided by the integrated dot (sfGFP) density of the same average. 34 A repeat: least-squares
    sinusoid amplitude vs period (20-70 A) on the detrended on-axis profile; reported relative to the on-axis level.
    Edge-on rings predict >= 2x duplex density per A and no 34 A repeat; a single duplex 1x; noise ~0.

  PY=/n/holylabs/hekstra_lab/Everyone/ianpoon/envs/dnaspring/bin/python
  sbatch -p test -c 8 --mem=32G -t 4:00:00 -o /n/netscratch/hekstra_lab/Lab/ian_poon/logs/A5_%j.log \
      --wrap "$PY -I scripts/43_A5_edge_on_ring_test.py run"
  $PY -I scripts/43_A5_edge_on_ring_test.py figure        (re-draw from WORK/A5_intermediate.npz + json)

Outputs: results/A5_edge_on.json, figures/qc/A5_edge_on.png; intermediates in WORK.
"""
import argparse
import importlib.util
import json
import os
import time
from multiprocessing import Pool

import mrcfile
import numpy as np
import pandas as pd
from scipy import ndimage as ndi
from scipy.spatial import cKDTree

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
_spec = importlib.util.spec_from_file_location("s17", os.path.join(REPO, "scripts", "17_rod_end_colocalisation.py"))
s17 = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(s17)

BASE = "/n/netscratch/hekstra_lab/Lab/ian_poon/sfgfp_chimera_screening_20261112"
EXPORTS = os.path.join(BASE, "exports")
DOTS_DIR = os.path.join(BASE, "rod_end_coloc", "dots")
CONV_JSON = os.path.join(BASE, "rod_end_coloc", "conventions.json")
CUR_NPZ = os.path.join(EXPORTS, "W2_curated__exposures.npz")
WORK = os.path.join(BASE, "A5_edge_on")
OUT_JSON = os.path.join(REPO, "results", "A5_edge_on.json")
FIG = os.path.join(REPO, "figures", "qc", "A5_edge_on.png")
STREAMS = {
    "rods": dict(npz=os.path.join(EXPORTS, "W2_2D_rods__particles.npz"),
                 mrc=os.path.join(EXPORTS, "W2_2D_rods", "J27_030_class_averages.mrc"),
                 labels=os.path.join(REPO, "results", "class2d_labels_W2_2D_rods.csv"),
                 dpt=[186, 101, 187, 188, 196]),
    "anchors": dict(npz=os.path.join(EXPORTS, "W2_2D_anchors__particles.npz"),
                    mrc=os.path.join(EXPORTS, "W2_2D_anchors", "J28_030_class_averages.mrc"),
                    labels=os.path.join(REPO, "results", "class2d_labels_W2_2D_anchors.csv"),
                    dpt=[5, 19, 33, 59, 101]),
}
A_MIC, A_P = 1.42, 2.84
MIC_N, BIN = 2048, 4
OFF = 4096.0
R_MATCH_A = 20.0
B_MIN_SEP_A = 35.0
EDGES = np.arange(40.0, 160.0 + 1e-6, 5.0)
CENT = 0.5 * (EDGES[1:] + EDGES[:-1])
RMAX_A = EDGES[-1]
WIN = (85.0, 110.0)
FLANKS = ((60.0, 85.0), (110.0, 135.0))
SHIFT_A = (150.0, 300.0)
NREP_SHIFT = 50
NREP_ROT = 20
NBOOT = 1000
SEED = 20261114 + 43
# thread / duplex density
DOT_EXCL_A = 22.0
STRIP_W_A = (15.0, 25.0)
BG_T_A = (25.0, 40.0)
SECTOR_DEG = 25.0
DOT_R_A = 18.0
DOT_BG_A = (25.0, 40.0)
PERIODS_A = np.arange(20.0, 70.0 + 1e-6, 1.0)


# ------------------------------------------------------------------------------------------- class-average geometry
def robust_z(img, apix):
    n = img.shape[0]
    yy, xx = np.indices(img.shape)
    rr = np.hypot(xx - n / 2, yy - n / 2) * apix
    win = rr < 0.45 * n * apix
    med = np.median(img[win])
    mad = 1.4826 * np.median(np.abs(img[win] - med)) + 1e-12
    return (img - med) / mad, rr


def two_dots(img, apix):
    """A = brightest local max (smoothed, robust z) within 60 A of the box centre; B = brightest other local max
    >= B_MIN_SEP_A from A inside the 2D window. Returns class-px coordinates (x = column, y = row) and z."""
    sm = ndi.gaussian_filter(img.astype(np.float64), 6.0 / apix)
    z, rr = robust_z(sm, apix)
    lm = (z == ndi.maximum_filter(z, size=5)) & (z > 3.0) & (rr < 115.0)
    ys, xs = np.nonzero(lm)
    zs = z[ys, xs]
    near = rr[ys, xs] < 60.0
    ia = int(np.argmax(np.where(near, zs, -np.inf)))
    ax_, ay_ = xs[ia], ys[ia]
    d = np.hypot(xs - ax_, ys - ay_) * apix
    cand = d >= B_MIN_SEP_A
    if not cand.any():
        return dict(ax=float(ax_), ay=float(ay_), az=float(zs[ia]), bx=None, by=None, bz=None, sep_A=None)
    ib = int(np.argmax(np.where(cand, zs, -np.inf)))
    # sub-pixel refinement (quadratic, 3x3)
    def refine(x, y):
        w = sm[y - 1:y + 2, x - 1:x + 2]
        if w.shape != (3, 3):
            return float(x), float(y)
        dx = 0.5 * (w[1, 2] - w[1, 0]) / (2 * w[1, 1] - w[1, 0] - w[1, 2] + 1e-12)
        dy = 0.5 * (w[2, 1] - w[0, 1]) / (2 * w[1, 1] - w[0, 1] - w[2, 1] + 1e-12)
        return float(x + np.clip(dx, -0.5, 0.5)), float(y + np.clip(dy, -0.5, 0.5))
    ax, ay = refine(ax_, ay_)
    bx, by = refine(xs[ib], ys[ib])
    return dict(ax=ax, ay=ay, az=float(zs[ia]), bx=bx, by=by, bz=float(zs[ib]),
                sep_A=float(np.hypot(bx - ax, by - ay) * apix))


def dot_and_axis(img, apix):
    """For dot + rod averages: A = brightest local max near centre; axis = principal axis of the z > 4 component
    containing A, oriented from A toward the component's far side; rod end = half-max point along the axis."""
    sm = ndi.gaussian_filter(img.astype(np.float64), 3.0 / apix)
    z, rr = robust_z(sm, apix)
    pk = np.unravel_index(np.argmax(np.where(rr < 40.0, z, -np.inf)), z.shape)
    lb, _ = ndi.label(z > 4.0, structure=np.ones((3, 3)))
    comp = lb == lb[pk]
    ys, xs = np.nonzero(comp)
    v = np.vstack([xs - xs.mean(), ys - ys.mean()])
    w, vec = np.linalg.eigh(np.cov(v))
    u = vec[:, 1]
    proj = (xs - pk[1]) * u[0] + (ys - pk[0]) * u[1]
    if np.sum(proj[proj > 0]) < -np.sum(proj[proj < 0]):
        u = -u
    s = np.arange(0.0, 140.0, apix / 4)
    prof = ndi.map_coordinates(z, [pk[0] + s / apix * u[1], pk[1] + s / apix * u[0]], order=1)
    rod = (s >= 25.0) & (s <= 55.0)
    lev = float(np.median(prof[rod]))
    i = int(np.searchsorted(s, 25.0))
    while i < len(s) - 1 and prof[i] >= 0.5 * lev:
        i += 1
    return dict(ax=float(pk[1]), ay=float(pk[0]), ux=float(u[0]), uy=float(u[1]), end_A=float(s[i]),
                rod_level_z=lev)


def rod_axis(img, apix):
    """rod_no_dot: component z > 4 containing the brightest pixel near centre; principal axis; half-max ends."""
    sm = ndi.gaussian_filter(img.astype(np.float64), 3.0 / apix)
    z, rr = robust_z(sm, apix)
    pk = np.unravel_index(np.argmax(np.where(rr < 40.0, z, -np.inf)), z.shape)
    lb, _ = ndi.label(z > 4.0, structure=np.ones((3, 3)))
    comp = lb == lb[pk]
    ys, xs = np.nonzero(comp)
    mx, my = xs.mean(), ys.mean()
    w, vec = np.linalg.eigh(np.cov(np.vstack([xs - mx, ys - my])))
    u = vec[:, 1]
    s = np.arange(-120.0, 120.0, apix / 4)
    prof = ndi.map_coordinates(z, [my + s / apix * u[1], mx + s / apix * u[0]], order=1)
    inside = prof > 4.0
    lev = float(np.median(prof[inside]))
    im = int(np.argmax(np.where(np.abs(s) < 30, prof, -np.inf)))
    i = im
    while i < len(s) - 1 and prof[i] >= 0.5 * lev:
        i += 1
    shi = s[i]
    i = im
    while i > 0 and prof[i] >= 0.5 * lev:
        i -= 1
    slo = s[i]
    c = 0.5 * (shi + slo)
    return dict(cx=float(mx + c / apix * u[0]), cy=float(my + c / apix * u[1]), ux=float(u[0]), uy=float(u[1]),
                h_A=float(0.5 * (shi - slo)), rod_level_z=lev, width_fwhm_A=float(np.sqrt(max(w[0], 0)) * 2.355 * apix))


def strip_profile(img, apix, x0, y0, ux, uy, s_A, w_A):
    """Integrated density across the line (|t| <= w_A) minus the per-s median of |t| in BG_T_A, per A of s.
    Returns D(s) (density x A) and the on-axis profile (mean |t| <= 1 px)."""
    dt = apix / 2
    tt = np.arange(-BG_T_A[1], BG_T_A[1] + 1e-6, dt)
    S, T = np.meshgrid(s_A, tt, indexing="ij")
    X = x0 + (S * ux - T * uy) / apix
    Y = y0 + (S * uy + T * ux) / apix
    v = ndi.map_coordinates(img.astype(np.float64), [Y, X], order=1, mode="nearest")
    bgm = (np.abs(tt) >= BG_T_A[0])
    bg = np.median(v[:, bgm], axis=1)
    inn = np.abs(tt) <= w_A
    D = ((v[:, inn] - bg[:, None]).sum(1)) * dt
    ax = v[:, np.abs(tt) <= apix].mean(1) - bg
    return D, ax


def polar_profile(img, apix, x0, y0, ux, uy, r_A):
    """Sector integral (+/- SECTOR_DEG about direction u) minus the same-radius median of the perpendicular sectors
    (60-120 deg on both sides), per A of radius (density x A)."""
    th0 = np.arctan2(uy, ux)
    dth = np.radians(1.0)
    th = np.arange(-180.0, 180.0, 1.0)
    R, TH = np.meshgrid(r_A, np.radians(th) + th0, indexing="ij")
    v = ndi.map_coordinates(img.astype(np.float64), [y0 + R * np.sin(TH) / apix, x0 + R * np.cos(TH) / apix],
                            order=1, mode="nearest")
    a = np.abs(th)
    bgm = (a >= 60) & (a <= 120)
    bg = np.median(v[:, bgm], axis=1)
    sec = a <= SECTOR_DEG
    return ((v[:, sec] - bg[:, None]).sum(1)) * dth * r_A


def dot_integral(img, apix, x0, y0):
    yy, xx = np.indices(img.shape)
    r = np.hypot(xx - x0, yy - y0) * apix
    bg = np.median(img[(r >= DOT_BG_A[0]) & (r <= DOT_BG_A[1])])
    return float((img[r <= DOT_R_A] - bg).sum() * apix * apix)


def periodogram(prof, s_A):
    """Least-squares sinusoid amplitude vs period on a quadratically detrended profile."""
    p = prof - np.polyval(np.polyfit(s_A, prof, 2), s_A)
    amps = []
    for P in PERIODS_A:
        M = np.c_[np.cos(2 * np.pi * s_A / P), np.sin(2 * np.pi * s_A / P)]
        c, *_ = np.linalg.lstsq(M, p, rcond=None)
        amps.append(float(np.hypot(*c)))
    return np.array(amps)


def density_tests(ca, apix, labels, dpt, counts):
    """Part (2) for one stream."""
    out = dict(apix_A=apix, thread={}, duplex_strip={}, duplex_polar={}, notes=[])
    # duplex references
    rod_k = [k for k in range(ca.shape[0]) if labels.get(k) == "rod_no_dot"]
    dotrod_k = [k for k in range(ca.shape[0]) if labels.get(k) in ("dot_long_rod", "dot_short_rod")]
    for k in rod_k:
        g = rod_axis(ca[k], apix)
        hs = g["h_A"] - 20.0
        if hs < 15:
            continue
        s = np.arange(-hs, hs + 1e-6, apix / 2)
        D15, axp = strip_profile(ca[k], apix, g["cx"], g["cy"], g["ux"], g["uy"], s, STRIP_W_A[0])
        D25, _ = strip_profile(ca[k], apix, g["cx"], g["cy"], g["ux"], g["uy"], s, STRIP_W_A[1])
        pg = periodogram(axp, s)
        out["duplex_strip"][int(k)] = dict(label="rod_no_dot", n=int(counts.get(k, 0)), L_A=2 * g["h_A"],
                                           D15=float(D15.mean()), D25=float(D25.mean()),
                                           axis_level=float(axp.mean()), width_fwhm_A=g["width_fwhm_A"],
                                           rel_amp_34=float(pg[np.abs(PERIODS_A - 34) <= 4].max() / abs(axp.mean())),
                                           best_period_A=float(PERIODS_A[np.argmax(pg)]),
                                           periodogram_rel=(pg / abs(axp.mean())).tolist())
    for k in dotrod_k:
        g = dot_and_axis(ca[k], apix)
        r_hi = g["end_A"] - 15.0
        if r_hi < 25.0 + 15.0:
            continue
        r = np.arange(25.0, r_hi + 1e-6, apix / 2)
        Dp = polar_profile(ca[k], apix, g["ax"], g["ay"], g["ux"], g["uy"], r)
        D15, axp = strip_profile(ca[k], apix, g["ax"], g["ay"], g["ux"], g["uy"], r, STRIP_W_A[0])
        Idot = dot_integral(ca[k], apix, g["ax"], g["ay"])
        out["duplex_polar"][int(k)] = dict(label=labels.get(k), n=int(counts.get(k, 0)), r_range_A=[25.0, r_hi],
                                           Dpolar=float(Dp.mean()), D15=float(D15.mean()), I_dot=Idot,
                                           Dpolar_per_dot=float(Dp.mean() / Idot) if Idot > 0 else None,
                                           Dpolar_r=Dp.tolist(), r_A=r.tolist())
    # dot_pair_thread classes
    for k in dpt:
        g = two_dots(ca[k], apix)
        if g["bx"] is None:
            out["thread"][int(k)] = dict(error="no second dot")
            continue
        ux, uy = (g["bx"] - g["ax"]) / (g["sep_A"] / apix), (g["by"] - g["ay"]) / (g["sep_A"] / apix)
        sfull = np.arange(-30.0, g["sep_A"] + 30.0 + 1e-6, apix / 2)
        Df15, axf = strip_profile(ca[k], apix, g["ax"], g["ay"], ux, uy, sfull, STRIP_W_A[0])
        mid = (sfull >= DOT_EXCL_A) & (sfull <= g["sep_A"] - DOT_EXCL_A)
        Df25, _ = strip_profile(ca[k], apix, g["ax"], g["ay"], ux, uy, sfull, STRIP_W_A[1])
        rmid = sfull[mid]
        Dp = polar_profile(ca[k], apix, g["ax"], g["ay"], ux, uy, rmid) if mid.sum() > 2 else np.array([np.nan])
        Idot = dot_integral(ca[k], apix, g["ax"], g["ay"])
        IdotB = dot_integral(ca[k], apix, g["bx"], g["by"])
        pg = periodogram(axf[mid], rmid) if mid.sum() >= 8 else np.full(len(PERIODS_A), np.nan)
        out["thread"][int(k)] = dict(n=int(counts.get(k, 0)), **g, mid_len_A=float(np.ptp(rmid)) if mid.any() else 0.0,
                                     D15=float(Df15[mid].mean()), D25=float(Df25[mid].mean()),
                                     Dpolar=float(np.nanmean(Dp)), I_dotA=Idot, I_dotB=IdotB,
                                     Dpolar_per_dot=float(np.nanmean(Dp) / Idot) if Idot > 0 else None,
                                     axis_level=float(axf[mid].mean()),
                                     rel_amp_34=float(np.nanmax(pg[np.abs(PERIODS_A - 34) <= 4]) / abs(axf[mid].mean())),
                                     best_period_A=float(PERIODS_A[np.nanargmax(pg)]) if np.isfinite(pg).any() else None,
                                     periodogram_rel=(pg / abs(axf[mid].mean())).tolist(),
                                     profile_s_A=sfull.tolist(), profile_D15=Df15.tolist())
    # references (particle-weighted medians) and ratios
    def wmed(vals, wts):
        vals, wts = np.asarray(vals, float), np.asarray(wts, float)
        o = np.argsort(vals)
        c = np.cumsum(wts[o])
        return float(vals[o][np.searchsorted(c, 0.5 * c[-1])])
    ref = {}
    if out["duplex_strip"]:
        d = out["duplex_strip"]
        for key in ("D15", "D25"):
            v = [x[key] for x in d.values()]
            ref[f"rod_no_dot_{key}"] = dict(wmedian=wmed(v, [x["n"] for x in d.values()]),
                                            p10=float(np.percentile(v, 10)), p90=float(np.percentile(v, 90)),
                                            n_classes=len(v))
    if out["duplex_polar"]:
        d = out["duplex_polar"]
        for key in ("Dpolar", "D15", "Dpolar_per_dot"):
            v = [x[key] for x in d.values() if x[key] is not None]
            ref[f"dot_rod_{key}"] = dict(wmedian=wmed(v, [x["n"] for x in d.values() if x[key] is not None]),
                                         p10=float(np.percentile(v, 10)), p90=float(np.percentile(v, 90)),
                                         n_classes=len(v))
    out["reference"] = ref
    for k, t in out["thread"].items():
        if "error" in t:
            continue
        r = {}
        if "rod_no_dot_D15" in ref:
            r["strip15_vs_rod_no_dot"] = t["D15"] / ref["rod_no_dot_D15"]["wmedian"]
            r["strip25_vs_rod_no_dot"] = t["D25"] / ref["rod_no_dot_D25"]["wmedian"]
        if "dot_rod_Dpolar" in ref:
            r["polar_vs_dot_rod"] = t["Dpolar"] / ref["dot_rod_Dpolar"]["wmedian"]
            r["strip15_vs_dot_rod"] = t["D15"] / ref["dot_rod_D15"]["wmedian"]
            if t["Dpolar_per_dot"] is not None:
                r["polar_per_dot_vs_dot_rod"] = t["Dpolar_per_dot"] / ref["dot_rod_Dpolar_per_dot"]["wmedian"]
        t["ratio_to_duplex"] = r
    if out["duplex_strip"]:
        a = np.array([x["periodogram_rel"] for x in out["duplex_strip"].values()])
        out["duplex_periodogram_rel_median"] = np.median(a, 0).tolist()
    return out


# ------------------------------------------------------------------------------------------- particles -> micrograph
def load_labels(path):
    lab = pd.read_csv(path).drop_duplicates("class")
    return dict(zip(lab["class"].astype(int), lab["consensus_label"]))


def map_points(P, idx, cx, cy, conv, ycfg):
    """class-frame points in particle px (2.84 A, from the box centre), one per particle -> micrograph px."""
    th = P["alignments2D__pose"][idx].astype(np.float64)
    sh = P["alignments2D__shift"][idx].astype(np.float64)
    px, py = s17.class_to_particle(cx, cy, th, sh[:, 0], sh[:, 1], conv)
    W = P["location__micrograph_shape"][idx, 1].astype(float)
    H = P["location__micrograph_shape"][idx, 0].astype(float)
    xf = P["location__center_x_frac"][idx].astype(np.float64)
    yf = P["location__center_y_frac"][idx].astype(np.float64)
    colc = xf * W if ycfg["x"] == "direct" else (1 - xf) * W
    rowc = yf * H if ycfg["y"] == "direct" else (1 - yf) * H
    fx = -1.0 if ycfg["fliplr"] else 1.0
    fy = -1.0 if ycfg["flipud"] else 1.0
    r = A_P / A_MIC
    return colc + r * fx * px, rowc + r * fy * py


def load_dots():
    xs, ys, zs, ms, paths, m4 = [], [], [], [], [], []
    off = 0
    for f in sorted(os.listdir(DOTS_DIR)):
        if not f.endswith(".npz"):
            continue
        p = np.load(os.path.join(DOTS_DIR, f))
        n = p["n"]
        xs.append(p["x"].astype(float))
        ys.append(p["y"].astype(float))
        zs.append(p["z"].astype(float))
        ms.append(np.repeat(np.arange(len(n)) + off, n))
        paths += list(p["mic_path"])
        m4.append(p["m4p"])
        off += len(n)
    return (np.array(paths), np.concatenate(xs), np.concatenate(ys), np.concatenate(zs), np.concatenate(ms),
            np.concatenate(m4))


# ------------------------------------------------------------------------------------------- per-micrograph g(r)
W_ = {}


def hist_d(c, d):
    if len(c) == 0 or len(d) == 0:
        return np.zeros(len(CENT))
    dd = np.hypot(c[:, None, 0] - d[None, :, 0], c[:, None, 1] - d[None, :, 1]).ravel() * A_MIC
    return np.histogram(dd, EDGES)[0].astype(float)


def mic_worker(m):
    """Observed and shift-null neighbour histograms for all groups on micrograph m."""
    lo, hi = W_["mstart"][m], W_["mstart"][m + 1]
    xy = np.c_[W_["dx"][lo:hi], W_["dy"][lo:hi]]
    strict = W_["dz"][lo:hi] >= 6.0
    valid = np.unpackbits(W_["m4p"][m], axis=1).astype(bool)
    edt = ndi.distance_transform_edt(np.pad(valid, 1))[1:-1, 1:-1] * (BIN * A_MIC) - 0.5 * BIN * A_MIC

    def clear(p):
        xi = np.clip((p[:, 0] // BIN).astype(int), 0, valid.shape[1] - 1)
        yi = np.clip((p[:, 1] // BIN).astype(int), 0, valid.shape[0] - 1)
        ok = (p[:, 0] >= 0) & (p[:, 0] < MIC_N) & (p[:, 1] >= 0) & (p[:, 1] < MIC_N)
        return ok & (edt[yi, xi] >= RMAX_A)
    elig = clear(xy)
    ng = len(W_["groups"])
    nb = len(CENT)
    obs = np.zeros((ng, nb))
    nobs = np.zeros(ng)
    nul = np.zeros((NREP_SHIFT, ng, nb))
    nnul = np.zeros((NREP_SHIFT, ng))
    gm = W_["gmask"][:, lo:hi] & elig[None, :]          # group membership of eligible centres
    rng = np.random.default_rng(SEED + 7919 * m)
    for gi, gname in enumerate(W_["groups"]):
        sel = gm[gi]
        if gname == "all_strict":
            sel = sel & strict
        if not sel.any():
            continue
        nbr = xy[strict] if gname == "all_strict" else xy
        c = xy[sel]
        # observed: exclude self (distance 0 < 40 A bin edge anyway)
        obs[gi] = hist_d(c, nbr)
        nobs[gi] = len(c)
        for rep in range(NREP_SHIFT):
            q = c.copy()
            todo = np.arange(len(c))
            for _ in range(50):
                if len(todo) == 0:
                    break
                rr = (SHIFT_A[0] + (SHIFT_A[1] - SHIFT_A[0]) * rng.random(len(todo))) / A_MIC
                ph = rng.random(len(todo)) * 2 * np.pi
                q[todo] = c[todo] + np.c_[rr * np.cos(ph), rr * np.sin(ph)]
                todo = todo[~clear(q[todo])]
            okq = np.ones(len(c), bool)
            okq[todo] = False
            nul[rep, gi] = hist_d(q[okq], nbr)
            nnul[rep, gi] = okq.sum()
    return m, obs, nobs, nul, nnul, int(elig.sum()), int(hi - lo)


# ------------------------------------------------------------------------------------------- main stage
def stage_run(a):
    t0 = time.time()
    os.makedirs(WORK, exist_ok=True)
    rng = np.random.default_rng(SEED)
    conv = json.load(open(CONV_JSON))
    ycfg = {k: conv["micrograph_best"][k] for k in ("x", "y", "flipud", "fliplr")}
    paths, dx, dy, dz, dm, m4p = load_dots()
    nmic = len(paths)
    assert np.all(np.diff(dm) >= 0)
    mstart = np.searchsorted(dm, np.arange(nmic + 1))
    lut = {p: i for i, p in enumerate(paths)}
    tree = cKDTree(np.c_[dx + dm * OFF, dy])
    print(f"dots: {len(dx)} on {nmic} micrographs ({time.time() - t0:.0f} s)", flush=True)

    res = dict(created=time.strftime("%Y-%m-%d %H:%M"), script="scripts/43_A5_edge_on_ring_test.py",
               hypothesis="H5: intact rings lie edge-on and were filed as dot_pair_thread classes",
               settings=dict(r_match_A=R_MATCH_A, b_min_sep_A=B_MIN_SEP_A, bins_A=EDGES.tolist(), window_A=WIN,
                             flanks_A=FLANKS, shift_null_A=SHIFT_A, n_rep_shift=NREP_SHIFT, n_rep_rot=NREP_ROT,
                             n_boot=NBOOT, centre_clearance_A=RMAX_A, dot_excl_A=DOT_EXCL_A, strip_w_A=STRIP_W_A,
                             bg_t_A=BG_T_A, sector_deg=SECTOR_DEG, seed=SEED),
               dots=dict(total=int(len(dx)), n_micrographs=int(nmic), strict_fraction=float(np.mean(dz >= 6))))
    groups_def = {}       # name -> (stream, classes, has_B)
    streams = {}
    for sname, sd in STREAMS.items():
        labels = load_labels(sd["labels"])
        ca = mrcfile.open(sd["mrc"], permissive=True).data.astype(np.float64)
        apix = float(mrcfile.open(sd["mrc"], permissive=True).voxel_size.x)
        Pz = np.load(sd["npz"])
        keys = ["alignments2D__class", "alignments2D__pose", "alignments2D__shift", "location__micrograph_path",
                "location__center_x_frac", "location__center_y_frac", "location__micrograph_shape"]
        P = {k: Pz[k] for k in keys}
        cls = P["alignments2D__class"].astype(int)
        counts = dict(zip(*np.unique(cls, return_counts=True)))
        streams[sname] = dict(labels=labels, ca=ca, apix=apix, P=P, cls=cls, counts=counts)
        dp = [k for k, v in labels.items() if v == "dot_pair"]
        groups_def[f"{sname}_dpt"] = (sname, sd["dpt"], True)
        groups_def[f"{sname}_dot_pair"] = (sname, dp, True)
        one = "dot_long_rod" if sname == "rods" else "single_dot"
        groups_def[f"{sname}_{one}"] = (sname, [k for k, v in labels.items() if v == one], False)
        for k in sd["dpt"]:
            groups_def[f"{sname}_c{k}"] = (sname, [k], True)

    # class two-dot geometry for every class used
    geo = {}
    for gname, (sname, ks, hasB) in groups_def.items():
        S = streams[sname]
        for k in ks:
            if (sname, k) not in geo:
                g = two_dots(S["ca"][k], S["apix"])
                geo[(sname, k)] = g
    res["class_two_dot_geometry"] = {f"{s}:{k}": dict(label=streams[s]["labels"].get(k), n=int(streams[s]["counts"].get(k, 0)), **g)
                                     for (s, k), g in geo.items() if k in STREAMS[s]["dpt"] or
                                     streams[s]["labels"].get(k) == "dot_pair"}

    # per-stream convention check: detected-dot hit rate at mapped A, all 16 candidates
    conv_check = {}
    best_conv = {}
    for sname in STREAMS:
        S = streams[sname]
        ks = groups_def[f"{sname}_dpt"][1] + groups_def[f"{sname}_dot_pair"][1]
        idx = np.nonzero(np.isin(S["cls"], ks))[0]
        idx = np.sort(rng.choice(idx, size=min(20000, len(idx)), replace=False))
        mi = np.array([lut.get(p, -1) for p in S["P"]["location__micrograph_path"][idx]])
        idx, mi = idx[mi >= 0], mi[mi >= 0]
        f = S["apix"] / A_P
        n0 = S["ca"].shape[1] / 2
        cx = np.array([(geo[(sname, k)]["ax"] - n0) * f for k in S["cls"][idx]])
        cy = np.array([(geo[(sname, k)]["ay"] - n0) * f for k in S["cls"][idx]])
        rates = []
        for c in s17.CANDS:
            X, Y = map_points(S["P"], idx, cx, cy, c, ycfg)
            d, _ = tree.query(np.c_[X + mi * OFF, Y], distance_upper_bound=R_MATCH_A / A_MIC)
            rates.append(float(np.mean(np.isfinite(d))))
        # control: random-angle A position (rotation about the box centre)
        ph = rng.random(len(idx)) * 2 * np.pi
        cxr, cyr = cx * np.cos(ph) - cy * np.sin(ph), cx * np.sin(ph) + cy * np.cos(ph)
        X, Y = map_points(S["P"], idx, cxr, cyr, s17.CANDS[int(np.argmax(rates))], ycfg)
        d, _ = tree.query(np.c_[X + mi * OFF, Y], distance_upper_bound=R_MATCH_A / A_MIC)
        b = int(np.argmax(rates))
        best_conv[sname] = s17.CANDS[b]
        conv_check[sname] = dict(n=int(len(idx)), rates={s17.cand_name(c): r for c, r in zip(s17.CANDS, rates)},
                                 best=s17.cand_name(s17.CANDS[b]), script17_best=conv["best"]["name"],
                                 rate_random_angle=float(np.mean(np.isfinite(d))),
                                 mean_A_radius_A=float(np.mean(np.hypot(cx, cy)) * A_P))
        print(f"{sname}: convention {conv_check[sname]['best']} A-hit {max(rates):.3f} (random-angle "
              f"{conv_check[sname]['rate_random_angle']:.3f}; script-17 best {conv['best']['name']})", flush=True)
    res["convention_check"] = conv_check

    # map A and B for each group, match to detected dots, rotation null for B
    ng_all = ["all", "all_strict"] + list(groups_def)
    gmask = np.zeros((len(ng_all), len(dx)), bool)
    gmask[0] = True
    gmask[1] = True
    pairs_out = {}
    sep_hists = {}
    for gname, (sname, ks, hasB) in groups_def.items():
        S = streams[sname]
        idx = np.nonzero(np.isin(S["cls"], ks))[0]
        mi = np.array([lut.get(p, -1) for p in S["P"]["location__micrograph_path"][idx]])
        idx, mi = idx[mi >= 0], mi[mi >= 0]
        f = S["apix"] / A_P
        n0 = S["ca"].shape[1] / 2
        G = [geo[(sname, k)] for k in S["cls"][idx]]
        cxa = np.array([(g["ax"] - n0) * f for g in G])
        cya = np.array([(g["ay"] - n0) * f for g in G])
        XA, YA = map_points(S["P"], idx, cxa, cya, best_conv[sname], ycfg)
        dA, iA = tree.query(np.c_[XA + mi * OFF, YA], distance_upper_bound=R_MATCH_A / A_MIC)
        okA = np.isfinite(dA)
        uA = np.unique(iA[okA])
        gmask[ng_all.index(gname), uA] = True
        rec = dict(n_particles=int(len(idx)), A_matched=float(okA.mean()), n_unique_A_dots=int(len(uA)),
                   A_offset_median_A=float(np.median(dA[okA]) * A_MIC) if okA.any() else None)
        if hasB:
            okg = np.array([g["bx"] is not None for g in G])
            cxb = np.array([(g["bx"] - n0) * f if g["bx"] is not None else 0.0 for g in G])
            cyb = np.array([(g["by"] - n0) * f if g["by"] is not None else 0.0 for g in G])
            XB, YB = map_points(S["P"], idx, cxb, cyb, best_conv[sname], ycfg)

            def match_B(xb, yb):
                kq = tree.query(np.c_[xb + mi * OFF, yb], k=2, distance_upper_bound=R_MATCH_A / A_MIC)
                d2, i2 = kq
                first = np.where((i2[:, 0] != iA) & np.isfinite(d2[:, 0]), i2[:, 0],
                                 np.where(np.isfinite(d2[:, 1]), i2[:, 1], -1))
                return np.where(okA & okg, first, -1)
            iB = match_B(XB, YB)
            okB = iB >= 0
            key = np.unique(np.c_[np.minimum(iA[okB], iB[okB]), np.maximum(iA[okB], iB[okB])], axis=0)
            sep = np.hypot(dx[key[:, 0]] - dx[key[:, 1]], dy[key[:, 0]] - dy[key[:, 1]]) * A_MIC
            sep_pred = np.array([g["sep_A"] if g["sep_A"] is not None else np.nan for g in G])
            h_obs = np.histogram(sep, EDGES)[0]
            # rotation null: B rotated about the mapped A position by a uniform angle
            rot_rates, rot_h = [], []
            vx, vy = XB - XA, YB - YA
            for rep in range(NREP_ROT):
                ph = rng.random(len(idx)) * 2 * np.pi
                xb = XA + vx * np.cos(ph) - vy * np.sin(ph)
                yb = YA + vx * np.sin(ph) + vy * np.cos(ph)
                ir = match_B(xb, yb)
                okr = ir >= 0
                kr = np.unique(np.c_[np.minimum(iA[okr], ir[okr]), np.maximum(iA[okr], ir[okr])], axis=0)
                rot_rates.append(float(np.sum(okr) / max(np.sum(okA & okg), 1)))
                rot_h.append(np.histogram(np.hypot(dx[kr[:, 0]] - dx[kr[:, 1]], dy[kr[:, 0]] - dy[kr[:, 1]]) * A_MIC,
                                          EDGES)[0])
            rot_h = np.array(rot_h, float)
            inwin = (sep >= WIN[0]) & (sep < WIN[1])
            rec.update(B_matched_given_A=float(okB.sum() / max(np.sum(okA & okg), 1)),
                       B_matched_rotation_null=dict(mean=float(np.mean(rot_rates)), sd=float(np.std(rot_rates))),
                       n_unique_pairs=int(len(key)), n_unique_pairs_in_window=int(inwin.sum()),
                       n_unique_pairs_in_window_rotation_null=float(rot_h[:, (CENT > WIN[0]) & (CENT < WIN[1])].sum(1).mean()),
                       sep_detected_median_A=float(np.median(sep)) if len(sep) else None,
                       sep_detected_p10_p90_A=[float(np.percentile(sep, 10)), float(np.percentile(sep, 90))] if len(sep) else None,
                       sep_class_pred_weighted_median_A=float(np.nanmedian(sep_pred)))
            sep_hists[gname] = dict(obs=h_obs.tolist(), rot_mean=rot_h.mean(0).tolist(),
                                    rot_lo=np.percentile(rot_h, 2.5, 0).tolist(),
                                    rot_hi=np.percentile(rot_h, 97.5, 0).tolist())
        pairs_out[gname] = rec
        print(f"  {gname}: {rec}", flush=True)
    res["class_pairs"] = pairs_out
    res["separation_histograms"] = sep_hists

    # pair correlation, per micrograph (Pool)
    W_.update(dx=dx, dy=dy, dz=dz, mstart=mstart, m4p=m4p, groups=ng_all, gmask=gmask)
    procs = int(os.environ.get("SLURM_CPUS_PER_TASK", 4))
    nb = len(CENT)
    OBS = np.zeros((nmic, len(ng_all), nb))
    NOBS = np.zeros((nmic, len(ng_all)))
    NULR = np.zeros((NREP_SHIFT, len(ng_all), nb))
    NNULR = np.zeros((NREP_SHIFT, len(ng_all)))
    NULM = np.zeros((nmic, len(ng_all), nb))
    NNULM = np.zeros((nmic, len(ng_all)))
    n_elig = np.zeros(nmic, int)
    with Pool(procs) as pool:
        for j, (m, o, no, nu, nn, ne, nt) in enumerate(pool.imap_unordered(mic_worker, range(nmic), chunksize=8)):
            OBS[m], NOBS[m] = o, no
            NULR += nu
            NNULR += nn
            NULM[m] = nu.sum(0)
            NNULM[m] = nn.sum(0)
            n_elig[m] = ne
            if (j + 1) % 1000 == 0:
                print(f"  g(r): {j + 1}/{nmic} micrographs ({time.time() - t0:.0f} s)", flush=True)
    print(f"g(r) done ({time.time() - t0:.0f} s); eligible centres {n_elig.sum()} of {len(dx)}", flush=True)

    # ice thickness per micrograph
    cur = np.load(CUR_NPZ)
    ice_l = dict(zip(cur["micrograph_blob__path"], cur["ctf_stats__ice_thickness_rel"]))
    ice = np.array([ice_l.get(p, np.nan) for p in paths], float)

    win_b = (CENT > WIN[0]) & (CENT < WIN[1])
    fl_b = np.zeros(nb, bool)
    for lo_, hi_ in FLANKS:
        fl_b |= (CENT > lo_) & (CENT < hi_)

    def gstats(mics_sel, gi):
        o = OBS[mics_sel, gi].sum(0)
        no = NOBS[mics_sel, gi].sum()
        nu = NULM[mics_sel, gi].sum(0)
        nn = NNULM[mics_sel, gi].sum()
        return o, no, nu, nn

    gres = {}
    boot_idx = rng.integers(0, nmic, size=(NBOOT, nmic))
    for gi, gname in enumerate(ng_all):
        o, no, nu, nn = gstats(slice(None), gi)
        if no == 0:
            continue
        g = (o / no) / (nu / nn)
        grep = (o / no)[None, :] / (NULR[:, gi] / NNULR[:, gi, None])       # obs vs each null replicate
        exc_win = float((o[win_b].sum() / no) - (nu[win_b].sum() / nn))     # excess neighbours per centre in window
        gw = float(o[win_b].sum() / no / (nu[win_b].sum() / nn))
        gf = float(o[fl_b].sum() / no / (nu[fl_b].sum() / nn))
        # bootstrap over micrographs
        bw, bf, be, bg = [], [], [], []
        Og, Ng, Ug, NNg = OBS[:, gi], NOBS[:, gi], NULM[:, gi], NNULM[:, gi]
        for b in range(NBOOT):
            ii = boot_idx[b]
            ob, nob, ub, nub = Og[ii].sum(0), Ng[ii].sum(), Ug[ii].sum(0), NNg[ii].sum()
            if nob == 0 or nub == 0:
                continue
            gb = (ob / nob) / (ub / nub)
            bg.append(gb)
            bw.append(ob[win_b].sum() / nob / (ub[win_b].sum() / nub))
            bf.append(ob[fl_b].sum() / nob / (ub[fl_b].sum() / nub))
            be.append(ob[win_b].sum() / nob - ub[win_b].sum() / nub)
        bg = np.array(bg)
        bump = np.array(bw) - np.array(bf)
        gres[gname] = dict(n_centres=int(no), g=g.tolist(), g_boot_lo=np.percentile(bg, 2.5, 0).tolist(),
                           g_boot_hi=np.percentile(bg, 97.5, 0).tolist(),
                           g_vs_null_reps_lo=np.percentile(grep, 2.5, 0).tolist(),
                           g_vs_null_reps_hi=np.percentile(grep, 97.5, 0).tolist(),
                           g_window=gw, g_window_ci95=[float(np.percentile(bw, 2.5)), float(np.percentile(bw, 97.5))],
                           g_flanks=gf, g_flanks_ci95=[float(np.percentile(bf, 2.5)), float(np.percentile(bf, 97.5))],
                           bump_window_minus_flanks=gw - gf,
                           bump_ci95=[float(np.percentile(bump, 2.5)), float(np.percentile(bump, 97.5))],
                           excess_nbrs_per_centre_window=exc_win,
                           excess_ci95=[float(np.percentile(be, 2.5)), float(np.percentile(be, 97.5))],
                           obs_nbrs_per_centre_window=float(o[win_b].sum() / no),
                           null_nbrs_per_centre_window=float(nu[win_b].sum() / nn))
    # ice-thickness terciles for all dots
    fin = np.isfinite(ice)
    qs = np.nanpercentile(ice[fin], [33.33, 66.67])
    ice_strata = {}
    for nm, sel in (("thin", fin & (ice < qs[0])), ("mid", fin & (ice >= qs[0]) & (ice < qs[1])),
                    ("thick", fin & (ice >= qs[1]))):
        o, no, nu, nn = gstats(sel, 0)
        ice_strata[nm] = dict(ice_rel_range=[float(np.nanmin(ice[sel])), float(np.nanmax(ice[sel]))],
                              n_mics=int(sel.sum()), g=((o / no) / (nu / nn)).tolist(),
                              g_window=float(o[win_b].sum() / no / (nu[win_b].sum() / nn)),
                              g_flanks=float(o[fl_b].sum() / no / (nu[fl_b].sum() / nn)))
    res["pair_correlation"] = gres
    res["pair_correlation_ice_terciles_all_dots"] = ice_strata

    # global excess pairs vs class pairs
    n_dots = len(dx)
    ex = gres["all"]["excess_nbrs_per_centre_window"]
    exci = gres["all"]["excess_ci95"]
    budget = dict(excess_pairs_all_dots_window=ex * n_dots / 2,
                  excess_pairs_all_dots_window_ci95=[exci[0] * n_dots / 2, exci[1] * n_dots / 2],
                  null_pairs_all_dots_window=gres["all"]["null_nbrs_per_centre_window"] * n_dots / 2,
                  note="unordered pairs at 85-110 A in the whole curated set beyond the shift null, extrapolated from "
                       "eligible centres to all dots")
    for gname in ("rods_dpt", "anchors_dpt", "rods_dot_pair", "anchors_dot_pair"):
        cp = pairs_out[gname]
        gg = gres.get(gname)
        budget[gname] = dict(
            unique_detected_pairs_window=cp["n_unique_pairs_in_window"],
            rotation_null_pairs_window=cp["n_unique_pairs_in_window_rotation_null"],
            class_excess_A_neighbours_window=None if gg is None else gg["excess_nbrs_per_centre_window"] * cp["n_unique_A_dots"],
            global_excess_over_class_pairs=(ex * n_dots / 2) / max(cp["n_unique_pairs_in_window"], 1))
    res["pair_budget"] = budget

    # part (2): thread density
    dens = {}
    for sname, sd in STREAMS.items():
        S = streams[sname]
        dens[sname] = density_tests(S["ca"], S["apix"], S["labels"], sd["dpt"], S["counts"])
    res["thread_density"] = dens
    res["runtime_s"] = round(time.time() - t0, 1)
    with open(OUT_JSON, "w") as fh:
        json.dump(res, fh, indent=1, default=float)
    np.savez(os.path.join(WORK, "A5_intermediate.npz"), OBS=OBS, NOBS=NOBS, NULM=NULM, NNULM=NNULM, NULR=NULR,
             NNULR=NNULR, groups=np.array(ng_all), ice=ice, n_elig=n_elig)
    print(f"wrote {OUT_JSON} ({time.time() - t0:.0f} s)", flush=True)
    make_figure()


# ------------------------------------------------------------------------------------------- figure
def make_figure():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    R = json.load(open(OUT_JSON))
    fig = plt.figure(figsize=(24, 21), dpi=85)
    sf = fig.subfigures(3, 1, height_ratios=[0.9, 1.05, 1.05])
    # row 1: class averages with A/B and the measured segment
    axs = sf[0].subplots(2, 9)
    shown = []
    for sname in STREAMS:
        ca = mrcfile.open(STREAMS[sname]["mrc"], permissive=True).data
        labels = load_labels(STREAMS[sname]["labels"])
        T = R["thread_density"][sname]["thread"]
        dsk = sorted(R["thread_density"][sname]["duplex_strip"].items(), key=lambda kv: -kv[1]["n"])[:2]
        dpk = sorted(R["thread_density"][sname]["duplex_polar"].items(), key=lambda kv: -kv[1]["n"])[:2]
        dp = [int(k.split(":")[1]) for k, v in R["class_two_dot_geometry"].items()
              if k.startswith(sname) and v["label"] == "dot_pair"]
        dp = sorted(dp, key=lambda k: -R["class_two_dot_geometry"][f"{sname}:{k}"]["n"])[:2 if sname == "anchors" else 1]
        items = [(int(k), "dpt") for k in STREAMS[sname]["dpt"]] + [(int(k), "dup") for k, _ in dsk + dpk][:2] + \
                [(k, "dp") for k in dp]
        row = [(sname, k, kind, ca[k], labels.get(k)) for (k, kind) in items[:9]]
        shown += row + [None] * (9 - len(row))
    for a_ in axs.ravel():
        a_.axis("off")
    for a_, it in zip(axs.ravel(), shown):
        if it is None:
            continue
        sname, k, kind, img, lab = it
        a_.imshow(img, cmap="gray")
        g = R["class_two_dot_geometry"].get(f"{sname}:{k}")
        if g is not None and g["bx"] is not None:
            a_.plot([g["ax"], g["bx"]], [g["ay"], g["by"]], "y-", lw=0.6)
            a_.plot(g["ax"], g["ay"], "r+", ms=7)
            a_.plot(g["bx"], g["by"], "c+", ms=7)
        t = R["thread_density"][sname]["thread"].get(str(k))
        extra = ""
        if t and "ratio_to_duplex" in t:
            rr = t["ratio_to_duplex"]
            extra = f"\nthread/duplex {rr.get('strip15_vs_rod_no_dot', rr.get('polar_vs_dot_rod', np.nan)):.2f}"
            if "polar_vs_dot_rod" in rr:
                extra += f" (polar {rr['polar_vs_dot_rod']:.2f})"
        sep = f" sep {g['sep_A']:.0f} A" if g is not None and g.get("sep_A") else ""
        a_.set_title(f"{sname[:-1]} {k} {lab}{sep}{extra}", fontsize=7)
        a_.set_xlim(8, 56)
        a_.set_ylim(56, 8)
    sf[0].suptitle("class averages: red + = dot A (alignment dot), cyan + = dot B, yellow = A-B segment "
                   "(thread measured on its middle, 22 A from each dot excluded); last columns = duplex references "
                   "and separated dot pairs", fontsize=10)
    # row 2: pair statistics
    ax = sf[1].subplots(1, 3)
    PC = R["pair_correlation"]
    for gname, col in (("all", "k"), ("all_strict", "0.5")):
        g = PC[gname]
        ax[0].plot(CENT, g["g"], color=col, lw=1.6, label=f"{gname} dots (n={g['n_centres']:,} centres)")
        ax[0].fill_between(CENT, g["g_boot_lo"], g["g_boot_hi"], color=col, alpha=0.25, lw=0)
    for nm, col in (("thin", "C0"), ("mid", "C1"), ("thick", "C3")):
        s_ = R["pair_correlation_ice_terciles_all_dots"][nm]
        ax[0].plot(CENT, s_["g"], color=col, lw=0.8, ls="--", label=f"all dots, ice {nm} ({s_['ice_rel_range'][0]:.2f}-"
                                                                     f"{s_['ice_rel_range'][1]:.2f})")
    ax[0].axvspan(*WIN, color="C2", alpha=0.12)
    ax[0].axhline(1, color="0.6", lw=0.6)
    a_ = PC["all"]
    ax[0].set_title(f"(a) pair correlation of ALL detected dots vs shift null (150-300 A)\n85-110 A: g = "
                    f"{a_['g_window']:.3f} [{a_['g_window_ci95'][0]:.3f}, {a_['g_window_ci95'][1]:.3f}]; window - flanks "
                    f"= {a_['bump_window_minus_flanks']:+.3f} [{a_['bump_ci95'][0]:+.3f}, {a_['bump_ci95'][1]:+.3f}]",
                    fontsize=9)
    ax[0].set_xlabel("dot-dot distance (A)")
    ax[0].set_ylabel("g(r) = obs / shift null")
    ax[0].legend(fontsize=7)
    for gname, col in (("all", "k"), ("rods_dpt", "C3"), ("anchors_dpt", "C1"), ("rods_dot_pair", "C4"),
                       ("anchors_dot_pair", "C0"), ("anchors_single_dot", "C2"), ("rods_dot_long_rod", "C5")):
        if gname not in PC:
            continue
        g = PC[gname]
        ax[1].plot(CENT, g["g"], color=col, lw=1.3, label=f"{gname} A dots (n={g['n_centres']:,})")
        ax[1].fill_between(CENT, g["g_boot_lo"], g["g_boot_hi"], color=col, alpha=0.18, lw=0)
    ax[1].axvspan(*WIN, color="C2", alpha=0.12)
    ax[1].axhline(1, color="0.6", lw=0.6)
    ax[1].set_title("(b) neighbour density around the class-mapped A dot vs shift null\n(class selection makes these "
                    "peaks; compare with (a) for whether the population carries an excess)", fontsize=9)
    ax[1].set_xlabel("distance from A dot (A)")
    ax[1].legend(fontsize=7)
    SH = R["separation_histograms"]
    for gname, col in (("rods_dpt", "C3"), ("anchors_dpt", "C1"), ("anchors_dot_pair", "C0")):
        h = SH[gname]
        tot = max(sum(h["obs"]), 1)
        ax[2].step(CENT, np.array(h["obs"]) / 5.0, where="mid", color=col, lw=1.5,
                   label=f"{gname}: detected A-B pairs (n={sum(h['obs']):,})")
        ax[2].step(CENT, np.array(h["rot_mean"]) / 5.0, where="mid", color=col, lw=1.0, ls=":",
                   label=f"{gname}: B rotated about A (null)")
    ax[2].axvspan(*WIN, color="C2", alpha=0.12)
    B = R["pair_budget"]
    ax[2].set_title(f"(c) unique detected A-B dot pairs per A of separation\nglobal excess pairs at 85-110 A: "
                    f"{B['excess_pairs_all_dots_window']:,.0f} [{B['excess_pairs_all_dots_window_ci95'][0]:,.0f}, "
                    f"{B['excess_pairs_all_dots_window_ci95'][1]:,.0f}] vs dot_pair_thread pairs there: rods "
                    f"{B['rods_dpt']['unique_detected_pairs_window']:,}, anchors "
                    f"{B['anchors_dpt']['unique_detected_pairs_window']:,}", fontsize=9)
    ax[2].set_xlabel("detected A-B separation (A)")
    ax[2].set_ylabel("pairs per A")
    ax[2].set_yscale("log")
    ax[2].legend(fontsize=7)
    # row 3: thread density
    ax = sf[2].subplots(1, 3)
    for sname, ls in (("rods", "-"), ("anchors", "--")):
        D = R["thread_density"][sname]
        ref = D["reference"]
        refv = ref.get("rod_no_dot_D15", ref.get("dot_rod_D15"))["wmedian"]
        for j, (k, t) in enumerate(D["thread"].items()):
            if "error" in t:
                continue
            s = np.array(t["profile_s_A"])
            ax[0].plot(s / t["sep_A"], np.array(t["profile_D15"]) / refv, ls=ls, lw=1.1, color=f"C{j}",
                       label=f"{sname} {k} (sep {t['sep_A']:.0f} A)")
    ax[0].axhline(1, color="k", lw=0.8)
    ax[0].axhline(2, color="C3", lw=0.8, ls=":")
    ax[0].axhline(0, color="0.6", lw=0.6)
    ax[0].axvspan(0, 1, color="0.9", zorder=0)
    ax[0].set_ylim(-0.5, 6)
    ax[0].set_xlabel("position along A->B (fraction of separation; 0 = A, 1 = B)")
    ax[0].set_ylabel("integrated density per A (|t|<=15 A) / single duplex")
    ax[0].set_title("(d) density along the A-B segment / duplex reference (rods: rod_no_dot; anchors: dot_short_rod)\n"
                    "black = 1x duplex, red dotted = 2x (edge-on ring prediction)", fontsize=9)
    ax[0].legend(fontsize=6, ncol=2)
    labs, vals = [], []
    for sname in STREAMS:
        D = R["thread_density"][sname]
        for k, t in D["thread"].items():
            if "error" in t:
                continue
            rr = t["ratio_to_duplex"]
            labs.append(f"{sname[0]}{k}")
            vals.append([rr.get("strip15_vs_rod_no_dot", np.nan), rr.get("strip25_vs_rod_no_dot", np.nan),
                         rr.get("strip15_vs_dot_rod", np.nan), rr.get("polar_vs_dot_rod", np.nan),
                         rr.get("polar_per_dot_vs_dot_rod", np.nan)])
    vals = np.array(vals)
    names = ["strip15 / rod_no_dot", "strip25 / rod_no_dot", "strip15 / dot+rod", "polar / dot+rod",
             "polar per dot / dot+rod per dot"]
    xx = np.arange(len(labs))
    for j in range(vals.shape[1]):
        ax[1].bar(xx + (j - 2) * 0.16, vals[:, j], 0.16, label=names[j])
    ax[1].axhline(1, color="k", lw=0.8)
    ax[1].axhline(2, color="C3", lw=0.8, ls=":")
    ax[1].set_xticks(xx)
    ax[1].set_xticklabels(labs, fontsize=8)
    ax[1].set_ylabel("thread density per A / single duplex")
    ax[1].set_title("(e) mid-segment thread density relative to a single duplex (same stream)", fontsize=9)
    ax[1].legend(fontsize=7)
    for sname, col in (("rods", "C3"), ("anchors", "C1")):
        D = R["thread_density"][sname]
        for k, t in D["thread"].items():
            if "error" in t or not np.isfinite(t["periodogram_rel"][0]):
                continue
            ax[2].plot(PERIODS_A, t["periodogram_rel"], color=col, lw=0.8, alpha=0.8)
        if "duplex_periodogram_rel_median" in D:
            ax[2].plot(PERIODS_A, D["duplex_periodogram_rel_median"], color="k", lw=2,
                       label=f"{sname} rod_no_dot (median of {len(D['duplex_strip'])})")
    ax[2].plot([], [], color="C3", label="rod-stream dot_pair_thread middles")
    ax[2].plot([], [], color="C1", label="anchor-stream dot_pair_thread middles")
    ax[2].axvspan(30, 38, color="C2", alpha=0.12)
    ax[2].set_xlabel("period (A)")
    ax[2].set_ylabel("sinusoid amplitude / on-axis level")
    ax[2].set_title("(f) on-axis periodicity (34 A B-DNA repeat band shaded); thread middles are only ~30-90 A long",
                    fontsize=9)
    ax[2].legend(fontsize=7)
    fig.savefig(FIG)
    plt.close(fig)
    print(f"wrote {FIG}", flush=True)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("stage", choices=["run", "figure", "density"])
    a = ap.parse_args()
    if a.stage == "run":
        stage_run(a)
    elif a.stage == "figure":
        make_figure()
    else:
        for sname, sd in STREAMS.items():
            labels = load_labels(sd["labels"])
            ca = mrcfile.open(sd["mrc"], permissive=True).data.astype(np.float64)
            apix = float(mrcfile.open(sd["mrc"], permissive=True).voxel_size.x)
            cls = np.load(sd["npz"])["alignments2D__class"].astype(int)
            counts = dict(zip(*np.unique(cls, return_counts=True)))
            d = density_tests(ca, apix, labels, sd["dpt"], counts)
            print(sname, json.dumps(d["reference"], indent=0))
            for k, t in d["thread"].items():
                print(k, {kk: (round(v, 3) if isinstance(v, float) else v) for kk, v in t.items()
                          if kk not in ("profile_s_A", "profile_D15", "periodogram_rel")})
