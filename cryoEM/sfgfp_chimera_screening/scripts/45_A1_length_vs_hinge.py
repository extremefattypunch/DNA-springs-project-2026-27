"""A1 follow-up: is the DNA attached to sfGFP genuinely SHORTER than free DNA (H1: yielded / incomplete loops), or
full length but BLURRED by a hinge at the attachment (H2/H3: linear sfGFP-dsDNA mono-adducts)?   CPU only.

Background. A1 (P6 J45: 223,536 re-centred rod particles, 454 A box, 400 A window, no recentring, 150 classes, 12 A)
gives dot-at-rod-end classes whose half-max axial length is ~117-149 A including the dot, i.e. DNA ~100-115 A, while
dot-free straight rods read ~123-170 A. 2D registers on the bright sfGFP dot; if the DNA pivots about the attachment,
the far end of the DNA is spread over an arc and a narrow-strip half-max length under-reads. Pure rotation about the
dot does NOT move the far end radially, so a hinge-insensitive length exists: the end of the RADIAL MASS profile
(density integrated across the rod over a wide angular sector about the dot).

Stages (run `all` as one Slurm CPU job on `test`; `classes` alone is quick and runs on class averages only):

  classes   (a) per-class polar profiles on the CryoSPARC class averages; (c) helical repeat; forward-blur calibration
  particles (b) pose-convention check, reading + class-frame mapping of particles, dot-free re-alignment
  report   (d) bootstrap length distributions, verdict numbers, JSON + figure

DEFINITIONS (fixed before looking at any of the A1 results beyond the first-pass table)
1. Class groups (classes with >= 500 particles; 105 and 12 = lone dot + haze are excluded by name):
   DOT  ("dot at one end of a rod"): round-1 dot-label fraction >= 0.5; brightest 12 A-smoothed peak / median rod
        ridge at 40-75 A from it >= 1.6; dot <= 100 A from the box centre (the re-centring put it ~90 A off; this
        drops ice-gradient class 58); one-sided (density in the +-30 deg sector opposite the rod, 35-100 A from the
        dot, < 25 % of that along the rod); no SEPARATE second blob (a local max >= 50 % of the dot, > 35 A away,
        with a gap < 35 % of its peak on the line back to the dot: drops dot-gap-dot class 53 but not duplex beads);
        rod body width (Gaussian sigma across the rod, 40-75 A) <= 25 A; ridge half-max reach >= 60 A.
        `strict_dot` flags the subset with ratio >= 2.3.
   FREE (dot-free straight rod, the no-hinge reference): round-1 dot-label fraction <= 0.35; the same dot/ridge ratio
        < 1.3; body width <= 15 A; total half-max length >= 80 A.
   CHANGES after the first trial run (justification): the ratio thresholds were first written as 2.3 / 2.0 by analogy
   with the first-pass `peak_body` column, but this script's ratio is defined differently (smoothed peak over the
   MEDIAN ridge 40-75 A) and 2.3 kept only 4 classes, including ice class 58 and excluding 8 of the 10 classes the
   task names; 1.6 / 1.3 separate the two populations (FREE <= 1.24 in the dot frame, the named DOT classes >= 1.69,
   intermediates 1.3-1.6 are excluded). A body-width cap of 15 A would have excluded class 130 for being WIDE, which
   is exactly the hinge signature under test, so it is 25 A for DOT (it only removes blobs: 117, 67 hit 40 A).
   Everything else (blobs, ice, dot pairs, rods with a dot in the middle, mixed classes) is listed as excluded.
2. Polar arm measurement (class-frame z map: (image - median) / MAD inside the 200 A radius 2D window). Origin = the
   dot (DOT) or the rod midpoint (FREE, two arms). Sampled every 1 A in r, 0.5 deg in angle (+-90 deg about the arm
   axis).
     mass m(r)    : z integrated across the rod along the arc |r * dphi| <= 30 A, but never narrower than +-30 deg
                    (so a rod swinging up to ~+-25 deg about the origin stays inside; hinge-insensitive length)
     ridge R(r)   : maximum of z within |t| <= 30 A of the axis (the "narrow strip" value; what a half-max reads)
     width s_t(r) : Gaussian sigma of the transverse profile (t = r * dphi, |t| <= 40 A) at each r
     body level   : median over r = 40-75 A (DOT) or 20-40 A from the midpoint (FREE; below 20 A the +-30 A arc
                    is truncated at the origin)
     L_ridge      : first r beyond the body where R < 0.5 * body (the first-pass style length)
     L_end, sig_end : fit m(r) = c + a * erfc((r - L_end) / (sqrt2 sig_end)) / 2 from the body start to r = 200 A
     hinge spread : fit s_t(r)^2 = s0^2 + (theta * r)^2 over the body (r from 45 A (DOT) / 20 A (FREE) out to
                    where the ridge falls below 30 % of its body level); theta is the rms in-plane angle about the
                    origin that the widening implies (0 for a rigid rod)
   DNA length of a DOT class = L_end - 17 A (dot centre to sfGFP surface, as in script 17); FREE length = sum of the
   two arms' L_end. B-DNA 50 bp = 170 A if straight and in-plane.
   A true end is sharp: sig_end ~ the point-spread sigma of the average, estimated per class as
   sqrt(s_t(body)^2 - 5^2) (a 20 A cylinder projects with transverse sigma 5 A).
3. Forward-blur calibration (made from REAL averages, no simulated particles): the four largest FREE averages are
   pivoted about a point 17 A beyond one end and averaged over Gaussian in-plane angles (sigma 0-25 deg), then
   measured exactly as the DOT classes. It maps the hinge angle theta to the L_ridge under-read, L_end and sig_end.
4. Re-alignment (b) on the particles of DOT classes 123, 130, 31, 147, 103, 114, 92, 145, 13, 52 (ALL particles, not a
   3-5k subset: the reading cost is per micrograph file and almost every file is touched anyway; the whole read
   takes minutes).
   - pose/shift convention re-verified on 3,000 A1 particles exactly as in script 17 (16 candidates; particle phase
     flipped with its own CTF, 15 A low-pass, contrast sign -1 as established there), required to be the same best
     convention as script 17 found for J27.
   - every particle is mapped into its class frame (3.55 A, 128 px, same grid as the class average); my plain
     re-average is checked against the CryoSPARC average (correlation).
   - dot masked (soft disc, radius 25 A, 5 A cosine edge); in-plane rotation about the dot from -30 to +30 deg in
     1 deg steps against a STRAIGHT-ROD reference: the transverse profile of the FREE-class particles (re-averaged the
     same way, central 60 A of each rod) laid along a line from the dot. The reference covers only r = 30-80 A, i.e.
     inside the rod that the largest DOT averages already show; a trial with 30-100 A showed the realigned half-max
     length snapping to the reference end (reference bias), so everything beyond 100 A ("far zone") is measured
     outside the reference and cannot be created by it.
   - realigned re-average measured as in item 2, plus the far-zone mass (100-170 A) / body.
   - per-particle angles from disjoint segments, S1 30-55 A, S2 55-80 A, S3 100-170 A (same reference profile).
     Independent noise does not covary between segments, so
       hinge spread sigma_h = sqrt(cov(phi_S1, phi_S2))  (true angle spread of the proximal rod, noise-free estimate)
       cov(phi_prox, phi_S3) ~ cov(S1, S2) if the particles carry a rigid rod at 100-170 A that turns with the
                            proximal rod (full-length duplex); ~0 if nothing is there. Reported as the difference
                            with a particle bootstrap CI, and as r(phi_prox, phi_S3) with a permutation null.
     DOT classes 97 and 87 (visible DNA ~90-110 A, reaching into S3) are added as positive controls.
   - NOTE (geometry): a sfGFP dot is near-round at 12 A, so the 2D in-plane rotation of a dot+rod particle is set by
     the rod; an in-plane hinge is therefore largely absorbed by the pose. What remains is rotational misalignment
     of weak rods (smear) and out-of-plane tilt (foreshortening, which this 2D test cannot separate from shortening).
   - linear density: transverse-integrated mass per A of the body in my re-averages (per-particle-sigma units), DOT
     vs FREE; out-of-plane tilt would RAISE it for foreshortened rods.
5. Helical repeat (c): ridge value along the body, detrended (minus a 20 A Gaussian smooth), Hann-windowed,
   zero-padded periodogram; peak period searched in 20-50 A (widened from 28-42 A so that an edge hit is visible
   as such); SNR = peak power / median power at periods 15-100 A outside +-25 % of the peak; snr34 = the same at
   34 A. Repeats along the body = body length / period (only meaningful where a 30-38 A peak with SNR >= 3 exists).
6. Length distributions (d): particle-weighted over classes; 95 % CIs by resampling classes with replacement
   (10,000 x), keeping each class's particle weight; difference FREE - DOT (DNA lengths).
7. Decision rule, written before the run:
   H2/H3 (full length, hinge-blurred) if the hinge-insensitive DOT DNA length (L_end - 17) is within 20 A of the FREE
     length, OR re-alignment lengthens L_end by >= 15 A, OR r(phi_prox, phi_S3) > 0 at p < 0.01 in most re-aligned
     classes;
   H1 (genuinely shorter) if DOT DNA L_end - 17 is >= 25 A below FREE with the class-bootstrap CI of the difference
     excluding 0, the DOT ends are about as sharp as the FREE ends (median sig_end ratio < 1.5 after PSF scaling),
     re-alignment changes L_end by < 10 A, and r(phi_prox, phi_S3) is not positive;
   otherwise mixed / inconclusive.

  PY=/n/holylabs/hekstra_lab/Everyone/ianpoon/envs/dnaspring/bin/python
  sbatch -p test -c 8 --mem=32G -t 2:00:00 -o /n/netscratch/hekstra_lab/Lab/ian_poon/logs/A1_hinge_%j.log \
         --wrap "$PY -I scripts/45_A1_length_vs_hinge.py all --nproc 8"

Outputs: results/A1_length_vs_hinge.json, figures/qc/A1_length_vs_hinge.png; intermediates under WORK below.
"""
import argparse
import glob
import importlib.util
import json
import os
import time
from multiprocessing import Pool

import mrcfile
import numpy as np
import pandas as pd
from scipy import ndimage as ndi
from scipy.optimize import curve_fit
from scipy.special import erfc

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
_spec = importlib.util.spec_from_file_location("rec17", os.path.join(REPO, "scripts", "17_rod_end_colocalisation.py"))
rec = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(rec)          # helpers only (pose convention, CTF phase flip); no side effects on import

PROJ = rec.PROJ
EXPORTS = rec.EXPORTS
A1_DIR = os.path.join(EXPORTS, "A1_2D")
PART_NPZ = os.path.join(EXPORTS, "A1_2D__particles.npz")
LAB_CSV = os.path.join(REPO, "results", "A1_class_by_round1_label.csv")
LEN_CSV = os.path.join(REPO, "results", "A1_class_lengths.csv")
CONV17 = os.path.join(rec.WORK, "conventions.json")
WORK = "/n/netscratch/hekstra_lab/Lab/ian_poon/sfgfp_chimera_screening_20261112/A1_length_vs_hinge"
OUT_JSON = os.path.join(REPO, "results", "A1_length_vs_hinge.json")
FIG = os.path.join(REPO, "figures", "qc", "A1_length_vs_hinge.png")

A_CA, N_CA = 3.55, 128            # class averages (CryoSPARC Fourier-cropped)
A_P, N_P = 2.84, 160              # particle stacks
WIN_A = 200.0                     # 2D window radius (400 A diameter)
DOT10 = [123, 130, 31, 147, 103, 114, 92, 145, 13, 52]
DOT_EXTRA = [97, 87]              # DOT classes whose visible DNA is long (L_end ~107-126 A): positive controls for S3
EXCLUDE = [105, 12]
MIN_N = 500
DOT_RATIO_MIN = 1.6               # dot (12 A smoothed peak) / median rod ridge 40-75 A from it
FREE_RATIO_MAX = 1.3              # the same ratio for a dot-free rod (beads of the duplex reach ~1.0-1.3)
R_DOT_A = 17.0
BODY_DOT = (40.0, 75.0)
BODY_FREE = (20.0, 40.0)        # from the midpoint; below 20 A the +-30 A arc is cut by the origin
WFIT0_DOT, WFIT0_FREE = 45.0, 20.0
RMAX = 200.0
DPHI = 0.5
MASS_T_A = 30.0
CAL_SIG = [0.0, 4.0, 8.0, 12.0, 16.0, 20.0, 25.0]
ANG = np.arange(-30, 31, 1.0)
SEGS = {"prox": (30.0, 80.0), "S1": (30.0, 55.0), "S2": (55.0, 80.0), "S3": (100.0, 170.0)}
FAR_A = (100.0, 170.0)            # far zone: beyond the reference, so re-alignment cannot put density there by bias
DOT_MASK_A = (20.0, 30.0)         # soft dot mask: 0 inside 20 A, cosine to 1 at 30 A (25 A radius, 5 A edge)
NCONV = 3000
NBOOT = 10000
SEED = 20261114

YY, XX = np.indices((N_CA, N_CA)).astype(np.float64)
C0 = N_CA / 2.0
RR = np.hypot(XX - C0, YY - C0)
WIN = RR < WIN_A / A_CA
S12 = 12.0 / 2.355 / A_CA


# ------------------------------------------------------------------------------------------------ image measurements
def zmap(im):
    v = im[WIN]
    med = np.median(v)
    mad = 1.4826 * np.median(np.abs(v - med)) + 1e-12
    return np.where(WIN, (im - med) / mad, 0.0)


def find_dot(z):
    sm = ndi.gaussian_filter(z, S12)
    cand = np.where(RR < 160 / A_CA, sm, -np.inf)
    py, px = np.unravel_index(np.argmax(cand), sm.shape)
    loc = (np.hypot(XX - px, YY - py) <= 2.0) & (sm > 0.5 * sm[py, px])
    w = sm * loc
    oy, ox = (w * YY).sum() / w.sum(), (w * XX).sum() / w.sum()
    # a SEPARATE second blob: a local maximum > 35 A away, >= 50 % of the dot, with a gap (< 35 % of its own peak)
    # on the straight line back to the dot (bright beads along an attached rod have no gap)
    mx = (ndi.maximum_filter(sm, size=int(round(20 / A_CA))) == sm) & WIN
    d = np.hypot(XX - ox, YY - oy) * A_CA
    second = 0.0
    for qy, qx in zip(*np.nonzero(mx & (d > 35) & (sm >= 0.5 * sm[py, px]))):
        tt = np.linspace(0, 1, 60)
        line = ndi.map_coordinates(sm, [oy + tt * (qy - oy), ox + tt * (qx - ox)], order=1)
        if line.min() < 0.35 * sm[qy, qx]:
            second = max(second, float(sm[qy, qx]))
    return float(ox), float(oy), float(sm[py, px]), second


def axis_from(z, ox, oy, r1=35.0, r2=110.0):
    d = np.hypot(XX - ox, YY - oy) * A_CA
    ann = (d > r1) & (d < r2) & WIN
    thr = 0.25 * z[ann].max()
    w = np.where(ann & (z > thr), z, 0.0)
    return float(np.arctan2((w * (YY - oy)).sum(), (w * (XX - ox)).sum()))


def polar(z, ox, oy, phi0):
    r = np.arange(0, RMAX + 1e-6, 1.0)
    ph = np.deg2rad(np.arange(-90, 90 + 1e-6, DPHI))
    R, PH = np.meshgrid(r, ph, indexing="ij")
    X = ox + R * np.cos(phi0 + PH) / A_CA
    Y = oy + R * np.sin(phi0 + PH) / A_CA
    return r, ph, ndi.map_coordinates(z, [Y, X], order=3, mode="constant", cval=0.0)


def _gauss(t, a, t0, s, c):
    return a * np.exp(-0.5 * ((t - t0) / s) ** 2) + c


def _erfend(r, c, a, L, s):
    return c + a * 0.5 * erfc((r - L) / (np.sqrt(2) * s))


def _first_below(r, y, start, level):
    i0 = int(np.searchsorted(r, start))
    for i in range(i0, len(r)):
        if y[i] < level:
            if i == 0:
                return float(r[0])
            f = (y[i - 1] - level) / max(y[i - 1] - y[i], 1e-12)
            return float(r[i - 1] + f * (r[i] - r[i - 1]))
    return float("nan")


def measure_arm(z, ox, oy, phi0, body, wfit0):
    r, ph, P = polar(z, ox, oy, phi0)
    dphi = np.deg2rad(DPHI)
    m = np.zeros(len(r))
    ridge = np.zeros(len(r))
    ridge_c = np.zeros(len(r))
    A = np.full(len(r), np.nan)
    T0 = np.full(len(r), np.nan)
    SG = np.full(len(r), np.nan)
    for i, ri in enumerate(r):
        if ri == 0:
            ridge[i] = ridge_c[i] = P[i, len(ph) // 2]
            continue
        Phi = np.clip(MASS_T_A / ri, np.deg2rad(30), np.pi / 2)
        m[i] = P[i, np.abs(ph) <= Phi].sum() * dphi * ri
        t = ri * ph
        s30 = np.abs(t) <= 30
        ridge[i] = P[i, s30].max()
        if ri >= 4 and ri % 2 == 0:
            s40 = np.abs(t) <= 40
            tt, yy = t[s40], P[i, s40]
            j = int(np.argmax(np.where(np.abs(tt) <= 30, yy, -np.inf)))
            try:
                p, _ = curve_fit(_gauss, tt, yy, p0=[yy[j] - np.median(yy), tt[j], 8.0, np.median(yy)],
                                 bounds=([0, -30, 2.0, -np.inf], [np.inf, 30, 40.0, np.inf]), maxfev=4000)
                A[i], T0[i], SG[i] = p[0], p[1], p[2]
            except Exception:
                pass
        # ridge-centre value (for the helical repeat): z on the ridge, i.e. at the narrow-strip maximum
        ridge_c[i] = ridge[i]
    bm = (r >= body[0]) & (r <= body[1])
    m_body = float(np.median(m[bm]))
    r_body = float(np.median(ridge[bm]))
    sg_body = float(np.nanmedian(SG[bm]))
    L_ridge = _first_below(r, ridge, body[1], 0.5 * r_body)
    msm = ndi.gaussian_filter1d(m, 3.0)
    L_mhalf = _first_below(r, msm, body[1], 0.5 * m_body)
    fit = dict(L_end=float("nan"), sig_end=float("nan"), L_end_err=float("nan"), sig_end_err=float("nan"),
               plateau=float("nan"), c=float("nan"))
    fr = r >= body[0]
    try:
        p0 = [0.0, m_body, L_mhalf if np.isfinite(L_mhalf) else body[1] + 20, 8.0]
        p, cov = curve_fit(_erfend, r[fr], m[fr], p0=p0,
                           bounds=([-np.inf, 0, body[0], 1.0], [np.inf, np.inf, RMAX, 80.0]), maxfev=20000)
        e = np.sqrt(np.diag(cov))
        fit = dict(L_end=float(p[2]), sig_end=float(p[3]), L_end_err=float(e[2]), sig_end_err=float(e[3]),
                   plateau=float(p[1]), c=float(p[0]))
    except Exception:
        pass
    # hinge spread from the widening of the rod with distance from the origin
    hinge = dict(theta_deg=float("nan"), s0=float("nan"), slope=float("nan"), r_range=None, n=0)
    if np.isfinite(fit["L_end"]):
        # body = from wfit0 out to where the ridge is still >= 30 % of its body level (the rod is clearly present)
        rhi = _first_below(r, ndi.gaussian_filter1d(ridge, 2.0), body[0], 0.3 * r_body)
        rhi = rhi if np.isfinite(rhi) else RMAX
        sel = (r >= wfit0) & (r <= rhi) & np.isfinite(SG)
        if sel.sum() >= 6:
            rr_, ss_ = r[sel], SG[sel]
            try:
                p, _ = curve_fit(lambda x, s0, th: np.sqrt(s0 ** 2 + (th * x) ** 2), rr_, ss_, p0=[8.0, 0.05],
                                 bounds=([0.5, 0.0], [40.0, 1.0]))
                sl = np.polyfit(rr_, ss_, 1)[0]
                hinge = dict(theta_deg=float(np.rad2deg(p[1])), s0=float(p[0]), slope=float(sl),
                             r_range=[float(rr_.min()), float(rr_.max())], n=int(sel.sum()))
            except Exception:
                pass
    psf = float(np.sqrt(max(sg_body ** 2 - 25.0, 1.0))) if np.isfinite(sg_body) else float("nan")
    return dict(r=r, m=m, ridge=ridge, ridge_c=ridge_c, A=A, T0=T0, SG=SG, m_body=m_body, ridge_body=r_body,
                width_body=sg_body, psf_sigma=psf, L_ridge=L_ridge, L_mhalf=L_mhalf, hinge=hinge, **fit)


def helical(prof_s, prof_v):
    """Periodogram of a detrended ridge profile sampled every 1 A."""
    if len(prof_v) < 40:
        return dict(period=float("nan"), snr=float("nan"), n_repeats=float("nan"), length=float(len(prof_v)))
    v = prof_v - ndi.gaussian_filter1d(prof_v, 20.0)
    v = (v - v.mean()) * np.hanning(len(v))
    nf = 4096
    f = np.fft.rfftfreq(nf, d=1.0)
    pw = np.abs(np.fft.rfft(v, nf)) ** 2
    per = np.where(f > 0, 1 / np.maximum(f, 1e-12), np.inf)
    sel = (per >= 20) & (per <= 50)
    k = np.nonzero(sel)[0][np.argmax(pw[sel])]
    pk = float(per[k])
    band = (per >= 15) & (per <= 100) & ((per < 0.75 * pk) | (per > 1.25 * pk))
    snr = float(pw[k] / np.median(pw[band]))
    L = float(prof_s[-1] - prof_s[0])
    k34 = int(np.argmin(np.abs(per - 34.0)))
    return dict(period=pk, snr=snr, n_repeats=L / pk, length=L, f=f, pw=pw,
                snr34=float(pw[k34] / np.median(pw[band])))


def arm_summary(a):
    keys = ["m_body", "ridge_body", "width_body", "psf_sigma", "L_ridge", "L_mhalf", "L_end", "sig_end", "L_end_err",
            "sig_end_err", "plateau"]
    out = {k: (None if not np.isfinite(a[k]) else round(float(a[k]), 2)) for k in keys}
    out["hinge"] = {k: (round(v, 3) if isinstance(v, float) and np.isfinite(v) else v) for k, v in a["hinge"].items()}
    return out


def analyse_dot(z, origin=None, phi0=None):
    ox, oy, peak, second = find_dot(z)
    if origin is not None:
        ox, oy = origin
    if phi0 is None:
        phi0 = axis_from(z, ox, oy)
    a = measure_arm(z, ox, oy, phi0, BODY_DOT, WFIT0_DOT)
    b = measure_arm(z, ox, oy, phi0 + np.pi, BODY_DOT, WFIT0_DOT)       # opposite side (one-sidedness)
    r = a["r"]
    sel = (r >= 35) & (r <= 100)
    one_sided = float(b["m"][sel].sum() / max(a["m"][sel].sum(), 1e-9))
    sm = ndi.gaussian_filter(z, S12)
    dot_ratio = float(sm[int(round(oy)), int(round(ox))] / max(a["ridge_body"], 1e-9))
    # helical profile along the body
    hb = (r >= 35) & (r <= (a["L_end"] - 10 if np.isfinite(a["L_end"]) else 80))
    hel = helical(r[hb], a["ridge_c"][hb])
    return dict(kind="dot", ox=ox, oy=oy, phi0=phi0, dot_peak=peak, second_peak=second, dot_ratio=dot_ratio,
                opp_ratio=one_sided, arm=a, hel=hel)


def analyse_free(z):
    sm = ndi.gaussian_filter(z, S12)
    cand = np.where(RR < 160 / A_CA, sm, -np.inf)
    py, px = np.unravel_index(np.argmax(cand), sm.shape)
    lb, _ = ndi.label(sm > 0.3 * sm[py, px])
    comp = lb == lb[py, px]
    ys, xs = np.nonzero(comp)
    w = sm[comp]
    mx, my = (w * xs).sum() / w.sum(), (w * ys).sum() / w.sum()
    ev, v = np.linalg.eigh(np.cov(np.vstack([xs - mx, ys - my])))
    u = v[:, 1]
    phi0 = float(np.arctan2(u[1], u[0]))
    # midpoint from a strip profile along the axis (half-max ends)
    s = np.arange(-120, 120.01, 1.0)
    prof = np.zeros_like(s)
    for t in np.arange(-6, 6.01, 2.0):
        prof += ndi.map_coordinates(z, [my + (s * u[1] + t * u[0]) / A_CA, mx + (s * u[0] - t * u[1]) / A_CA],
                                    order=1) / 7
    lev = 0.5 * np.median(prof[prof > 0.3 * prof.max()])
    i0 = int(np.argmax(prof))
    hi = i0
    while hi < len(s) - 1 and prof[hi] >= lev:
        hi += 1
    lo = i0
    while lo > 0 and prof[lo] >= lev:
        lo -= 1
    sc = 0.5 * (s[lo] + s[hi])
    ox, oy = mx + sc * u[0] / A_CA, my + sc * u[1] / A_CA
    a = measure_arm(z, ox, oy, phi0, BODY_FREE, WFIT0_FREE)
    b = measure_arm(z, ox, oy, phi0 + np.pi, BODY_FREE, WFIT0_FREE)
    body_ridge = 0.5 * (a["ridge_body"] + b["ridge_body"])
    dot_ratio = float(sm[comp].max() / max(body_ridge, 1e-9))
    rs = np.concatenate([-b["r"][::-1], a["r"][1:]])
    rv = np.concatenate([b["ridge_c"][::-1], a["ridge_c"][1:]])
    lo_, hi_ = -(b["L_end"] - 10), a["L_end"] - 10
    hb = (rs >= lo_) & (rs <= hi_) if np.isfinite(lo_ + hi_) else (np.abs(rs) <= 50)
    hel = helical(rs[hb], rv[hb])
    return dict(kind="free", ox=float(ox), oy=float(oy), phi0=phi0, dot_ratio=dot_ratio, arm=a, arm2=b, hel=hel,
                L_total=float(a["L_end"] + b["L_end"]), L_ridge_total=float(a["L_ridge"] + b["L_ridge"]),
                sig_end=float(np.nanmean([a["sig_end"], b["sig_end"]])),
                width_body=float(np.nanmean([a["width_body"], b["width_body"]])),
                psf_sigma=float(np.nanmean([a["psf_sigma"], b["psf_sigma"]])),
                theta_deg=float(np.nanmean([a["hinge"]["theta_deg"], b["hinge"]["theta_deg"]])))


def rotate_about(img, ox, oy, deg):
    a = np.deg2rad(deg)
    dx, dy = XX - ox, YY - oy
    sx = ox + np.cos(a) * dx - np.sin(a) * dy
    sy = oy + np.sin(a) * dx + np.cos(a) * dy
    return ndi.map_coordinates(img, [sy, sx], order=1, mode="constant", cval=0.0)


def hinge_blur(z, ox, oy, sig):
    if sig <= 0:
        return z.copy()
    angs = np.arange(-np.ceil(3 * sig), np.ceil(3 * sig) + 0.5, 1.0)
    w = np.exp(-0.5 * (angs / sig) ** 2)
    out = np.zeros_like(z)
    for a, wi in zip(angs, w):
        out += wi * rotate_about(z, ox, oy, a)
    return out / w.sum()


def fnum(x, nd=2):
    return None if x is None or not np.isfinite(x) else round(float(x), nd)


# ------------------------------------------------------------------------------------------------ stage: classes
def load_classes():
    mrcs = sorted(glob.glob(os.path.join(A1_DIR, "*class_averages*.mrc")))
    with mrcfile.open(mrcs[-1], permissive=True) as m:
        ca = np.asarray(m.data, np.float64)
        apix = float(m.voxel_size.x)
    assert abs(apix - A_CA) < 0.01 and ca.shape[1] == N_CA, (apix, ca.shape)
    P = np.load(PART_NPZ)
    n = np.bincount(P["alignments2D__class"].astype(int), minlength=ca.shape[0])
    lab = pd.read_csv(LAB_CSV).set_index("a1")
    frac = (lab["dot_long_rod"] + lab["dot_short_rod"]) / lab["n"]
    return ca, n, frac, os.path.basename(mrcs[-1])


def classify(ca, n, frac):
    rows, dot, free = [], [], []
    for k in np.argsort(-n):
        k = int(k)
        if n[k] < MIN_N:
            continue
        z = zmap(ca[k])
        fk = float(frac.get(k, np.nan))
        row = dict(cls=k, n=int(n[k]), frac_dotlabel=round(fk, 3))
        if k in EXCLUDE:
            row.update(group="excluded", why="lone dot + haze (rod not averaged)")
            rows.append(row)
            continue
        D = analyse_dot(z)
        F = analyse_free(z)
        a = D["arm"]
        reach = a["L_ridge"] if np.isfinite(a["L_ridge"]) else 0.0
        dot_off = float(np.hypot(D["ox"] - C0, D["oy"] - C0) * A_CA)
        crit_dot = dict(frac=fk >= 0.5, dot_ratio=D["dot_ratio"] >= DOT_RATIO_MIN, one_sided=D["opp_ratio"] < 0.25,
                        single=D["second_peak"] < 0.5 * D["dot_peak"], dot_near_centre=dot_off <= 100.0,
                        width=np.isfinite(a["width_body"]) and a["width_body"] <= 25.0, reach=reach >= 60.0)
        crit_free = dict(frac=fk <= 0.35, no_dot=D["dot_ratio"] < FREE_RATIO_MAX,
                         width=np.isfinite(F["width_body"]) and F["width_body"] <= 15.0,
                         length=np.isfinite(F["L_ridge_total"]) and F["L_ridge_total"] >= 80.0)
        row.update(dot_ratio_dotframe=fnum(D["dot_ratio"]), dot_offset_from_box_centre_A=fnum(dot_off, 1), opp_ratio=fnum(D["opp_ratio"], 3),
                   second_over_dot=fnum(D["second_peak"] / max(D["dot_peak"], 1e-9), 3),
                   width_body_dotframe=fnum(a["width_body"]), ridge_reach_dotframe=fnum(reach),
                   dot_ratio_freeframe=fnum(F["dot_ratio"]), width_body_freeframe=fnum(F["width_body"]),
                   L_ridge_total_freeframe=fnum(F["L_ridge_total"]))
        if all(crit_dot.values()):
            row.update(group="DOT", strict_dot=bool(D["dot_ratio"] >= 2.3))
            dot.append((k, D))
        elif all(crit_free.values()):
            row.update(group="FREE")
            free.append((k, F))
        else:
            row.update(group="excluded", why="DOT fails: " + ",".join(c for c, v in crit_dot.items() if not v)
                       + "; FREE fails: " + ",".join(c for c, v in crit_free.items() if not v))
        rows.append(row)
    return rows, dot, free


def calibration(ca, free, n):
    """Forward hinge blur of the largest FREE averages (see docstring item 3)."""
    big = sorted(free, key=lambda kf: -n[kf[0]])[:4]
    out, curves = [], {}
    for k, F in big:
        z = zmap(ca[k])
        L2 = F["arm2"]["L_end"]
        ux, uy = np.cos(F["phi0"]), np.sin(F["phi0"])
        px = F["ox"] - (L2 + R_DOT_A) * ux / A_CA
        py = F["oy"] - (L2 + R_DOT_A) * uy / A_CA
        for sg in CAL_SIG:
            zb = hinge_blur(z, px, py, sg)
            a = measure_arm(zb, px, py, F["phi0"], BODY_DOT, WFIT0_DOT)
            if k == big[0][0]:
                curves[f"cal{int(sg)}"] = np.vstack([a["r"], a["m"], a["ridge"], a["SG"]])
            out.append(dict(cls=k, sigma_deg=sg, L_ridge=fnum(a["L_ridge"]), L_end=fnum(a["L_end"]),
                            sig_end=fnum(a["sig_end"]), theta_fit_deg=fnum(a["hinge"]["theta_deg"]),
                            width_body=fnum(a["width_body"]), psf_sigma=fnum(a["psf_sigma"])))
    tab = pd.DataFrame(out)
    base = tab[tab.sigma_deg == 0].set_index("cls")
    tab["dL_ridge"] = tab.apply(lambda r: r.L_ridge - base.loc[r.cls, "L_ridge"], axis=1)
    tab["dL_end"] = tab.apply(lambda r: r.L_end - base.loc[r.cls, "L_end"], axis=1)
    tab["ridge_over_end"] = (tab.L_ridge - R_DOT_A) / (tab.L_end - R_DOT_A)
    summ = tab.groupby("sigma_deg")[["dL_ridge", "dL_end", "sig_end", "theta_fit_deg", "ridge_over_end"]].mean()
    return tab, summ, curves


def stage_classes(a):
    os.makedirs(WORK, exist_ok=True)
    ca, n, frac, mrcname = load_classes()
    rows, dot, free = classify(ca, n, frac)
    print(f"DOT {len(dot)}: {[k for k, _ in dot]}\nFREE {len(free)}: {[k for k, _ in free]}", flush=True)
    per = []
    for k, D in dot:
        a_ = D["arm"]
        per.append(dict(cls=k, group="DOT", n=int(n[k]), **arm_summary(a_),
                        DNA_L_end=fnum(a_["L_end"] - R_DOT_A), DNA_L_ridge=fnum(a_["L_ridge"] - R_DOT_A),
                        theta_deg=fnum(a_["hinge"]["theta_deg"]),
                        sig_end_over_psf=fnum(a_["sig_end"] / a_["psf_sigma"]),
                        hel_period=fnum(D["hel"]["period"]), hel_snr=fnum(D["hel"]["snr"]), hel_snr34=fnum(D["hel"].get("snr34")),
                        hel_repeats=fnum(D["hel"]["n_repeats"]), hel_len=fnum(D["hel"]["length"]),
                        dot_xy_A=[fnum((D["ox"] - C0) * A_CA, 1), fnum((D["oy"] - C0) * A_CA, 1)]))
    for k, F in free:
        per.append(dict(cls=k, group="FREE", n=int(n[k]), arm_plus=arm_summary(F["arm"]),
                        arm_minus=arm_summary(F["arm2"]), L_total=fnum(F["L_total"]),
                        L_ridge_total=fnum(F["L_ridge_total"]), sig_end=fnum(F["sig_end"]),
                        width_body=fnum(F["width_body"]), psf_sigma=fnum(F["psf_sigma"]),
                        sig_end_over_psf=fnum(F["sig_end"] / F["psf_sigma"]), theta_deg=fnum(F["theta_deg"]),
                        hel_period=fnum(F["hel"]["period"]), hel_snr=fnum(F["hel"]["snr"]), hel_snr34=fnum(F["hel"].get("snr34")),
                        hel_repeats=fnum(F["hel"]["n_repeats"]), hel_len=fnum(F["hel"]["length"])))
    tab, summ, curves = calibration(ca, free, n)
    print(pd.DataFrame([{k: v for k, v in p.items() if not isinstance(v, dict)} for p in per]).to_string(), flush=True)
    print(summ.round(2).to_string(), flush=True)
    # profiles for the figure
    prof = dict(curves)
    for k, D in dot:
        prof[f"dot{k}"] = np.vstack([D["arm"]["r"], D["arm"]["m"], D["arm"]["ridge"], D["arm"]["SG"]])
        prof[f"dotgeo{k}"] = np.array([D["ox"], D["oy"], D["phi0"]])
    for k, F in free:
        prof[f"free{k}"] = np.vstack([F["arm"]["r"], F["arm"]["m"], F["arm"]["ridge"], F["arm"]["SG"],
                                      F["arm2"]["m"], F["arm2"]["ridge"], F["arm2"]["SG"]])
        prof[f"freegeo{k}"] = np.array([F["ox"], F["oy"], F["phi0"]])
    for k, D in dot + free:
        if np.isfinite(D["hel"]["period"]):
            prof[f"hel{k}"] = np.vstack([D["hel"]["f"], D["hel"]["pw"]])
    np.savez_compressed(os.path.join(WORK, "class_profiles.npz"), **prof)
    res = dict(created=time.strftime("%Y-%m-%d %H:%M"), class_average_file=mrcname, selection=rows,
               per_class=per, calibration=dict(table=tab.round(2).to_dict("records"),
                                               mean_by_sigma=summ.round(2).reset_index().to_dict("records")))
    with open(os.path.join(WORK, "classes.json"), "w") as fh:
        json.dump(res, fh, indent=1)
    print("classes stage done", flush=True)


# ------------------------------------------------------------------------------------------------ stage: particles
G = {}                                    # read-only state shared with forked workers

PKEYS = ["alignments2D__class", "alignments2D__pose", "alignments2D__shift", "blob__path", "blob__idx",
         "ctf__df1_A", "ctf__df2_A", "ctf__df_angle_rad", "uid"]
JP = np.arange(N_P) - N_P / 2.0
YP, XP = np.meshgrid(JP, JP, indexing="ij")                      # particle-grid offsets (2.84 A px)
CIRC_P = np.hypot(XP, YP) < WIN_A / A_P - 1
XC = (XX - C0) * A_CA / A_P                                       # class-frame 3.55 A grid in 2.84 A px units
YC = (YY - C0) * A_CA / A_P


def _task(args):
    path, idx = args
    P, mode = G["P"], G["mode"]
    out = []
    with mrcfile.mmap(os.path.join(PROJ, path), permissive=True, mode="r") as m:
        data = m.data if m.data.ndim == 3 else m.data[None]          # a one-particle stack is stored 2D
        for i in idx:
            img = np.array(data[int(P["blob__idx"][i])], np.float64)
            pimg = G["psign"] * rec.flip_lp(img, P["ctf__df1_A"][i], P["ctf__df2_A"][i], P["ctf__df_angle_rad"][i],
                                            A_P)
            th = float(P["alignments2D__pose"][i])
            sx, sy = (float(v) for v in P["alignments2D__shift"][i])
            if mode == "conv":
                ref = G["ca160"][int(P["alignments2D__class"][i])][CIRC_P]
                row = []
                for cand in rec.CANDS:
                    px, py = rec.class_to_particle(XP, YP, th, sx, sy, cand)
                    smp = ndi.map_coordinates(pimg, [N_P / 2 + py, N_P / 2 + px], order=1, mode="constant")
                    row.append(np.corrcoef(smp[CIRC_P], ref)[0, 1])
                row.append(np.corrcoef(pimg[CIRC_P], ref)[0, 1])
                out.append(row)
            else:
                px, py = rec.class_to_particle(XC, YC, th, sx, sy, G["conv"])
                smp = ndi.map_coordinates(pimg, [N_P / 2 + py, N_P / 2 + px], order=1, mode="constant")
                v = smp[WIN]
                smp = np.where(WIN, (smp - v.mean()) / (v.std() + 1e-9), 0.0)
                out.append(smp.astype(np.float16))
    return idx, np.array(out)


def run_pool(idx_all, nproc):
    P = G["P"]
    by = {}
    for i in idx_all:
        by.setdefault(str(P["blob__path"][i]), []).append(int(i))
    tasks = list(by.items())
    res = {}
    t0 = time.time()
    with Pool(nproc) as pool:
        for j, (idx, arr) in enumerate(pool.imap_unordered(_task, tasks, chunksize=4)):
            for i, a in zip(idx, arr):
                res[i] = a
            if j % 1000 == 0:
                print(f"  {j}/{len(tasks)} files, {len(res)} particles, {time.time() - t0:.0f} s", flush=True)
    return res, len(tasks)


def soft_disc(ox, oy, r_in, r_out):
    d = np.hypot(XX - ox, YY - oy) * A_CA
    return np.clip((d - r_in) / (r_out - r_in), 0, 1) ** 2 * (3 - 2 * np.clip((d - r_in) / (r_out - r_in), 0, 1))


def transverse_profile(z, F):
    """Mean transverse profile of a FREE re-average over its central body (|s| <= 30 A, |t| <= 40 A, 1 A steps)."""
    u = np.array([np.cos(F["phi0"]), np.sin(F["phi0"])])
    s = np.arange(-30, 30.001, 1.0)             # central body only (every FREE class is >= 90 A long)
    t = np.arange(-40, 40.001, 1.0)
    S, T = np.meshgrid(s, t, indexing="ij")
    X = F["ox"] + (S * u[0] - T * u[1]) / A_CA
    Y = F["oy"] + (S * u[1] + T * u[0]) / A_CA
    return t, ndi.map_coordinates(z, [Y, X], order=1).mean(0)


def rod_refs(ox, oy, phi0, tprof, seg):
    """Straight-rod references (transverse profile tprof laid from the dot) for every angle in ANG; (nang, npix)."""
    t_ax, T = tprof
    lo, hi = seg
    out = []
    dx, dy = (XX - ox) * A_CA, (YY - oy) * A_CA
    for a in ANG:
        ps = phi0 + np.deg2rad(a)
        s = dx * np.cos(ps) + dy * np.sin(ps)
        t = -dx * np.sin(ps) + dy * np.cos(ps)
        tap = np.clip((s - lo) / 5.0, 0, 1) * np.clip((hi - s) / 5.0, 0, 1)
        tap = np.sin(0.5 * np.pi * tap) ** 2
        R = np.interp(t, t_ax, T, left=0.0, right=0.0) * tap * (np.abs(t) <= 40)
        out.append(R.ravel())
    return np.array(out, np.float32)


def best_angle(S):
    """argmax over ANG with parabolic sub-degree refinement; S (n, nang)."""
    k = np.argmax(S, 1)
    a = ANG[k].astype(np.float64)
    inner = (k > 0) & (k < len(ANG) - 1)
    ii = np.nonzero(inner)[0]
    y0, y1, y2 = S[ii, k[ii] - 1], S[ii, k[ii]], S[ii, k[ii] + 1]
    den = y0 - 2 * y1 + y2
    off = np.where(np.abs(den) > 1e-12, 0.5 * (y0 - y2) / den, 0.0)
    a[ii] += np.clip(off, -0.5, 0.5)
    return a, ~inner


def corr_stats(x, y, rng, nperm=2000, nboot=1000):
    r = float(np.corrcoef(x, y)[0, 1])
    perm = np.array([np.corrcoef(x, rng.permutation(y))[0, 1] for _ in range(nperm)])
    n = len(x)
    boot = []
    for _ in range(nboot):
        j = rng.integers(0, n, n)
        boot.append(np.corrcoef(x[j], y[j])[0, 1])
    return dict(r=round(r, 4), p_perm=float((1 + np.sum(perm >= r)) / (nperm + 1)),
                null_95=round(float(np.quantile(perm, 0.95)), 4),
                ci95=[round(float(np.quantile(boot, 0.025)), 4), round(float(np.quantile(boot, 0.975)), 4)], n=n)


def stage_particles(a):
    t0 = time.time()
    rng = np.random.default_rng(SEED)
    C = json.load(open(os.path.join(WORK, "classes.json")))
    free_cls = [p["cls"] for p in C["per_class"] if p["group"] == "FREE"]
    d = np.load(PART_NPZ)
    P = {k: d[k] for k in PKEYS}
    cls = P["alignments2D__class"].astype(int)
    ca, n, frac, _ = load_classes()
    conv17 = json.load(open(CONV17))
    psign = float(conv17["particle_contrast_sign_vs_micrograph"])
    G.update(P=P, psign=psign)

    # 1. pose / shift convention, re-verified on A1 particles (script 17 recipe)
    dot_cls = DOT10 + DOT_EXTRA
    pool_cls = dot_cls + free_cls
    cand_idx = np.nonzero(np.isin(cls, pool_cls))[0]
    sub = np.sort(rng.choice(cand_idx, size=min(a.nconv, len(cand_idx)), replace=False))
    ca160 = {}
    for k in pool_cls:
        up = ndi.map_coordinates(ca[k], [C0 + YP * A_P / A_CA, C0 + XP * A_P / A_CA], order=3, mode="constant")
        ca160[k] = rec.lp_only(up, A_P)
    G.update(mode="conv", ca160=ca160)
    print(f"convention test on {len(sub)} particles", flush=True)
    res, nf = run_pool(sub, a.nproc)
    corr = np.array([res[i] for i in sub])
    mean = corr.mean(0)
    sem = corr.std(0) / np.sqrt(len(sub))
    best = int(np.argmax(mean[:-1]))
    order = np.argsort(-mean[:-1])
    conv_tab = [dict(name=rec.cand_name(rec.CANDS[c]), mean_corr=round(float(mean[c]), 4), sem=round(float(sem[c]), 4),
                     frac_best=round(float(np.mean(np.argmax(corr[:, :-1], 1) == c)), 3)) for c in range(len(rec.CANDS))]
    conv_tab.append(dict(name="identity_no_alignment", mean_corr=round(float(mean[-1]), 4), sem=round(float(sem[-1]), 4)))
    conv = rec.CANDS[best]
    same17 = rec.cand_name(conv) == conv17["best"]["name"]
    print(f"best {rec.cand_name(conv)} r={mean[best]:.4f} (2nd {rec.cand_name(rec.CANDS[order[1]])} "
          f"{mean[order[1]]:.4f}); script 17 best {conv17['best']['name']}; same={same17}", flush=True)
    G.pop("ca160")

    # 2. map every particle of DOT10 + FREE classes into its class frame
    G.update(mode="map", conv=conv)
    idx_all = np.nonzero(np.isin(cls, pool_cls))[0]
    if a.maxper:                                    # debugging only: cap the particles per class
        idx_all = np.sort(np.concatenate([rng.permutation(idx_all[cls[idx_all] == k])[:a.maxper] for k in pool_cls]))
    print(f"mapping {len(idx_all)} particles into class frames", flush=True)
    res, nfiles = run_pool(idx_all, a.nproc)
    print(f"read {len(res)} particles from {nfiles} files in {time.time() - t0:.0f} s", flush=True)

    out = dict(convention=dict(table=conv_tab, best=rec.cand_name(conv), best_mean_corr=round(float(mean[best]), 4),
                               second=rec.cand_name(rec.CANDS[order[1]]),
                               second_mean_corr=round(float(mean[order[1]]), 4), n_particles=int(len(sub)),
                               same_as_script17=bool(same17), particle_contrast_sign=psign),
               n_files_read=int(nfiles), classes={}, free_reaverages={})
    # 3. plain re-averages; FREE transverse profile
    stacks, avgs = {}, {}
    for k in pool_cls:
        ii = idx_all[cls[idx_all] == k]
        n[k] = len(ii)
        stacks[k] = np.stack([res[i] for i in ii]).reshape(len(ii), -1)
        avgs[k] = stacks[k].astype(np.float64).mean(0).reshape(N_CA, N_CA)
    del res
    tps, wts = [], []
    for k in free_cls:
        z = zmap(avgs[k])
        F = analyse_free(z)
        t_ax, tp = transverse_profile(z, F)
        tps.append(tp / max(tp.max(), 1e-9))
        wts.append(n[k])
        Fr = measure_arm(avgs[k], F["ox"], F["oy"], F["phi0"], BODY_FREE, WFIT0_FREE)
        out["free_reaverages"][int(k)] = dict(n=int(n[k]), L_total=fnum(F["L_total"]),
                                             lin_density_raw=fnum(Fr["m_body"] / 1.0, 4),
                                             L_ridge_total=fnum(F["L_ridge_total"]), sig_end=fnum(F["sig_end"]),
                                             width_body=fnum(F["width_body"]), theta_deg=fnum(F["theta_deg"]),
                                             corr_with_cs_average=fnum(np.corrcoef(avgs[k][WIN], ca[k][WIN])[0, 1], 3))
    T = np.average(np.array(tps), axis=0, weights=wts)
    T = T - T.mean()
    T = T / np.sqrt((T ** 2).sum())
    tprof = (t_ax, T)
    figdat = dict(tprof=np.vstack([t_ax, T]))

    # 4. dot-free re-alignment of the rod body about the dot
    for k in dot_cls:
        X = stacks[k]
        zp = zmap(avgs[k])
        D = analyse_dot(zp)
        ox, oy, phi0 = D["ox"], D["oy"], D["phi0"]
        msk = soft_disc(ox, oy, *DOT_MASK_A).ravel().astype(np.float32)
        Xm = X.astype(np.float32) * msk
        ang = {}
        edge = {}
        for sname, seg in SEGS.items():
            S = Xm @ rod_refs(ox, oy, phi0, tprof, seg).T
            ang[sname], edge[sname] = best_angle(S)
        # realign on the proximal angle
        acc = np.zeros((N_CA, N_CA))
        dx, dy = XX - ox, YY - oy
        for j in range(len(X)):
            th = np.deg2rad(ang["prox"][j])
            sx = ox + np.cos(th) * dx - np.sin(th) * dy
            sy = oy + np.sin(th) * dx + np.cos(th) * dy
            acc += ndi.map_coordinates(X[j].astype(np.float64).reshape(N_CA, N_CA), [sy, sx], order=1)
        ra = acc / len(X)
        zr = zmap(ra)
        before = measure_arm(zp, ox, oy, phi0, BODY_DOT, WFIT0_DOT)
        after = measure_arm(zr, ox, oy, phi0, BODY_DOT, WFIT0_DOT)
        p1, p2, pp, p3 = ang["S1"], ang["S2"], ang["prox"], ang["S3"]
        c12 = float(np.cov(p1, p2)[0, 1])
        cp3 = float(np.cov(pp, p3)[0, 1])
        sh_boot, rat_boot, dif_boot = [], [], []
        for _ in range(1000):
            j = rng.integers(0, len(p1), len(p1))
            b12 = np.cov(p1[j], p2[j])[0, 1]
            bp3 = np.cov(pp[j], p3[j])[0, 1]
            sh_boot.append(np.sqrt(max(b12, 0.0)))
            dif_boot.append(bp3 - b12)
            rat_boot.append(bp3 / b12 if b12 > 0 else np.nan)
        # far-zone mass (100-170 A) relative to the body, before / after; raw linear density of the body
        rawb = measure_arm(avgs[k], ox, oy, phi0, BODY_DOT, WFIT0_DOT)
        def far_frac(arm):
            sel = (arm["r"] >= FAR_A[0]) & (arm["r"] <= FAR_A[1])
            return float(arm["m"][sel].mean() / max(arm["m_body"], 1e-9))
        out["classes"][int(k)] = dict(
            n=int(len(X)), corr_plain_vs_cs_average=fnum(np.corrcoef(avgs[k][WIN], ca[k][WIN])[0, 1], 3),
            dot_xy_A=[fnum((ox - C0) * A_CA, 1), fnum((oy - C0) * A_CA, 1)], phi0_deg=fnum(np.rad2deg(phi0), 1),
            before=arm_summary(before), after=arm_summary(after),
            dL_end=fnum(after["L_end"] - before["L_end"]), dL_ridge=fnum(after["L_ridge"] - before["L_ridge"]),
            far_mass_frac_before=fnum(far_frac(before), 4), far_mass_frac_after=fnum(far_frac(after), 4),
            lin_density_raw=fnum(rawb["m_body"], 4), extra_control=bool(k in DOT_EXTRA),
            cov_test=dict(cov_S1_S2=fnum(c12, 3), cov_prox_S3=fnum(cp3, 3),
                          diff=fnum(cp3 - c12, 3), diff_ci95=[fnum(np.nanquantile(dif_boot, 0.025), 3),
                                                              fnum(np.nanquantile(dif_boot, 0.975), 3)],
                          ratio=fnum(cp3 / c12 if c12 > 0 else np.nan, 3),
                          ratio_ci95=[fnum(np.nanquantile(rat_boot, 0.025), 3),
                                      fnum(np.nanquantile(rat_boot, 0.975), 3)]),
            angles=dict(prox_mean=fnum(pp.mean()), prox_sd=fnum(pp.std()),
                        prox_frac_edge=fnum(edge["prox"].mean(), 3), S3_frac_edge=fnum(edge["S3"].mean(), 3),
                        S3_sd=fnum(p3.std()),
                        sigma_hinge_cov12_deg=fnum(np.sqrt(max(c12, 0.0))),
                        sigma_hinge_ci95=[fnum(np.quantile(sh_boot, 0.025)), fnum(np.quantile(sh_boot, 0.975))],
                        cov12=fnum(c12, 3)),
            corr_S1_S2=corr_stats(p1, p2, rng), corr_prox_S3=corr_stats(pp, p3, rng),
            corr_prox_S3_nonedge=corr_stats(pp[~edge["prox"] & ~edge["S3"]], p3[~edge["prox"] & ~edge["S3"]], rng))
        c = out["classes"][int(k)]
        print(f"class {k}: n={len(X)} corr(cs)={c['corr_plain_vs_cs_average']} L_end {before['L_end']:.1f}->"
              f"{after['L_end']:.1f}  L_ridge {before['L_ridge']:.1f}->{after['L_ridge']:.1f}  "
              f"theta {before['hinge']['theta_deg']:.1f}->{after['hinge']['theta_deg']:.1f}  "
              f"sd(prox) {pp.std():.1f}  sig_h {np.sqrt(max(c12, 0)):.1f}  r(S1,S2) {c['corr_S1_S2']['r']}  "
              f"r(prox,S3) {c['corr_prox_S3']['r']} p={c['corr_prox_S3']['p_perm']:.4f}", flush=True)
        figdat[f"plain{k}"] = zp.astype(np.float32)
        figdat[f"realigned{k}"] = zr.astype(np.float32)
        figdat[f"angles{k}"] = np.vstack([pp, p1, p2, p3]).astype(np.float32)
        figdat[f"prof{k}"] = np.vstack([before["r"], before["m"], before["ridge"], before["SG"],
                                        after["m"], after["ridge"], after["SG"]])
        figdat[f"geo{k}"] = np.array([ox, oy, phi0])
    np.savez_compressed(os.path.join(WORK, "particle_profiles.npz"), **figdat)
    out["created"] = time.strftime("%Y-%m-%d %H:%M")
    out["runtime_s"] = round(time.time() - t0)
    with open(os.path.join(WORK, "particles.json"), "w") as fh:
        json.dump(out, fh, indent=1)
    print(f"particles stage done in {time.time() - t0:.0f} s", flush=True)


# ------------------------------------------------------------------------------------------------ stage: report
def wquant(x, w, q):
    o = np.argsort(x)
    x, w = np.asarray(x)[o], np.asarray(w, float)[o]
    c = (np.cumsum(w) - 0.5 * w) / w.sum()
    return float(np.interp(q, c, x))


def boot_groups(xd, wd, xf, wf, rng, nb=NBOOT):
    """Particle-weighted mean / median per group, class-resampling bootstrap; difference FREE - DOT."""
    xd, wd, xf, wf = map(np.asarray, (xd, wd, xf, wf))
    est = dict(dot_mean=np.average(xd, weights=wd), free_mean=np.average(xf, weights=wf),
               dot_median=wquant(xd, wd, 0.5), free_median=wquant(xf, wf, 0.5))
    est["diff_mean"] = est["free_mean"] - est["dot_mean"]
    est["diff_median"] = est["free_median"] - est["dot_median"]
    B = {k: [] for k in est}
    for _ in range(nb):
        i = rng.integers(0, len(xd), len(xd))
        j = rng.integers(0, len(xf), len(xf))
        dm, fm = np.average(xd[i], weights=wd[i]), np.average(xf[j], weights=wf[j])
        dmed, fmed = wquant(xd[i], wd[i], 0.5), wquant(xf[j], wf[j], 0.5)
        for k, v in zip(B, (dm, fm, dmed, fmed, fm - dm, fmed - dmed)):
            B[k].append(v)
    return {k: dict(est=round(float(v), 2), ci95=[round(float(np.quantile(B[k], 0.025)), 2),
                                                   round(float(np.quantile(B[k], 0.975)), 2)]) for k, v in est.items()}


def stage_report(a):
    rng = np.random.default_rng(SEED + 7)
    C = json.load(open(os.path.join(WORK, "classes.json")))
    Pj = json.load(open(os.path.join(WORK, "particles.json")))
    cp = np.load(os.path.join(WORK, "class_profiles.npz"))
    pp = np.load(os.path.join(WORK, "particle_profiles.npz"))
    per = {p["cls"]: p for p in C["per_class"]}
    D = [p for p in C["per_class"] if p["group"] == "DOT" and p["L_end"] is not None]
    F = [p for p in C["per_class"] if p["group"] == "FREE" and p["L_total"] is not None]
    cal = pd.DataFrame(C["calibration"]["mean_by_sigma"])

    # (d) length distributions
    xd = np.array([p["L_end"] - R_DOT_A for p in D])
    xd0 = np.array([p["L_end"] for p in D])
    xdr = np.array([p["L_ridge"] - R_DOT_A for p in D])
    wd = np.array([p["n"] for p in D], float)
    xf = np.array([p["L_total"] for p in F])
    xfr = np.array([p["L_ridge_total"] for p in F])
    wf = np.array([p["n"] for p in F], float)
    lengths = dict(
        hinge_insensitive_L_end=boot_groups(xd, wd, xf, wf, rng),
        junction_offset_0A=boot_groups(xd0, wd, xf, wf, rng),
        half_max_ridge=boot_groups(xdr, wd, xfr, wf, rng),
        DOT10_only=boot_groups(np.array([p["L_end"] - R_DOT_A for p in D if p["cls"] in DOT10]),
                               np.array([p["n"] for p in D if p["cls"] in DOT10], float), xf, wf, rng))
    free_med = lengths["hinge_insensitive_L_end"]["free_median"]["est"]
    long_frac = float(wd[xd >= free_med - 15].sum() / wd.sum())
    # ends and widening
    sd_ = np.array([p["sig_end"] for p in D])
    sf_ = np.array([p["sig_end"] for p in F])
    rd_ = np.array([p["sig_end_over_psf"] for p in D])
    rf_ = np.array([p["sig_end_over_psf"] for p in F])
    td_ = np.array([p["theta_deg"] if p["theta_deg"] is not None else np.nan for p in D])
    tf_ = np.array([p["theta_deg"] if p["theta_deg"] is not None else np.nan for p in F])
    okd, okf = np.isfinite(td_), np.isfinite(tf_)
    ends = dict(sig_end=boot_groups(sd_, wd, sf_, wf, rng),
                sig_end_over_psf=boot_groups(rd_, wd, rf_, wf, rng),
                theta_widening_deg=boot_groups(td_[okd], wd[okd], tf_[okf], wf[okf], rng))
    th_dot = ends["theta_widening_deg"]["dot_median"]["est"]
    th_hi = max(ends["theta_widening_deg"]["dot_median"]["ci95"][1], float(np.nanmax(td_)))
    pred = dict(dL_ridge_at_dot_median_theta=fnum(np.interp(th_dot, cal.sigma_deg, cal.dL_ridge)),
                dL_ridge_at_max_dot_theta=fnum(np.interp(th_hi, cal.sigma_deg, cal.dL_ridge)),
                dL_end_at_max_dot_theta=fnum(np.interp(th_hi, cal.sigma_deg, cal.dL_end)),
                theta_needed_for_observed_gap_deg=fnum(np.interp(
                    -lengths["hinge_insensitive_L_end"]["diff_median"]["est"], cal.dL_ridge[::-1],
                    cal.sigma_deg[::-1])) if lengths["hinge_insensitive_L_end"]["diff_median"]["est"] <=
                -cal.dL_ridge.min() else "> 25 deg (beyond calibration range)")
    # (b) re-alignment summary
    R = {int(k): v for k, v in Pj["classes"].items()}
    main = [k for k in DOT10 if k in R]
    w10 = np.array([R[k]["n"] for k in main], float)
    dLe = np.array([R[k]["dL_end"] if R[k]["dL_end"] is not None else np.nan for k in main])
    dLr = np.array([R[k]["dL_ridge"] if R[k]["dL_ridge"] is not None else np.nan for k in main])
    farb = np.array([R[k]["far_mass_frac_before"] for k in main])
    fara = np.array([R[k]["far_mass_frac_after"] for k in main])
    sigh = np.array([R[k]["angles"]["sigma_hinge_cov12_deg"] for k in main])
    pos_S3 = [k for k in main if R[k]["corr_prox_S3"]["p_perm"] < 0.01 and R[k]["corr_prox_S3"]["r"] > 0]
    pos_12 = [k for k in main if R[k]["corr_S1_S2"]["p_perm"] < 0.01 and R[k]["corr_S1_S2"]["r"] > 0]
    okL = np.isfinite(dLe)
    realign = dict(
        classes=main, n_particles=int(w10.sum()),
        dL_end_weighted_mean=fnum(np.average(dLe[okL], weights=w10[okL])),
        dL_ridge_weighted_mean=fnum(np.average(dLr[np.isfinite(dLr)], weights=w10[np.isfinite(dLr)])),
        far_mass_frac_before_wmean=fnum(np.average(farb, weights=w10), 4),
        far_mass_frac_after_wmean=fnum(np.average(fara, weights=w10), 4),
        sigma_hinge_deg_weighted_mean=fnum(np.average(sigh, weights=w10)),
        sigma_hinge_deg_range=[fnum(sigh.min()), fnum(sigh.max())],
        classes_with_r_S1_S2_positive_p01=pos_12, classes_with_r_prox_S3_positive_p01=pos_S3,
        convention=Pj["convention"])
    # decision rule (docstring item 7)
    gap = lengths["hinge_insensitive_L_end"]["diff_median"]
    gap_mean = lengths["hinge_insensitive_L_end"]["diff_mean"]
    endratio = ends["sig_end_over_psf"]["dot_median"]["est"] / ends["sig_end_over_psf"]["free_median"]["est"]
    h2 = (abs(gap["est"]) <= 20 or (realign["dL_end_weighted_mean"] or 0) >= 15 or len(pos_S3) > len(main) / 2)
    h1 = (gap["est"] >= 25 and gap["ci95"][0] > 0 and endratio < 1.5
          and abs(realign["dL_end_weighted_mean"] or 0) < 10 and len(pos_S3) == 0)
    rule = dict(gap_median_A=gap, gap_mean_A=gap_mean, sig_end_ratio_dot_over_free_psf_scaled=fnum(endratio, 3),
                realign_dL_end=realign["dL_end_weighted_mean"], n_classes_S3_positive=len(pos_S3),
                H2_criteria_met=bool(h2), H1_criteria_met=bool(h1),
                outcome="H1" if (h1 and not h2) else ("H2/H3" if (h2 and not h1) else "mixed/inconclusive"))
    # helical repeat
    hel = dict(
        DOT=[dict(cls=p["cls"], period=p["hel_period"], snr=p["hel_snr"], snr34=p.get("hel_snr34"),
                  repeats=p["hel_repeats"], body_len=p["hel_len"]) for p in D],
        FREE=[dict(cls=p["cls"], period=p["hel_period"], snr=p["hel_snr"], snr34=p.get("hel_snr34"),
                   repeats=p["hel_repeats"], body_len=p["hel_len"]) for p in F])
    for g in ("DOT", "FREE"):
        s34 = [h["snr34"] for h in hel[g] if h["snr34"] is not None]
        near = [h for h in hel[g] if h["period"] is not None and 30 <= h["period"] <= 38 and (h["snr"] or 0) >= 3]
        hel[g + "_summary"] = dict(median_snr_at_34A=fnum(np.median(s34)) if s34 else None,
                                   n_classes_peak_30_38A_snr3=len(near), n_classes=len(hel[g]),
                                   repeats_where_peak=[h["repeats"] for h in near])
    # linear density (per-particle-sigma units, my phase-flipped re-averages)
    lin = dict(DOT={int(k): R[k]["lin_density_raw"] for k in R},
               FREE={int(k): v["lin_density_raw"] for k, v in Pj["free_reaverages"].items()})
    lin["DOT_wmean"] = fnum(np.average([R[k]["lin_density_raw"] for k in main], weights=w10), 3)
    lin["FREE_wmean"] = fnum(np.average([v["lin_density_raw"] for v in Pj["free_reaverages"].values()],
                                        weights=[v["n"] for v in Pj["free_reaverages"].values()]), 3)

    out = dict(created=time.strftime("%Y-%m-%d %H:%M"), script="scripts/45_A1_length_vs_hinge.py",
               inputs=dict(class_averages=os.path.join(A1_DIR, C["class_average_file"]), particles=PART_NPZ,
                           labels=LAB_CSV, first_pass=LEN_CSV),
               definitions=dict(junction_offset_A=R_DOT_A, body_dot_A=BODY_DOT, body_free_A=BODY_FREE,
                                mass_arc_A=MASS_T_A, dot_ratio_min=DOT_RATIO_MIN, free_ratio_max=FREE_RATIO_MAX,
                                realign_segments_A=SEGS, far_zone_A=FAR_A, angles_deg=[float(ANG[0]), float(ANG[-1])],
                                dot_mask_A=DOT_MASK_A, B_DNA_50bp_A=170.0),
               selection=C["selection"], per_class=C["per_class"], calibration=C["calibration"],
               lengths=lengths, dot_frac_particles_at_free_length=fnum(long_frac, 3), ends=ends,
               hinge_prediction_from_calibration=pred, realign=realign, realign_per_class=Pj["classes"],
               free_reaverages=Pj["free_reaverages"], helical=hel, linear_density=lin, decision=rule)
    with open(OUT_JSON, "w") as fh:
        json.dump(out, fh, indent=1)
    print(json.dumps(dict(lengths=lengths, ends=ends, pred=pred, realign={k: v for k, v in realign.items()
                                                                           if k != "convention"},
                          long_frac=long_frac, rule=rule, lin={k: v for k, v in lin.items() if "wmean" in k},
                          hel={k: v for k, v in hel.items() if "summary" in k}), indent=1), flush=True)
    make_figure(C, Pj, cp, pp, out)


def make_figure(C, Pj, cp, pp, out):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    ca, n, frac, _ = load_classes()
    per = {p["cls"]: p for p in C["per_class"]}
    D = [p for p in C["per_class"] if p["group"] == "DOT" and p["L_end"] is not None]
    F = [p for p in C["per_class"] if p["group"] == "FREE" and p["L_total"] is not None]
    cD, cF = "#c0392b", "#2c6fbb"
    fig, ax = plt.subplots(3, 5, figsize=(26, 15.5))
    lim = (C0 - 48, C0 + 48)
    # row 0: example averages with measured ends
    for j, k in enumerate([123, 130, 13, F[0]["cls"]]):
        A = ax[0, j]
        A.imshow(ca[k], cmap="gray", origin="upper")
        p = per[k]
        if p["group"] == "DOT":
            ox, oy, ph = cp[f"dotgeo{k}"]
            for L, col, lab in ((p["L_end"], "r", "L_end (mass erf)"), (p["L_ridge"], "c", "L_ridge (half-max)")):
                A.plot([ox, ox + L * np.cos(ph) / A_CA], [oy, oy + L * np.sin(ph) / A_CA], col, lw=1.2, label=lab)
            for s in (-30, 30):
                A.plot([ox, ox + 170 * np.cos(ph + np.deg2rad(s)) / A_CA],
                       [oy, oy + 170 * np.sin(ph + np.deg2rad(s)) / A_CA], "y:", lw=0.7)
            A.set_title(f"DOT class {k} (n={p['n']:,})\nDNA L_end-17 = {p['DNA_L_end']:.0f} A, "
                        f"sig_end {p['sig_end']:.0f} A, widening {p['theta_deg'] or 0:.1f} deg", fontsize=9)
        else:
            ox, oy, ph = cp[f"freegeo{k}"]
            L1, L2 = p["arm_plus"]["L_end"], p["arm_minus"]["L_end"]
            A.plot([ox - L2 * np.cos(ph) / A_CA, ox + L1 * np.cos(ph) / A_CA],
                   [oy - L2 * np.sin(ph) / A_CA, oy + L1 * np.sin(ph) / A_CA], "r", lw=1.2)
            A.set_title(f"FREE class {k} (n={p['n']:,})\nL_end total = {p['L_total']:.0f} A, "
                        f"sig_end {p['sig_end']:.0f} A", fontsize=9)
        A.plot([lim[0] + 4, lim[0] + 4 + 100 / A_CA], [lim[1] - 5] * 2, color="gold", lw=3)
        A.set_xlim(*lim)
        A.set_ylim(lim[1], lim[0])
        A.set_xticks([])
        A.set_yticks([])
        if j == 0:
            A.legend(fontsize=7, loc="upper left")
    ax[0, 0].text(0.01, 0.01, "(a) yellow bar 100 A; dotted = +-30 deg mass sector", transform=ax[0, 0].transAxes,
                  color="w", fontsize=7)
    # (c) helical periodograms
    A = ax[0, 4]
    for grp, col in (("DOT", cD), ("FREE", cF)):
        for p in (D if grp == "DOT" else F):
            key = f"hel{p['cls']}"
            if key in cp.files:
                f, pw = cp[key]
                sel = (f > 1 / 100) & (f < 1 / 15)
                A.plot(1 / f[sel], pw[sel] / np.median(pw[sel]), color=col, alpha=0.35, lw=0.8)
    A.axvline(34, color="k", ls="--", lw=0.8)
    A.set_xlim(15, 60)
    A.set_ylim(0, 12)
    A.set_xlabel("period (A)")
    A.set_ylabel("power / median")
    A.set_title("(c) ridge-profile periodograms (red DOT, blue FREE);\n34 A = B-DNA helical repeat", fontsize=9)
    # (b) normalised mass profiles along the DNA
    A = ax[1, 0]
    for p in D:
        r, m, rd, sg = cp[f"dot{p['cls']}"]
        b = (r >= BODY_DOT[0]) & (r <= BODY_DOT[1])
        A.plot(r - R_DOT_A, m / np.median(m[b]), color=cD, alpha=0.25 + 0.6 * p["n"] / 14371, lw=1)
    for p in F:
        r, m, rd, sg, m2, rd2, sg2 = cp[f"free{p['cls']}"]
        b = (r >= BODY_FREE[0]) & (r <= BODY_FREE[1])
        L2 = p["arm_minus"]["L_end"]
        A.plot(r + L2, m / np.median(m[b]), color=cF, alpha=0.4, lw=1)
    A.axhline(0.5, color="k", ls=":", lw=0.7)
    A.axvline(170, color="g", ls="--", lw=0.8)
    A.text(171, 1.1, "50 bp B-DNA", color="g", fontsize=7, rotation=90)
    A.set_xlim(0, 200)
    A.set_ylim(-0.3, 2.0)
    A.set_xlabel("distance along the DNA (A): DOT r - 17; FREE from the far end")
    A.set_ylabel("transverse-integrated mass / body")
    A.set_title("(b) hinge-insensitive mass profiles (class averages)", fontsize=9)
    # (a) width profiles + calibration
    A = ax[1, 1]
    for p in D:
        r, m, rd, sg = cp[f"dot{p['cls']}"]
        ok = (r >= WFIT0_DOT) & (r <= (p["L_ridge"] or 100))
        A.plot(r[ok], sg[ok], ".", color=cD, ms=2, alpha=0.5)
    for p in F:
        r, m, rd, sg, m2, rd2, sg2 = cp[f"free{p['cls']}"]
        ok = (r >= WFIT0_FREE) & (r <= min(p["arm_plus"]["L_ridge"] or 50, 60))
        A.plot(r[ok], sg[ok], ".", color=cF, ms=2, alpha=0.5)
    for sg_, ls in ((0, "-"), (8, "--"), (16, ":"), (25, "-.")):
        key = f"cal{sg_}"
        if key in cp.files:
            r, m, rd, sg = cp[key]
            ok = (r >= WFIT0_DOT) & (r <= 125) & np.isfinite(sg)
            A.plot(r[ok], ndi.median_filter(sg[ok], 5), "k", ls=ls, lw=1, label=f"FREE avg pivoted, {sg_} deg")
    A.set_ylim(0, 25)
    A.set_xlabel("distance from origin (dot / FREE midpoint), A")
    A.set_ylabel("transverse Gaussian sigma (A)")
    A.legend(fontsize=7)
    A.set_title("(a) rod width vs distance (hinge blur widens linearly)", fontsize=9)
    # (a) length vs end sharpness
    A = ax[1, 2]
    for grp, col in (("DOT", cD), ("FREE", cF)):
        P_ = D if grp == "DOT" else F
        x = [p["L_end"] - R_DOT_A if grp == "DOT" else p["L_total"] for p in P_]
        y = [p["sig_end"] for p in P_]
        A.scatter(x, y, s=[8 + p["n"] / 60 for p in P_], color=col, alpha=0.6, label=grp)
        for p, xi, yi in zip(P_, x, y):
            if p["n"] > 2500:
                A.annotate(str(p["cls"]), (xi, yi), fontsize=7)
    A.set_xlabel("DNA length (A): DOT L_end - 17, FREE end-to-end")
    A.set_ylabel("end width sigma_end (A)")
    A.legend(fontsize=8)
    A.set_title("(a) hinge-insensitive length vs end sharpness (size ~ n)", fontsize=9)
    # (a) widening angle and calibration
    A = ax[1, 3]
    cal = pd.DataFrame(C["calibration"]["mean_by_sigma"])
    A.plot(cal.sigma_deg, -cal.dL_ridge, "k-o", ms=3, label="half-max under-read")
    A.plot(cal.sigma_deg, cal.dL_end, "g-s", ms=3, label="L_end change")
    A.plot(cal.sigma_deg, cal.sig_end, "m-^", ms=3, label="sigma_end")
    td = [p["theta_deg"] for p in D if p["theta_deg"] is not None]
    tf = [p["theta_deg"] for p in F if p["theta_deg"] is not None]
    A.plot(td, np.full(len(td), -2), "|", color=cD, ms=12, label="DOT classes: widening angle")
    A.plot(tf, np.full(len(tf), -4), "|", color=cF, ms=12, label="FREE classes")
    gap = out["lengths"]["hinge_insensitive_L_end"]["diff_median"]["est"]
    A.axhline(gap, color=cD, ls="--", lw=0.8)
    A.text(0.5, gap + 0.5, f"observed DOT shortfall {gap:.0f} A", color=cD, fontsize=7)
    A.set_xlabel("in-plane hinge sigma (deg)")
    A.set_ylabel("A")
    A.legend(fontsize=7, loc="upper left")
    A.set_title("(a) forward hinge blur of real FREE averages", fontsize=9)
    # (d) particle-weighted ECDF
    A = ax[1, 4]
    for grp, col in (("DOT", cD), ("FREE", cF)):
        P_ = D if grp == "DOT" else F
        x = np.array([p["L_end"] - R_DOT_A if grp == "DOT" else p["L_total"] for p in P_])
        w = np.array([p["n"] for p in P_], float)
        o = np.argsort(x)
        A.step(x[o], np.cumsum(w[o]) / w.sum(), where="post", color=col, label=grp)
        x2 = np.array([p["L_ridge"] - R_DOT_A if grp == "DOT" else p["L_ridge_total"] for p in P_])
        o = np.argsort(x2)
        A.step(x2[o], np.cumsum(w[o]) / w.sum(), where="post", color=col, ls=":", lw=0.8)
    L = out["lengths"]["hinge_insensitive_L_end"]
    A.set_title(f"(d) particle-weighted length ECDF (solid L_end, dotted half-max)\nmedian DOT "
                f"{L['dot_median']['est']:.0f} {L['dot_median']['ci95']}, FREE {L['free_median']['est']:.0f} "
                f"{L['free_median']['ci95']}; diff {L['diff_median']['est']:.0f} {L['diff_median']['ci95']} A",
                fontsize=8)
    A.axvline(170, color="g", ls="--", lw=0.8)
    A.set_xlabel("DNA length (A)")
    A.legend(fontsize=8)
    # row 2: re-alignment
    A = ax[2, 0]
    k1, k2 = 123, 130
    tiles = []
    for k in (k1, k2):
        row = []
        for key in ("plain", "realigned"):
            im = pp[f"{key}{k}"]
            ox, oy, ph = pp[f"geo{k}"]
            row.append(im[int(oy) - 48:int(oy) + 16, int(ox) - 32:int(ox) + 32])
        tiles.append(np.hstack([t / np.percentile(t, 99.5) for t in row]))
    A.imshow(np.vstack(tiles), cmap="gray", vmin=-0.3, vmax=1)
    A.set_xticks([])
    A.set_yticks([])
    A.set_title(f"(b) my re-averages: left plain, right rod re-aligned\n(top {k1}, bottom {k2}; dot at lower "
                f"centre)", fontsize=9)
    A = ax[2, 1]
    R = {int(k): v for k, v in Pj["classes"].items()}
    ks = list(R)
    for i, k in enumerate(ks):
        b, a_ = R[k]["before"]["L_end"], R[k]["after"]["L_end"]
        col = "0.5" if R[k]["extra_control"] else cD
        if b is not None and a_ is not None:
            A.plot([0, 1], [b - R_DOT_A, a_ - R_DOT_A], "-o", color=col, ms=3, alpha=0.7)
        b, a_ = R[k]["before"]["L_ridge"], R[k]["after"]["L_ridge"]
        if b is not None and a_ is not None:
            A.plot([2, 3], [b - R_DOT_A, a_ - R_DOT_A], "-o", color=col, ms=3, alpha=0.7)
    A.set_xticks([0, 1, 2, 3])
    A.set_xticklabels(["L_end\nplain", "L_end\nrealigned", "half-max\nplain", "half-max\nrealigned"], fontsize=8)
    A.set_ylabel("DNA length (A)")
    A.set_title("(b) re-alignment about the dot (grey = long-DNA controls 97, 87)", fontsize=9)
    A = ax[2, 2]
    for k in ks:
        p = pp[f"prof{k}"]
        r = p[0]
        b = (r >= BODY_DOT[0]) & (r <= BODY_DOT[1])
        col = "0.5" if R[k]["extra_control"] else cD
        A.plot(r - R_DOT_A, p[4] / np.median(p[4][b]), color=col, lw=0.8, alpha=0.7)
        A.plot(r - R_DOT_A, p[1] / np.median(p[1][b]), color=col, lw=0.8, alpha=0.7, ls=":")
    A.axvspan(SEGS["prox"][0] - R_DOT_A, SEGS["prox"][1] - R_DOT_A, color="y", alpha=0.15, label="reference")
    A.axvspan(FAR_A[0] - R_DOT_A, FAR_A[1] - R_DOT_A, color="g", alpha=0.08, label="far zone (bias-free)")
    A.axhline(0.5, color="k", ls=":", lw=0.7)
    A.set_xlim(20, 190)
    A.set_ylim(-0.6, 1.8)
    A.legend(fontsize=7)
    A.set_xlabel("r - 17 (A)")
    A.set_ylabel("mass / body (solid realigned, dotted plain)")
    A.set_title("(b) mass profiles of my re-averages before / after", fontsize=9)
    A = ax[2, 3]
    for i, k in enumerate(ks):
        ct = R[k]["cov_test"]
        col = "0.5" if R[k]["extra_control"] else cD
        A.errorbar(i, ct["diff"], yerr=[[ct["diff"] - ct["diff_ci95"][0]], [ct["diff_ci95"][1] - ct["diff"]]],
                   fmt="o", color=col, ms=4)
        A.plot(i, -ct["cov_S1_S2"], "kx", ms=5)
    A.axhline(0, color="k", lw=0.7)
    A.set_xticks(range(len(ks)))
    A.set_xticklabels([str(k) for k in ks], fontsize=7, rotation=90)
    A.set_ylabel("cov(prox, S3) - cov(S1, S2)  (deg^2)")
    A.set_title("(b) far-segment angle test: 0 = rod present at 100-170 A and rigid;\n"
                "x = expected value if NO rod at 100-170 A (= -cov(S1,S2))", fontsize=8)
    A = ax[2, 4]
    for k in [kk for kk in ks if not R[kk]["extra_control"]][:6]:
        ang = pp[f"angles{k}"]
        A.hist(ang[0], bins=np.arange(-30.5, 31, 1), histtype="step", density=True, label=f"{k}")
    A.set_xlabel("proximal-rod angle about the dot (deg)")
    A.set_ylabel("density")
    A.legend(fontsize=7, ncol=2)
    sh = out["realign"]["sigma_hinge_deg_weighted_mean"]
    A.set_title(f"(b) per-particle rod angles (noise-dominated);\nnoise-free spread sqrt(cov(S1,S2)) ~ {sh} deg "
                f"(weighted mean)", fontsize=9)
    dec = out["decision"]
    fig.suptitle(f"A1 attached-DNA length vs hinge blur.  DOT {len(D)} classes / FREE {len(F)} classes.  Pre-set "
                 f"rule outcome: {dec['outcome']}  (median gap {dec['gap_median_A']['est']} A "
                 f"{dec['gap_median_A']['ci95']}; re-align dL_end {dec['realign_dL_end']} A; S3-positive classes "
                 f"{dec['n_classes_S3_positive']})", fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.97])
    os.makedirs(os.path.dirname(FIG), exist_ok=True)
    fig.savefig(FIG, dpi=110)
    print("wrote", FIG, flush=True)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("stage", choices=["classes", "particles", "report", "all"])
    ap.add_argument("--nproc", type=int, default=8)
    ap.add_argument("--nconv", type=int, default=NCONV)
    ap.add_argument("--maxper", type=int, default=0, help="debug: cap particles per class (0 = all)")
    a = ap.parse_args()
    if a.stage in ("classes", "all"):
        stage_classes(a)
    if a.stage in ("particles", "all"):
        stage_particles(a)
    if a.stage in ("report", "all"):
        stage_report(a)
