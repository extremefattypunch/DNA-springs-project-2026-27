"""Are sfGFP dots attached to DNA-rod ENDS more often than chance?  (CPU only; real coordinates only.)

The round-1 rod 2D (W2_2D_rods, J27) shows many "dot at one end of a rod" averages, but the ice holds ~4,500 sfGFP
dots per um^2, so a dot lands within a few nm of any point by chance and chance neighbours also make sharp class
averages. This script tests end-colocalisation against nulls built from the SAME real rod and dot coordinates.

Stages (each run as a Slurm CPU job on `test`; see the bottom of this docstring):

  conventions   class-average rod geometry + empirical pose/shift and micrograph y-axis conventions + QC figure
  detect        G0 dot detector on every micrograph holding a rod-set particle (job array; --task/--ntasks)
  analyse       rod mapping, dedup, statistic T, nulls, angular / along-axis profiles, strata, outputs

CHOICES (all fixed before looking at the T results)
1. Rod set R = every particle in a class whose consensus label (results/class2d_labels_W2_2D_rods.csv) is
   rod_no_dot, dot_long_rod or dot_short_rod -- the union, so membership does not depend on whether a dot is present.
   Excluded: dot_pair_thread, dot_pair, ice_contaminant, noise_junk, unclear.
2. Rod geometry per class, measured on the 64 px / 5.68 A class average: robust z map (median / MAD of the 8-19 px
   ring, i.e. inside CryoSPARC's soft mask, where the rod covers ~12% of the pixels); 8-connected component of z > 4
   that contains the brightest pixel within 12 px of the box centre; rod axis = principal axis of that binary
   component (unweighted, so a bright dot cannot dominate it); along-axis profile = mean z over a 3 px strip
   (|t| <= 1 px); "rod level" = median profile inside the component; the rod's two ends = the half-max points
   (profile = 0.5 * rod level) walking outward from the profile maximum. For dot_long_rod / dot_short_rod classes
   the averaged dot is the profile maximum at s_dot (always at one end in these averages); +axis is oriented toward
   it, and the DNA end on that side is put at s_dot - 17 A (one sfGFP radius, so the averaged dot sits at
   L/2 + 17 A: the place a covalently attached sfGFP should sit). For rod_no_dot classes both ends are the half-max
   points and the +axis sign is arbitrary. Rod centre = midpoint of the two ends, L = their separation. The
   dot_short_rod averages are blobby, so their axes are the least reliable (10k particles, ~4% of R).
3. Pose/shift convention, determined EMPIRICALLY and without dots: 3,000 random R particles; each extracted 128 px
   particle (2.84 A/px) is phase-flipped with its own CTF and low-passed (Gaussian, half-max 1/15 A^-1); for each of
   16 candidate conventions -- rotate by +theta or -theta; add or subtract the shift; shift applied in the particle
   frame (p = R c +/- s) or in the class frame (p = R (c +/- s)); shift components (x, y) or swapped (y, x) -- the
   particle is resampled into the class frame and correlated (Pearson, r < 56 px) with its class average resampled
   to 2.84 A/px (and low-passed the same way). The candidate with the clearly highest mean correlation is used.
   Micrograph orientation, same idea: 400 particles on 40 micrographs; the 15 A low-passed micrograph is binned 2x
   (2.84 A/px), a 128 px patch is cut at (x_frac*W or (1-x_frac)*W, y_frac*H or (1-y_frac)*H), optionally flipped
   up-down / left-right (16 options), and correlated with the (unflipped, low-passed) particle image. The option
   with the largest |mean r| wins; its SIGN is the contrast of the stored particle images relative to the
   micrograph (found to be inverted: r < 0) and is applied before the pose/shift test above, which is run after
   this one so that the phase-flipped particles are protein-positive like the class averages.
   The QC figure (figures/qc/rod_end_coloc_conventions.png) draws the mapped rod axes on real micrograph crops.
4. Dots: the validated G0 detector of scripts/12_dot_spatial_stats.py (imported, unmodified): own-CTF phase flip,
   15 A low-pass, -sigma^2 LoG (sigma 12 A), threshold median + 4.5 * 1.4826 * MAD, roundness >= 0.2, 20 A NMS,
   52 px edge + contamination mask (bin-4 valid mask kept for the nulls). CTF from W2_curated exposures.
   A strict dot set (LoG robust z >= 6) is carried along as a sensitivity check.
5. Rods are kept only if their centre and both points at L/2 + 15 A lie in the micrograph's valid mask (dots
   cannot be detected outside it). Duplicates (the same physical rod picked twice, e.g. once on the dot and once on
   the rod middle; picks were only 80 A apart) are removed: rods on the same micrograph whose centres are within
   40 A and axes within 30 deg are merged, keeping one at random (seeded; independent of dots).
6. Statistic T(d, end offset) = fraction of rods with >= 1 ELIGIBLE dot centre within d of either end point, where an
   end point is the centre +/- (L/2 + offset) along the axis, d in {25, 35, 50} A, offset in {0, 15} A (15 A allows
   for the sfGFP radius; offset 15 A is the primary readout). ELIGIBLE = not inside the rod's middle region
   (|t| <= 30 A from the axis line AND |s| <= L/2 - 10 A, s along the axis from the centre) -- the rod's own density
   can make LoG peaks there -- and not inside the middle region of ANY real R rod on that micrograph (removed once,
   globally, so that the nulls cannot pick up rod-body false peaks either). The same rule is applied to the virtual
   rods of the nulls. T_both = fraction with an eligible dot at BOTH ends.
   For dot_* classes the dot-side end (+axis) is the end the classifier selected on; the FAR end (-axis) is the
   selection-free test within those classes, reported separately ("near" / "far").
7. Nulls (100 replicates each, real coordinates only):
   (a) rotation: each rod is rotated about its own centre by a uniform random angle;
   (b) shift: each rod keeps its axis and is moved 150-300 A (uniform) in a uniform random direction, same micrograph.
   A proposal whose centre or L/2 + 15 A end points fall outside the valid mask is redrawn (up to 50 times; the last
   proposal is kept otherwise, and counted). Ratio = T_obs / mean(T_null); one-sided p = (1 + #{T_null >= T_obs}) /
   (1 + n_rep) (min 0.0099 with 100 replicates) plus z = (T_obs - mean) / SD.
8. Angular profile: dots (eligible set) within 120 A of a rod centre; angle phi between the centre->dot vector and the
   +axis (0 = along +axis = dot side for dot_* classes, 180 = the other end), 10 deg bins; also restricted to the
   "end shell" r - L/2 in [10, 45] A, where the rod body excludes no dot positions at any angle. Expected = null (b)
   (same rods, same micrographs, no rod at the virtual position) mean and 95% band. Along-axis profile: ALL dots
   (no exclusion, so rod-end false peaks are visible) with |t| <= 25 A, histogram of |s| - L/2 from -100 to +120 A
   (5 A bins), observed vs null (b). Covalent end attachment predicts excess at |s| - L/2 ~ +15-25 A and at
   phi ~ 0 / 180 deg; chance predicts the null shape; a detector artefact at rod ends predicts excess at ~ -10 A.
   False-peak check: dots ON rods (|t| <= 15 A and |s| <= L/2), observed vs null (b).
9. Strata: consensus label (rod_no_dot / dot_long_rod / dot_short_rod) and, for dot_* classes, near vs far end.

Run (logs to /n/netscratch/hekstra_lab/Lab/ian_poon/logs/):
  PY=/n/holylabs/hekstra_lab/Everyone/ianpoon/envs/dnaspring/bin/python
  sbatch -p test -c 8  --mem=32G  -t 1:00:00 --wrap "$PY -I scripts/17_rod_end_colocalisation.py conventions"
  sbatch -p test -c 24 --mem=64G  -t 3:00:00 --array=0-3 \
      --wrap "$PY -I scripts/17_rod_end_colocalisation.py detect --task \$SLURM_ARRAY_TASK_ID --ntasks 4"
  sbatch -p test -c 16 --mem=96G  -t 4:00:00 --wrap "$PY -I scripts/17_rod_end_colocalisation.py analyse"

Outputs: results/rod_end_coloc.json, figures/qc/rod_end_coloc.png, figures/qc/rod_end_coloc_conventions.png;
intermediates (dots, masks, convention table) under WORK below.
"""
import argparse
import importlib.util
import itertools
import json
import os
import time
from multiprocessing import Pool

import numpy as np
import mrcfile
import pandas as pd
from scipy import ndimage as ndi
from scipy.spatial import cKDTree

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
_spec = importlib.util.spec_from_file_location("g0det", os.path.join(REPO, "scripts", "12_dot_spatial_stats.py"))
g0 = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(g0)

PROJ = "/n/netscratch/hekstra_lab/Lab/ian_poon/cryosparc_projects/CS-cs-sfgfp-chimera-screening"
EXPORTS = "/n/netscratch/hekstra_lab/Lab/ian_poon/sfgfp_chimera_screening_20261112/exports"
WORK = "/n/netscratch/hekstra_lab/Lab/ian_poon/sfgfp_chimera_screening_20261112/rod_end_coloc"
PART_NPZ = os.path.join(EXPORTS, "W2_2D_rods__particles.npz")
CA_MRC = os.path.join(EXPORTS, "W2_2D_rods", "J27_030_class_averages.mrc")
CUR_NPZ = os.path.join(EXPORTS, "W2_curated__exposures.npz")
LABELS = os.path.join(REPO, "results", "class2d_labels_W2_2D_rods.csv")
OUT_JSON = os.path.join(REPO, "results", "rod_end_coloc.json")
FIG = os.path.join(REPO, "figures", "qc", "rod_end_coloc.png")
FIG_CONV = os.path.join(REPO, "figures", "qc", "rod_end_coloc_conventions.png")

A_MIC, A_P, A_CA = 1.42, 2.84, 5.68
BOX_P, BOX_CA, MIC_N = 128, 64, 2048
R_LABELS = ["rod_no_dot", "dot_long_rod", "dot_short_rod"]
DOT_LABELS = ["dot_long_rod", "dot_short_rod"]
GEOM_ZT = 4.0
R_DOT_A = 17.0
D_LIST = [25.0, 35.0, 50.0]
ENDOFF = [0.0, 15.0]
MID_T_A, MID_INSET_A = 30.0, 10.0
ONROD_T_A = 15.0
STRICT_Z = 6.0
SHIFT_A = (150.0, 300.0)
ANG_R_A = 120.0
SHELL_A = (10.0, 45.0)
AX_EDGES = np.arange(-100.0, 120.0 + 1e-6, 5.0)
ANG_EDGES = np.arange(0.0, 180.0 + 1e-6, 10.0)
MAP_T_EDGES = np.arange(0.0, 120.0 + 1e-6, 5.0)
DEDUP_R_A, DEDUP_ANG = 40.0, 30.0
OFF = 4096.0          # x offset per micrograph index in the global dot tree (keeps micrographs apart)
NREP = 100
SEED = 20261114


# ----------------------------------------------------------------------------------------------- common helpers
def lowpass_gauss(k2):
    s = (1.0 / g0.LOWPASS_A) / np.sqrt(2 * np.log(2))
    G = np.exp(-k2 / (2 * s * s))
    G.flat[0] = 0.0
    return G


def flip_lp(img, df1, df2, ang, apix, flip=True):
    a = np.asarray(img, np.float64)
    a = a - a.mean()
    a /= a.std() + 1e-12
    F = np.fft.rfft2(a)
    T, k2 = g0.ctf_2d(a.shape, df1, df2, ang, apix=apix)
    f = np.sign(T) if flip else 1.0
    return np.fft.irfft2(F * f * lowpass_gauss(k2), s=a.shape)


def lp_only(img, apix):
    a = np.asarray(img, np.float64)
    a = a - a.mean()
    a /= a.std() + 1e-12
    ky = np.fft.fftfreq(a.shape[0], d=apix)[:, None]
    kx = np.fft.rfftfreq(a.shape[1], d=apix)[None, :]
    return np.fft.irfft2(np.fft.rfft2(a) * lowpass_gauss(kx ** 2 + ky ** 2), s=a.shape)


def class_labels():
    lab = pd.read_csv(LABELS)
    arr = np.full(512, "none", dtype=object)
    arr[lab["class"].values] = lab["consensus_label"].values
    return arr


def class_geometry(img, has_dot):
    """See docstring item 2. Returns class-frame geometry in 2.84 A particle pixels relative to the box centre."""
    n = img.shape[0]
    c0 = n / 2.0
    yy, xx = np.indices(img.shape)
    rr = np.hypot(xx - c0, yy - c0)
    bg = img[(rr > 8) & (rr < 19)]
    med = np.median(bg)
    mad = 1.4826 * np.median(np.abs(bg - med))
    z = (img - med) / mad
    lb, _ = ndi.label(z > GEOM_ZT, structure=np.ones((3, 3)))
    pk = np.unravel_index(np.argmax(np.where(rr < 12, z, -np.inf)), z.shape)
    comp = lb == lb[pk]
    ys, xs = np.nonzero(comp)
    mx, my = xs.mean(), ys.mean()
    w, v = np.linalg.eigh(np.cov(np.vstack([xs - mx, ys - my])))
    u = v[:, 1].copy()
    s = np.arange(-40, 40.001, 0.25)

    def profile(uu):
        p = np.zeros_like(s)
        for t in (-1, 0, 1):
            p += ndi.map_coordinates(z, [my + s * uu[1] + t * uu[0], mx + s * uu[0] - t * uu[1]], order=1,
                                     mode="constant") / 3
        return p
    prof = profile(u)
    inside = prof > GEOM_ZT
    rodlev = float(np.median(prof[inside]))
    im = int(np.argmax(np.where(inside, prof, -np.inf)))
    i = im
    while i < len(s) - 1 and prof[i] >= 0.5 * rodlev:
        i += 1
    shi = s[i]
    i = im
    while i > 0 and prof[i] >= 0.5 * rodlev:
        i -= 1
    slo = s[i]
    sdot = None
    if has_dot:
        sdot = s[im]
        if sdot - slo < shi - sdot:          # orient +axis toward the averaged dot
            u, sdot, slo, shi = -u, -sdot, -shi, -slo
            prof = prof[::-1]
        e_plus = sdot - R_DOT_A / A_CA       # DNA end on the dot side
    else:
        e_plus = shi
    e_minus = slo
    cs = 0.5 * (e_plus + e_minus)
    # convert to 2.84 A px relative to the box centre (64 px class box centre = 32 -> particle frame x2)
    f = A_CA / A_P
    cen = np.array([(mx - c0 + cs * u[0]) * f, (my - c0 + cs * u[1]) * f])
    # where the 2D box centre (= pick / alignment centre) sits in rod coordinates: s along +axis, t across
    s_box = float(-(cen[0] * u[0] + cen[1] * u[1]) * A_P)
    t_box = float(abs(-cen[0] * u[1] + cen[1] * u[0]) * A_P)
    return dict(cx=float(cen[0]), cy=float(cen[1]), ux=float(u[0]), uy=float(u[1]),
                box_centre_s_A=s_box, box_centre_t_A=t_box,
                h_A=float(0.5 * (e_plus - e_minus) * A_CA), L_halfmax_A=float((shi - slo) * A_CA),
                sdot_from_centre_A=None if sdot is None else float((sdot - cs) * A_CA),
                rod_level_z=rodlev, dot_ratio=float(prof.max() / rodlev), has_dot=bool(has_dot),
                width_fwhm_A=float(np.sqrt(max(w[0], 0)) * 2.355 * A_CA),
                comp=comp, mx=mx, my=my, slo=slo, shi=shi, sdot=sdot, e_plus=e_plus, e_minus=e_minus, u=u)


CANDS = [dict(rot=rs, shift=ss, order=o, swap=sw) for rs in (1, -1) for ss in (1, -1) for o in ("particle", "class")
         for sw in (False, True)]


def cand_name(c):
    return (f"rot{'+' if c['rot'] > 0 else '-'}th_shift{'+' if c['shift'] > 0 else '-'}_{c['order']}frame"
            f"{'_swapxy' if c['swap'] else ''}")


def class_to_particle(cx, cy, th, sx, sy, c):
    """class-frame point(s) (2.84 A px from box centre, x = column, y = row) -> particle-image frame."""
    if c["swap"]:
        sx, sy = sy, sx
    co, si = np.cos(c["rot"] * th), np.sin(c["rot"] * th)
    if c["order"] == "class":
        ax, ay = cx + c["shift"] * sx, cy + c["shift"] * sy
        return co * ax - si * ay, si * ax + co * ay
    return co * cx - si * cy + c["shift"] * sx, si * cx + co * cy + c["shift"] * sy


def load_particles():
    d = np.load(PART_NPZ)
    keys = ["alignments2D__class", "alignments2D__pose", "alignments2D__shift", "blob__path", "blob__idx",
            "location__micrograph_path", "location__center_x_frac", "location__center_y_frac",
            "location__micrograph_shape", "ctf__df1_A", "ctf__df2_A", "ctf__df_angle_rad",
            "alignments2D__class_posterior"]
    return {k: d[k] for k in keys}


def all_geometries():
    ca = mrcfile.open(CA_MRC, permissive=True).data.astype(np.float64)
    labs = class_labels()
    geo = {}
    for k in range(ca.shape[0]):
        if labs[k] in R_LABELS:
            geo[k] = class_geometry(ca[k], labs[k] in DOT_LABELS)
            geo[k]["label"] = labs[k]
    return ca, labs, geo


def map_rods(P, idx, geo, conv, ycfg):
    """Map class geometry of particles idx into micrograph px (1.42 A). Returns X, Y, ux, uy, h_px."""
    cls = P["alignments2D__class"][idx].astype(int)
    th = P["alignments2D__pose"][idx].astype(np.float64)
    sh = P["alignments2D__shift"][idx].astype(np.float64)
    gcx = np.array([geo[c]["cx"] for c in cls])
    gcy = np.array([geo[c]["cy"] for c in cls])
    gux = np.array([geo[c]["ux"] for c in cls])
    guy = np.array([geo[c]["uy"] for c in cls])
    h = np.array([geo[c]["h_A"] for c in cls]) / A_MIC
    px, py = class_to_particle(gcx, gcy, th, sh[:, 0], sh[:, 1], conv)
    qx, qy = class_to_particle(gcx + gux, gcy + guy, th, sh[:, 0], sh[:, 1], conv)
    vx, vy = qx - px, qy - py
    W = P["location__micrograph_shape"][idx, 1].astype(float)
    H = P["location__micrograph_shape"][idx, 0].astype(float)
    xf = P["location__center_x_frac"][idx].astype(np.float64)
    yf = P["location__center_y_frac"][idx].astype(np.float64)
    colc = xf * W if ycfg["x"] == "direct" else (1 - xf) * W
    rowc = yf * H if ycfg["y"] == "direct" else (1 - yf) * H
    fx = -1.0 if ycfg["fliplr"] else 1.0
    fy = -1.0 if ycfg["flipud"] else 1.0
    r = A_P / A_MIC
    X, Y = colc + r * fx * px, rowc + r * fy * py
    ux, uy = fx * vx, fy * vy
    nrm = np.hypot(ux, uy)
    return X, Y, ux / nrm, uy / nrm, h


# ----------------------------------------------------------------------------------------------- stage 1
def stage_conventions(a):
    t0 = time.time()
    os.makedirs(WORK, exist_ok=True)
    rng = np.random.default_rng(SEED)
    ca, labs, geo = all_geometries()
    P = load_particles()
    inR = np.isin(labs[P["alignments2D__class"].astype(int)], R_LABELS)
    ridx = np.nonzero(inR)[0]
    sub = np.sort(rng.choice(ridx, size=a.nconv, replace=False))
    jj = np.arange(BOX_P) - BOX_P / 2
    YY, XX = np.meshgrid(jj, jj, indexing="ij")
    circ = np.hypot(XX, YY) < 56
    # micrograph orientation: particle image vs patch cut from the binned, low-passed micrograph
    cur = np.load(CUR_NPZ)
    mpaths = np.unique(P["location__micrograph_path"][ridx])
    mic_sel = rng.choice(mpaths, size=a.nymic, replace=False)
    YOPTS = [dict(x=xo, y=yo, flipud=fu, fliplr=fl) for xo in ("direct", "inverted") for yo in ("direct", "inverted")
             for fu in (False, True) for fl in (False, True)]
    ycorr = []
    for mp in mic_sel:
        with mrcfile.open(os.path.join(PROJ, mp), permissive=True) as m:
            mic = np.asarray(m.data, np.float64)
        b2 = mic.reshape(MIC_N // 2, 2, MIC_N // 2, 2).mean(axis=(1, 3))
        b2 = lp_only(b2, A_P)
        on = ridx[P["location__micrograph_path"][ridx] == mp]
        on = rng.choice(on, size=min(len(on), a.nyper), replace=False)
        for i in on:
            with mrcfile.mmap(os.path.join(PROJ, str(P["blob__path"][i])), permissive=True, mode="r") as m:
                pimg = lp_only(np.array(m.data[int(P["blob__idx"][i])], np.float64), A_P)
            xf, yf = float(P["location__center_x_frac"][i]), float(P["location__center_y_frac"][i])
            row = []
            for o in YOPTS:
                cc = (xf if o["x"] == "direct" else 1 - xf) * MIC_N / 2
                rc = (yf if o["y"] == "direct" else 1 - yf) * MIC_N / 2
                x0, y0 = int(round(cc)) - BOX_P // 2, int(round(rc)) - BOX_P // 2
                if x0 < 0 or y0 < 0 or x0 + BOX_P > MIC_N // 2 or y0 + BOX_P > MIC_N // 2:
                    row.append(np.nan)
                    continue
                patch = b2[y0:y0 + BOX_P, x0:x0 + BOX_P]
                if o["flipud"]:
                    patch = patch[::-1]
                if o["fliplr"]:
                    patch = patch[:, ::-1]
                row.append(np.corrcoef(patch[circ], pimg[circ])[0, 1])
            ycorr.append(row)
    ycorr = np.array(ycorr)
    ymean = np.nanmean(ycorr, 0)
    ybest = int(np.nanargmax(np.abs(ymean)))         # contrast of the stored particles may be inverted
    psign = float(np.sign(ymean[ybest]))             # particle = psign * micrograph (low-pass, same location)
    ytab = [dict(**YOPTS[k], mean_corr=float(ymean[k]), n=int(np.isfinite(ycorr[:, k]).sum()),
                 sem=float(np.nanstd(ycorr[:, k]) / np.sqrt(max(np.isfinite(ycorr[:, k]).sum(), 1))))
            for k in range(len(YOPTS))]
    print("micrograph orientation:", YOPTS[ybest], f"{ymean[ybest]:.3f}", "next best |r|",
          sorted(np.abs(ymean[np.isfinite(ymean)]))[-2], f"particle contrast sign {psign:+.0f}", flush=True)

    # class averages on the 2.84 A grid, low-passed like the particles
    jj = np.arange(BOX_P) - BOX_P / 2
    YY, XX = np.meshgrid(jj, jj, indexing="ij")
    ca128 = {}
    for k in geo:
        up = ndi.map_coordinates(ca[k], [BOX_CA / 2 + YY / 2, BOX_CA / 2 + XX / 2], order=3, mode="constant")
        ca128[k] = lp_only(up, A_P)
    circ = np.hypot(XX, YY) < 56
    corr = np.zeros((len(sub), len(CANDS) + 1))
    t_read = time.time()
    by_file = {}
    for j, i in enumerate(sub):
        by_file.setdefault(str(P["blob__path"][i]), []).append((j, i))
    for f, lst in by_file.items():
        with mrcfile.mmap(os.path.join(PROJ, f), permissive=True, mode="r") as m:
            for j, i in lst:
                img = np.array(m.data[int(P["blob__idx"][i])], np.float64)
                pimg = psign * flip_lp(img, P["ctf__df1_A"][i], P["ctf__df2_A"][i], P["ctf__df_angle_rad"][i],
                                       A_P)                      # protein-positive, like the class average
                ref = ca128[int(P["alignments2D__class"][i])][circ]
                th = float(P["alignments2D__pose"][i])
                sx, sy = (float(v) for v in P["alignments2D__shift"][i])
                for c, cand in enumerate(CANDS):
                    px, py = class_to_particle(XX, YY, th, sx, sy, cand)
                    smp = ndi.map_coordinates(pimg, [BOX_P / 2 + py, BOX_P / 2 + px], order=1, mode="constant")
                    corr[j, c] = np.corrcoef(smp[circ], ref)[0, 1]
                corr[j, -1] = np.corrcoef(pimg[circ], ref)[0, 1]       # identity (no alignment) baseline
    print(f"convention correlations on {len(sub)} particles from {len(by_file)} files ({time.time() - t_read:.0f} s)",
          flush=True)
    mean = corr.mean(0)
    sem = corr.std(0) / np.sqrt(len(sub))
    best = int(np.argmax(mean[:-1]))
    order = np.argsort(-mean[:-1])
    conv_tab = [dict(name=cand_name(CANDS[c]), **CANDS[c], mean_corr=float(mean[c]), sem=float(sem[c]),
                     frac_particles_best=float(np.mean(np.argmax(corr[:, :-1], 1) == c))) for c in range(len(CANDS))]
    conv_tab.append(dict(name="identity_no_alignment", mean_corr=float(mean[-1]), sem=float(sem[-1])))
    for c in order[:4]:
        print(f"  {cand_name(CANDS[c]):45s} {mean[c]:.4f} +/- {sem[c]:.4f}", flush=True)


    geo_out = {int(k): {kk: vv for kk, vv in g.items() if kk in (
        "label", "cx", "cy", "ux", "uy", "h_A", "L_halfmax_A", "sdot_from_centre_A", "rod_level_z", "dot_ratio",
        "box_centre_s_A", "box_centre_t_A",
        "has_dot", "width_fwhm_A")} for k, g in geo.items()}
    res = dict(created=time.strftime("%Y-%m-%d %H:%M"), n_particles=int(len(sub)), n_files=len(by_file),
               candidates=conv_tab, best=conv_tab[best], best_index=best,
               second_best=conv_tab[int(order[1])], micrograph_orientation=ytab, micrograph_best=ytab[ybest],
               particle_contrast_sign_vs_micrograph=psign,
               n_orientation_particles=int(len(ycorr)), n_orientation_micrographs=int(len(mic_sel)),
               class_geometry=geo_out)
    with open(os.path.join(WORK, "conventions.json"), "w") as fh:
        json.dump(res, fh, indent=1)
    fig_conventions(ca, geo, P, ridx, cur, CANDS[best], CANDS[int(order[1])], YOPTS[ybest], conv_tab, ytab, rng)
    print(f"conventions done in {time.time() - t0:.0f} s", flush=True)


def mic_ctf(cur, path):
    k = int(np.nonzero(cur["micrograph_blob__path"] == path)[0][0])
    return float(cur["ctf__df1_A"][k]), float(cur["ctf__df2_A"][k]), float(cur["ctf__df_angle_rad"][k])


def fig_conventions(ca, geo, P, ridx, cur, conv, conv2, ycfg, conv_tab, ytab, rng):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    ks = sorted(geo, key=lambda k: (R_LABELS.index(geo[k]["label"]), -np.sum(P["alignments2D__class"] == k)))
    fig = plt.figure(figsize=(22, 17), dpi=90)
    sf = fig.subfigures(3, 1, height_ratios=[1.0, 1.25, 1.6])
    ax = sf[0].subplots(1, 2, gridspec_kw=dict(width_ratios=[1.3, 1]))
    m = [c["mean_corr"] for c in conv_tab]
    e = [c.get("sem", 0) for c in conv_tab]
    ax[0].barh(range(len(m)), m, xerr=e, color=["C2" if v == max(m[:-1]) else "0.6" for v in m[:-1]] + ["C3"])
    ax[0].set_yticks(range(len(m)))
    ax[0].set_yticklabels([c["name"] for c in conv_tab], fontsize=7)
    ax[0].invert_yaxis()
    ax[0].set_xlabel("mean Pearson r, particle -> class frame vs class average")
    ax[0].set_title("pose / shift convention (3,000 rod particles, no dots involved)", fontsize=10)
    ym = [y["mean_corr"] for y in ytab]
    ax[1].barh(range(len(ym)), ym, color=["C2" if abs(v) == np.nanmax(np.abs(ym)) else "0.6" for v in ym])
    ax[1].set_yticks(range(len(ym)))
    ax[1].set_yticklabels([f"x {y['x']}, y {y['y']}, ud {int(y['flipud'])}, lr {int(y['fliplr'])}" for y in ytab],
                          fontsize=7)
    ax[1].invert_yaxis()
    ax[1].set_xlabel("mean Pearson r, particle image vs micrograph patch")
    ax[1].set_title("micrograph orientation of location fractions", fontsize=10)
    # class geometry montage
    n = len(ks)
    axs = sf[1].subplots(3, int(np.ceil(n / 3)))
    for a_ in axs.ravel():
        a_.axis("off")
    for a_, k in zip(axs.ravel(), ks):
        g = geo[k]
        a_.imshow(ca[k], cmap="gray")
        a_.contour(g["comp"], [0.5], colors="y", linewidths=0.4)
        for sv, c in ((g["e_minus"], "c"), (g["e_plus"], "m")):
            a_.plot(g["mx"] + sv * g["u"][0], g["my"] + sv * g["u"][1], "o", mfc="none", mec=c, ms=4)
        a_.plot([g["mx"] + g["e_minus"] * g["u"][0], g["mx"] + g["e_plus"] * g["u"][0]],
                [g["my"] + g["e_minus"] * g["u"][1], g["my"] + g["e_plus"] * g["u"][1]], "g-", lw=0.6)
        if g["sdot"] is not None:
            a_.plot(g["mx"] + g["sdot"] * g["u"][0], g["my"] + g["sdot"] * g["u"][1], "r+", ms=5)
        a_.set_xlim(12, 52)
        a_.set_ylim(52, 12)
        a_.set_title(f"{k} {g['label'].replace('_rod', '')}\nL {2 * g['h_A']:.0f} A", fontsize=6)
    sf[1].suptitle("class-average rod geometry: yellow = z>4 component, green = rod axis between DNA ends (cyan = -axis "
                   "end, magenta = +axis end), red + = averaged dot (DNA end placed 17 A inside it)", fontsize=9)
    # overlay on micrographs
    axs = sf[2].subplots(1, 3)
    mp_all = P["location__micrograph_path"][ridx]
    for a_, mp in zip(axs, rng.choice(np.unique(mp_all), 3, replace=False)):
        df1, df2, an = mic_ctf(cur, mp)
        pos, _ = g0.load_and_filter(os.path.join(PROJ, mp), df1, df2, an)
        on = ridx[mp_all == mp]
        X, Y, ux, uy, h = map_rods(P, on, geo, conv, ycfg)
        x0 = int(np.clip(np.median(X) - 300, 0, MIC_N - 600))
        y0 = int(np.clip(np.median(Y) - 300, 0, MIC_N - 600))
        show = -pos[y0:y0 + 600, x0:x0 + 600]
        lo, hi = np.percentile(show, [1, 99])
        a_.imshow(show, cmap="gray", vmin=lo, vmax=hi, extent=(x0 - .5, x0 + 599.5, y0 + 599.5, y0 - .5))
        alt = dict(ycfg, y="inverted" if ycfg["y"] == "direct" else "direct")
        for cfg, cv, col, ls in ((ycfg, conv2, "orange", ":"), (alt, conv, "red", "--"), (ycfg, conv, "lime", "-")):
            X2, Y2, ux2, uy2, h2 = map_rods(P, on, geo, cv, cfg)
            for xx, yy, vx, vy, hh in zip(X2, Y2, ux2, uy2, h2):
                a_.plot([xx - hh * vx, xx + hh * vx], [yy - hh * vy, yy + hh * vy], color=col, ls=ls, lw=1.2)
        a_.set_xlim(x0, x0 + 600)
        a_.set_ylim(y0 + 600, y0)
        a_.set_title(os.path.basename(mp)[22:52] + "...", fontsize=8)
    sf[2].suptitle(f"rod axes on phase-flipped 15 A micrograph crops (protein dark, 852 A wide). lime = chosen "
                   f"({cand_name(conv)}; y {ycfg['y']}); orange dotted = 2nd-best pose convention; red dashed = "
                   f"chosen pose with y inverted", fontsize=9)
    fig.savefig(FIG_CONV)
    plt.close(fig)


# ----------------------------------------------------------------------------------------------- stage 2
def detect_one(args):
    j, path, df1, df2, ang = args
    pos, _ = g0.load_and_filter(path, df1, df2, ang)
    m4, frac_bad = g0.valid_mask(pos)
    mfull = np.repeat(np.repeat(m4, g0.BIN, axis=0), g0.BIN, axis=1)
    x, y, zz, resp, n_mirror, _, rnd = g0.detect_dots(pos, mfull, 4.5)
    return j, x.astype(np.int16), y.astype(np.int16), zz.astype(np.float32), rnd.astype(np.float32), \
        np.packbits(m4, axis=1), float(frac_bad), int(n_mirror)


def mic_list():
    P = np.load(PART_NPZ)
    labs = class_labels()
    inR = np.isin(labs[P["alignments2D__class"].astype(int)], R_LABELS)
    return np.unique(P["location__micrograph_path"][inR])


def stage_detect(a):
    t0 = time.time()
    os.makedirs(os.path.join(WORK, "dots"), exist_ok=True)
    mics = mic_list()
    mine = np.array_split(np.arange(len(mics)), a.ntasks)[a.task]
    cur = np.load(CUR_NPZ)
    cpath = cur["micrograph_blob__path"]
    lut = {p: k for k, p in enumerate(cpath)}
    jobs = []
    for j in mine:
        k = lut[mics[j]]
        jobs.append((int(j), os.path.join(PROJ, mics[j]), float(cur["ctf__df1_A"][k]), float(cur["ctf__df2_A"][k]),
                     float(cur["ctf__df_angle_rad"][k])))
    procs = int(os.environ.get("SLURM_CPUS_PER_TASK", 4))
    print(f"task {a.task}/{a.ntasks}: {len(jobs)} micrographs, {procs} procs", flush=True)
    out = []
    with Pool(procs) as pool:
        for r in pool.imap_unordered(detect_one, jobs, chunksize=2):
            out.append(r)
            if len(out) % 200 == 0:
                print(f"  {len(out)}/{len(jobs)} ({time.time() - t0:.0f} s)", flush=True)
    out.sort(key=lambda r: r[0])
    n = np.array([len(r[1]) for r in out])
    np.savez(os.path.join(WORK, "dots", f"task_{a.task:03d}.npz"),
             mic_index=np.array([r[0] for r in out]), mic_path=np.array([mics[r[0]] for r in out]), n=n,
             x=np.concatenate([r[1] for r in out]), y=np.concatenate([r[2] for r in out]),
             z=np.concatenate([r[3] for r in out]), rnd=np.concatenate([r[4] for r in out]),
             m4p=np.stack([r[5] for r in out]), frac_bad=np.array([r[6] for r in out]),
             n_mirror=np.array([r[7] for r in out]))
    print(f"task {a.task} done: {n.sum()} dots on {len(out)} micrographs, median {np.median(n):.0f}/mic "
          f"({time.time() - t0:.0f} s)", flush=True)


# ----------------------------------------------------------------------------------------------- stage 3
G = {}


def in_valid(mic, X, Y):
    ok = (X >= 0) & (X < MIC_N) & (Y >= 0) & (Y < MIC_N)
    xi = np.clip((X // g0.BIN).astype(int), 0, MIC_N // g0.BIN - 1)
    yi = np.clip((Y // g0.BIN).astype(int), 0, MIC_N // g0.BIN - 1)
    byte = G["m4p"][mic, yi, xi >> 3]
    return ok & (((byte >> (7 - (xi & 7))) & 1) == 1)


def rods_valid(mic, X, Y, ux, uy, h):
    e = h + ENDOFF[-1] / A_MIC
    return in_valid(mic, X, Y) & in_valid(mic, X + e * ux, Y + e * uy) & in_valid(mic, X - e * ux, Y - e * uy)


def pairs(mic, X, Y, ux, uy, h, rmax_px):
    """All (rod, dot) pairs within rmax_px of each rod centre (same micrograph). Returns rod, dot, s, t (px)."""
    lists = G["tree"].query_ball_point(np.c_[X + mic * OFF, Y], r=rmax_px)
    lens = np.fromiter((len(v) for v in lists), int, len(lists))
    dot = np.fromiter(itertools.chain.from_iterable(lists), int, int(lens.sum()))
    rod = np.repeat(np.arange(len(X)), lens)
    rx = G["dx"][dot] - X[rod]
    ry = G["dy"][dot] - Y[rod]
    s = rx * ux[rod] + ry * uy[rod]
    t = -rx * uy[rod] + ry * ux[rod]
    return rod, dot, s, t


def rod_stats(mic, X, Y, ux, uy, h, profiles):
    """Per-rod end hits (+ profile histograms) for one geometry. h in px. Returns dict of arrays."""
    n = len(X)
    rod, dot, s, t = pairs(mic, X, Y, ux, uy, h, h + ANG_R_A / A_MIC + 1)
    hr = h[rod]
    mid = (np.abs(t) <= MID_T_A / A_MIC) & (np.abs(s) <= hr - MID_INSET_A / A_MIC)
    elig = ~mid & ~G["onrod"][dot]
    out = {}
    for eo in ENDOFF:
        for sg, nm in ((1, "plus"), (-1, "minus")):
            de = np.hypot(s - sg * (hr + eo / A_MIC), t) * A_MIC
            for d in D_LIST:
                for dset, dm in (("all", elig), ("strict", elig & G["strict"][dot])):
                    sel = dm & (de <= d)
                    out[(dset, eo, d, nm)] = np.bincount(rod[sel], minlength=n) > 0
    sA, tA = s * A_MIC, np.abs(t) * A_MIC
    beyond = np.abs(s) * A_MIC - hr * A_MIC
    out["onrod_count"] = np.bincount(rod[(tA <= ONROD_T_A) & (beyond <= 0)], minlength=n)
    if profiles:
        rA = np.hypot(sA, tA)
        phi = np.degrees(np.arctan2(tA, sA))
        within = elig & (rA <= ANG_R_A)
        shell = elig & (rA - hr * A_MIC >= SHELL_A[0]) & (rA - hr * A_MIC <= SHELL_A[1])
        ax = (tA <= 25.0)
        nb_ang, nb_ax, nb_t = len(ANG_EDGES) - 1, len(AX_EDGES) - 1, len(MAP_T_EDGES) - 1
        ia = np.clip(np.digitize(phi, ANG_EDGES) - 1, 0, nb_ang - 1)
        ix = np.digitize(beyond, AX_EDGES) - 1
        it = np.digitize(tA, MAP_T_EDGES) - 1
        okx = (ix >= 0) & (ix < nb_ax)
        okt = (it >= 0) & (it < nb_t)
        side = np.where(s >= 0, 0, 1)            # 0 = +axis end (dot side for dot_* classes), 1 = -axis end
        out["prof"] = dict(rod=rod, ia=ia, within=within, shell=shell, ix=ix, ax=ax & okx, side=side,
                           it=it, map=okx & okt)
    return out


def strata_masks(lab):
    return {"all": np.ones(len(lab), bool), "rod_no_dot": lab == "rod_no_dot", "dot_long_rod": lab == "dot_long_rod",
            "dot_short_rod": lab == "dot_short_rod", "dot_classes": np.isin(lab, DOT_LABELS)}


def summarise(st, lab):
    """Scalar statistics per stratum from per-rod hits."""
    S = {}
    for sname, m in strata_masks(lab).items():
        if m.sum() == 0:
            continue
        for dset in ("all", "strict"):
            for eo in ENDOFF:
                for d in D_LIST:
                    p, q = st[(dset, eo, d, "plus")][m], st[(dset, eo, d, "minus")][m]
                    key = f"{sname}|{dset}|off{int(eo)}|d{int(d)}"
                    S[key + "|either"] = float(np.mean(p | q))
                    S[key + "|both"] = float(np.mean(p & q))
                    S[key + "|plus"] = float(np.mean(p))
                    S[key + "|minus"] = float(np.mean(q))
        S[f"{sname}|onrod_per_rod"] = float(st["onrod_count"][m].mean())
    return S


def profile_hists(prof, lab):
    """Histogram counts per stratum group from a profile pair set."""
    rod = prof["rod"]
    groups = {"all": np.ones(len(lab), bool), "rod_no_dot": lab == "rod_no_dot",
              "dot_classes": np.isin(lab, DOT_LABELS)}
    H = {}
    nb_ang, nb_ax, nb_t = len(ANG_EDGES) - 1, len(AX_EDGES) - 1, len(MAP_T_EDGES) - 1
    for gname, gm in groups.items():
        sel = gm[rod]
        H[f"{gname}|ang_within"] = np.bincount(prof["ia"][sel & prof["within"]], minlength=nb_ang)
        H[f"{gname}|ang_shell"] = np.bincount(prof["ia"][sel & prof["shell"]], minlength=nb_ang)
        for sd, snm in ((0, "plus"), (1, "minus")):
            m = sel & prof["ax"] & (prof["side"] == sd)
            H[f"{gname}|ax_{snm}"] = np.bincount(prof["ix"][m], minlength=nb_ax)
        m = sel & prof["map"]
        H[f"{gname}|map"] = np.bincount(prof["ix"][m] * nb_t + prof["it"][m], minlength=nb_ax * nb_t)
    return H


def null_geometry(kind, rng):
    X0, Y0, ux0, uy0, h, mic = (G[k] for k in ("X", "Y", "ux", "uy", "h", "mic"))
    n = len(X0)
    X, Y, ux, uy = X0.copy(), Y0.copy(), ux0.copy(), uy0.copy()
    todo = np.arange(n)
    for _ in range(50):
        if len(todo) == 0:
            break
        if kind == "rotation":
            ph = rng.random(len(todo)) * 2 * np.pi
            c, s = np.cos(ph), np.sin(ph)
            ux[todo] = c * ux0[todo] - s * uy0[todo]
            uy[todo] = s * ux0[todo] + c * uy0[todo]
        else:
            r = (SHIFT_A[0] + (SHIFT_A[1] - SHIFT_A[0]) * rng.random(len(todo))) / A_MIC
            ps = rng.random(len(todo)) * 2 * np.pi
            X[todo] = X0[todo] + r * np.cos(ps)
            Y[todo] = Y0[todo] + r * np.sin(ps)
        ok = rods_valid(mic[todo], X[todo], Y[todo], ux[todo], uy[todo], h[todo])
        todo = todo[~ok]
    return X, Y, ux, uy, len(todo)


def null_worker(args):
    kind, rep, seed = args
    rng = np.random.default_rng(seed)
    X, Y, ux, uy, nfail = null_geometry(kind, rng)
    st = rod_stats(G["mic"], X, Y, ux, uy, G["h"], profiles=True)
    S = summarise(st, G["lab"])
    H = profile_hists(st["prof"], G["lab"])
    return kind, rep, S, H, nfail


def load_dots():
    files = sorted(f for f in os.listdir(os.path.join(WORK, "dots")) if f.endswith(".npz"))
    parts = [np.load(os.path.join(WORK, "dots", f)) for f in files]
    mics = mic_list()
    nm = len(mics)
    m4p = np.zeros((nm, MIC_N // g0.BIN, MIC_N // g0.BIN // 8), np.uint8)
    have = np.zeros(nm, bool)
    xs, ys, zs, ms = [], [], [], []
    nd = np.zeros(nm, int)
    frac_bad = np.full(nm, np.nan)
    n_mirror = np.zeros(nm, int)
    for p in parts:
        mi = p["mic_index"]
        assert np.all(mics[mi] == p["mic_path"])
        m4p[mi] = p["m4p"]
        have[mi] = True
        nd[mi] = p["n"]
        frac_bad[mi] = p["frac_bad"]
        n_mirror[mi] = p["n_mirror"]
        xs.append(p["x"].astype(float))
        ys.append(p["y"].astype(float))
        zs.append(p["z"])
        ms.append(np.repeat(mi, p["n"]))
    return mics, have, m4p, nd, frac_bad, n_mirror, np.concatenate(xs), np.concatenate(ys), np.concatenate(zs), \
        np.concatenate(ms)


def band(v):
    v = np.asarray(v, float)
    return dict(mean=float(v.mean()), sd=float(v.std()), lo95=float(np.percentile(v, 2.5)),
                hi95=float(np.percentile(v, 97.5)))


def test_vs(obs, null):
    null = np.asarray(null, float)
    b = band(null)
    return dict(obs=obs, null=b, ratio=obs / b["mean"] if b["mean"] > 0 else None,
                z=(obs - b["mean"]) / b["sd"] if b["sd"] > 0 else None,
                p_greater=float((1 + np.sum(null >= obs)) / (1 + len(null))),
                p_less=float((1 + np.sum(null <= obs)) / (1 + len(null))))


def stage_analyse(a):
    t0 = time.time()
    rng = np.random.default_rng(SEED + 1)
    conv = json.load(open(os.path.join(WORK, "conventions.json")))
    cbest = {k: conv["best"][k] for k in ("rot", "shift", "order", "swap")}
    ycfg = {k: conv["micrograph_best"][k] for k in ("x", "y", "flipud", "fliplr")}
    ca, labs, geo = all_geometries()
    P = load_particles()
    plab = labs[P["alignments2D__class"].astype(int)]
    ridx = np.nonzero(np.isin(plab, R_LABELS))[0]
    n_R = len(ridx)
    mics, have, m4p, nd, frac_bad, n_mirror, dx, dy, dz, dm = load_dots()
    assert have.all(), f"dot detection missing for {np.sum(~have)} micrographs"
    lut = {p: k for k, p in enumerate(mics)}
    mic = np.array([lut[p] for p in P["location__micrograph_path"][ridx]])
    X, Y, ux, uy, h = map_rods(P, ridx, geo, cbest, ycfg)
    lab = plab[ridx].astype(str)
    G.update(m4p=m4p)
    valid = rods_valid(mic, X, Y, ux, uy, h)
    print(f"R = {n_R} particles; {valid.sum()} with centre and L/2+15 A ends in the valid mask", flush=True)
    # dedup
    keep_v = np.nonzero(valid)[0]
    pri = rng.random(len(keep_v))
    xy = np.c_[X[keep_v] + mic[keep_v] * OFF, Y[keep_v]]
    nb = cKDTree(xy).query_ball_point(xy, DEDUP_R_A / A_MIC)
    ang = np.arctan2(uy[keep_v], ux[keep_v])
    blocked = np.zeros(len(keep_v), bool)
    kept = []
    for i in np.argsort(-pri):
        if blocked[i]:
            continue
        kept.append(i)
        js = np.array(nb[i])
        da = np.abs((ang[js] - ang[i] + np.pi / 2) % np.pi - np.pi / 2)     # axis (undirected) difference
        blocked[js[np.degrees(da) < DEDUP_ANG]] = True
    sel = keep_v[np.sort(np.array(kept))]
    print(f"after dedup: {len(sel)} rods ({len(keep_v) - len(sel)} duplicates removed)", flush=True)
    X, Y, ux, uy, h, mic, lab = X[sel], Y[sel], ux[sel], uy[sel], h[sel], mic[sel], lab[sel]
    pidx = ridx[sel]

    # dots: global tree; flag dots in ANY real rod's middle region
    G.update(tree=cKDTree(np.c_[dx + dm * OFF, dy]), dx=dx, dy=dy, strict=dz >= STRICT_Z,
             onrod=np.zeros(len(dx), bool))
    rod, dot, s, t = pairs(mic, X, Y, ux, uy, h, h + 1)
    onrod = np.zeros(len(dx), bool)
    onrod[dot[(np.abs(t) <= MID_T_A / A_MIC) & (np.abs(s) <= h[rod] - MID_INSET_A / A_MIC)]] = True
    G.update(onrod=onrod, X=X, Y=Y, ux=ux, uy=uy, h=h, mic=mic, lab=lab)
    print(f"dots: {len(dx)} on {len(mics)} micrographs; {onrod.sum()} ({onrod.mean():.3%}) in some rod's middle "
          f"region removed from the eligible set", flush=True)

    st_obs = rod_stats(mic, X, Y, ux, uy, h, profiles=True)
    S_obs = summarise(st_obs, lab)
    H_obs = profile_hists(st_obs["prof"], lab)
    for k in ("all|all|off15|d35|either", "all|all|off0|d35|either"):
        print(f"  observed {k}: {S_obs[k]:.4f}", flush=True)

    jobs = [(kind, r, SEED * 13 + 1000 * j + r) for j, kind in enumerate(("rotation", "shift")) for r in range(a.nrep)]
    procs = int(os.environ.get("SLURM_CPUS_PER_TASK", 4))
    with Pool(procs) as pool:
        res = pool.map(null_worker, jobs, chunksize=1)
    print(f"nulls done ({time.time() - t0:.0f} s)", flush=True)
    NS = {kind: [r[2] for r in res if r[0] == kind] for kind in ("rotation", "shift")}
    NH = {kind: [r[3] for r in res if r[0] == kind] for kind in ("rotation", "shift")}
    nfail = {kind: int(np.mean([r[4] for r in res if r[0] == kind])) for kind in ("rotation", "shift")}

    tests = {}
    for key, ov in S_obs.items():
        tests[key] = {kind: test_vs(ov, [s_[key] for s_ in NS[kind]]) for kind in NS}
    hists = {}
    for key, ov in H_obs.items():
        hists[key] = dict(obs=ov.tolist())
        for kind in NH:
            arr = np.array([hh[key] for hh in NH[kind]], float)
            hists[key][kind] = dict(mean=arr.mean(0).tolist(), lo95=np.percentile(arr, 2.5, 0).tolist(),
                                    hi95=np.percentile(arr, 97.5, 0).tolist())

    # angular verdict numbers: end-shell, phi within 20 deg of the axis (bins 0-20 + 160-180) vs null (b)
    def ang_excess(key, kind="shift"):
        o = np.array(hists[key]["obs"], float)
        nl = np.array([hh[key] for hh in NH[kind]], float)
        end_b = np.r_[0, 1, len(o) - 2, len(o) - 1]
        oo, nn = o[end_b].sum(), nl[:, end_b].sum(1)
        frac_o = oo / o.sum()
        frac_n = nn / nl.sum(1)
        return dict(obs_count_0_20_160_180=float(oo), null_count=band(nn), ratio_count=float(oo / nn.mean()),
                    p_greater=float((1 + np.sum(nn >= oo)) / (1 + len(nn))),
                    obs_fraction=float(frac_o), null_fraction=band(frac_n),
                    total_obs=float(o.sum()), total_null_mean=float(nl.sum(1).mean()))

    def ax_window(key, lo, hi, kind="shift"):
        c = 0.5 * (AX_EDGES[1:] + AX_EDGES[:-1])
        w = (c > lo) & (c < hi)
        o = float(np.array(hists[key]["obs"])[w].sum())
        nl = np.array([hh[key] for hh in NH[kind]], float)[:, w].sum(1)
        return test_vs(o, nl)

    angular = {g: dict(shell=ang_excess(f"{g}|ang_shell"), within120=ang_excess(f"{g}|ang_within"))
               for g in ("all", "rod_no_dot", "dot_classes")}
    along = {}
    for g in ("all", "rod_no_dot", "dot_classes"):
        for sd in ("plus", "minus"):
            k = f"{g}|ax_{sd}"
            along[k] = dict(beyond_end_10_30A=ax_window(k, 10, 30), beyond_end_0_40A=ax_window(k, 0, 40),
                            rod_end_cap_minus20_0A=ax_window(k, -20, 0), rod_body_minus100_minus20A=ax_window(k, -100, -20))

    # class geometry (recomputed here) and where the pick / 2D box centre sits relative to the +axis DNA end
    keep_keys = ("label", "cx", "cy", "ux", "uy", "h_A", "L_halfmax_A", "sdot_from_centre_A", "rod_level_z",
                 "dot_ratio", "box_centre_s_A", "box_centre_t_A", "width_fwhm_A")
    geo_json = {int(k): {kk: g[kk] for kk in keep_keys} for k, g in geo.items()}
    box_summary = {}
    for l_ in R_LABELS:
        ks = [k for k in geo if geo[k]["label"] == l_]
        w = np.array([np.sum(P["alignments2D__class"] == k) for k in ks], float)
        dsb = np.array([geo[k]["box_centre_s_A"] - geo[k]["h_A"] for k in ks])
        tb = np.array([geo[k]["box_centre_t_A"] for k in ks])
        box_summary[l_] = dict(particle_weighted_mean_s_minus_halfL_A=float(np.average(dsb, weights=w)),
                               range_s_minus_halfL_A=[float(dsb.min()), float(dsb.max())],
                               particle_weighted_mean_t_A=float(np.average(tb, weights=w)),
                               note="box centre position along +axis relative to the +axis DNA end (negative = "
                                    "inside the rod); the pick/alignment centre is where the selection acts")
    print("box centre vs +end:", {k: round(v["particle_weighted_mean_s_minus_halfL_A"], 1)
                                  for k, v in box_summary.items()}, flush=True)
    # rods per label, geometry summaries
    counts = {s_: int(m.sum()) for s_, m in strata_masks(lab).items()}
    hA = h * A_MIC
    summary = dict(
        created=time.strftime("%Y-%m-%d %H:%M"), script="scripts/17_rod_end_colocalisation.py",
        question="Are sfGFP dots attached to DNA rod ends more often than chance?",
        inputs=dict(particles=PART_NPZ, class_averages=CA_MRC, labels=LABELS, exposures=CUR_NPZ,
                    detector="scripts/12_dot_spatial_stats.py (imported unmodified)"),
        settings=dict(rod_labels=R_LABELS, d_A=D_LIST, end_offsets_A=ENDOFF, middle_exclusion=dict(
            t_A=MID_T_A, s_inset_from_end_A=MID_INSET_A), dot_side_end_inset_A=R_DOT_A, strict_dot_z=STRICT_Z,
            shift_null_A=SHIFT_A, n_rep=a.nrep, dedup=dict(r_A=DEDUP_R_A, angle_deg=DEDUP_ANG),
            angular_radius_A=ANG_R_A, end_shell_A=SHELL_A, onrod_t_A=ONROD_T_A, seed=SEED),
        conventions=dict(pose_shift=conv["candidates"], best=conv["best"], second_best=conv["second_best"],
                         micrograph_orientation=conv["micrograph_orientation"], micrograph_best=conv["micrograph_best"],
                         n_particles=conv["n_particles"], n_orientation_particles=conv["n_orientation_particles"]),
        class_geometry=geo_json,
        box_centre_vs_rod_end=box_summary,
        rod_set=dict(R_particles=int(n_R),
                     R_by_label={l_: int(np.sum(plab[ridx] == l_)) for l_ in R_LABELS},
                     valid_mask_pass=int(valid.sum()), after_dedup=int(len(sel)), analysed_by_label=counts,
                     n_micrographs=int(len(np.unique(mic))),
                     half_length_A=dict(median=float(np.median(hA)), p10=float(np.percentile(hA, 10)),
                                        p90=float(np.percentile(hA, 90)))),
        dots=dict(total=int(len(dx)), per_mic_median=float(np.median(nd)),
                  per_um2_median=float(np.median(nd / (np.array([np.unpackbits(m4p[i], axis=1).sum()
                                                                for i in range(len(mics))]) * 16 * A_MIC ** 2 / 1e8))),
                  mirror_false_peak_ratio=float(n_mirror.sum() / max(len(dx), 1)),
                  strict_fraction=float(np.mean(dz >= STRICT_Z)), removed_on_rod_middles=int(onrod.sum())),
        null_unplaced_rods_mean=nfail,
        false_peak_check=dict(
            onrod_dots_per_rod=tests["all|onrod_per_rod"],
            definition="dots with |t| <= 15 A and |s| <= L/2 (ON the rod between its ends), all dots, vs null (b)"),
        T=tests, angular=angular, along_axis=along, hists=hists,
        bins=dict(angle_deg=ANG_EDGES.tolist(), along_axis_A=AX_EDGES.tolist(), map_t_A=MAP_T_EDGES.tolist()),
        runtime_s=round(time.time() - t0, 1))
    with open(OUT_JSON, "w") as fh:
        json.dump(summary, fh, indent=1)
    np.savez(os.path.join(WORK, "rods_analysed.npz"), pidx=pidx, X=X, Y=Y, ux=ux, uy=uy, h=h, mic=mic, lab=lab,
             hit35=st_obs[("all", 15.0, 35.0, "plus")] | st_obs[("all", 15.0, 35.0, "minus")])
    print(f"wrote {OUT_JSON} ({time.time() - t0:.0f} s)", flush=True)
    fig_main(summary, conv, mics, dx, dy, dm, onrod, st_obs, rng)
    print(f"done ({time.time() - t0:.0f} s)", flush=True)


def fig_main(S, conv, mics, dx, dy, dm, onrod, st_obs, rng):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import Circle
    T = S["T"]
    fig = plt.figure(figsize=(24, 19), dpi=90)
    sf = fig.subfigures(3, 1, height_ratios=[1, 1, 1.05])
    top = sf[0].subplots(1, 4, gridspec_kw=dict(width_ratios=[1.25, 1, 1, 1.1]))
    # (a) convention bars
    ax = top[0]
    cands = S["conventions"]["pose_shift"]
    m = [c["mean_corr"] for c in cands]
    ax.barh(range(len(m)), m, color=["C2" if v == max(m[:-1]) else "0.6" for v in m[:-1]] + ["C3"])
    ax.set_yticks(range(len(m)))
    ax.set_yticklabels([c["name"] for c in cands], fontsize=6)
    ax.invert_yaxis()
    ax.set_xlabel("mean r (particle in class frame vs class avg)", fontsize=8)
    yb = S["conventions"]["micrograph_best"]
    ax.set_title(f"(a) pose/shift convention, n={S['conventions']['n_particles']}\nmicrograph: x {yb['x']}, "
                 f"y {yb['y']}, flipud {yb['flipud']} (r={yb['mean_corr']:.2f})", fontsize=9)
    # (b,c) T vs nulls
    for j, eo in enumerate(ENDOFF):
        ax = top[1 + j]
        for k, d in enumerate(D_LIST):
            key = f"all|all|off{int(eo)}|d{int(d)}|either"
            o = T[key]["rotation"]["obs"]
            ax.plot(k, o, "ko", ms=7, zorder=3)
            for dk, (kind, col) in enumerate((("rotation", "C0"), ("shift", "C1"))):
                b = T[key][kind]["null"]
                xx = k + (dk - 0.5) * 0.25
                ax.errorbar(xx, b["mean"], yerr=[[b["mean"] - b["lo95"]], [b["hi95"] - b["mean"]]], fmt="s",
                            color=col, ms=4, capsize=3, label=f"null {kind}" if k == 0 else None)
            ax.text(k, o + 0.03, f"x{T[key]['rotation']['ratio']:.2f} / x{T[key]['shift']['ratio']:.2f}",
                    ha="center", fontsize=8)
            kb = f"all|all|off{int(eo)}|d{int(d)}|both"
            ax.plot(k, T[kb]["rotation"]["obs"], "kD", mfc="none", ms=6)
            ax.plot(k + 0.12, T[kb]["shift"]["null"]["mean"], "D", color="C1", mfc="none", ms=4)
        ax.set_xticks(range(len(D_LIST)))
        ax.set_xticklabels([f"d = {int(d)} A" for d in D_LIST])
        ax.set_ylim(0, 1)
        ax.set_ylabel("fraction of rods with a dot at an end")
        ax.legend(fontsize=7, loc="upper left")
        ax.set_title(f"({'bc'[j]}) T, end = L/2 + {int(eo)} A, n={S['rod_set']['after_dedup']} rods\n"
                     f"filled = either end, open diamond = both ends; ratio obs/null rot / shift", fontsize=9)
    # (d) per-label ratios
    ax = top[3]
    strata = ["all", "rod_no_dot", "dot_long_rod", "dot_short_rod"]
    for k, sname in enumerate(strata):
        for dk, (kind, col) in enumerate((("rotation", "C0"), ("shift", "C1"))):
            key = f"{sname}|all|off15|d35|either"
            if key not in T:
                continue
            tt = T[key][kind]
            b = tt["null"]
            r, lo, hi = tt["obs"] / b["mean"], tt["obs"] / b["hi95"], tt["obs"] / b["lo95"]
            ax.errorbar(k + (dk - 0.5) * 0.25, r, yerr=[[r - lo], [hi - r]], fmt="o", color=col, capsize=3,
                        label=f"either end vs {kind}" if k == 0 else None)
    for k, sname in enumerate(["dot_classes"]):
        for end, mk in (("plus", "^"), ("minus", "v")):
            key = f"dot_classes|all|off15|d35|{end}"
            tt = T[key]["shift"]
            ax.plot(4 + (0.15 if end == "minus" else -0.15), tt["obs"] / tt["null"]["mean"], mk, color="C3", ms=8,
                    label=f"dot classes {'near (dot side)' if end == 'plus' else 'FAR'} end vs shift")
    ax.axhline(1, color="0.5", lw=0.8)
    ax.set_xticks(range(5))
    ax.set_xticklabels([s_.replace("_rod", "") + f"\nn={S['rod_set']['analysed_by_label'].get(s_, '')}"
                        for s_ in strata] + ["dot classes\nnear/far"], fontsize=7)
    ax.set_ylabel("T_obs / null mean (bar: null 95% band)")
    ax.legend(fontsize=6)
    ax.set_title("(d) per consensus label, d = 35 A, end = L/2 + 15 A", fontsize=9)

    mid = sf[1].subplots(1, 4)
    H = S["hists"]
    ca_ = 0.5 * (ANG_EDGES[1:] + ANG_EDGES[:-1])
    for j, (key, ttl) in enumerate((("ang_within", f"dots within {int(ANG_R_A)} A of rod centre (rod middle excluded)"),
                                    ("ang_shell", f"end shell: r - L/2 in {SHELL_A[0]:.0f}-{SHELL_A[1]:.0f} A"))):
        ax = mid[j]
        for g, col in (("rod_no_dot", "C0"), ("dot_classes", "C3")):
            o = np.array(H[f"{g}|{key}"]["obs"], float)
            nl = H[f"{g}|{key}"]["shift"]
            mn = np.array(nl["mean"])
            ax.step(ca_, o / mn, where="mid", color=col, lw=1.6, label=f"{g} obs / null(b)")
            ax.fill_between(ca_, np.array(nl["lo95"]) / mn, np.array(nl["hi95"]) / mn, step="mid", color=col,
                            alpha=0.2, lw=0)
        ax.axhline(1, color="0.5", lw=0.8)
        ax.set_xlabel("angle to +axis (deg); 0 = +end (dot side for dot classes), 180 = -end")
        ax.set_ylabel("observed / shift-null dots per 10 deg")
        ax.set_xticks(range(0, 181, 30))
        ax.legend(fontsize=7)
        ax.set_title(f"({'ef'[j]}) angular profile: {ttl}", fontsize=9)
    ax = mid[2]
    cx_ = 0.5 * (AX_EDGES[1:] + AX_EDGES[:-1])
    for g, sd, col, ls in (("rod_no_dot", "plus", "C0", "-"), ("dot_classes", "plus", "C3", "-"),
                           ("dot_classes", "minus", "C3", "--")):
        k = f"{g}|ax_{sd}"
        o = np.array(H[k]["obs"], float)
        mn = np.array(H[k]["shift"]["mean"])
        mn = np.where(mn > 0, mn, np.nan)
        ax.plot(cx_, o / mn, color=col, ls=ls, lw=1.5,
                label=f"{g} {'+' if sd == 'plus' else '-'} end" + (" (dot side)" if g == "dot_classes" and sd == "plus"
                                                                   else " (FAR)" if g == "dot_classes" else ""))
        ax.fill_between(cx_, np.array(H[k]["shift"]["lo95"]) / mn, np.array(H[k]["shift"]["hi95"]) / mn,
                        color=col, alpha=0.12, lw=0)
    ax.axhline(1, color="0.5", lw=0.8)
    ax.axvline(0, color="k", lw=0.8)
    ax.axvspan(15, 25, color="g", alpha=0.12)
    ax.set_xlabel("|s| - L/2 (A): distance beyond the rod end along the axis (|t| <= 25 A); all dots")
    ax.set_ylabel("observed / shift-null")
    ax.legend(fontsize=7)
    ax.set_title("(g) along-axis profile (green = where an attached sfGFP centre should sit)", fontsize=9)
    ax = mid[3]
    nb_t = len(MAP_T_EDGES) - 1
    o = np.array(H["all|map"]["obs"], float).reshape(-1, nb_t)
    mn = np.array(H["all|map"]["shift"]["mean"]).reshape(-1, nb_t)
    im = ax.imshow((o / np.where(mn > 0, mn, np.nan)).T, origin="lower", cmap="RdBu_r", vmin=0, vmax=2,
                   extent=(AX_EDGES[0], AX_EDGES[-1], MAP_T_EDGES[0], MAP_T_EDGES[-1]), aspect="auto")
    plt.colorbar(im, ax=ax, label="observed / shift-null dot density")
    ax.axvline(0, color="k", lw=0.8)
    ax.set_xlabel("|s| - L/2 (A)  (negative = alongside the rod)")
    ax.set_ylabel("|t| distance from the axis line (A)")
    ax.set_title("(h) dot density around rods (all R, all dots), folded over both ends", fontsize=9)

    # examples
    bot = sf[2].subplots(1, 6)
    R = np.load(os.path.join(WORK, "rods_analysed.npz"))
    hit = R["hit35"]
    cur = np.load(CUR_NPZ)
    picks = list(rng.choice(np.nonzero(hit)[0], 3, replace=False)) + list(rng.choice(np.nonzero(~hit)[0], 3,
                                                                                     replace=False))
    half = 180
    for ax, i in zip(bot, picks):
        mi = int(R["mic"][i])
        df1, df2, an = mic_ctf(cur, mics[mi])
        pos, _ = g0.load_and_filter(os.path.join(PROJ, mics[mi]), df1, df2, an)
        x0 = int(np.clip(R["X"][i] - half, 0, MIC_N - 2 * half))
        y0 = int(np.clip(R["Y"][i] - half, 0, MIC_N - 2 * half))
        show = -pos[y0:y0 + 2 * half, x0:x0 + 2 * half]
        lo, hi = np.percentile(show, [1, 99])
        ax.imshow(show, cmap="gray", vmin=lo, vmax=hi, extent=(x0 - .5, x0 + 2 * half - .5, y0 + 2 * half - .5,
                                                                y0 - .5))
        on = np.nonzero((R["mic"] == mi) & (np.abs(R["X"] - R["X"][i]) < half + 150) &
                        (np.abs(R["Y"] - R["Y"][i]) < half + 150))[0]
        for j in on:
            xx, yy, vx, vy, hh = R["X"][j], R["Y"][j], R["ux"][j], R["uy"][j], R["h"][j]
            lw = 2.0 if j == i else 1.0
            ax.plot([xx - hh * vx, xx + hh * vx], [yy - hh * vy, yy + hh * vy], color="lime", lw=lw)
            for sg in (1, -1):
                e = hh + 15 / A_MIC
                ax.add_patch(Circle((xx + sg * e * vx, yy + sg * e * vy), 35 / A_MIC, fill=False,
                                    ec="m" if sg > 0 else "c", lw=0.9, ls="--"))
        dsel = np.nonzero((dm == mi) & (dx >= x0) & (dx < x0 + 2 * half) & (dy >= y0) & (dy < y0 + 2 * half))[0]
        for k in dsel:
            ax.add_patch(Circle((dx[k], dy[k]), 10, fill=False, ec="orange" if onrod[k] else "red", lw=1.0))
        ax.set_xlim(x0, x0 + 2 * half)
        ax.set_ylim(y0 + 2 * half, y0)
        ax.set_xticks([])
        ax.set_yticks([])
        ax.plot([x0 + 15, x0 + 15 + 100 / A_MIC], [y0 + 2 * half - 15] * 2, "w-", lw=3)
        ax.set_title(f"{R['lab'][i]} ({'end dot' if hit[i] else 'no end dot'} at d=35)", fontsize=8)
    sf[2].suptitle("(i) random rods (3 with, 3 without an end dot): phase-flipped 15 A crops, protein dark, 511 A wide, "
                   "bar 10 nm. lime = rod axis between DNA ends (thick = the example); dashed circles = d 35 A around "
                   "L/2+15 A (magenta +end, cyan -end); red = detected dots, orange = dots on a rod middle (excluded)",
                   fontsize=9)
    fig.savefig(FIG)
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("stage", choices=["conventions", "detect", "analyse"])
    ap.add_argument("--task", type=int, default=0)
    ap.add_argument("--ntasks", type=int, default=1)
    ap.add_argument("--nconv", type=int, default=3000)
    ap.add_argument("--nymic", type=int, default=40)
    ap.add_argument("--nyper", type=int, default=10)
    ap.add_argument("--nrep", type=int, default=NREP)
    a = ap.parse_args()
    {"conventions": stage_conventions, "detect": stage_detect, "analyse": stage_analyse}[a.stage](a)


if __name__ == "__main__":
    main()
