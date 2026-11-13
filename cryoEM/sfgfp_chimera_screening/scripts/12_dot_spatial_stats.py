"""G0 dot spatial statistics on the 300 pilot micrographs (CPU only; no CryoSPARC, no synthetic images).

Question: the pilot 2D (P0_2D) shows many "two dots ~120 A apart" classes and dot pairs touching at ~40 A. Are
sfGFP-like dots really paired at those distances on the micrographs, or is this what a random arrangement of dots
looks like once a 2D classifier centres on one dot? We answer with coordinates only:

  1. per micrograph: CTF phase flip (own astigmatic CTF, 200 kV, Cs 2.7 mm, amp 0.1) -> protein-POSITIVE image
     -> Gaussian low-pass (half-max at 1/15 A^-1);
  2. dots = local maxima of the scale-normalised LoG response (-sigma^2 LoG, sigma = 12 A) with a robust
     per-micrograph threshold (median + k * 1.4826 * MAD over the valid area), 20 A non-max suppression,
     50 px edge exclusion and a simple contamination mask (low-frequency mean / RMS outliers, bin 4);
  3. NN-distance histogram + pair-correlation g(r) (0-300 A, 5 A bins) vs two nulls built only from real
     coordinates: (a) permutation (pooled coordinates re-assigned to each micrograph's valid mask, counts
     preserved, same 20 A hard core) and (b) toroidal shift (another micrograph of the same defocus bin, randomly
     shifted with wrap-around, used as an independent partner pattern for cross-NN / cross-pair counts);
  4. a Hessian ridge ("DNA thread") score sampled on rings 25-120 A around each dot -> bare / one-arm /
     two-arm-open / two-arm-closed (loop candidate), with the same scoring at random background positions of the
     same micrographs as the false-arm baseline;
  5. dot densities per micrograph and per defocus bin.

  /n/holylabs/hekstra_lab/Everyone/ianpoon/envs/dnaspring/bin/python -I scripts/12_dot_spatial_stats.py [--n 5]

Outputs: results/G0_spatial.json, results/G0_dots_pilot.csv.gz, figures/qc/G0_nn_hist.png,
figures/qc/G0_dot_detection_examples.png
"""
import argparse
import gzip
import json
import os
import sys
import time
from multiprocessing import Pool

import numpy as np
import mrcfile
from scipy import ndimage as ndi
from scipy.spatial import cKDTree

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PROJ = "/n/netscratch/hekstra_lab/Lab/ian_poon/cryosparc_projects/CS-cs-sfgfp-chimera-screening"
EXPORTS = "/n/netscratch/hekstra_lab/Lab/ian_poon/sfgfp_chimera_screening_20261112/exports"
CTF_NPZ = os.path.join(EXPORTS, "P0_ctf__exposures.npz")

APIX = 1.42
KV, CS_MM, AMP = 200.0, 2.7, 0.10
LOWPASS_A = 15.0
DOT_SIGMA_A = 12.0          # LoG scale: blob radius ~ sigma*sqrt2 ~ 17 A, i.e. a ~3.4 nm dot
NMS_A = 20.0                # minimum dot separation
MIN_ROUND = 0.2             # l2/l1 Hessian roundness at the dot scale: drops thread segments (rare at k=4.5)
EDGE_PX = 52                # 50 px requested, rounded up to a multiple of the bin-4 mask grid
BIN = 4                     # mask grid (5.68 A)
RIDGE_SIGMA_A = 7.0
ARM_R_A = (25.0, 120.0)
ARM_DR_A = 2.5
ARM_NTH = 72                # 5 deg angular bins
ARM_SHORT_A = 60.0          # an arm = a ridge path from r = 25 A (touching the dot) out to 60 A
ARM_PEN = 0.3               # DP penalty (ridge noise units) per 5-deg lateral step (paths may curve <= 2 deg/A)
ARM_NMS_DEG = 45.0          # distinct arms leave the dot >= 45 deg apart and their 25-60 A paths never meet
ARM_CLOSE_MIN_R_A = 60.0    # two arms whose full (25-120 A) paths meet (<= 10 deg) at r >= 60 A -> converging
ARM_MIN_SEP_DEG = 60.0      # ... and which leave the dot >= 60 deg apart -> "two_arm_closed" (loop candidate)
DOT_EXCL_A = 22.0           # ridge zeroed within 22 A of every detected dot (neighbour-dot halos are not arms)
JITTER_A = 200.0            # null (c): each dot displaced uniformly within a 200 A disk (keeps >~0.5 um clustering)
R_MAX_A, R_BIN_A = 300.0, 5.0
DF_BINS_UM = [0.0, 1.0, 1.5, 2.0, 10.0]
WINDOWS = {"touching_35_50A": (35.0, 50.0), "two_dot_110_130A": (110.0, 130.0),
           "posthoc_contact_30_40A": (30.0, 40.0)}   # added after the 8-micrograph test showed a 30-40 A bump


# ----------------------------------------------------------------------------------------------- image processing
def wavelength_A(kv):
    v = kv * 1e3
    return 12.2643247 / np.sqrt(v * (1 + 0.978466e-6 * v))


def ctf_2d(shape, df1, df2, ang, apix=APIX):
    """Transfer T(k) = -sin(chi + asin(w)); chi = pi*lam*df*k^2 - pi/2*Cs*lam^3*k^4 (df > 0 underfocus), so the
    low-frequency lobe is negative: protein (positive potential) images dark. rfft2 layout."""
    ny, nx = shape
    ky = np.fft.fftfreq(ny, d=apix)[:, None]
    kx = np.fft.rfftfreq(nx, d=apix)[None, :]
    k2 = kx ** 2 + ky ** 2
    th = np.arctan2(ky, kx)
    df = 0.5 * (df1 + df2) + 0.5 * (df1 - df2) * np.cos(2 * (th - ang))
    lam = wavelength_A(KV)
    cs = CS_MM * 1e7
    chi = np.pi * lam * df * k2 - 0.5 * np.pi * cs * lam ** 3 * k2 ** 2
    return -np.sin(chi + np.arcsin(AMP)), k2


def load_and_filter(path, df1, df2, ang):
    with mrcfile.open(path, permissive=True) as m:
        a = np.asarray(m.data, dtype=np.float32)
    a = a - a.mean()
    a /= a.std()
    F = np.fft.rfft2(a)
    T, k2 = ctf_2d(a.shape, df1, df2, ang)
    s = (1.0 / LOWPASS_A) / np.sqrt(2 * np.log(2))       # Gaussian, half-max at 1/15 A^-1
    G = np.exp(-k2 / (2 * s * s))
    G[0, 0] = 0.0
    pos = np.fft.irfft2(F * np.sign(T) * G, s=a.shape).astype(np.float32)   # phase-flipped, protein POSITIVE
    raw_lp = np.fft.irfft2(F * G, s=a.shape).astype(np.float32)             # not flipped (protein dark)
    return pos, raw_lp


def robust_z(x, sel):
    v = x[sel]
    med = np.median(v)
    mad = 1.4826 * np.median(np.abs(v - med))
    return (x - med) / (mad if mad > 0 else 1.0), float(med), float(mad)


def valid_mask(pos):
    """Bin-4 mask: edge exclusion + contamination = local (68 A Gaussian) mean or RMS of the protein-positive
    15 A image beyond a robust 5-sigma (mean two-sided, RMS high side), dilated by ~57 A. Returns (mask4, frac
    flagged as contamination, n_components)."""
    ny, nx = pos.shape
    p4 = pos.reshape(ny // BIN, BIN, nx // BIN, BIN).mean(axis=(1, 3))
    m4 = np.zeros(p4.shape, bool)
    e = EDGE_PX // BIN
    m4[e:-e, e:-e] = True
    loc_mean = ndi.gaussian_filter(p4, 12)
    loc_rms = np.sqrt(ndi.gaussian_filter(p4 ** 2, 12))
    zm, _, _ = robust_z(loc_mean, m4)
    zs, _, _ = robust_z(loc_rms, m4)
    bad = (np.abs(zm) > 5) | (zs > 5)
    bad = ndi.binary_dilation(bad, iterations=10)
    frac_bad = float((bad & m4).sum() / m4.sum())
    return m4 & ~bad, frac_bad


def detect_dots(pos, mask_full, k):
    sig = DOT_SIGMA_A / APIX
    R = -(sig ** 2) * ndi.gaussian_laplace(pos, sig)
    z, med, mad = robust_z(R, mask_full)
    half = int(np.floor(NMS_A / APIX / np.sqrt(2)))      # square inscribed in the 20 A circle
    mx = ndi.maximum_filter(z, size=2 * half + 1, mode="nearest")
    pk = (z == mx) & mask_full
    yy, xx = np.nonzero(pk & (z > k))
    zz = z[yy, xx]
    # mirror count: same detector on the negated response (dark-on-positive = noise-only false-peak estimate)
    mn = ndi.minimum_filter(z, size=2 * half + 1, mode="nearest")
    neg = (z == mn) & mask_full & (z < -k)
    n_mirror = int(neg.sum())
    # roundness at the dot scale: Hessian eigenvalue ratio l2/l1 (1 = round blob, ~0 = thread segment)
    hyy = ndi.gaussian_filter(pos, sig, order=(2, 0))[yy, xx]
    hxx = ndi.gaussian_filter(pos, sig, order=(0, 2))[yy, xx]
    hxy = ndi.gaussian_filter(pos, sig, order=(1, 1))[yy, xx]
    tr, dd = 0.5 * (hxx + hyy), np.sqrt((0.5 * (hxx - hyy)) ** 2 + hxy ** 2)
    rnd = np.clip((tr + dd) / np.minimum(tr - dd, -1e-12), -1, 1)
    ok = rnd >= MIN_ROUND
    xx, yy, zz, rnd = xx[ok], yy[ok], zz[ok], rnd[ok]
    keep = greedy_nms(np.c_[xx, yy].astype(float), zz, NMS_A / APIX)
    return xx[keep], yy[keep], zz[keep], R[yy[keep], xx[keep]], n_mirror, z, rnd[keep]


def greedy_nms(xy, strength, r):
    if len(xy) == 0:
        return np.zeros(0, int)
    order = np.argsort(-strength)
    tree = cKDTree(xy)
    nb = tree.query_ball_point(xy, r)
    blocked = np.zeros(len(xy), bool)
    keep = []
    for i in order:
        if blocked[i]:
            continue
        keep.append(i)
        blocked[nb[i]] = True
    return np.array(sorted(keep), int)


def ridge_map(pos, mask_full):
    """Bright-line (DNA) Hessian ridge measure sigma^2 * max(0, -l1 - |l2|), l1 <= l2: ~ -l1 on a line, ~0 on a
    round blob (l1 ~ l2) -> robust z over the valid area."""
    s = RIDGE_SIGMA_A / APIX
    hyy = ndi.gaussian_filter(pos, s, order=(2, 0))
    hxx = ndi.gaussian_filter(pos, s, order=(0, 2))
    hxy = ndi.gaussian_filter(pos, s, order=(1, 1))
    tr = 0.5 * (hxx + hyy)
    d = np.sqrt((0.5 * (hxx - hyy)) ** 2 + hxy ** 2)
    l1, l2 = tr - d, tr + d
    rid = np.maximum(0.0, -l1 - np.abs(l2))
    # noise unit = robust SD of the second-derivative filter output (hxx) over the valid area
    v = hxx[mask_full][::7]
    unit = 1.4826 * np.median(np.abs(v - np.median(v)))
    return (rid / unit).astype(np.float32)


RADII_A = np.arange(ARM_R_A[0], ARM_R_A[1] + 1e-6, ARM_DR_A)
THETA = np.arange(ARM_NTH) * 2 * np.pi / ARM_NTH


def polar_samples(rz, x, y):
    r = RADII_A / APIX
    X = x[:, None, None] + r[None, :, None] * np.cos(THETA)[None, None, :]
    Y = y[:, None, None] + r[None, :, None] * np.sin(THETA)[None, None, :]
    v = ndi.map_coordinates(rz, [Y.ravel(), X.ravel()], order=1, mode="nearest")
    return v.reshape(len(x), len(RADII_A), ARM_NTH)


N_SHORT = int(np.sum(RADII_A <= ARM_SHORT_A))
CLOSE_R_IDX = int(np.searchsorted(RADII_A, ARM_CLOSE_MIN_R_A))


def dp_paths(P, r1, pen=ARM_PEN):
    """Backward dynamic programming over rings 0..r1-1: for every start angle at r = 25 A, the outward path (<= 1
    angular bin per radial step, penalty pen per lateral step) maximising the summed ridge score. Returns the
    mean score per start angle and the successor table for tracing."""
    nt = P.shape[1]
    B = P[r1 - 1].astype(float).copy()
    succ = np.zeros((r1, nt), int)
    steps = np.array([-1, 0, 1])
    for r in range(r1 - 2, -1, -1):
        c = np.stack([np.roll(B, 1) - pen, B, np.roll(B, -1) - pen])   # value of moving to theta-1, theta, theta+1
        k = np.argmax(c, axis=0)
        succ[r] = (np.arange(nt) + steps[k]) % nt
        B = P[r] + c[k, np.arange(nt)]
    return B / r1, succ


def trace(succ, t0, r1):
    path = [t0]
    for r in range(r1 - 1):
        path.append(succ[r, path[-1]])
    return np.array(path)


def _adist(a, b, n=ARM_NTH):
    d = np.abs(np.asarray(a) - np.asarray(b)) % n
    return np.minimum(d, n - d)


def classify_arms(P, thr):
    """P: (n_rings, n_theta) ridge score on rings 25-120 A around one position.
    arms: start angles whose best 25-60 A path has mean ridge > thr, greedy by score, >= 45 deg apart at the start,
    and whose 25-60 A paths never come within 1 bin of an already accepted arm's path (else the same thread).
    two_arm_closed: a pair of arms leaving >= 60 deg apart whose best full-range (25-120 A) paths come within
    2 bins (10 deg) of each other at r >= 60 A (two strands converging on the far side, as a ~43 A-radius loop
    with the dot on its rim would: far point at ~86 A)."""
    ms, ss = dp_paths(P, N_SHORT)
    nms = int(round(ARM_NMS_DEG / (360.0 / ARM_NTH)))
    arms, paths = [], []
    for t0 in np.argsort(-ms):
        if ms[t0] <= thr:
            break
        if any(_adist(t0, a) < nms for a in arms):
            continue
        pth = trace(ss, t0, N_SHORT)
        if any((_adist(pth, q) <= 1).any() for q in paths):
            continue
        arms.append(int(t0))
        paths.append(pth)
    n = len(arms)
    closed = False
    if n >= 2:
        _, sf = dp_paths(P, len(RADII_A))
        full = {a: trace(sf, a, len(RADII_A)) for a in arms}
        sep_min = ARM_MIN_SEP_DEG / (360.0 / ARM_NTH)
        for ia in range(n):
            for ib in range(ia + 1, n):
                a, b = arms[ia], arms[ib]
                if _adist(a, b) < sep_min:
                    continue
                if (_adist(full[a][CLOSE_R_IDX:], full[b][CLOSE_R_IDX:]) <= 2).any():
                    closed = True
    if n == 0:
        return 0, "bare"
    if n == 1:
        return 1, "one_arm"
    return n, ("two_arm_closed" if closed else "two_arm_open")


CLASSES = ["bare", "one_arm", "two_arm_open", "two_arm_closed"]


def set_covariance(m4, r_px_full):
    """Window set covariance gamma(r) (rotation-averaged area of W & (W+r)), in full-res px^2."""
    ny, nx = m4.shape
    f = np.fft.rfft2(m4.astype(np.float64), s=(2 * ny, 2 * nx))
    ac = np.fft.irfft2(f * np.conj(f), s=(2 * ny, 2 * nx))
    ac = np.fft.fftshift(ac)
    yy, xx = np.indices(ac.shape)
    rr = np.hypot(yy - ny, xx - nx)
    ri = np.rint(rr).astype(int)
    nmax = int(np.ceil(r_px_full.max() / BIN)) + 2
    sel = ri <= nmax
    prof = np.bincount(ri[sel], weights=ac[sel], minlength=nmax + 1) / np.maximum(
        np.bincount(ri[sel], minlength=nmax + 1), 1)
    return np.interp(r_px_full / BIN, np.arange(nmax + 1), prof[: nmax + 1]) * BIN * BIN


EDGES_A = np.arange(0.0, R_MAX_A + R_BIN_A / 2, R_BIN_A)
CENT_A = 0.5 * (EDGES_A[1:] + EDGES_A[:-1])


# ----------------------------------------------------------------------------------------------- per micrograph
def process(args):
    i, path, df1, df2, ang, k, ridge_thr, want_images, seed = args
    t0 = time.time()
    pos, raw_lp = load_and_filter(path, df1, df2, ang)
    m4, frac_bad = valid_mask(pos)
    mfull = np.repeat(np.repeat(m4, BIN, axis=0), BIN, axis=1)
    x, y, zz, resp, n_mirror, zmap, rnd = detect_dots(pos, mfull, k)
    rz = ridge_map(pos, mfull)
    if len(x):   # zero the ridge map within 22 A of every dot: halos of neighbouring dots must not count as arms
        r_ex = int(np.ceil(DOT_EXCL_A / APIX))
        yy_, xx_ = np.mgrid[-r_ex:r_ex + 1, -r_ex:r_ex + 1]
        disk = (xx_ ** 2 + yy_ ** 2) <= (DOT_EXCL_A / APIX) ** 2
        dm = np.zeros(rz.shape, bool)
        dm[y, x] = True
        dm = ndi.binary_dilation(dm, structure=disk)
        rz[dm] = 0.0
    arms = np.zeros(len(x), int)
    cls = np.empty(len(x), object)
    if len(x):
        P = polar_samples(rz, x.astype(float), y.astype(float))
        for j in range(len(x)):
            arms[j], cls[j] = classify_arms(P[j], ridge_thr)
    # control: random valid background positions >= 30 A from every dot (same count, same micrograph)
    rng = np.random.default_rng(seed)
    vy, vx = np.nonzero(mfull)
    cand = rng.integers(0, len(vx), size=4 * max(len(x), 50))
    cx, cy = vx[cand].astype(float), vy[cand].astype(float)
    if len(x):
        dmin, _ = cKDTree(np.c_[x, y]).query(np.c_[cx, cy])
        ok = dmin * APIX >= 30.0
        cx, cy = cx[ok], cy[ok]
    cx, cy = cx[: len(x)], cy[: len(x)]
    ctrl_cls = []
    if len(cx):
        Pc = polar_samples(rz, cx, cy)
        ctrl_cls = [classify_arms(Pc[j], ridge_thr)[1] for j in range(len(cx))]
    gamma = set_covariance(m4, CENT_A / APIX)
    # sign check: phase-flipped protein-positive image should anti-correlate with the unflipped image
    sl = (slice(200, -200), slice(200, -200))
    cc = float(np.corrcoef(pos[sl].ravel()[::7], raw_lp[sl].ravel()[::7])[0, 1])
    out = dict(i=i, x=x.astype(np.int32), y=y.astype(np.int32), z=zz.astype(np.float32), rnd=rnd.astype(np.float32),
               resp=resp.astype(np.float32), arms=arms, cls=cls, ctrl_cls=np.array(ctrl_cls, object),
               m4=m4, area_px=float(m4.sum() * BIN * BIN), frac_contam=frac_bad, n_mirror=n_mirror,
               gamma=gamma, cc_flip_vs_raw=cc, sec=time.time() - t0)
    if want_images:
        out["pos"] = pos
        out["raw_lp"] = raw_lp
        out["rz"] = rz
    return out


# ----------------------------------------------------------------------------------------------- nulls
G = {}   # globals for forked workers


def pair_counts(tree_a, tree_b, self_pairs):
    c = tree_a.count_neighbors(tree_b, EDGES_A / APIX, cumulative=True).astype(float)
    if self_pairs:
        c = c - tree_a.n          # remove the zero-distance self pairs (every edge >= 0 includes them)
        c[0] = 0
    return np.diff(c)


def nn_hist(d_A):
    h, _ = np.histogram(d_A, bins=EDGES_A)
    return h.astype(float)


def obs_stats(xy):
    if len(xy) < 2:
        return np.zeros(len(CENT_A)), np.zeros(len(CENT_A)), np.zeros(0)
    t = cKDTree(xy)
    d = t.query(xy, k=2)[0][:, 1] * APIX
    return nn_hist(d), pair_counts(t, t, True), d


def in_mask(m4, xy):
    ix = np.clip((xy[:, 0] // BIN).astype(int), 0, m4.shape[1] - 1)
    iy = np.clip((xy[:, 1] // BIN).astype(int), 0, m4.shape[0] - 1)
    return m4[iy, ix]


def perm_null_one(i, rng):
    """Permutation null for micrograph i: pooled real coordinates (other micrographs) that fall inside i's valid
    mask, taken in random order with the detector's 20 A hard core, until N_i dots."""
    n = G["N"][i]
    pool_xy, pool_mic = G["pool_xy"], G["pool_mic"]
    m4 = G["masks"][i]
    acc = np.zeros((0, 2))
    for _ in range(200):
        idx = rng.integers(0, len(pool_xy), size=int(2.5 * n) + 50)
        idx = idx[pool_mic[idx] != i]
        c = pool_xy[idx]
        c = c[in_mask(m4, c)]
        c = np.vstack([acc, c])
        keep = greedy_nms(c, rng.random(len(c)) + 1e9 * (np.arange(len(c)) < len(acc)), G["r_hc_px"])
        acc = c[keep]
        if len(acc) >= n:
            acc = acc[rng.permutation(len(acc))[:n]]
            break
    return obs_stats(acc)[:2]


def toro_null_one(i, rng):
    """Toroidal-shift null for micrograph i: dots of other micrographs in the same defocus bin, each randomly
    shifted with wrap-around in the 2048^2 frame and cut to i's valid mask, until N_i-1 partner dots; statistics are
    CROSS (each real dot of i -> partner pattern), with the same 20 A hard core on the NN."""
    n = G["N"][i]
    xy_i = G["xy"][i]
    partners = G["bin_members"][G["dfbin"][i]]
    partners = partners[partners != i]
    if len(partners) == 0:
        partners = np.array([j for j in range(len(G["N"])) if j != i])
    m4 = G["masks"][i]
    parts = [np.zeros((0, 2))]
    tot = 0
    for _ in range(500):
        if tot >= n - 1:
            break
        j = partners[rng.integers(len(partners))]
        if G["N"][j] == 0:
            continue
        sh = rng.random(2) * 2048
        c = np.mod(G["xy"][j] + sh, 2048.0)
        c = c[in_mask(m4, c)]
        parts.append(c)
        tot += len(c)
    B = np.vstack(parts)
    B = B[rng.permutation(len(B))[: n - 1]]
    if len(B) == 0:
        return np.zeros(len(CENT_A)), np.zeros(len(CENT_A))
    tb = cKDTree(B)
    ta = cKDTree(xy_i)
    dd, _ = tb.query(xy_i, k=min(6, len(B)))
    dd = np.asarray(dd).reshape(len(xy_i), -1) * APIX
    dd[dd < G["r_hc_px"] * APIX] = np.inf
    d = dd.min(axis=1)
    d = d[np.isfinite(d)]
    return nn_hist(d), pair_counts(ta, tb, False)


def jitter_null_one(i, rng):
    """Local-displacement null for micrograph i: each real dot moved uniformly within a 200 A disk (rejected moves
    outside the valid mask or inside the 20 A hard core are redrawn), so density variations on scales >~ 0.5 um
    (thickness gradients, aggregates' surroundings) are kept but structure at the 40-120 A scale is destroyed."""
    xy0 = G["xy"][i]
    m4 = G["masks"][i]
    n = len(xy0)
    rmax = JITTER_A / APIX
    cur = np.full((n, 2), np.nan)
    todo = np.arange(n)
    for _ in range(20):
        if len(todo) == 0:
            break
        rr = rmax * np.sqrt(rng.random(len(todo)))
        th = rng.random(len(todo)) * 2 * np.pi
        prop = xy0[todo] + np.c_[rr * np.cos(th), rr * np.sin(th)]
        okb = (prop[:, 0] >= 0) & (prop[:, 0] < 2048) & (prop[:, 1] >= 0) & (prop[:, 1] < 2048)
        okb[okb] = in_mask(m4, prop[okb])
        cur[todo[okb]] = prop[okb]
        have = np.nonzero(np.isfinite(cur[:, 0]))[0]
        pri = rng.random(len(have)) + 2.0 * ~np.isin(have, todo)       # previously placed dots win conflicts
        keep = have[greedy_nms(cur[have], pri, G["r_hc_px"])]
        drop = np.setdiff1d(have, keep)
        cur[drop] = np.nan
        todo = np.nonzero(~np.isfinite(cur[:, 0]))[0]
    return obs_stats(cur[np.isfinite(cur[:, 0])])[:2]


def null_worker(args):
    rep, seed = args
    rng = np.random.default_rng(seed)
    nb = len(CENT_A)
    out = np.zeros((6, nb))
    for i in G["active"]:
        for k, fn in enumerate((perm_null_one, toro_null_one, jitter_null_one)):
            h, pcs = fn(i, rng)
            out[2 * k] += h
            out[2 * k + 1] += pcs
    return out


# ----------------------------------------------------------------------------------------------- figures
def fig_examples(examples, rows, path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import Circle
    col = {"bare": "#1f77b4", "one_arm": "#ff7f0e", "two_arm_open": "#d62728", "two_arm_closed": "#ff00ff"}
    n = len(examples)
    fig, axs = plt.subplots(2, n, figsize=(5.2 * n, 10.8), dpi=120)
    axs = np.atleast_2d(axs).reshape(2, n)
    for c, (ex, row) in enumerate(zip(examples, rows)):
        x0, y0 = ex["crop"]
        show = -ex["pos"][y0:y0 + 600, x0:x0 + 600]          # phase-flipped, displayed protein-DARK
        lo, hi = np.percentile(show, [1, 99])
        for r in range(2):
            ax = axs[r, c]
            ax.imshow(show, cmap="gray", vmin=lo, vmax=hi, interpolation="nearest")
            ax.set_xticks([])
            ax.set_yticks([])
            m = ex["m4"][y0 // BIN:(y0 + 600) // BIN, x0 // BIN:(x0 + 600) // BIN]
            if (~m).any():
                ax.contour(np.kron(~m, np.ones((BIN, BIN))), levels=[0.5], colors="yellow", linewidths=0.8)
        sel = (ex["x"] >= x0) & (ex["x"] < x0 + 600) & (ex["y"] >= y0) & (ex["y"] < y0 + 600)
        for xx, yy, cl in zip(ex["x"][sel], ex["y"][sel], ex["cls"][sel]):
            axs[1, c].add_patch(Circle((xx - x0, yy - y0), 14, fill=False, ec=col[cl], lw=1.3))
        axs[0, c].set_title(f"{row['short']}\ndf {row['df_um']:.2f} um, fit {row['fit_A']:.1f} A, "
                            f"{row['n_dots']} dots/mic", fontsize=9)
        axs[1, c].set_title(f"crop x{x0}-{x0 + 600}, y{y0}-{y0 + 600}: {int(sel.sum())} dots", fontsize=9)
        for r in range(2):
            axs[r, c].set_xlim(-0.5, 599.5)
            axs[r, c].set_ylim(599.5, -0.5)
        axs[0, c].plot([20, 20 + 100 / APIX], [580, 580], "w-", lw=3)
        axs[0, c].text(20, 565, "10 nm", color="w", fontsize=8)
    handles = [plt.Line2D([], [], marker="o", ls="", mfc="none", mec=v, label=k) for k, v in col.items()]
    fig.legend(handles=handles, loc="lower center", ncol=4, fontsize=10, frameon=False)
    fig.suptitle("G0 dot detector on pilot micrographs: CTF phase-flipped, 15 A low-pass, protein dark; "
                 "circles = detected dots (r 20 A) coloured by arm class; yellow = masked contamination", fontsize=11)
    fig.tight_layout(rect=(0, 0.03, 1, 0.97))
    fig.savefig(path)
    plt.close(fig)


def fig_nn(res, path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, axs = plt.subplots(3, 1, figsize=(9, 11), dpi=120, sharex=True)
    x = CENT_A
    nn, g = res["nn"], res["g"]
    ax = axs[0]
    ax.step(x, nn["obs"], where="mid", color="k", lw=1.5, label="observed")
    for key, c, lab in [("perm", "#1f77b4", "permutation null"), ("toro", "#d62728", "toroidal-shift null"),
                        ("jitter", "#2ca02c", "200 A local-jitter null")]:
        ax.plot(x, nn[key]["mean"], color=c, lw=1, label=f"{lab} mean")
        ax.fill_between(x, nn[key]["lo"], nn[key]["hi"], color=c, alpha=0.25, lw=0, label=f"{lab} 95% band")
    ax.set_ylabel("dots per 5 A bin (NN distance)")
    ax.legend(fontsize=8)
    ax = axs[1]
    ax.plot(x, g["obs"], "k-", lw=1.5, label="observed g(r)")
    for key, c in [("perm", "#1f77b4"), ("toro", "#d62728"), ("jitter", "#2ca02c")]:
        ax.plot(x, g[key]["mean"], color=c, lw=1)
        ax.fill_between(x, g[key]["lo"], g[key]["hi"], color=c, alpha=0.25, lw=0)
    ax.axhline(1, color="0.5", lw=0.5)
    ax.set_ylabel("pair correlation g(r)")
    ax.set_ylim(0, max(1.6, np.nanmax(g["obs"][x > 25]) * 1.1))
    ax.legend(fontsize=8)
    ax = axs[2]
    for key, c in [("perm", "#1f77b4"), ("toro", "#d62728"), ("jitter", "#2ca02c")]:
        m = np.where(g[key]["mean"] > 0, g[key]["mean"], np.nan)
        ax.plot(x, g["obs"] / m, color=c, lw=1.2, label=f"g_obs / g_{key} null")
        ax.fill_between(x, g[key]["lo"] / m, g[key]["hi"] / m, color=c, alpha=0.2, lw=0)
    ax.axhline(1, color="0.5", lw=0.5)
    ax.set_ylim(0.5, 1.6)
    ax.set_ylabel("observed / null")
    ax.set_xlabel("distance (A)")
    ax.legend(fontsize=8)
    for a in axs:
        for lo, hi in WINDOWS.values():
            a.axvspan(lo, hi, color="0.6", alpha=0.15)
        a.axvline(res["r_hc_A"], color="0.4", ls=":", lw=0.8)
    axs[0].set_title(f"G0 pilot: {res['n_dots']} dots on {res['n_mics']} micrographs; nulls from real coordinates "
                     f"only ({res['n_rep']} reps)\ndotted = 20 A hard core; grey = 35-50 / 110-130 A pre-set windows "
                     f"(+ 30-40 A post hoc)", fontsize=9)
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)


# ----------------------------------------------------------------------------------------------- main
def band(arr):
    return dict(mean=arr.mean(0), lo=np.percentile(arr, 2.5, axis=0), hi=np.percentile(arr, 97.5, axis=0),
                sd=arr.std(0))


def window_test(obs, null, lo, hi, n_dots):
    sel = (CENT_A > lo) & (CENT_A < hi)
    o = float(obs[sel].sum())
    ns = null[:, sel].sum(1)
    mu, sd = float(ns.mean()), float(ns.std())
    return dict(window_A=[lo, hi], obs=o, null_mean=mu, null_lo95=float(np.percentile(ns, 2.5)),
                null_hi95=float(np.percentile(ns, 97.5)), z=(o - mu) / sd if sd > 0 else None,
                p_excess=float((1 + np.sum(ns >= o)) / (1 + len(ns))),
                ratio=o / mu if mu > 0 else None, excess=o - mu,
                excess_fraction_of_dots=(o - mu) / n_dots if n_dots else None)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, default=0, help="first N micrographs only (test)")
    ap.add_argument("--k", type=float, default=4.5, help="dot threshold: median + k*MAD of the LoG response")
    ap.add_argument("--ridge-thr", type=float, default=2.5,
                    help="arm threshold: mean ridge (noise units) along the best 25-60 A path")
    ap.add_argument("--nrep", type=int, default=200)
    ap.add_argument("--procs", type=int, default=int(os.environ.get("SLURM_CPUS_PER_TASK", 4)))
    ap.add_argument("--results", default=os.path.join(REPO, "results"))
    ap.add_argument("--figs", default=os.path.join(REPO, "figures", "qc"))
    ap.add_argument("--seed", type=int, default=20261113)
    a = ap.parse_args()
    t0 = time.time()
    d = np.load(CTF_NPZ, allow_pickle=False)
    rel = [str(p) for p in d["micrograph_blob__path"]]
    df1, df2, ang = d["ctf__df1_A"].astype(float), d["ctf__df2_A"].astype(float), d["ctf__df_angle_rad"].astype(float)
    fit, ice = d["ctf__ctf_fit_to_A"].astype(float), d["ctf_stats__ice_thickness_rel"].astype(float)
    nm = len(rel) if a.n == 0 else a.n
    dfm = 0.5 * (df1 + df2) / 1e4
    dfbin = np.digitize(dfm, DF_BINS_UM[1:-1])
    # example micrographs for the figure: defocus quartile medians among micrographs with CTF fit <= 6 A
    good = [i for i in range(nm) if fit[i] <= 6.0] or list(range(nm))
    q = np.argsort(dfm[good])
    ex_ids = sorted({good[q[int(f * (len(q) - 1))]] for f in (0.15, 0.4, 0.65, 0.9)})
    jobs = [(i, os.path.join(PROJ, rel[i]), df1[i], df2[i], ang[i], a.k, a.ridge_thr, i in ex_ids, a.seed + i)
            for i in range(nm)]
    print(f"{nm} micrographs, k={a.k}, ridge_thr={a.ridge_thr}, procs={a.procs}, examples={ex_ids}", flush=True)
    with Pool(a.procs) as pool:
        res = []
        for r in pool.imap(process, jobs, chunksize=1):
            res.append(r)
            if len(res) % 25 == 0 or len(res) == nm:
                print(f"  processed {len(res)}/{nm} ({time.time() - t0:.0f} s)", flush=True)
    res.sort(key=lambda r: r["i"])
    N = np.array([len(r["x"]) for r in res])
    xy = [np.c_[r["x"], r["y"]].astype(float) for r in res]
    area_um2 = np.array([r["area_px"] * APIX ** 2 / 1e8 for r in res])
    n_dots = int(N.sum())
    print(f"dots: {n_dots} total, per micrograph median {np.median(N):.0f} (IQR {np.percentile(N, 25):.0f}-"
          f"{np.percentile(N, 75):.0f}); mirror false-peak ratio {sum(r['n_mirror'] for r in res) / max(n_dots, 1):.3f}; "
          f"flip-vs-raw corr median {np.median([r['cc_flip_vs_raw'] for r in res]):.2f}", flush=True)

    # observed statistics
    nn_obs = np.zeros(len(CENT_A))
    pc_obs = np.zeros(len(CENT_A))
    nn_all = []
    for i in range(nm):
        h, p, dd = obs_stats(xy[i])
        nn_obs += h
        pc_obs += p
        nn_all.append(dd)
    nn_all = np.concatenate(nn_all) if nn_all else np.zeros(0)
    # g(r) denominator: sum_i N_i(N_i-1)/A_i^2 * gamma_i(r) * 2 pi r dr   (ordered pairs, px units)
    dr = R_BIN_A / APIX
    den = np.zeros(len(CENT_A))
    for i in range(nm):
        if N[i] > 1:
            A = res[i]["area_px"]
            den += N[i] * (N[i] - 1) / A ** 2 * res[i]["gamma"] * 2 * np.pi * (CENT_A / APIX) * dr
    den[den <= 0] = np.nan

    # nulls
    G.update(N=N, xy=xy, masks=[r["m4"] for r in res], dfbin=dfbin[:nm],
             bin_members={b: np.nonzero(dfbin[:nm] == b)[0] for b in range(len(DF_BINS_UM) - 1)},
             pool_xy=np.vstack([v for v in xy if len(v)]),
             pool_mic=np.concatenate([np.full(len(v), i) for i, v in enumerate(xy) if len(v)]),
             r_hc_px=NMS_A / APIX, active=[i for i in range(nm) if N[i] >= 2])
    print(f"running {a.nrep} null replicates ...", flush=True)
    with Pool(a.procs) as pool:
        nulls = np.array(pool.map(null_worker, [(r, a.seed * 7 + r) for r in range(a.nrep)], chunksize=1))
    print(f"  nulls done ({time.time() - t0:.0f} s)", flush=True)
    nn_perm, pc_perm, nn_toro, pc_toro, nn_jit, pc_jit = (nulls[:, k] for k in range(6))
    g_obs = pc_obs / den
    out_nn = dict(obs=nn_obs, perm=band(nn_perm), toro=band(nn_toro), jitter=band(nn_jit))
    out_g = dict(obs=g_obs, perm=band(pc_perm / den), toro=band(pc_toro / den), jitter=band(pc_jit / den))

    tests = {}
    for name, (lo, hi) in WINDOWS.items():
        tests[name] = dict(
            nn_vs_perm=window_test(nn_obs, nn_perm, lo, hi, n_dots),
            nn_vs_toro=window_test(nn_obs, nn_toro, lo, hi, n_dots),
            pairs_vs_perm=window_test(pc_obs, pc_perm, lo, hi, n_dots),
            pairs_vs_toro=window_test(pc_obs, pc_toro, lo, hi, n_dots),
            nn_vs_jitter=window_test(nn_obs, nn_jit, lo, hi, n_dots),
            pairs_vs_jitter=window_test(pc_obs, pc_jit, lo, hi, n_dots))
        # ordered pair counts: excess ordered pairs / N = fraction of dots that have an extra partner there
    # local-peak check: is the window a bump (above both flanks) rather than part of a broad trend?
    for name, (lo, hi) in WINDOWS.items():
        w = (CENT_A > lo) & (CENT_A < hi)
        fl = ((CENT_A > lo - 20) & (CENT_A <= lo)) | ((CENT_A >= hi) & (CENT_A < hi + 20))
        for key in ("perm", "toro", "jitter"):
            ratio = g_obs / out_g[key]["mean"]
            tests[name][f"g_ratio_window_vs_flanks_{key}"] = [float(np.nanmean(ratio[w])), float(np.nanmean(ratio[fl]))]

    per_bin = []
    for b, c in enumerate(CENT_A):
        if c < 20 or c > 160:
            continue
        row = dict(r_A=float(c), nn_obs=float(nn_obs[b]), g_obs=float(g_obs[b]))
        for key in ("perm", "toro", "jitter"):
            row[f"nn_{key}"] = [float(out_nn[key][q][b]) for q in ("mean", "lo", "hi")]
            row[f"g_{key}"] = [float(out_g[key][q][b]) for q in ("mean", "lo", "hi")]
        per_bin.append(row)

    # arms
    cls_all = np.concatenate([r["cls"] for r in res]) if n_dots else np.zeros(0, object)
    ctrl_all = np.concatenate([r["ctrl_cls"] for r in res])
    arm_frac = {c: float(np.mean(cls_all == c)) for c in CLASSES}
    ctrl_frac = {c: float(np.mean(ctrl_all == c)) if len(ctrl_all) else None for c in CLASSES}
    arm_by_bin = {}
    for b in range(len(DF_BINS_UM) - 1):
        cc = np.concatenate([res[i]["cls"] for i in range(nm) if dfbin[i] == b] or [np.zeros(0, object)])
        kk = np.concatenate([res[i]["ctrl_cls"] for i in range(nm) if dfbin[i] == b] or [np.zeros(0, object)])
        lab = f"{DF_BINS_UM[b]:.1f}-{DF_BINS_UM[b + 1]:.1f}um"
        arm_by_bin[lab] = dict(n_dots=int(len(cc)),
                               dots={c: float(np.mean(cc == c)) if len(cc) else None for c in CLASSES},
                               random_positions={c: float(np.mean(kk == c)) if len(kk) else None for c in CLASSES})

    # density by defocus bin
    per_mic = []
    for i in range(nm):
        per_mic.append(dict(path=rel[i], df_um=round(float(dfm[i]), 4), df_bin=int(dfbin[i]),
                            ctf_fit_A=round(float(fit[i]), 2), ice_rel=round(float(ice[i]), 4), n_dots=int(N[i]),
                            valid_area_um2=round(float(area_um2[i]), 5),
                            dots_per_um2=round(float(N[i] / area_um2[i]), 1) if area_um2[i] > 0 else None,
                            frac_contam_masked=round(res[i]["frac_contam"], 4), n_mirror_peaks=res[i]["n_mirror"],
                            flip_vs_raw_corr=round(res[i]["cc_flip_vs_raw"], 3)))
    by_bin = {}
    for b in range(len(DF_BINS_UM) - 1):
        s = [p for p in per_mic if p["df_bin"] == b]
        lab = f"{DF_BINS_UM[b]:.1f}-{DF_BINS_UM[b + 1]:.1f}um"
        if not s:
            by_bin[lab] = dict(n_mics=0)
            continue
        nn_ = np.array([p["n_dots"] for p in s])
        dd_ = np.array([p["dots_per_um2"] for p in s], float)
        by_bin[lab] = dict(n_mics=len(s), dots_per_mic_median=float(np.median(nn_)),
                           dots_per_mic_iqr=[float(np.percentile(nn_, 25)), float(np.percentile(nn_, 75))],
                           dots_per_um2_median=float(np.nanmedian(dd_)),
                           median_ctf_fit_A=float(np.median([p["ctf_fit_A"] for p in s])),
                           mirror_false_peak_ratio=float(sum(p["n_mirror_peaks"] for p in s) / max(nn_.sum(), 1)))

    # outputs
    os.makedirs(a.results, exist_ok=True)
    os.makedirs(a.figs, exist_ok=True)
    rows = ["micrograph_path,df_um,x_px,y_px,x_A,y_A,log_z,log_resp,nn_A,arm_count,arm_class"]
    for i in range(nm):
        r = res[i]
        if len(r["x"]) == 0:
            continue
        dd = cKDTree(xy[i]).query(xy[i], k=2)[0][:, 1] * APIX if len(r["x"]) > 1 else np.full(len(r["x"]), np.nan)
        for j in range(len(r["x"])):
            rows.append(f"{rel[i]},{dfm[i]:.4f},{r['x'][j]},{r['y'][j]},{r['x'][j] * APIX:.1f},"
                        f"{r['y'][j] * APIX:.1f},{r['z'][j]:.3f},{r['resp'][j]:.4f},{dd[j]:.1f},{r['arms'][j]},"
                        f"{r['cls'][j]}")
    csv_path = os.path.join(a.results, "G0_dots_pilot.csv.gz")
    with gzip.open(csv_path, "wt") as fh:
        fh.write("\n".join(rows) + "\n")
    if os.path.getsize(csv_path) > 20e6:
        alt = os.path.join(EXPORTS, "G0_dots_pilot.csv.gz")
        os.replace(csv_path, alt)
        csv_path = alt

    def tolist(o):
        if isinstance(o, dict):
            return {k: tolist(v) for k, v in o.items()}
        if isinstance(o, np.ndarray):
            return [None if not np.isfinite(v) else round(float(v), 5) for v in o.ravel()]
        if isinstance(o, (np.floating, np.integer)):
            return o.item()
        return o

    summary = dict(
        created=time.strftime("%Y-%m-%d %H:%M"), script="scripts/12_dot_spatial_stats.py",
        inputs=dict(micrographs=f"{PROJ}/J3/motioncorrected/*_patch_aligned_doseweighted.mrc", ctf=CTF_NPZ,
                    n_micrographs=nm),
        settings=dict(apix_A=APIX, kv=KV, cs_mm=CS_MM, amp_contrast=AMP,
                      phase_flip="sign(-sin(chi+asin(0.1))), own astigmatic CTF; output protein-positive",
                      lowpass="Gaussian, half-max at 1/15 A^-1", dot_detector="-sigma^2 LoG, sigma 12 A (8.45 px)",
                      threshold_k_MAD=a.k, threshold="median + k*1.4826*MAD of LoG response over valid area",
                      nms_A=NMS_A, edge_px=EDGE_PX,
                      contamination_mask="bin-4 local (68 A Gaussian) mean |robust z|>5 or RMS robust z>5, dilated 57 A",
                      ridge="max(0,-l1-|l2|) Hessian, sigma 7 A, in units of the robust SD of d2/dx2", arm_thr=a.ridge_thr,
                      arm_rings_A=[ARM_R_A[0], ARM_R_A[1], ARM_DR_A], arm_theta_bins=ARM_NTH,
                      arm_rules=dict(method="DP ridge path from r=25 A, <=1 theta bin per 2.5 A ring, lateral penalty",
                                     short_path_to_A=ARM_SHORT_A, lateral_penalty=ARM_PEN, arm_nms_deg=ARM_NMS_DEG,
                                     closed_meet_beyond_A=ARM_CLOSE_MIN_R_A, closed_meet_tol_deg=10,
                                     closed_min_start_sep_deg=ARM_MIN_SEP_DEG, dot_exclusion_A=DOT_EXCL_A),
                      jitter_null_radius_A=JITTER_A,
                      n_null_reps=a.nrep, df_bins_um=DF_BINS_UM, hard_core_in_nulls_A=NMS_A),
        dots=dict(total=n_dots, per_mic_median=float(np.median(N)), per_mic_iqr=[float(np.percentile(N, 25)),
                  float(np.percentile(N, 75))], per_mic_min_max=[int(N.min()), int(N.max())],
                  per_um2_median=float(np.median(N / area_um2)),
                  mirror_false_peak_ratio=float(sum(r["n_mirror"] for r in res) / max(n_dots, 1)),
                  median_masked_contam_frac=float(np.median([r["frac_contam"] for r in res])),
                  nn_median_A=float(np.median(nn_all)) if len(nn_all) else None,
                  flip_vs_raw_corr_median=float(np.median([r["cc_flip_vs_raw"] for r in res]))),
        density_by_defocus=by_bin,
        bins_A=dict(edges=EDGES_A, centres=CENT_A),
        nn_hist=out_nn, g_r=out_g, pair_counts=dict(obs=pc_obs, perm_mean=pc_perm.mean(0), toro_mean=pc_toro.mean(0),
                                                    jitter_mean=pc_jit.mean(0)),
        window_tests=tests, per_bin_20_160A=per_bin,
        arms=dict(dots=arm_frac, random_background_positions=ctrl_frac, n_dots=int(len(cls_all)),
                  n_random=int(len(ctrl_all)), by_defocus=arm_by_bin),
        per_micrograph=per_mic, dots_csv=csv_path, example_micrographs=[rel[i] for i in ex_ids],
        runtime_s=round(time.time() - t0, 1))
    with open(os.path.join(a.results, "G0_spatial.json"), "w") as fh:
        json.dump(tolist(summary), fh, indent=1)

    fig_nn(dict(nn=out_nn, g=out_g, n_dots=n_dots, n_mics=nm, n_rep=a.nrep, r_hc_A=NMS_A),
           os.path.join(a.figs, "G0_nn_hist.png"))
    exs, rws = [], []
    for i in ex_ids:
        r = res[i]
        # crop: 600x600 inside the edge, centred on the densest valid region
        best, bxy = -1, (EDGE_PX, EDGE_PX)
        for y0 in range(EDGE_PX, 2048 - 600 - EDGE_PX + 1, 64):
            for x0 in range(EDGE_PX, 2048 - 600 - EDGE_PX + 1, 64):
                s = int(((r["x"] >= x0) & (r["x"] < x0 + 600) & (r["y"] >= y0) & (r["y"] < y0 + 600)).sum())
                if s > best:
                    best, bxy = s, (x0, y0)
        exs.append(dict(pos=r["pos"], m4=r["m4"], x=r["x"], y=r["y"], cls=r["cls"], crop=bxy))
        rws.append(dict(short="FoilHole_" + os.path.basename(rel[i]).split("FoilHole_")[-1].split("_Fractions")[0], df_um=dfm[i], fit_A=fit[i],
                        n_dots=int(N[i])))
    fig_examples(exs, rws, os.path.join(a.figs, "G0_dot_detection_examples.png"))
    print(json.dumps(tolist(dict(dots=summary["dots"], window_tests=tests, arms=summary["arms"]["dots"],
                                 arms_random=ctrl_frac)), indent=1))
    print(f"done in {time.time() - t0:.0f} s", flush=True)


if __name__ == "__main__":
    main()
