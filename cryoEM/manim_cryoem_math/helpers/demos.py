"""Precomputed, deterministic data for each clip (wrapped by cache.cached).

Every number a clip displays (SNR points, loss curves, PCA scores, FSC
curves, local-resolution values) is computed here from the synthetic
phantom; nothing is a hand-drawn curve. Each function takes only plain
parameters so the cache key is stable.
"""
from __future__ import annotations

import numpy as np
from scipy import ndimage

from . import fourier as fo
from . import optics as op
from . import projections as pj
from . import synthetic_phantom as sp

SEED = sp.SEED

# a canonical, asymmetric viewing direction used whenever "the" view is needed
VIEW = dict(rot=30.0, tilt=55.0, psi=15.0)
# the side-on display direction for 3D renders
DISPLAY = dict(rot=-20.0, tilt=62.0, psi=0.0)


def view_R():
    return pj.rotation_from_euler(**VIEW)


def display_R():
    return pj.rotation_from_euler(**DISPLAY)


# ------------------------------------------------------------------ clip 01
def image_formation(n: int = 128, voxel: float = 1.5,
                    defoci=(0.6, 1.0, 1.6, 2.4), snrs=(4.0, 1.0, 0.5, 0.25, 0.1, 0.05)) -> dict:
    """Clean projection, CTF-modulated versions, a CTF curve per defocus and a
    noise ladder (one fixed white-noise pattern scaled to each SNR)."""
    rng = np.random.default_rng(SEED + 1)
    R = view_R()
    clean = pj.project_analytic(R, "A", n, voxel)
    k = np.linspace(0, 1 / (2 * voxel), 600)
    out = dict(clean=clean, k=k, defoci=np.array(defoci), snrs=np.array(snrs), voxel=voxel)
    out["ctf_curves"] = np.stack([op.ctf_1d(k, d) for d in defoci])
    out["ctf_images"] = np.stack([op.apply_ctf(clean, voxel, defocus_um=d) for d in defoci])
    out["ctf2d"] = np.stack([op.ctf_2d(n, voxel, shifted=True, defocus_um=d) for d in defoci])
    out["first_zero"] = np.array([op.first_zero(d) for d in defoci])
    white = rng.standard_normal((n, n))
    ref = out["ctf_images"][1]                       # 1.0 um is the working defocus
    sig = np.array([op.sigma_for_snr(ref, s) for s in snrs])
    out["noise_sigma"] = sig
    out["noisy"] = np.stack([ref + s * white for s in sig])
    out["signal_rms"] = np.std(ref)
    return out


# ------------------------------------------------------------------ clip 02
def alignment(n: int = 96, voxel: float = 1.6, n_particles: int = 256, snr: float = 0.3,
              shift_sd: float = 6.0, ns=(1, 4, 16, 64, 256), jitter_deg: float = 1.5,
              jitter_A: float = 0.4) -> dict:
    """2D alignment and averaging of one view with random in-plane rotations,
    shifts and i.i.d. Gaussian noise (identical underlying 2D signal)."""
    rng = np.random.default_rng(SEED + 2)
    S = pj.project_analytic(view_R(), "A", n, voxel)
    sigma = op.sigma_for_snr(S, snr)
    phi = rng.uniform(0, 360, n_particles)
    t = np.clip(rng.normal(0, shift_sd, (n_particles, 2)), -15, 15)
    phi_est = phi + rng.normal(0, jitter_deg, n_particles)
    t_est = t + rng.normal(0, jitter_A, (n_particles, 2))
    imgs, aligned, aligned_noise = [], [], []
    for i in range(n_particles):
        clean_i = pj.shift_image(pj.rotate_image_exact(S, phi[i]), t[i], voxel)
        noise_i = sigma * rng.standard_normal((n, n))
        img = clean_i + noise_i
        undo = lambda x: pj.rotate_image_exact(pj.shift_image(x, -t_est[i], voxel), -phi_est[i])
        imgs.append(img)
        aligned.append(undo(img))
        aligned_noise.append(undo(noise_i))
    imgs, aligned, aligned_noise = map(np.array, (imgs, aligned, aligned_noise))
    ns = np.array(ns)
    out = dict(signal=S, sigma=sigma, phi=phi, t=t, phi_est=phi_est, t_est=t_est, ns=ns,
               images=imgs[:12], aligned_examples=aligned[:12], voxel=voxel, snr1=snr)
    out["unaligned_mean"] = np.stack([imgs[:m].mean(0) for m in ns])
    out["aligned_mean"] = np.stack([aligned[:m].mean(0) for m in ns])
    # noise std of the mean, measured on the transformed noise alone (inside the
    # circular field of view, where the rotation does not pad with zeros)
    yy, xx = np.indices((n, n)) - n // 2
    fov = (xx ** 2 + yy ** 2) < (0.42 * n) ** 2
    out["sigma_mean_measured"] = np.array([aligned_noise[:m].mean(0)[fov].std() for m in ns])
    out["sigma_mean_theory"] = out["sigma_mean_measured"][0] / np.sqrt(ns)
    out["snr_measured"] = np.std(S) / out["sigma_mean_measured"]
    out["snr_theory"] = out["snr_measured"][0] * np.sqrt(ns)
    return out


# ------------------------------------------------------------------ clip 03
def reconstruction(n: int = 64, voxel: float = 2.0, counts=(2, 6, 20, 80, 300),
                   biased_count: int = 300, biased_tilt: float = 12.0) -> dict:
    """Ideal projections with KNOWN poses -> direct Fourier reconstructions."""
    rng = np.random.default_rng(SEED + 3)
    truth = sp.ground_truth(n, voxel, atoms=False)
    show_R = np.stack([pj.rotation_from_euler(*a) for a in
                       [(0, 0, 0), (40, 70, 20), (-60, 110, 75), (120, 35, -40)]])
    out = dict(truth=truth.astype(np.float32), show_R=show_R, voxel=voxel,
               counts=np.array(counts))
    out["show_proj"] = np.stack([pj.project_analytic(R, "A", n, voxel, atoms=False)
                                 for R in show_R])
    out["show_ft"] = np.stack([fo.log_magnitude(fo.fft2c(p)) for p in out["show_proj"]])
    # Fourier-slice check for display: FT2(projection) vs central slice of FT3
    F3 = fo.fft3c(fo.pad_volume(truth, 2))
    R = show_R[1]
    p = pj.project_analytic(R, "A", n, voxel, atoms=False) / voxel   # line integral in voxels
    P2 = fo.fft2c(np.pad(p, n // 2))
    out["slice_ft2"] = fo.log_magnitude(P2)
    out["slice_central"] = fo.log_magnitude(fo.central_slice(F3, R))
    all_R = pj.random_rotations(max(counts), rng)
    projs = np.stack([pj.project_analytic(R, "A", n, voxel, atoms=False) / voxel for R in all_R])
    out["rotations"] = all_R
    vols = [fo.reconstruct_direct_fourier(projs[:c], all_R[:c]) for c in counts]
    out["recons"] = np.stack(vols).astype(np.float32)
    bR = pj.biased_rotations(biased_count, rng, biased_tilt)
    bprojs = np.stack([pj.project_analytic(R, "A", n, voxel, atoms=False) / voxel for R in bR])
    out["biased_rotations"] = bR
    out["biased_recon"] = fo.reconstruct_direct_fourier(bprojs, bR).astype(np.float32)
    out["recon_corr"] = np.array([np.corrcoef(v.ravel(), truth.ravel())[0, 1] for v in vols])
    out["biased_corr"] = np.corrcoef(out["biased_recon"].ravel(), truth.ravel())[0, 1]
    # side-on projections (perpendicular to the biased axis) of every map
    side = pj.rotation_from_euler(0, 90, 0)
    out["side_truth"] = pj.project_voxel(truth, side)
    out["side_recons"] = np.stack([pj.project_voxel(v, side) for v in vols])
    out["side_biased"] = pj.project_voxel(out["biased_recon"].astype(float), side)
    cov_u = fo.sampling_coverage(all_R[:max(counts)], 40)
    cov_b = fo.sampling_coverage(bR, 40)
    out["coverage_counts"] = np.array([fo.sampling_coverage(all_R[:c], 40).mean() for c in counts])
    out["coverage_uniform"] = cov_u.mean()
    out["coverage_biased"] = cov_b.mean()
    return out


# ------------------------------------------------------------------ clip 04
def _small_rotations(step_deg: float) -> list[np.ndarray]:
    d = np.radians(step_deg)
    rs = [np.eye(3)]
    for f in (pj.rot_x, pj.rot_y, pj.rot_z):
        rs += [f(d), f(-d)]
    for f, g in ((pj.rot_x, pj.rot_y), (pj.rot_y, pj.rot_z), (pj.rot_x, pj.rot_z)):
        rs += [f(d) @ g(d), f(-d) @ g(-d), f(d) @ g(-d)]
    return rs


def _angle_between(Ra, Rb) -> float:
    c = (np.trace(Ra @ Rb.T) - 1) / 2
    return float(np.degrees(np.arccos(np.clip(c, -1, 1))))


def _random_axis_rotation(angle_deg, rng):
    ax = rng.standard_normal(3)
    ax /= np.linalg.norm(ax)
    a = np.radians(angle_deg)
    K = np.array([[0, -ax[2], ax[1]], [ax[2], 0, -ax[0]], [-ax[1], ax[0], 0]])
    return np.eye(3) + np.sin(a) * K + (1 - np.cos(a)) * K @ K


def refinement(n: int = 48, voxel: float = 2.5, n_particles: int = 120, snr: float = 1.0,
               init_err_deg: float = 20.0, iters: int = 8, smooth_sigma: float = 0.8,
               lam: float = 2.0, steps=(8.0, 6.0, 4.0, 3.0, 2.0, 1.5, 1.0, 1.0)) -> dict:
    """Toy projection matching: alternate pose search (rotation + shift) and a
    regularised direct-Fourier volume update. Pedagogical, not CryoSPARC.

    Loss L = sum_i ||I_i - T_{t_i} P_{R_i} rho||^2 + lam * ||grad rho||^2
    (real space; equivalent to the Fourier form by Parseval, CTF omitted).
    """
    rng = np.random.default_rng(SEED + 4)
    truth = sp.ground_truth(n, voxel, atoms=False)
    R_true = pj.random_rotations(n_particles, rng)
    t_true = np.clip(rng.normal(0, 2.5, (n_particles, 2)), -6, 6)
    # images are line integrals in voxel units, matching project_fourier / reconstruction
    clean = np.stack([pj.project_analytic(R, "A", n, voxel, atoms=False, shift=t) / voxel
                      for R, t in zip(R_true, t_true)])
    sigma = op.sigma_for_snr(clean[0], snr)
    images = clean + sigma * rng.standard_normal(clean.shape)
    R_est = np.stack([R @ _random_axis_rotation(init_err_deg, rng) for R in R_true])
    t_est = np.zeros_like(t_true)

    def grad_penalty(v):
        return float(sum(np.sum(np.diff(v, axis=a) ** 2) for a in range(3)))

    def update_volume(R_set, t_set):
        unshifted = np.stack([pj.shift_image(im, -t, voxel) for im, t in zip(images, t_set)])
        v = fo.reconstruct_direct_fourier(unshifted, R_set)
        return ndimage.gaussian_filter(v, smooth_sigma)   # smoothness prior (regulariser)

    def residuals(vol, R_set, t_set):
        F3 = fo.fft3c(fo.pad_volume(vol, 2))
        res = []
        for im, R, t in zip(images, R_set, t_set):
            pr = pj.shift_image(fo.project_fourier(F3, R, n), t, voxel)
            res.append(im - pr)
        return np.array(res)

    def loss(vol, R_set, t_set):
        r = residuals(vol, R_set, t_set)
        return float(np.sum(r ** 2)), lam * grad_penalty(vol)

    vol = ndimage.gaussian_filter(update_volume(R_est, t_est), 2.5)   # crude start
    track = 0
    hist = dict(loss_data=[], loss_reg=[], label=[], pose_err=[], map_corr=[])
    track_proj, track_resid, vols = [], [], [vol.copy()]

    def record(tag):
        d, r = loss(vol, R_est, t_est)
        hist["loss_data"].append(d)
        hist["loss_reg"].append(r)
        hist["label"].append(tag)
        hist["pose_err"].append(np.mean([_angle_between(a, b) for a, b in zip(R_est, R_true)]))
        hist["map_corr"].append(np.corrcoef(vol.ravel(), truth.ravel())[0, 1])
        F3 = fo.fft3c(fo.pad_volume(vol, 2))
        pr = pj.shift_image(fo.project_fourier(F3, R_est[track], n), t_est[track], voxel)
        track_proj.append(pr)
        track_resid.append(images[track] - pr)

    record(0)   # tag 0 = start, 1 = after pose step, 2 = after volume step
    for it in range(iters):
        F3 = fo.fft3c(fo.pad_volume(vol, 2))
        step = steps[min(it, len(steps) - 1)]
        cands = _small_rotations(step)
        for i in range(n_particles):
            best = (np.inf, None, None)
            FI = np.fft.fft2(images[i])
            for dR in cands:
                Rc = R_est[i] @ dR
                pr = fo.project_fourier(F3, Rc, n)
                cc = np.real(np.fft.ifft2(FI * np.conj(np.fft.fft2(pr))))
                ix, iy = np.unravel_index(np.argmax(cc), cc.shape)
                sx = ((ix + n // 2) % n - n // 2) * voxel
                sy = ((iy + n // 2) % n - n // 2) * voxel
                if abs(sx) > 8 or abs(sy) > 8:
                    sx, sy = t_est[i]
                r = np.sum((images[i] - pj.shift_image(pr, (sx, sy), voxel)) ** 2)
                if r < best[0]:
                    best = (r, Rc, (sx, sy))
            R_est[i], t_est[i] = best[1], best[2]
        record(1)
        vol = update_volume(R_est, t_est)
        vols.append(vol.copy())
        record(2)
    out = {k: np.array(v) for k, v in hist.items()}
    out.update(truth=truth.astype(np.float32), vol_init=vols[0].astype(np.float32),
               vol_final=vols[-1].astype(np.float32), track_image=images[track],
               track_clean=clean[track], track_proj=np.array(track_proj),
               track_resid=np.array(track_resid), R_true=R_true, t_true=t_true,
               R_final=R_est, t_final=t_est, sigma=sigma, voxel=voxel, lam=lam)
    out["loss_total"] = out["loss_data"] + out["loss_reg"]
    # expected data term if the model were perfect: the noise alone
    out["loss_noise_floor"] = float(images.size * sigma ** 2)
    # an illustrative spatial regularisation-strength map: strong where the
    # (smoothed) final map carries little signal, weak in the reliable core
    sig = ndimage.gaussian_filter(np.abs(vols[-1]), 2.0)
    out["reg_strength"] = (1 - sig / sig.max()).astype(np.float32)
    return out


# ------------------------------------------------------------------ clip 05
def _gmm2(X, iters=40, record_every=5):
    """2-component full-covariance Gaussian mixture by EM (deterministic init)."""
    o = np.argsort(X[:, 0])
    mu = np.stack([X[o[: len(X) // 5]].mean(0), X[o[-len(X) // 5:]].mean(0)])
    cov = np.stack([np.cov(X.T)] * 2)
    pi = np.array([0.5, 0.5])
    snaps = []
    for it in range(iters):
        logp = []
        for c in range(2):
            d = X - mu[c]
            inv = np.linalg.inv(cov[c])
            logp.append(np.log(pi[c]) - 0.5 * np.einsum("ij,jk,ik->i", d, inv, d)
                        - 0.5 * np.log(np.linalg.det(cov[c])))
        logp = np.stack(logp, 1)
        logp -= logp.max(1, keepdims=True)
        r = np.exp(logp)
        r /= r.sum(1, keepdims=True)
        if it % record_every == 0:
            snaps.append(r[:, 0].copy())
        Nk = r.sum(0)
        pi = Nk / len(X)
        mu = (r.T @ X) / Nk[:, None]
        for c in range(2):
            d = X - mu[c]
            cov[c] = (r[:, c, None] * d).T @ d / Nk[c] + 1e-6 * np.eye(2)
    snaps.append(r[:, 0].copy())
    return r, mu, cov, pi, np.array(snaps)


def classification(n: int = 64, voxel: float = 2.0, n_particles: int = 400,
                   frac_A: float = 0.55, snr: float = 0.12, lowpass_px: float = 1.5) -> dict:
    """Two synthetic conformers (appendage at A or B), one shared view for an
    illustrative PCA feature plot, an EM Gaussian mixture for soft class
    assignment and the resulting (idealised) class densities."""
    rng = np.random.default_rng(SEED + 5)
    states = np.where(rng.random(n_particles) < frac_A, 0, 1)    # 0 = A, 1 = B
    R = pj.rotation_from_euler(-35, 25, 10)                        # appendage clearly visible
    base = {s: pj.project_analytic(R, s, n, voxel, atoms=False) for s in "AB"}
    clean = np.stack([base["AB"[s]] for s in states])
    sigma = op.sigma_for_snr(base["A"], snr)
    imgs = clean + sigma * rng.standard_normal(clean.shape)
    def pca2(F):
        F = F - F.mean(0)
        U, Sv, Vt = np.linalg.svd(F, full_matrices=False)
        pcs = U[:, :2] * Sv[:2]
        if pcs[states == 0, 0].mean() > 0:      # fix the arbitrary PCA sign: A on the left
            pcs[:, 0] *= -1
            Vt[0] *= -1
        return pcs / pcs.std(0), Vt
    # raw pixels vs low-pass-filtered pixels: the chosen feature representation
    pcs_raw, _ = pca2(imgs.reshape(n_particles, -1))
    feats = np.stack([ndimage.gaussian_filter(im, lowpass_px) for im in imgs])
    pcs, Vt = pca2(feats.reshape(n_particles, -1))
    r, mu, cov, pi, snaps = _gmm2(pcs)
    if mu[0, 0] > mu[1, 0]:                  # class 1 = the left (A-like) cluster
        r, mu, cov, pi = r[:, ::-1], mu[::-1], cov[::-1], pi[::-1]
        snaps = 1 - snaps
    wA = np.array([(r[:, c] * (states == 0)).sum() / r[:, c].sum() for c in range(2)])
    consensus_wA = float((states == 0).mean())
    maskA = sp.appendage_mask("A", n, voxel)
    vA = sp.conformer("A", n, voxel, atoms=False)
    vB = sp.conformer("B", n, voxel, atoms=False)

    def app_contrast(w):
        v = w * vA + (1 - w) * vB
        return float(v[maskA].mean())

    out = dict(states=states, pcs=pcs, pcs_raw=pcs_raw, resp=r[:, 0], resp_snaps=snaps, gmm_mu=mu, gmm_cov=cov,
               gmm_pi=pi, class_wA=wA, consensus_wA=consensus_wA, sigma=sigma,
               examples=imgs[:10], example_states=states[:10], base_A=base["A"],
               base_B=base["B"], pc_images=Vt[:2].reshape(2, n, n), voxel=voxel)
    out["contrast_consensus"] = app_contrast(consensus_wA)
    out["contrast_classes"] = np.array([app_contrast(w) for w in wA])
    out["contrast_pure"] = np.array([app_contrast(1.0), app_contrast(0.0)])
    out["ambiguous"] = np.nonzero((r[:, 0] > 0.2) & (r[:, 0] < 0.8))[0]
    out["accuracy"] = float(((r[:, 0] > 0.5) == (states == 0)).mean())
    return out


# ------------------------------------------------------------------ clip 06
def _soft_sphere(n, voxel, center, radius, edge=4.0):
    ax = sp.grid_coords(n, voxel)
    X, Y, Z = np.meshgrid(ax, ax, ax, indexing="ij")
    d = np.sqrt((X - center[0]) ** 2 + (Y - center[1]) ** 2 + (Z - center[2]) ** 2)
    return np.clip(0.5 + 0.5 * np.cos(np.pi * np.clip((d - radius) / edge, 0, 1)), 0, 1) \
        * (d < radius + edge)


def _partial_volume(n, voxel, indices, state="A"):
    """Density of a subset of blobs (domains/appendage) including their pseudo-atoms."""
    ax = sp.grid_coords(n, voxel)
    X, Y, Z = np.meshgrid(ax, ax, ax, indexing="ij")
    vol = np.zeros((n, n, n))
    for b in sp.blob_params(state):
        if b.index in indices:
            vol += sp._eval_blob(b, X, Y, Z)
    pos, owner = sp.pseudo_atoms(state)
    sp._stamp_atoms(vol, pos[np.isin(owner, list(indices))], voxel, 0.35)
    return vol


def band_edges(n: int, shells_per_band: int = 2):
    """Integer shell edges so every band holds the same number of whole shells."""
    return np.arange(2, n // 2 + 1, shells_per_band)


def band_inner_products(signal, noise_a, noise_b, weights, shells_per_band=2):
    """Per Fourier band j and window w: weighted inner products of the
    band-pass-filtered signal (s) and the two noise maps (a, b).

    Each map is band-filtered on the whole box FIRST and windowed afterwards,
    which avoids the spurious high-frequency correlation a shared mask would
    impose. Everything downstream is linear in the noise gain g, so
    corr_w,j(g) follows from these six sums without refiltering.
    Returns ip[w, j, 6] with columns (ss, sa, sb, aa, bb, ab) and band centres s_j
    in index units.
    """
    n = signal.shape[0]
    FS, FA, FB = fo.fft3c(signal), fo.fft3c(noise_a), fo.fft3c(noise_b)
    r = fo.shell_index(signal.shape)
    edges = band_edges(n, shells_per_band)
    n_bands = len(edges) - 1
    ip = np.zeros((len(weights), n_bands, 6))
    for j in range(n_bands):
        band = (r >= edges[j]) & (r < edges[j + 1])
        s_ = np.real(fo.ifft3c(FS * band))
        a_ = np.real(fo.ifft3c(FA * band))
        b_ = np.real(fo.ifft3c(FB * band))
        for wi, w in enumerate(weights):
            ip[wi, j] = [np.sum(w * s_ * s_), np.sum(w * s_ * a_), np.sum(w * s_ * b_),
                         np.sum(w * a_ * a_), np.sum(w * b_ * b_), np.sum(w * a_ * b_)]
    centres = 0.5 * (edges[1:] + edges[:-1])
    return ip, centres


def band_corr(ip, g):
    """Correlation of (s + g a) with (s + g b) per window and band."""
    ss, sa, sb, aa, bb, ab = np.moveaxis(ip, -1, 0)
    num = ss + g * (sa + sb) + g * g * ab
    da = ss + 2 * g * sa + g * g * aa
    db = ss + 2 * g * sb + g * g * bb
    den = np.sqrt(np.maximum(da * db, 1e-300))
    return num / den


def molecular_mask(vol, voxel, frac=0.08, dilate_A=3.0, soft_A=2.0):
    m = ndimage.gaussian_filter(vol, 2.0 / voxel) > frac * vol.max()
    m = ndimage.binary_dilation(m, iterations=max(1, int(round(dilate_A / voxel))))
    return np.clip(ndimage.gaussian_filter(m.astype(float), soft_A / voxel), 0, 1)


FLEX = {  # blob index -> extra blur (angstrom): rigid core, mobile tips, flexible appendage
    **{i: 0.0 for i in range(3, 12)}, 0: 0.9, 1: 0.6, 2: 0.4, 12: 0.4, 13: 0.6, 14: 0.9, 15: 4.5}


def fsc_resolution(n: int = 128, voxel: float = 1.0, target_A: float = 3.2,
                   coord_err_A: float = 0.4, local_radius_A: float = 5.0,
                   shells_per_band: int = 2) -> dict:
    """Illustrative gold-standard resolution analysis on synthetic half-maps.

    truth   = phantom whose mobile parts (tips, appendage) are extra blurred
    half_X  = truth + g * independent noise_X  (white noise, same spectrum shape)
    The single noise gain g is solved so the masked half-map FSC crosses 0.143
    at target_A (the "3.2 A" value is chosen to echo the paper; it is not a
    measurement). Local resolution = the same band-limited correlation inside a
    small soft sphere around each blob. Map-model FSC = full map vs a map made
    from a toy atomic model with coordinate error, read at 0.5.
    """
    rng = np.random.default_rng(SEED + 6)
    truth = np.zeros((n, n, n))
    for blur in sorted(set(FLEX.values())):
        part = _partial_volume(n, voxel, [i for i, b in FLEX.items() if b == blur])
        truth += ndimage.gaussian_filter(part, blur / voxel) if blur > 0 else part
    noise_a = np.random.default_rng(SEED + 61).standard_normal(truth.shape)
    noise_b = np.random.default_rng(SEED + 62).standard_normal(truth.shape)
    mask = molecular_mask(truth, voxel)
    centers = np.array([b.center for b in sp.blob_params("A")])
    windows = [_soft_sphere(n, voxel, c, local_radius_A, edge=2.5) for c in centers]
    ip, cen = band_inner_products(truth, noise_a, noise_b, [mask] + windows,
                                  shells_per_band)
    s_b = cen / (n * voxel)
    # solve the noise gain so the masked global FSC crosses 0.143 at target_A
    lo, hi = 1e-6, 1e3
    for _ in range(80):
        g = np.sqrt(lo * hi)
        x = 1 / fo.crossing(s_b, band_corr(ip[0], g), 0.143)
        lo, hi = (g, hi) if x < target_A else (lo, g)
    g = np.sqrt(lo * hi)
    curves = band_corr(ip, g)
    half_A, half_B = truth + g * noise_a, truth + g * noise_b
    res = dict(s=s_b, fsc_half=curves[0], voxel=voxel, target_A=target_A, noise_gain=g,
               mask_fraction=mask.mean())
    res["half_crossing_s"] = fo.crossing(s_b, curves[0], 0.143)
    s_u, f_u, _ = fo.fsc(half_A, half_B, voxel)          # for reference: no mask
    res["s_unmasked"], res["fsc_unmasked"] = s_u, f_u
    res["local_curves"] = curves[1:]
    res["local_res_A"] = np.array([1 / fo.crossing(s_b, c, 0.143) for c in curves[1:]])
    res["blob_centers"] = centers
    res["flex_blur_A"] = np.array([FLEX[i] for i in range(16)])
    # map-model FSC: toy atomic model with coordinate error vs the full (averaged) map
    pos, owner = sp.pseudo_atoms("A")
    pos_model = pos + rng.normal(0, coord_err_A, pos.shape)
    model_map = sp.phantom_volume(n, voxel, state="A", atoms=False)   # smooth domains
    sp._stamp_atoms(model_map, pos_model, voxel, 0.35)
    full_map = 0.5 * (half_A + half_B)
    ipm, _ = band_inner_products(model_map, full_map - model_map, np.zeros_like(full_map),
                                 [mask], shells_per_band)
    # corr(model, full) = corr(s, s + 1 * d) with d = full - model; reuse band_corr
    ss, sa, _, aa, _, _ = np.moveaxis(ipm[0], -1, 0)
    res["fsc_model"] = (ss + sa) / np.sqrt(np.maximum(ss * (ss + 2 * sa + aa), 1e-300))
    res["model_crossing_s"] = fo.crossing(s_b, res["fsc_model"], 0.5)
    res["atoms"] = pos_model
    res["atom_owner"] = owner
    # 2D displays: central xy sections of both half-maps and their transforms
    c = n // 2
    res["sec_A"] = half_A[:, :, c]
    res["sec_B"] = half_B[:, :, c]
    res["sec_truth"] = truth[:, :, c]
    res["ft_A"] = fo.log_magnitude(fo.fft2c(res["sec_A"]))
    res["ft_B"] = fo.log_magnitude(fo.fft2c(res["sec_B"]))
    return res


# ------------------------------------------------------------------ renders
def _render(rot, tilt, psi, state, wA, n, level, extent):
    R = pj.rotation_from_euler(rot, tilt, psi)
    fn = None if wA is None else pj.mixture_density_fn({"A": wA, "B": 1 - wA})
    r = pj.render_density(R, state=state, n=n, level=level, extent=extent, nz=n, density_fn=fn)
    return dict(shade=r["shade"], mask=r["mask"], label=r["label"])


def cached_render(rot=DISPLAY["rot"], tilt=DISPLAY["tilt"], psi=DISPLAY["psi"], state="A",
                  wA=None, n=220, level=0.35, extent=64.0) -> dict:
    """Shaded isodensity render (~7 s uncached), cached on disk.

    wA=None renders the pure `state`; wA in [0, 1] renders the mixture
    wA*rho_A + (1-wA)*rho_B (labels refer to `state`'s blobs; index 15 = appendage).
    Feed the result to styling.density_image(...).
    """
    from .cache import cached
    return cached("render", _render, rot=float(rot), tilt=float(tilt), psi=float(psi),
                  state=state, wA=None if wA is None else float(wA), n=int(n),
                  level=float(level), extent=float(extent))
