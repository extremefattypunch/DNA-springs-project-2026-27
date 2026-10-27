"""Reproducible numbers behind every graph, heatmap and moving dot in the clips.

Each function returns a dict of NumPy arrays/scalars and is deterministic
(seeded from point_constellation.SEED). Scenes load them through
cache.cached(name, fn, **params) so nothing heavy runs inside a frame.
All lengths are dimensionless (same units as the dot coordinates).
"""
from __future__ import annotations

import numpy as np
from scipy import optimize

from . import fourier as fo
from . import noise_and_ctf as nc
from . import point_constellation as pc
from . import projections as pj

# A fixed reference view used where one generic projection is needed.
VIEW = (35.0, 60.0, 20.0)       # (alpha, theta, psi) in degrees, ZYZ


def view_R(angles=VIEW) -> np.ndarray:
    return pj.rotation_from_euler(*angles)


def _np(d: dict) -> dict:
    return {k: np.asarray(v) for k, v in d.items()}


# ======================================================================= clip 01
def forward_model(n1d: int = 256, extent1d: float = 2.0, snr_1d: float = 1.0,
                  snr_2d: float = 0.6, n_noise_frames: int = 6, band_width: float = 0.08
                  ) -> dict:
    """1D illustration and 2D inset of projection -> CTF -> noise.

    1D: s(x) is the line integral of the projected dot image along y (so it is
    the same source, seen as a profile). CTF is applied in Fourier space,
    then i.i.d. Gaussian noise is added. Band components show one frequency
    band passed with negative sign (lobe 1), one at a zero, one positive.
    """
    R = view_R()
    P = pj.project_points(R, pc.points())
    # compute on a 4x wider periodic interval (finer frequency spacing), show [-extent1d, extent1d)
    pad = 4
    N = n1d * pad
    dx = 2 * extent1d / n1d
    xb = pc.grid_coords(N, extent1d * pad)
    crop = slice(N // 2 - n1d // 2, N // 2 + n1d // 2)
    x = xb[crop]
    sb = np.zeros(N)
    for (px, _), a, sg in zip(P, pc.AMPLITUDES, pc.SIGMAS):
        sb += a * 2 * np.pi * sg ** 2 * np.exp(-0.5 * (xb - px) ** 2 / sg ** 2)
    s, s_ctf = sb[crop], nc.apply_ctf_1d(sb, dx)[crop]
    g = nc.rng(1)
    sigma1 = nc.sigma_for_snr(s_ctf, snr_1d)
    eps = sigma1 * g.standard_normal(n1d)
    I = s_ctf + eps

    kmax = 8.0
    k = np.linspace(0, kmax, 801)
    ctf = nc.ctf_1d(k)
    zs = nc.zeros(kmax)
    # bands: middle of lobe 1 (negative), first zero, middle of lobe 2 (positive)
    k_neg = 0.5 * zs[0]
    k_zero = zs[0]
    k_pos = 0.5 * (zs[0] + zs[1])
    kk = np.fft.fftfreq(N, d=dx)
    S = np.fft.fft(sb)
    bands_before, bands_after = [], []
    for k0 in (k_neg, k_zero, k_pos):
        w = np.exp(-0.5 * ((np.abs(kk) - k0) / band_width) ** 2)
        bands_before.append(np.real(np.fft.ifft(S * w))[crop])
        bands_after.append(np.real(np.fft.ifft(S * w * nc.ctf_1d(np.abs(kk))))[crop])

    # 2D inset
    img = pj.project_image(R)
    img_ctf = nc.apply_ctf(img)
    sigma2 = nc.sigma_for_snr(img_ctf, snr_2d)
    noise2 = sigma2 * nc.rng(2).standard_normal(img.shape)
    alphas = np.linspace(0, 1, n_noise_frames)
    frames = np.stack([img_ctf + a * noise2 for a in alphas])
    return _np(dict(
        x=x, s=s, s_ctf=s_ctf, eps=eps, I=I, sigma_1d=sigma1, snr_1d=snr_1d,
        k=k, ctf=ctf, ctf_zeros=zs, k_bands=np.array([k_neg, k_zero, k_pos]),
        ctf_at_bands=nc.ctf_1d(np.array([k_neg, k_zero, k_pos])),
        bands_before=np.stack(bands_before), bands_after=np.stack(bands_after),
        img=img, img_ctf=img_ctf, noise_frames=frames, noise_alphas=alphas, sigma_2d=sigma2,
        ctf2d=nc.ctf_2d(shifted=True), points3d=pc.points(), proj2d=P, R=R))


# ======================================================================= clip 02
def alignment(n: int = 64, extent: float = 2.4, n_obs: int = 256, snr: float = 0.5,
              shift_sd: float = 0.22, ns=(1, 4, 16, 64, 256)) -> dict:
    """Identical 4-dot signal under random in-plane rotation + translation + noise.

    Alignment uses the TRUE (phi_i, t_i) (alignment assumed accurate), undoing
    them exactly with a Fourier phase ramp and three-shear rotation, which
    keep the noise i.i.d. with unchanged sigma. Measured amplitude SNR of the
    aligned mean is std(S) / std(mean - S).
    """
    g = nc.rng(20)
    S = pj.dots_image(pc.PATTERN_2D, pc.PATTERN_AMPLITUDES, pc.PATTERN_SIGMAS, n, extent)
    sigma = nc.sigma_for_snr(S, snr)
    phi = g.uniform(-180, 180, n_obs)
    t = np.clip(g.normal(0, shift_sd, (n_obs, 2)), -0.45, 0.45)
    obs, aligned, clean = [], [], []
    for i in range(n_obs):
        P = pj.transform_2d(pc.PATTERN_2D, phi[i], t[i])
        c = pj.dots_image(P, pc.PATTERN_AMPLITUDES, pc.PATTERN_SIGMAS, n, extent)
        o = c + sigma * g.standard_normal((n, n))
        a = pj.rotate_image(pj.shift_image(o, -t[i], extent), -phi[i])
        obs.append(o)
        aligned.append(a)
        clean.append(c)
    obs, aligned, clean = map(np.stack, (obs, aligned, clean))
    unaligned_mean = np.stack([obs[:N].mean(0) for N in ns])
    aligned_mean = np.stack([aligned[:N].mean(0) for N in ns])
    snr_meas = np.array([np.std(S) / np.std(m - S) for m in aligned_mean])
    sd_meas = np.array([np.std(m - S) for m in aligned_mean])
    corr_un = np.array([np.corrcoef(m.ravel(), S.ravel())[0, 1] for m in unaligned_mean])
    corr_al = np.array([np.corrcoef(m.ravel(), S.ravel())[0, 1] for m in aligned_mean])
    return _np(dict(
        S=S, sigma=sigma, snr1=snr, phi=phi, t=t, obs12=obs[:12], clean12=clean[:12],
        aligned12=aligned[:12], ns=np.array(ns), unaligned_mean=unaligned_mean,
        aligned_mean=aligned_mean, residual=aligned_mean - S, snr_measured=snr_meas,
        snr_theory=snr * np.sqrt(np.array(ns)), sd_measured=sd_meas,
        sd_theory=sigma / np.sqrt(np.array(ns)), corr_unaligned=corr_un,
        corr_aligned=corr_al, pattern=pc.PATTERN_2D, extent=extent))


# ======================================================================= clip 03
VIEWS_03 = np.array([(0, 0, 0), (40, 70, 0), (120, 50, 0), (-70, 100, 0)], float)


def orientation_fourier(n_cov: int = 60, max_theta_biased: float = 15.0,
                        ks=(1, 2, 3, 5, 10, 20), meas_noise: float = 0.03,
                        n_planes_draw: int = 12) -> dict:
    """Projections at several views, Fourier-slice check, coverage and a
    least-squares 'inversion' of the dot coordinates from known poses."""
    X = pc.points()
    Rs = np.stack([pj.rotation_from_euler(*a) for a in VIEWS_03])
    proj = np.stack([pj.project_points(R, X) for R in Rs])
    imgs = np.stack([pj.project_image(R) for R in Rs])
    ft_mag = np.stack([fo.log_magnitude(fo.fft2c(im)) for im in imgs])
    ft_re = np.stack([np.real(fo.fft2c(im)) for im in imgs])
    ft_im = np.stack([np.imag(fo.fft2c(im)) for im in imgs])
    vdirs = np.stack([pj.viewing_direction(R) for R in Rs])

    # numerical Fourier-slice check (same grid for field and images)
    n, ext = pc.FIELD_N, pc.FIELD_EXTENT
    vol = pc.field(n=n, extent=ext)
    F3 = fo.fft3c(fo.pad(vol, 2))
    h = pc.spacing(n, ext)
    m = fo.shell_index((2 * n, 2 * n)) < n // 2
    corr_RT, corr_R = [], []
    for R in Rs[1:]:
        p = pj.project_image(R, n=n, extent=ext) / h
        P2 = fo.fft2c(fo.pad(p, 2))
        a, b = fo.central_slice(F3, R)[m], fo.central_slice(F3, R.T)[m]
        corr_RT.append(np.abs(np.vdot(P2[m], a)) / (np.linalg.norm(P2[m]) * np.linalg.norm(a)))
        corr_R.append(np.abs(np.vdot(P2[m], b)) / (np.linalg.norm(P2[m]) * np.linalg.norm(b)))

    # coverage: uniform vs biased orientation sets
    g = nc.rng(30)
    Ru = pj.random_rotations(n_cov, g)
    Rb = pj.biased_rotations(n_cov, g, max_theta_biased)
    cov_u = fo.ball_fraction(fo.sampling_coverage(Ru, 32))
    cov_b = fo.ball_fraction(fo.sampling_coverage(Rb, 32))
    normals_u = np.stack([fo.plane_normal(R) for R in Ru[:n_planes_draw]])
    normals_b = np.stack([fo.plane_normal(R) for R in Rb[:n_planes_draw]])

    # inversion of dot coordinates from k known views (ridge toward a neutral start)
    x0 = np.zeros_like(X)
    def solve(R_set, k):
        A = np.concatenate([R[:2] for R in R_set[:k]])            # (2k, 3)
        out = []
        for j in range(len(X)):
            b = np.concatenate([(R @ X[j])[:2] + meas_noise * gn.standard_normal(2)
                                for R in R_set[:k]])
            mu = 1e-3
            out.append(np.linalg.solve(A.T @ A + mu * np.eye(3), A.T @ b + mu * x0[j]))
        return np.array(out)
    gn = nc.rng(31)
    est_u = np.stack([solve(Ru, k) for k in ks])
    gn = nc.rng(31)
    est_b = np.stack([solve(Rb, k) for k in ks])
    err_u = np.linalg.norm(est_u - X, axis=-1).mean(-1)
    err_b = np.linalg.norm(est_b - X, axis=-1).mean(-1)
    return _np(dict(
        X=X, views=VIEWS_03, Rs=Rs, proj=proj, imgs=imgs, ft_mag=ft_mag, ft_re=ft_re,
        ft_im=ft_im, vdirs=vdirs, corr_RT=np.array(corr_RT), corr_R=np.array(corr_R),
        cov_uniform=cov_u, cov_biased=cov_b, normals_uniform=normals_u,
        normals_biased=normals_b, n_cov=n_cov, max_theta_biased=max_theta_biased,
        ks=np.array(ks), est_uniform=est_u, est_biased=est_b, err_uniform=err_u,
        err_biased=err_b, meas_noise=meas_noise))


# ======================================================================= clip 04
def _axis_angle(w: np.ndarray) -> np.ndarray:
    th = np.linalg.norm(w)
    if th < 1e-12:
        return np.eye(3)
    k = w / th
    K = np.array([[0, -k[2], k[1]], [k[2], 0, -k[0]], [-k[1], k[0], 0]])
    return np.eye(3) + np.sin(th) * K + (1 - np.cos(th)) * K @ K


def _angle_between(Ra: np.ndarray, Rb: np.ndarray) -> float:
    c = (np.trace(Ra.T @ Rb) - 1) / 2
    return float(np.degrees(np.arccos(np.clip(c, -1, 1))))


def refinement(n: int = 40, extent: float = 2.0, n_obs: int = 24, snr: float = 1.0,
               init_pose_deg: float = 15.0, init_shift: float = 0.15,
               init_point_sd: float = 0.25, lam: float = 2.0, n_iter: int = 6,
               track: int = 0) -> dict:
    """Toy alternating minimisation of the displayed loss

        L(X, {R_i, t_i}) = sum_i ||I_i - P_{R_i, t_i}[X]||^2 + lam * sum_j |x_j|^2

    where P_{R,t}[X] is the exact Gaussian-kernel projection of the dot set X
    (CTF = 1 here, to keep the toy small). Pose step: each (R_i, t_i) is
    refined by local L-BFGS with X fixed. Source step: X is refined with all
    poses fixed. Every stored loss is evaluated from the stored estimates.
    hist_X / pose_err / shift_err / track_R / track_t are expressed in the true
    frame after removing the global rigid-motion ambiguity (Kabsch alignment);
    hist_X_raw keeps the optimiser's own frame.
    """
    g = nc.rng(40)
    X_true = pc.points()
    R_true = pj.random_rotations(n_obs, g)
    t_true = g.normal(0, 0.1, (n_obs, 2))
    clean = np.stack([pj.dots_image(pj.project_points(R, X_true, t), n=n, extent=extent)
                      for R, t in zip(R_true, t_true)])
    sigma = nc.sigma_for_snr(clean[0], snr)
    imgs = clean + sigma * g.standard_normal(clean.shape)

    def render(X, R, t):
        return pj.dots_image(pj.project_points(R, X, t), n=n, extent=extent)

    def loss_terms(X, Rs, ts):
        data = sum(np.sum((I - render(X, R, t)) ** 2) for I, R, t in zip(imgs, Rs, ts))
        return float(data), float(lam * np.sum(X ** 2))

    # initial estimates
    R_est = np.stack([_axis_angle(np.radians(init_pose_deg) * _unit(g.standard_normal(3))) @ R
                      for R in R_true])
    t_est = t_true + init_shift * np.stack([_unit(g.standard_normal(2)) for _ in range(n_obs)])
    X_est = X_true + init_point_sd * g.standard_normal(X_true.shape)

    hist_X, hist_loss, hist_data, hist_reg, hist_step = [], [], [], [], []
    hist_pose_err, hist_shift_err, hist_track_R, hist_track_t = [], [], [], []
    hist_track_pred, hist_track_res = [], []

    hist_X_raw = []

    def record(step):
        d, r = loss_terms(X_est, R_est, t_est)
        # the loss is unchanged by a global rigid motion of X compensated in every
        # pose (gauge freedom); express estimates in the true frame for display
        Q, me, mt = _kabsch(X_est, X_true)
        X_al = (X_est - me) @ Q.T + mt
        R_al = R_est @ Q.T
        t_al = t_est + (me @ R_est.transpose(0, 2, 1) - mt @ R_al.transpose(0, 2, 1))[:, :2]
        hist_X.append(X_al); hist_X_raw.append(X_est.copy())
        hist_loss.append(d + r); hist_data.append(d); hist_reg.append(r); hist_step.append(step)
        hist_pose_err.append([_angle_between(a, b) for a, b in zip(R_al, R_true)])
        hist_shift_err.append(np.linalg.norm(t_al - t_true, axis=1))
        hist_track_R.append(R_al[track].copy()); hist_track_t.append(t_al[track].copy())
        pred = render(X_est, R_est[track], t_est[track])
        hist_track_pred.append(pred); hist_track_res.append(imgs[track] - pred)

    record(0)                                     # 0 = initial
    for _ in range(n_iter):
        # pose step
        for i in range(n_obs):
            R0 = R_est[i]
            def f(p):
                return np.sum((imgs[i] - render(X_est, _axis_angle(p[:3]) @ R0, p[3:])) ** 2)
            p = optimize.minimize(f, np.r_[0, 0, 0, t_est[i]], method="L-BFGS-B").x
            R_est[i] = _axis_angle(p[:3]) @ R0
            t_est[i] = p[3:]
        record(1)                                 # 1 = after pose step
        # source step
        def fX(v):
            Xv = v.reshape(-1, 3)
            d, r = loss_terms(Xv, R_est, t_est)
            return d + r
        X_est = optimize.minimize(fX, X_est.ravel(), method="L-BFGS-B").x.reshape(-1, 3)
        record(2)                                 # 2 = after source step

    # regularization example (1D): noisy samples of the source profile
    reg = _regularization_example()
    return _np(dict(
        imgs=imgs, clean=clean, sigma=sigma, lam=lam, R_true=R_true, t_true=t_true,
        X_true=X_true, hist_X=np.stack(hist_X), hist_X_raw=np.stack(hist_X_raw), loss=np.array(hist_loss),
        loss_data=np.array(hist_data), loss_reg=np.array(hist_reg), step=np.array(hist_step),
        pose_err=np.array(hist_pose_err), shift_err=np.array(hist_shift_err),
        track=track, track_R=np.stack(hist_track_R), track_t=np.stack(hist_track_t),
        track_pred=np.stack(hist_track_pred), track_res=np.stack(hist_track_res),
        track_res_norm=np.array([np.linalg.norm(r) for r in hist_track_res]),
        noise_floor=float(imgs.size * sigma ** 2), **{f"reg_{k}": v for k, v in reg.items()}))


def _kabsch(A: np.ndarray, B: np.ndarray):
    """Rotation Q and centroids so that (A - ma) Q^T + mb best matches B."""
    ma, mb = A.mean(0), B.mean(0)
    U, _, Vt = np.linalg.svd((A - ma).T @ (B - mb))
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    return Vt.T @ np.diag([1, 1, d]) @ U.T, ma, mb


def _unit(v: np.ndarray) -> np.ndarray:
    return v / np.linalg.norm(v)


def _regularization_example(n_samp: int = 40, n_basis: int = 22, lam_smooth: float = 3e-4,
                            noise: float = 0.12) -> dict:
    """Least squares with and without a smoothness penalty on a 1D profile.

    f(x) = 1D projection profile of the dots; y = f(x_i) + noise. Both fits use
    the same cosine/sine basis; the regularized one adds lam * ||D2 f||^2
    (second differences on a dense grid). Reports the data residual and the
    error against f for both, showing that a smaller data residual is not a
    better estimate.
    """
    fm = forward_model()
    xd, f = fm["x"], fm["s"]
    g = nc.rng(41)
    xs = np.sort(g.uniform(-1.9, 1.9, n_samp))
    fs = np.interp(xs, xd, f)
    y = fs + noise * g.standard_normal(n_samp)
    L = 4.0

    def basis(x):
        cols = [np.ones_like(x)]
        for m in range(1, (n_basis + 1) // 2 + 1):
            cols += [np.cos(2 * np.pi * m * x / L), np.sin(2 * np.pi * m * x / L)]
        return np.stack(cols[:n_basis], -1)
    A, Ad = basis(xs), basis(xd)
    D2 = np.diff(np.eye(len(xd)), 2, axis=0) @ Ad / (xd[1] - xd[0]) ** 2
    c_ls = np.linalg.lstsq(A, y, rcond=None)[0]
    c_rg = np.linalg.solve(A.T @ A + lam_smooth * D2.T @ D2 * (xd[1] - xd[0]), A.T @ y)
    fit_ls, fit_rg = Ad @ c_ls, Ad @ c_rg
    out = dict(x=xd, f=f, xs=xs, y=y, fit_ls=fit_ls, fit_rg=fit_rg,
               data_res_ls=float(np.sum((A @ c_ls - y) ** 2)),
               data_res_rg=float(np.sum((A @ c_rg - y) ** 2)),
               true_err_ls=float(np.sqrt(np.mean((fit_ls - f) ** 2))),
               true_err_rg=float(np.sqrt(np.mean((fit_rg - f) ** 2))),
               penalty_ls=float(np.sum((D2 @ c_ls) ** 2) * (xd[1] - xd[0])),
               penalty_rg=float(np.sum((D2 @ c_rg) ** 2) * (xd[1] - xd[0])),
               lam_smooth=lam_smooth)
    return out


# ======================================================================= clip 05
def classification(n: int = 48, extent: float = 2.0, n_obs: int = 240, frac_A: float = 0.6,
                   snr: float = 0.15, view=(10.0, 20.0, 0.0), em_iters: int = 60) -> dict:
    """Two states (dot 5 moved) seen at one matched orientation, mixed.

    PCA of the noisy images gives PC1/PC2 scores; a 2-component Gaussian
    mixture (EM, full covariances) on those scores gives soft responsibilities
    p(c | I_i). Class means are responsibility-weighted averages of the images.
    """
    g = nc.rng(50)
    R = pj.rotation_from_euler(*view)
    cA, cB = pj.project_image(R, "A", n=n, extent=extent), pj.project_image(R, "B", n=n, extent=extent)
    labels = (g.uniform(size=n_obs) >= frac_A).astype(int)      # 0 = A, 1 = B
    sigma = nc.sigma_for_snr(cA, snr)
    clean = np.where(labels[:, None, None] == 0, cA, cB)
    imgs = clean + sigma * g.standard_normal(clean.shape)
    Xf = imgs.reshape(n_obs, -1)
    mu = Xf.mean(0)
    U, Sv, Vt = np.linalg.svd(Xf - mu, full_matrices=False)
    scores = (Xf - mu) @ Vt[:2].T
    expl = Sv ** 2 / np.sum(Sv ** 2)

    # EM for a 2-component Gaussian mixture in (PC1, PC2)
    order = np.argsort(scores[:, 0])
    means = np.stack([scores[order[: n_obs // 4]].mean(0), scores[order[-n_obs // 4:]].mean(0)])
    covs = np.stack([np.cov(scores.T)] * 2)
    pis = np.array([0.5, 0.5])
    resp_hist = []
    for _ in range(em_iters):
        lik = np.stack([pis[c] * _gauss2(scores, means[c], covs[c]) for c in range(2)], -1)
        resp = lik / lik.sum(1, keepdims=True)
        resp_hist.append(resp[:, 0].copy())
        Nk = resp.sum(0)
        pis = Nk / n_obs
        means = (resp.T @ scores) / Nk[:, None]
        covs = np.stack([((resp[:, c, None] * (scores - means[c])).T @ (scores - means[c])) / Nk[c]
                         + 1e-6 * np.eye(2) for c in range(2)])
    # name component 0 as the one closer to state A
    mA = (resp[:, 0] @ Xf) / resp[:, 0].sum()
    if np.linalg.norm(mA - cA.ravel()) > np.linalg.norm(mA - cB.ravel()):
        resp = resp[:, ::-1]; means = means[::-1]; covs = covs[::-1]; pis = pis[::-1]
        resp_hist = [1 - r for r in resp_hist]
    pA = resp[:, 0]
    class_mean = np.stack([(resp[:, c] @ Xf / resp[:, c].sum()).reshape(n, n) for c in range(2)])
    mixed = imgs.mean(0)
    assign = (pA < 0.5).astype(int)
    amb = np.nonzero((pA > 0.25) & (pA < 0.75))[0]
    # visibility of the moved dot at its two sites in each mean
    p5A = pj.project_points(R, pc.points("A"))[pc.MOVED]
    p5B = pj.project_points(R, pc.points("B"))[pc.MOVED]
    def at(img, p):
        ax = pc.grid_coords(n, extent)
        return float(img[np.argmin(np.abs(ax - p[0])), np.argmin(np.abs(ax - p[1]))])
    peaks = np.array([[at(m, p5A), at(m, p5B)] for m in (cA, cB, mixed, class_mean[0], class_mean[1])])
    return _np(dict(
        cleanA=cA, cleanB=cB, imgs12=imgs[:12], labels=labels, sigma=sigma, mixed=mixed,
        mean_true_A=imgs[labels == 0].mean(0), mean_true_B=imgs[labels == 1].mean(0),
        scores=scores, explained=expl[:5], pc_images=Vt[:2].reshape(2, n, n), pA=pA,
        resp_hist=np.stack(resp_hist), em_means=means, em_covs=covs, em_pis=pis,
        class_mean=class_mean, assign=assign, ambiguous=amb,
        accuracy=float(np.mean(assign == labels)), frac_A=frac_A,
        site_A=p5A, site_B=p5B, peaks=peaks, points_A=pc.points("A"),
        points_B=pc.points("B"), proj_A=pj.project_points(R, pc.points("A")),
        proj_B=pj.project_points(R, pc.points("B"))))


def _gauss2(x: np.ndarray, m: np.ndarray, C: np.ndarray) -> np.ndarray:
    d = x - m
    q = np.einsum("ni,ij,nj->n", d, np.linalg.inv(C), d)
    return np.exp(-0.5 * q) / (2 * np.pi * np.sqrt(np.linalg.det(C)))


# ======================================================================= clip 06
REGION_1 = (0, 2, 3)      # dots 1, 3, 4: consistent estimates
REGION_2 = (1, 4)         # dots 2, 5: less consistent estimates


def fsc_resolution(n: int = 64, extent: float = 3.0, n_copies: int = 24,
                   jitter=(0.03, 0.18), noise_rel: float = 0.05, mask_radius: float = 0.55
                   ) -> dict:
    """Independent half-maps from independently estimated dot positions.

    Each half-map = average of n_copies fields whose dot positions are jittered
    (sd jitter[0] in region 1, jitter[1] in region 2) + independent white noise.
    Averaging jittered copies attenuates high frequencies more in region 2.
    FSC is computed on the 3D arrays (masked by soft spheres around the dots);
    region-wise FSC uses masks around each region; the map-model curve compares
    the full map with the field of the ideal point model.
    """
    X = pc.points()
    sd = np.array([jitter[0] if j in REGION_1 else jitter[1] for j in range(len(X))])
    h = pc.spacing(n, extent)

    def half(gen):
        est = X[None] + sd[None, :, None] * gen.standard_normal((n_copies,) + X.shape)
        m = np.mean([pc.field(n=n, extent=extent, pts_override=e) for e in est], axis=0)
        return m, est

    hA, estA = half(nc.rng(60))
    hB, estB = half(nc.rng(61))
    model = pc.field(n=n, extent=extent)
    s_noise = noise_rel * model.max()
    hA = hA + s_noise * nc.rng(62).standard_normal(hA.shape)
    hB = hB + s_noise * nc.rng(63).standard_normal(hB.shape)
    full = 0.5 * (hA + hB)

    def mask_for(idx):
        return np.max([fo.soft_mask(n, extent, X[j], mask_radius) for j in idx], axis=0)
    m_all = mask_for(range(len(X)))
    m1, m2 = mask_for(REGION_1), mask_for(REGION_2)
    s, f_glob, cnt = fo.fsc(hA * m_all, hB * m_all, h)
    _, f_r1, _ = fo.fsc(hA * m1, hB * m1, h)
    _, f_r2, _ = fo.fsc(hA * m2, hB * m2, h)
    _, f_mm, _ = fo.fsc(full * m_all, model * m_all, h)
    _, f_unmasked, _ = fo.fsc(hA, hB, h)
    s_half = fo.crossing(s, f_glob, 0.143)
    s_r1, s_r2 = fo.crossing(s, f_r1, 0.143), fo.crossing(s, f_r2, 0.143)
    s_model = fo.crossing(s, f_mm, 0.5)
    mid = n // 2
    return _np(dict(
        s=s, fsc_half=f_glob, fsc_region1=f_r1, fsc_region2=f_r2, fsc_model=f_mm,
        fsc_unmasked=f_unmasked, counts=cnt, s_half=s_half, d_half=1 / s_half,
        s_region1=s_r1, d_region1=1 / s_r1, s_region2=s_r2, d_region2=1 / s_r2,
        s_model=s_model, d_model=1 / s_model, X=X, est_A=estA, est_B=estB,
        jitter=np.array(jitter), region1=np.array(REGION_1), region2=np.array(REGION_2),
        noise_sd=s_noise, nyquist=1 / (2 * h),
        sec_A=hA[:, :, mid], sec_B=hB[:, :, mid], sec_model=model[:, :, mid],
        proj_A=hA.sum(2) * h, proj_B=hB.sum(2) * h, proj_model=model.sum(2) * h,
        ft2_A=fo.log_magnitude(fo.fft2c(hA.sum(2))), ft2_B=fo.log_magnitude(fo.fft2c(hB.sum(2)))))
