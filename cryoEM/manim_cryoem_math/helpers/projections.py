"""Rotations, ideal parallel projections and 2D image transforms.

Convention used everywhere in this project (verified in tests):

    rotated object      rho_R(r) = rho(R^T r)
    projection          (P_R rho)(x, y) = integral of rho(R^T (x, y, z)) dz
    Fourier slice       FT2[P_R rho](kx, ky) = FT3[rho](R^T (kx, ky, 0))

i.e. R rotates the object into the microscope frame and the beam runs along z.
The viewing direction expressed in the object frame is v = R^T e_z.
"""
from __future__ import annotations

import numpy as np
from scipy import ndimage

from . import synthetic_phantom as sp


# ---------------------------------------------------------------- rotations
def rot_x(a: float) -> np.ndarray:
    c, s = np.cos(a), np.sin(a)
    return np.array([[1, 0, 0], [0, c, -s], [0, s, c]])


def rot_y(a: float) -> np.ndarray:
    c, s = np.cos(a), np.sin(a)
    return np.array([[c, 0, s], [0, 1, 0], [-s, 0, c]])


def rot_z(a: float) -> np.ndarray:
    c, s = np.cos(a), np.sin(a)
    return np.array([[c, -s, 0], [s, c, 0], [0, 0, 1]])


def rotation_from_euler(rot: float, tilt: float, psi: float, degrees: bool = True) -> np.ndarray:
    """ZYZ Euler angles: R = Rz(psi) Ry(tilt) Rz(rot).

    (rot, tilt) choose the viewing direction (2 DOF); psi is the in-plane
    rotation about the beam axis (the 3rd DOF).
    """
    if degrees:
        rot, tilt, psi = np.radians([rot, tilt, psi])
    return rot_z(psi) @ rot_y(tilt) @ rot_z(rot)


def viewing_direction(R: np.ndarray) -> np.ndarray:
    """Beam direction expressed in the object frame: v = R^T e_z."""
    return R.T @ np.array([0.0, 0.0, 1.0])


def random_rotations(n: int, rng: np.random.Generator) -> np.ndarray:
    """Uniform random rotations on SO(3) via unit quaternions."""
    q = rng.standard_normal((n, 4))
    q /= np.linalg.norm(q, axis=1, keepdims=True)
    w, x, y, z = q.T
    return np.stack([
        np.stack([1 - 2 * (y * y + z * z), 2 * (x * y - z * w), 2 * (x * z + y * w)], -1),
        np.stack([2 * (x * y + z * w), 1 - 2 * (x * x + z * z), 2 * (y * z - x * w)], -1),
        np.stack([2 * (x * z - y * w), 2 * (y * z + x * w), 1 - 2 * (x * x + y * y)], -1),
    ], axis=1)


def biased_rotations(n: int, rng: np.random.Generator, max_tilt_deg: float = 15.0) -> np.ndarray:
    """Preferred orientation: views clustered around one direction (tilt <= max_tilt)."""
    rots = rng.uniform(0, 360, n)
    tilts = rng.uniform(0, max_tilt_deg, n)
    psis = rng.uniform(0, 360, n)
    return np.stack([rotation_from_euler(a, b, c) for a, b, c in zip(rots, tilts, psis)])


# ---------------------------------------------------------------- projections
def image_coords(n: int = sp.BOX_N, voxel: float = sp.VOXEL_A):
    ax = sp.grid_coords(n, voxel)
    return np.meshgrid(ax, ax, indexing="ij")     # X[i, j] = x_i, Y[i, j] = y_j


def project_analytic(R: np.ndarray, state: str | None = "A", n: int = sp.BOX_N,
                     voxel: float = sp.VOXEL_A, atoms: bool = True,
                     atom_weight: float = 0.35, shift=(0.0, 0.0)) -> np.ndarray:
    """Exact line integral of the Gaussian phantom after rotation by R.

    A 3D Gaussian (mass m, mean mu, covariance S) rotated by R has mean R mu and
    covariance R S R^T; integrating over z gives a 2D Gaussian with the xy
    marginal covariance and the same mass. Returned in units of density*angstrom.
    shift = (tx, ty) in angstrom translates the image (the t_i of the model).
    """
    X, Y = image_coords(n, voxel)
    img = np.zeros((n, n))
    for b in sp.blob_params(state):
        mu = R @ b.center
        S = R @ b.cov @ R.T
        S2 = S[:2, :2]
        mass = b.amplitude * (2 * np.pi) ** 1.5 * np.prod(b.sigmas)
        inv = np.linalg.inv(S2)
        dx, dy = X - mu[0] - shift[0], Y - mu[1] - shift[1]
        q = inv[0, 0] * dx * dx + 2 * inv[0, 1] * dx * dy + inv[1, 1] * dy * dy
        img += mass / (2 * np.pi * np.sqrt(np.linalg.det(S2))) * np.exp(-0.5 * q)
    if atoms:
        pos, _ = sp.pseudo_atoms(state)
        p2 = pos @ R.T
        s2 = sp.ATOM_SIGMA_A ** 2
        mass = atom_weight * (2 * np.pi * s2) ** 1.5
        for p in p2:
            r2 = (X - p[0] - shift[0]) ** 2 + (Y - p[1] - shift[1]) ** 2
            img += mass / (2 * np.pi * s2) * np.exp(-0.5 * r2 / s2)
    return img


def rotate_volume(vol: np.ndarray, R: np.ndarray, order: int = 1) -> np.ndarray:
    """rho_R(r) = rho(R^T r) on the same grid (rotation about the box centre)."""
    n = vol.shape[0]
    c = n // 2
    idx = np.indices(vol.shape).reshape(3, -1) - c
    src = R.T @ idx + c
    return ndimage.map_coordinates(vol, src, order=order, mode="constant").reshape(vol.shape)


def project_voxel(vol: np.ndarray, R: np.ndarray, voxel: float = sp.VOXEL_A) -> np.ndarray:
    """Numerical line integral of a voxel volume after rotation (sum along z)."""
    return rotate_volume(vol, R).sum(axis=2) * voxel


# ---------------------------------------------------------------- 2D transforms
def shift_image(img: np.ndarray, t, voxel: float = sp.VOXEL_A) -> np.ndarray:
    """Translate by t = (tx, ty) angstrom via the Fourier phase ramp e^{-2 pi i k.t}."""
    n = img.shape[0]
    k = np.fft.fftfreq(n, d=voxel)
    KX, KY = np.meshgrid(k, k, indexing="ij")
    ramp = np.exp(-2j * np.pi * (KX * t[0] + KY * t[1]))
    return np.real(np.fft.ifft2(np.fft.fft2(img) * ramp))


def rotate_image(img: np.ndarray, phi_deg: float, order: int = 1) -> np.ndarray:
    """In-plane rotation by phi (counter-clockwise in the displayed x-right/y-up frame)."""
    n = img.shape[0]
    c = n // 2
    a = np.radians(phi_deg)
    Rm = np.array([[np.cos(a), -np.sin(a)], [np.sin(a), np.cos(a)]])
    idx = np.indices(img.shape).reshape(2, -1) - c
    src = Rm.T @ idx + c
    return ndimage.map_coordinates(img, src, order=order, mode="constant").reshape(img.shape)


def _shear_axis(img: np.ndarray, factor: float, axis: int) -> np.ndarray:
    """Shift every line along `axis` by factor * (coordinate along the other axis),
    using exact Fourier phase ramps (unitary: preserves white-noise variance)."""
    n = img.shape[0]
    c = np.arange(n) - n // 2
    k = np.fft.fftfreq(n)
    if axis == 0:    # shift along x by factor * y   (lines are columns j)
        ramp = np.exp(-2j * np.pi * k[:, None] * (factor * c)[None, :])
        return np.real(np.fft.ifft(np.fft.fft(img, axis=0) * ramp, axis=0))
    ramp = np.exp(-2j * np.pi * (factor * c)[:, None] * k[None, :])
    return np.real(np.fft.ifft(np.fft.fft(img, axis=1) * ramp, axis=1))


def rotate_image_exact(img: np.ndarray, phi_deg: float) -> np.ndarray:
    """Counter-clockwise in-plane rotation by three Fourier shears (+ exact
    90-degree turns). Unlike interpolation it does not smooth the noise, so
    averaging statistics stay i.i.d. with unchanged sigma."""
    q = int(np.round(phi_deg / 90.0))
    rem = np.radians(phi_deg - 90.0 * q)
    out = img
    for _ in range(q % 4):
        out = _rot90_ccw(out)
    a, b = -np.tan(rem / 2), np.sin(rem)
    out = _shear_axis(out, a, 0)
    out = _shear_axis(out, b, 1)
    out = _shear_axis(out, a, 0)
    return out


def _rot90_ccw(img: np.ndarray) -> np.ndarray:
    """Exact 90-degree CCW rotation of an (x, y)-indexed array about index n//2."""
    n = img.shape[0]
    out = np.zeros_like(img)
    c = n // 2
    i, j = np.indices(img.shape)
    # new(x, y) = old(R^T (x, y)) with R = 90 deg: old = (y, -x)
    xs, ys = i - c, j - c
    oi, oj = ys + c, -xs + c
    ok = (oi >= 0) & (oi < n) & (oj >= 0) & (oj < n)
    out[ok] = img[oi[ok], oj[ok]]
    return out


# ---------------------------------------------------------------- isodensity render
def mixture_density_fn(weights: dict):
    """rho(pts) = sum_s w_s rho_s(pts) for smooth conformers, e.g. {"A": .55, "B": .45}."""
    def fn(pts):
        return sum(w * sp.density(pts, s, atoms=False) for s, w in weights.items())
    return fn


def voxel_density_fn(vol: np.ndarray, voxel: float):
    """Trilinear lookup into a voxel map (object frame, angstrom coordinates)."""
    n = vol.shape[0]

    def fn(pts):
        idx = (pts / voxel + n // 2).T
        return ndimage.map_coordinates(vol, idx, order=1, mode="constant")
    return fn


def render_density(R: np.ndarray, state: str | None = "A", n: int = 220,
                   extent: float = 64.0, level: float = 0.35, nz: int = 220,
                   light=(-0.45, 0.55, 0.70), density_fn=None) -> dict:
    """Shaded pseudo-3D isodensity view of the smooth phantom (no pseudo-atoms).

    density_fn(pts[..., 3]) -> values overrides the phantom (e.g. a mixture of
    conformers or a reconstructed voxel map); labels still refer to the blobs
    of `state` (pass state=None to label scaffold domains only).

    Ray-marches rho(R^T r) along -z from the viewer and Lambert-shades the first
    crossing of `level`. Returns dict(shade[n, n] in 0..1, mask[n, n] bool,
    label[n, n] = index of the blob that dominates the hit point (-1 outside),
    hit[n, n, 3] object-frame hit points). Arrays are in (x, y) index order.
    """
    ax = np.linspace(-extent, extent, n)
    zs = np.linspace(extent, -extent, nz)
    X, Y = np.meshgrid(ax, ax, indexing="ij")
    blobs = sp.blob_params(state if state is not None else None)
    depth = np.full((n, n), np.nan)
    remaining = np.ones((n, n), bool)
    for z in zs:
        if not remaining.any():
            break
        pts = np.stack([X[remaining], Y[remaining], np.full(remaining.sum(), z)], -1) @ R
        val = (density_fn(pts) if density_fn is not None
               else sp.density(pts, state, atoms=False))
        hit = val > level
        if hit.any():
            ii, jj = np.nonzero(remaining)
            depth[ii[hit], jj[hit]] = z
            remaining[ii[hit], jj[hit]] = False
    mask = ~np.isnan(depth)
    d = np.where(mask, depth, -extent)
    dzdx = np.gradient(d, ax, axis=0)
    dzdy = np.gradient(d, ax, axis=1)
    nrm = np.stack([-dzdx, -dzdy, np.ones_like(d)], -1)
    nrm /= np.linalg.norm(nrm, axis=-1, keepdims=True)
    L = np.asarray(light, float)
    L /= np.linalg.norm(L)
    lam = np.clip(nrm @ L, 0, 1)
    H = (L + np.array([0, 0, 1.0]))
    H /= np.linalg.norm(H)
    spec = np.clip(nrm @ H, 0, 1) ** 24
    shade = np.where(mask, 0.25 + 0.65 * lam + 0.25 * spec, 0.0)
    shade = np.clip(shade, 0, 1)
    hit = np.stack([X, Y, d], -1) @ R          # back to the object frame
    label = np.full((n, n), -1)
    if mask.any():
        contrib = np.stack([sp._eval_blob(b, hit[..., 0], hit[..., 1], hit[..., 2])
                            for b in blobs], -1)
        label = np.where(mask, np.array([b.index for b in blobs])[contrib.argmax(-1)], -1)
    return dict(shade=shade, mask=mask, label=label, hit=hit)
