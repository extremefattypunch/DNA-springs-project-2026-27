"""Rotations, translations and ideal parallel projections.

Convention used everywhere in this project (verified in tests/test_math.py):

    rotated source      x -> R x            (fields: rho_R(r) = rho(R^T r))
    projection          (P_R rho)(x, y) = integral of rho(R^T (x, y, z)) dz
    projected dots      p_j = (R x_j)_{x,y}
    Fourier slice       FT2[P_R rho](kx, ky) = FT3[rho](R^T (kx, ky, 0))

R rotates the source into the detector frame; projection runs along z.
The viewing direction in the source frame is v = R^T e_z (2 DOF); the third
DOF, psi, spins the image in the detector plane.
"""
from __future__ import annotations

import numpy as np
from scipy import ndimage

from . import point_constellation as pc


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


def rot2(phi_deg: float) -> np.ndarray:
    """2D counter-clockwise rotation (x right, y up)."""
    a = np.radians(phi_deg)
    return np.array([[np.cos(a), -np.sin(a)], [np.sin(a), np.cos(a)]])


def rotation_from_euler(alpha: float, theta: float, psi: float, degrees: bool = True
                        ) -> np.ndarray:
    """ZYZ Euler angles: R = Rz(psi) Ry(theta) Rz(alpha).

    (alpha, theta) fix the viewing direction v = R^T e_z (2 DOF);
    psi is the in-plane rotation about the projection axis (3rd DOF).
    """
    if degrees:
        alpha, theta, psi = np.radians([alpha, theta, psi])
    return rot_z(psi) @ rot_y(theta) @ rot_z(alpha)


def viewing_direction(R: np.ndarray) -> np.ndarray:
    """Projection direction expressed in the source frame: v = R^T e_z."""
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


def biased_rotations(n: int, rng: np.random.Generator, max_theta_deg: float = 15.0
                     ) -> np.ndarray:
    """Preferred orientation: viewing directions within max_theta of +z."""
    al = rng.uniform(0, 360, n)
    th = rng.uniform(0, max_theta_deg, n)
    ps = rng.uniform(0, 360, n)
    return np.stack([rotation_from_euler(a, b, c) for a, b, c in zip(al, th, ps)])


# ---------------------------------------------------------------- dots
def project_points(R: np.ndarray, X: np.ndarray, t=(0.0, 0.0)) -> np.ndarray:
    """Projected 2D dot positions p_j = (R x_j)_{xy} + t, shape (m, 2)."""
    return (np.asarray(X) @ R.T)[:, :2] + np.asarray(t)


def transform_2d(P: np.ndarray, phi_deg: float = 0.0, t=(0.0, 0.0)) -> np.ndarray:
    """In-plane rotation by phi about the origin, then translation by t."""
    return np.asarray(P) @ rot2(phi_deg).T + np.asarray(t)


# ---------------------------------------------------------------- images
def image_coords(n: int = pc.IMAGE_N, extent: float = pc.IMAGE_EXTENT):
    ax = pc.grid_coords(n, extent)
    return np.meshgrid(ax, ax, indexing="ij")     # X[i, j] = x_i, Y[i, j] = y_j


def dots_image(P: np.ndarray, amplitudes=None, sigmas=None, n: int = pc.IMAGE_N,
               extent: float = pc.IMAGE_EXTENT) -> np.ndarray:
    """2D image of Gaussian kernels at 2D positions P (m, 2), indexed [x, y].

    With the 3D kernel weights (a_j, sigma_j) this is the exact line integral of
    rho: a 3D isotropic Gaussian projects to a 2D Gaussian of mass
    a_j (2 pi)^{3/2} sigma_j^3. The constant sqrt(2 pi) sigma_j is kept.
    """
    m = len(P)
    a = pc.AMPLITUDES[:m] if amplitudes is None else np.asarray(amplitudes)
    s = pc.SIGMAS[:m] if sigmas is None else np.asarray(sigmas)
    X, Y = image_coords(n, extent)
    img = np.zeros((n, n))
    for (px, py), aj, sj in zip(P, a, s):
        img += aj * np.sqrt(2 * np.pi) * sj * np.exp(-0.5 * ((X - px) ** 2 + (Y - py) ** 2) / sj ** 2)
    return img


def project_image(R: np.ndarray, state: str = "A", t=(0.0, 0.0), n: int = pc.IMAGE_N,
                  extent: float = pc.IMAGE_EXTENT, pts_override=None) -> np.ndarray:
    """Exact projection P_R rho (translated by t) of the Gaussian-kernel field."""
    X = pc.points(state) if pts_override is None else pts_override
    return dots_image(project_points(R, X, t), n=n, extent=extent)


def rotate_field(vol: np.ndarray, R: np.ndarray, order: int = 1) -> np.ndarray:
    """rho_R(r) = rho(R^T r) on the same grid (rotation about the grid centre)."""
    n = vol.shape[0]
    c = n // 2
    idx = np.indices(vol.shape).reshape(3, -1) - c
    src = R.T @ idx + c
    return ndimage.map_coordinates(vol, src, order=order, mode="constant").reshape(vol.shape)


def project_field(vol: np.ndarray, R: np.ndarray, extent: float = pc.FIELD_EXTENT) -> np.ndarray:
    """Numerical line integral of a sampled field after rotation (sum along z)."""
    return rotate_field(vol, R).sum(axis=2) * pc.spacing(vol.shape[0], extent)


# ---------------------------------------------------------------- 2D transforms
def shift_image(img: np.ndarray, t, extent: float = pc.IMAGE_EXTENT) -> np.ndarray:
    """Translate by t = (tx, ty) via the Fourier phase ramp e^{-2 pi i k.t}."""
    n = img.shape[0]
    k = np.fft.fftfreq(n, d=pc.spacing(n, extent))
    KX, KY = np.meshgrid(k, k, indexing="ij")
    ramp = np.exp(-2j * np.pi * (KX * t[0] + KY * t[1]))
    return np.real(np.fft.ifft2(np.fft.fft2(img) * ramp))


def _shear_axis(img: np.ndarray, factor: float, axis: int) -> np.ndarray:
    """Shift every line along `axis` by factor * (coordinate along the other axis)
    using exact Fourier phase ramps (unitary: preserves white-noise variance)."""
    n = img.shape[0]
    c = np.arange(n) - n // 2
    k = np.fft.fftfreq(n)
    if axis == 0:
        ramp = np.exp(-2j * np.pi * k[:, None] * (factor * c)[None, :])
        return np.real(np.fft.ifft(np.fft.fft(img, axis=0) * ramp, axis=0))
    ramp = np.exp(-2j * np.pi * (factor * c)[:, None] * k[None, :])
    return np.real(np.fft.ifft(np.fft.fft(img, axis=1) * ramp, axis=1))


def _rot90_ccw(img: np.ndarray) -> np.ndarray:
    """Exact 90-degree CCW rotation of an [x, y]-indexed array about index n//2."""
    n = img.shape[0]
    out = np.zeros_like(img)
    c = n // 2
    i, j = np.indices(img.shape)
    oi, oj = (j - c) + c, -(i - c) + c
    ok = (oi >= 0) & (oi < n) & (oj >= 0) & (oj < n)
    out[ok] = img[oi[ok], oj[ok]]
    return out


def rotate_image(img: np.ndarray, phi_deg: float) -> np.ndarray:
    """Counter-clockwise in-plane rotation by three Fourier shears (+ exact
    90-degree turns). Unlike interpolation it does not smooth the noise, so
    averages of rotated noisy images keep i.i.d. noise with unchanged sigma."""
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
