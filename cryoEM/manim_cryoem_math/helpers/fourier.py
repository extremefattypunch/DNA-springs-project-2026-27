"""Centred FFTs, central slices, direct Fourier inversion and FSC.

"Centred" means the array origin sits at index n//2 in real and Fourier space
(ifftshift -> fft -> fftshift), matching point_constellation.grid_coords.
"""
from __future__ import annotations

import numpy as np
from scipy import ndimage


def fft2c(x: np.ndarray) -> np.ndarray:
    return np.fft.fftshift(np.fft.fft2(np.fft.ifftshift(x)))


def ifft2c(F: np.ndarray) -> np.ndarray:
    return np.fft.fftshift(np.fft.ifft2(np.fft.ifftshift(F)))


def fft3c(x: np.ndarray) -> np.ndarray:
    return np.fft.fftshift(np.fft.fftn(np.fft.ifftshift(x)))


def ifft3c(F: np.ndarray) -> np.ndarray:
    return np.fft.fftshift(np.fft.ifftn(np.fft.ifftshift(F)))


def log_magnitude(F: np.ndarray) -> np.ndarray:
    """Display convention for Fourier panels: log(1 + |F|) scaled to 0..1.
    A magnitude plot is a visualization; the theorem is about complex values."""
    m = np.log1p(np.abs(F))
    return m / m.max() if m.max() > 0 else m


def pad(a: np.ndarray, factor: int = 2) -> np.ndarray:
    n = a.shape[0]
    N = n * factor
    out = np.zeros((N,) * a.ndim, dtype=a.dtype)
    o = N // 2 - n // 2
    out[(slice(o, o + n),) * a.ndim] = a
    return out


def crop_center(a: np.ndarray, n: int) -> np.ndarray:
    N = a.shape[0]
    o = N // 2 - n // 2
    return a[(slice(o, o + n),) * a.ndim]


def central_slice(F3: np.ndarray, R: np.ndarray, n_out: int | None = None,
                  order: int = 1) -> np.ndarray:
    """Sample a centred 3D transform on the plane {R^T (kx, ky, 0)}."""
    N = F3.shape[0]
    n_out = n_out or N
    k = np.arange(n_out) - n_out // 2
    KX, KY = np.meshgrid(k, k, indexing="ij")
    kpl = np.stack([KX.ravel(), KY.ravel(), np.zeros(KX.size)])
    src = R.T @ kpl + N // 2
    re = ndimage.map_coordinates(F3.real, src, order=order, mode="constant")
    im = ndimage.map_coordinates(F3.imag, src, order=order, mode="constant")
    return (re + 1j * im).reshape(n_out, n_out)


def plane_normal(R: np.ndarray) -> np.ndarray:
    """Normal of the central plane R^T (kx, ky, 0): the viewing direction R^T e_z."""
    return R.T @ np.array([0.0, 0.0, 1.0])


def reconstruct_direct_fourier(images: np.ndarray, rotations: np.ndarray, factor: int = 2,
                               floor: float = 1e-3) -> np.ndarray:
    """Idealised direct Fourier inversion with KNOWN poses (no CTF, no noise model).

    Each image's padded 2D FFT is trilinearly gridded onto its central plane
    R^T (kx, ky, 0); the sum is divided by the sampling weight. Unsampled
    Fourier voxels stay zero, so missing views leave missing information.
    Images are line integrals (sum * spacing) -> pass images / spacing to get
    the field back on its own scale.
    """
    n = images.shape[1]
    N = n * factor
    acc = np.zeros((N, N, N), complex)
    wsum = np.zeros((N, N, N))
    k = np.arange(N) - N // 2
    KX, KY = np.meshgrid(k, k, indexing="ij")
    inside = (KX ** 2 + KY ** 2) <= (N // 2 - 1) ** 2
    kpl = np.stack([KX[inside], KY[inside], np.zeros(inside.sum())])
    for img, R in zip(images, rotations):
        F = fft2c(pad(img, factor))[inside]
        c = R.T @ kpl + N // 2
        f = np.floor(c).astype(int)
        d = c - f
        for dx in (0, 1):
            for dy in (0, 1):
                for dz in (0, 1):
                    w = ((d[0] if dx else 1 - d[0]) * (d[1] if dy else 1 - d[1])
                         * (d[2] if dz else 1 - d[2]))
                    ix, iy, iz = f[0] + dx, f[1] + dy, f[2] + dz
                    ok = (ix >= 0) & (ix < N) & (iy >= 0) & (iy < N) & (iz >= 0) & (iz < N)
                    np.add.at(acc, (ix[ok], iy[ok], iz[ok]), w[ok] * F[ok])
                    np.add.at(wsum, (ix[ok], iy[ok], iz[ok]), w[ok])
    F3 = np.where(wsum > floor, acc / np.maximum(wsum, floor), 0)
    return crop_center(np.real(ifft3c(F3)), n)


def sampling_coverage(rotations: np.ndarray, N: int = 32) -> np.ndarray:
    """Boolean N^3 map of Fourier voxels touched by the central planes."""
    cov = np.zeros((N, N, N), bool)
    k = np.arange(N) - N // 2
    KX, KY = np.meshgrid(k, k, indexing="ij")
    inside = (KX ** 2 + KY ** 2) <= (N // 2 - 1) ** 2
    kpl = np.stack([KX[inside], KY[inside], np.zeros(inside.sum())])
    for R in rotations:
        c = np.round(R.T @ kpl + N // 2).astype(int)
        ok = np.all((c >= 0) & (c < N), axis=0)
        cov[c[0, ok], c[1, ok], c[2, ok]] = True
    return cov


def ball_fraction(cov: np.ndarray) -> float:
    """Fraction of the Fourier ball |k| < N/2 - 1 that is sampled."""
    N = cov.shape[0]
    ball = shell_index(cov.shape) < N // 2 - 1
    return float(cov[ball].mean())


# ---------------------------------------------------------------- FSC
def shell_index(shape) -> np.ndarray:
    """Integer Fourier-shell radius for a centred array of the given shape."""
    grids = np.meshgrid(*[np.arange(s) - s // 2 for s in shape], indexing="ij")
    return np.round(np.sqrt(sum(g.astype(float) ** 2 for g in grids))).astype(int)


def fsc(A: np.ndarray, B: np.ndarray, spacing: float = 1.0):
    """Fourier shell (ring, in 2D) correlation of two real arrays of equal shape.

        FSC(s) = Re sum F_A F_B^* / sqrt(sum |F_A|^2 sum |F_B|^2)   over shell s

    Returns (s, fsc, counts) for shells 1..n/2-1; s is in cycles per unit length.
    Shells with zero power in either map return 0 instead of dividing by zero.
    """
    FA = fft3c(A) if A.ndim == 3 else fft2c(A)
    FB = fft3c(B) if B.ndim == 3 else fft2c(B)
    n = A.shape[0]
    sh = shell_index(A.shape).ravel()
    nshell = n // 2
    keep = sh < nshell
    sh = sh[keep]
    fa, fb = FA.ravel()[keep], FB.ravel()[keep]
    num = np.bincount(sh, weights=np.real(fa * np.conj(fb)), minlength=nshell)
    pa = np.bincount(sh, weights=np.abs(fa) ** 2, minlength=nshell)
    pb = np.bincount(sh, weights=np.abs(fb) ** 2, minlength=nshell)
    cnt = np.bincount(sh, minlength=nshell)
    den = np.sqrt(pa * pb)
    out = np.divide(num, den, out=np.zeros_like(num), where=den > 1e-300)
    s = np.arange(nshell) / (n * spacing)
    return s[1:], out[1:], cnt[1:]


def crossing(s: np.ndarray, curve: np.ndarray, threshold: float) -> float:
    """First frequency where the curve falls below threshold (linear interpolation)."""
    below = np.nonzero(curve < threshold)[0]
    if len(below) == 0:
        return float(s[-1])
    i = below[0]
    if i == 0:
        return float(s[0])
    s0, s1, c0, c1 = s[i - 1], s[i], curve[i - 1], curve[i]
    return float(s0 + (threshold - c0) * (s1 - s0) / (c1 - c0))


def soft_mask(n: int, extent: float, center, radius: float, edge: float = 0.3) -> np.ndarray:
    """Spherical mask with a raised-cosine edge (for masked / region-wise FSC)."""
    from .point_constellation import grid_coords
    ax = grid_coords(n, extent)
    X, Y, Z = np.meshgrid(ax, ax, ax, indexing="ij")
    r = np.sqrt((X - center[0]) ** 2 + (Y - center[1]) ** 2 + (Z - center[2]) ** 2)
    m = np.clip((radius + edge - r) / edge, 0, 1)
    return 0.5 - 0.5 * np.cos(np.pi * m)
