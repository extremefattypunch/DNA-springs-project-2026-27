"""Centred FFTs, central slices, direct Fourier reconstruction and FSC.

"Centred" means the array origin is at index n//2 in both real and Fourier
space (ifftshift -> fft -> fftshift), matching synthetic_phantom.grid_coords.
"""
from __future__ import annotations

import numpy as np
from scipy import ndimage

from . import synthetic_phantom as sp


def fft2c(x: np.ndarray) -> np.ndarray:
    return np.fft.fftshift(np.fft.fft2(np.fft.ifftshift(x)))


def ifft2c(F: np.ndarray) -> np.ndarray:
    return np.fft.fftshift(np.fft.ifft2(np.fft.ifftshift(F)))


def fft3c(x: np.ndarray) -> np.ndarray:
    return np.fft.fftshift(np.fft.fftn(np.fft.ifftshift(x)))


def ifft3c(F: np.ndarray) -> np.ndarray:
    return np.fft.fftshift(np.fft.ifftn(np.fft.ifftshift(F)))


def log_magnitude(F: np.ndarray) -> np.ndarray:
    """Display convention for Fourier panels: log(1 + |F|), scaled to 0..1."""
    m = np.log1p(np.abs(F))
    return m / m.max() if m.max() > 0 else m


def pad_volume(vol: np.ndarray, factor: int = 2) -> np.ndarray:
    n = vol.shape[0]
    N = n * factor
    out = np.zeros((N, N, N), dtype=vol.dtype)
    o = N // 2 - n // 2
    out[o:o + n, o:o + n, o:o + n] = vol
    return out


def crop_center(a: np.ndarray, n: int) -> np.ndarray:
    N = a.shape[0]
    o = N // 2 - n // 2
    sl = slice(o, o + n)
    return a[(sl,) * a.ndim]


def central_slice(F3: np.ndarray, R: np.ndarray, n_out: int | None = None,
                  order: int = 1) -> np.ndarray:
    """Sample F3 on the plane {R^T (kx, ky, 0)} (Fourier slice theorem).

    F3 is a centred 3D transform of size N. The returned 2D slice has size
    n_out (default N) on the same frequency spacing.
    """
    N = F3.shape[0]
    n_out = n_out or N
    k = np.arange(n_out) - n_out // 2
    KX, KY = np.meshgrid(k, k, indexing="ij")
    kpl = np.stack([KX.ravel(), KY.ravel(), np.zeros(KX.size)])
    src = R.T @ kpl + N // 2
    re = ndimage.map_coordinates(F3.real, src, order=order, mode="constant")
    im = ndimage.map_coordinates(F3.imag, src, order=order, mode="constant")
    return (re + 1j * im).reshape(n_out, n_out)


def project_fourier(F3_padded: np.ndarray, R: np.ndarray, n: int) -> np.ndarray:
    """Projection of an (oversampled) volume through its central slice.

    F3_padded = fft3c(pad_volume(vol, p)); returns the n x n real-space projection.
    """
    N = F3_padded.shape[0]
    sl = central_slice(F3_padded, R)
    img = np.real(ifft2c(sl))
    return crop_center(img, n)


def reconstruct_direct_fourier(images: np.ndarray, rotations: np.ndarray,
                               pad: int = 2, weights_floor: float = 1e-3,
                               weights: np.ndarray | None = None) -> np.ndarray:
    """Idealised direct Fourier inversion with KNOWN poses (no CTF, no noise model).

    Each image's (padded) 2D FFT is trilinearly gridded onto the central plane
    R^T (kx, ky, 0) of a 3D Fourier grid; the sum is divided by the sampling
    weight; the inverse 3D FFT is cropped back to the original box. Fourier
    voxels never sampled stay zero, so missing views leave missing information.
    weights (optional, per image) supports soft class assignments.
    """
    n = images.shape[1]
    N = n * pad
    acc = np.zeros((N, N, N), complex)
    wsum = np.zeros((N, N, N))
    k = np.arange(N) - N // 2
    KX, KY = np.meshgrid(k, k, indexing="ij")
    inside = (KX ** 2 + KY ** 2) <= (N // 2 - 1) ** 2
    kpl = np.stack([KX[inside], KY[inside], np.zeros(inside.sum())])
    o = N // 2 - n // 2
    for idx, (img, R) in enumerate(zip(images, rotations)):
        w_img = 1.0 if weights is None else float(weights[idx])
        if w_img <= 0:
            continue
        p = np.zeros((N, N))
        p[o:o + n, o:o + n] = img
        F = fft2c(p)[inside]
        c = R.T @ kpl + N // 2
        f = np.floor(c).astype(int)
        d = c - f
        for dx in (0, 1):
            for dy in (0, 1):
                for dz in (0, 1):
                    w = ((d[0] if dx else 1 - d[0]) * (d[1] if dy else 1 - d[1])
                         * (d[2] if dz else 1 - d[2])) * w_img
                    ix, iy, iz = f[0] + dx, f[1] + dy, f[2] + dz
                    ok = (ix >= 0) & (ix < N) & (iy >= 0) & (iy < N) & (iz >= 0) & (iz < N)
                    np.add.at(acc, (ix[ok], iy[ok], iz[ok]), w[ok] * F[ok])
                    np.add.at(wsum, (ix[ok], iy[ok], iz[ok]), w[ok])
    F3 = np.where(wsum > weights_floor, acc / np.maximum(wsum, weights_floor), 0)
    # numpy's unnormalised forward FFTs make FT2(sum_z rho) == FT3(rho)|_{kz=0},
    # so the normalised inverse returns density on the original voxel scale
    return crop_center(np.real(ifft3c(F3)), n)


def sampling_coverage(rotations: np.ndarray, N: int = 48) -> np.ndarray:
    """Boolean N^3 map of Fourier voxels touched by the central planes (for display)."""
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


# ---------------------------------------------------------------- FSC
def shell_index(shape) -> np.ndarray:
    """Integer Fourier-shell radius for a centred array of the given shape."""
    grids = np.meshgrid(*[np.arange(s) - s // 2 for s in shape], indexing="ij")
    return np.round(np.sqrt(sum(g.astype(float) ** 2 for g in grids))).astype(int)


def fsc(A: np.ndarray, B: np.ndarray, voxel: float = sp.VOXEL_A):
    """Fourier shell (ring, in 2D) correlation of two real maps of equal shape.

    Returns (s, fsc, counts): s = shell spatial frequency in 1/angstrom for
    shells 1..n/2-1, the correlation, and the number of Fourier voxels per shell.
    Shells with zero power in either map return 0 rather than dividing by zero.
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
    out = np.divide(num, den, out=np.zeros_like(num), where=den > 0)
    s = np.arange(nshell) / (n * voxel)
    return s[1:], out[1:], cnt[1:]


def crossing(s: np.ndarray, curve: np.ndarray, threshold: float) -> float:
    """First spatial frequency where the curve falls below threshold (linear interp)."""
    below = np.nonzero(curve < threshold)[0]
    if len(below) == 0:
        return float(s[-1])
    i = below[0]
    if i == 0:
        return float(s[0])
    s0, s1, c0, c1 = s[i - 1], s[i], curve[i - 1], curve[i]
    return float(s0 + (threshold - c0) * (s1 - s0) / (c1 - c0))


def target_fsc(s: np.ndarray, s_half: float, width: float = 0.025) -> np.ndarray:
    """Smooth logistic FSC-shaped curve passing 0.143 at s_half."""
    # logistic L(s) = 1/(1+exp((s-s0)/w)), solve L(s_half)=0.143
    s0 = s_half - width * np.log(1 / 0.143 - 1)
    return 1.0 / (1.0 + np.exp((s - s0) / width))


def shaped_noise(signal: np.ndarray, fsc_target_fn, rng: np.random.Generator,
                 voxel: float = sp.VOXEL_A) -> np.ndarray:
    """Real-space noise whose per-shell power N(s) = S(s) (1/FSC(s) - 1).

    For two maps signal + independent noise of this spectrum, the expected FSC
    is S/(S+N) = FSC_target(s). Used to make illustrative half-maps.
    """
    n = signal.shape[0]
    F = fft3c(signal) if signal.ndim == 3 else fft2c(signal)
    sh = shell_index(signal.shape)
    nsh = sh.max() + 1
    S = np.bincount(sh.ravel(), weights=np.abs(F.ravel()) ** 2, minlength=nsh)
    cnt = np.maximum(np.bincount(sh.ravel(), minlength=nsh), 1)
    S_per = S / cnt
    s_of_shell = np.arange(nsh) / (n * voxel)
    t = np.clip(fsc_target_fn(s_of_shell), 1e-4, 0.9999)
    Npow = S_per * (1 / t - 1)
    white = rng.standard_normal(signal.shape)
    W = fft3c(white) if signal.ndim == 3 else fft2c(white)
    # white noise has E|W|^2 = number of voxels per Fourier sample
    W *= np.sqrt(Npow[sh] / signal.size)
    # the radial filter keeps Hermitian symmetry, so the inverse is (numerically) real
    return np.real(ifft3c(W) if signal.ndim == 3 else ifft2c(W))
