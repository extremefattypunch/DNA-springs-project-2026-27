"""Simplified weak-phase CTF and additive-noise model.

CTF(k) = -sin(chi(k)),   chi(k) = pi*lambda*df*k^2 - (pi/2)*Cs*lambda^3*k^4 + w,
with w = arcsin(amplitude contrast) and df > 0 meaning underfocus. This is an
oscillating, sign-reversing contrast transfer (with zeros), not a blur.
Envelope functions, astigmatism and the detector MTF are deliberately omitted.
"""
from __future__ import annotations

import numpy as np

from . import synthetic_phantom as sp

KV = 300.0
CS_MM = 2.7
AMP_CONTRAST = 0.07


def electron_wavelength(kv: float = KV) -> float:
    """Relativistic electron wavelength in angstrom."""
    v = kv * 1e3
    return 12.2643 / np.sqrt(v * (1 + 0.978466e-6 * v))


def ctf_1d(k: np.ndarray, defocus_um: float = 1.2, kv: float = KV, cs_mm: float = CS_MM,
           amp: float = AMP_CONTRAST, bfactor: float = 0.0) -> np.ndarray:
    """CTF as a function of spatial frequency |k| in 1/angstrom.

    bfactor > 0 multiplies by exp(-B k^2 / 4) (optional, illustrative envelope).
    """
    lam = electron_wavelength(kv)
    df = defocus_um * 1e4
    cs = cs_mm * 1e7
    chi = np.pi * lam * df * k ** 2 - 0.5 * np.pi * cs * lam ** 3 * k ** 4 + np.arcsin(amp)
    out = -np.sin(chi)
    if bfactor:
        out = out * np.exp(-bfactor * k ** 2 / 4)
    return out


def freq_grid(n: int = sp.BOX_N, voxel: float = sp.VOXEL_A, shifted: bool = False):
    k = np.fft.fftfreq(n, d=voxel)
    if shifted:
        k = np.fft.fftshift(k)
    KX, KY = np.meshgrid(k, k, indexing="ij")
    return KX, KY


def ctf_2d(n: int = sp.BOX_N, voxel: float = sp.VOXEL_A, shifted: bool = False, **kw) -> np.ndarray:
    KX, KY = freq_grid(n, voxel, shifted)
    return ctf_1d(np.hypot(KX, KY), **kw)


def apply_ctf(img: np.ndarray, voxel: float = sp.VOXEL_A, **kw) -> np.ndarray:
    """Multiply the image's Fourier transform by the CTF (unshifted FFT layout)."""
    c = ctf_2d(img.shape[0], voxel, shifted=False, **kw)
    return np.real(np.fft.ifft2(np.fft.fft2(img) * c))


def first_zero(defocus_um: float = 1.2, kmax: float = 0.5, **kw) -> float:
    """Spatial frequency (1/angstrom) of the first CTF zero crossing."""
    k = np.linspace(0, kmax, 200001)
    c = ctf_1d(k, defocus_um, **kw)
    i = np.nonzero(np.diff(np.sign(c)) != 0)[0]
    return float(k[i[0]]) if len(i) else np.nan


def add_noise(img: np.ndarray, sigma: float, rng: np.random.Generator) -> np.ndarray:
    """Independent zero-mean Gaussian noise (an approximation to detector/shot noise)."""
    return img + sigma * rng.standard_normal(img.shape)


def snr_amplitude(signal: np.ndarray, sigma: float) -> float:
    """SNR convention used in the clips: signal RMS (std over the box) / noise std."""
    return float(np.std(signal) / sigma) if sigma > 0 else np.inf


def sigma_for_snr(signal: np.ndarray, snr: float) -> float:
    return float(np.std(signal) / snr)
