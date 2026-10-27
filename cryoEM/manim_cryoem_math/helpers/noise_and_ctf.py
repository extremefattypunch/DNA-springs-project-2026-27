"""Illustrative, dimensionless CTF and deterministic additive noise.

    CTF(k) = -sin(chi(k)),   chi(k) = pi * D * k^2 - (pi/2) * C * k^4 + w

k is a dimensionless spatial frequency (cycles per unit length of the dot
coordinates). D plays the role of defocus, C of spherical aberration and
w = arcsin(0.07) of a small amplitude-contrast phase. The parameters are
chosen so that zeros and sign reversals fall inside the frequency band of the
64-pixel images; they are not fitted to any instrument. This is an
oscillating, sign-changing transfer function, not a blur. Envelopes,
astigmatism and detector response are omitted.
"""
from __future__ import annotations

import numpy as np

from . import point_constellation as pc

D_DEFAULT = 0.12
C_DEFAULT = 0.0006
W_AMP = float(np.arcsin(0.07))


def ctf_1d(k: np.ndarray, D: float = D_DEFAULT, C: float = C_DEFAULT, w: float = W_AMP,
           envelope: float = 0.0) -> np.ndarray:
    """CTF as a function of |k|; envelope > 0 multiplies by exp(-envelope k^2)."""
    k = np.asarray(k, float)
    chi = np.pi * D * k ** 2 - 0.5 * np.pi * C * k ** 4 + w
    out = -np.sin(chi)
    if envelope:
        out = out * np.exp(-envelope * k ** 2)
    return out


def freq_grid(n: int = pc.IMAGE_N, extent: float = pc.IMAGE_EXTENT, shifted: bool = False):
    k = np.fft.fftfreq(n, d=pc.spacing(n, extent))
    if shifted:
        k = np.fft.fftshift(k)
    return np.meshgrid(k, k, indexing="ij")


def ctf_2d(n: int = pc.IMAGE_N, extent: float = pc.IMAGE_EXTENT, shifted: bool = False,
           **kw) -> np.ndarray:
    KX, KY = freq_grid(n, extent, shifted)
    return ctf_1d(np.hypot(KX, KY), **kw)


def apply_ctf(img: np.ndarray, extent: float = pc.IMAGE_EXTENT, **kw) -> np.ndarray:
    """Multiply the image's Fourier transform by the CTF."""
    c = ctf_2d(img.shape[0], extent, shifted=False, **kw)
    return np.real(np.fft.ifft2(np.fft.fft2(img) * c))


def apply_ctf_1d(signal: np.ndarray, dx: float, **kw) -> np.ndarray:
    """1D analogue: multiply the FFT of a sampled profile by CTF(|k|)."""
    k = np.fft.fftfreq(len(signal), d=dx)
    return np.real(np.fft.ifft(np.fft.fft(signal) * ctf_1d(np.abs(k), **kw)))


def zeros(kmax: float, **kw) -> np.ndarray:
    """Frequencies of CTF zero crossings in (0, kmax]."""
    k = np.linspace(0, kmax, 200001)
    c = ctf_1d(k, **kw)
    i = np.nonzero(np.diff(np.sign(c)) != 0)[0]
    return k[i] - c[i] * (k[i + 1] - k[i]) / (c[i + 1] - c[i])


def rng(offset: int = 0) -> np.random.Generator:
    """Deterministic generator: default_rng(2027 + offset)."""
    return np.random.default_rng(pc.SEED + offset)


def add_noise(x: np.ndarray, sigma: float, gen: np.random.Generator) -> np.ndarray:
    """Independent zero-mean Gaussian noise. Real detector/counting noise is
    more complicated; this is the textbook assumption behind the sqrt(N) law."""
    return x + sigma * gen.standard_normal(np.shape(x))


def snr_amplitude(signal: np.ndarray, sigma: float) -> float:
    """Amplitude SNR used throughout: RMS of the (mean-removed) signal / noise std."""
    return float(np.std(signal) / sigma) if sigma > 0 else np.inf


def sigma_for_snr(signal: np.ndarray, snr: float) -> float:
    return float(np.std(signal) / snr)
