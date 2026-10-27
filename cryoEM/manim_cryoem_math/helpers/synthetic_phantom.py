"""Deterministic PP2A-inspired synthetic phantom.

The phantom is an abstract density: 15 elongated Gaussian "domains" on a
curved, horseshoe-like scaffold, plus one B-region appendage that can sit in
two places (state "A" or "B"). Each domain also carries a cloud of small
pseudo-atoms so the density has high-resolution content for the FSC demos.
It is NOT an atomic model, a measured PP2A conformation or a cryo-EM map.

All lengths are in angstroms. Volumes are sampled on a cubic grid whose voxel
(i, j, k) sits at r = ((i, j, k) - n//2) * voxel_size, with array axes (x, y, z).
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

SEED = 2027
BOX_N = 128           # default grid size (voxels)
VOXEL_A = 1.0         # default voxel size (angstrom)
ATOMS_PER_DOMAIN = 24
ATOM_SIGMA_A = 1.1


@dataclass(frozen=True)
class Blob:
    """One anisotropic Gaussian: a * exp(-1/2 (r-mu)^T Sigma^-1 (r-mu))."""
    center: np.ndarray        # (3,) angstrom
    axes: np.ndarray          # (3, 3) columns = principal directions (orthonormal)
    sigmas: np.ndarray        # (3,) angstrom, along each principal direction
    amplitude: float
    kind: str                 # "domain" or "appendage"
    index: int                # 0..14 for domains, 15 for the appendage

    @property
    def cov(self) -> np.ndarray:
        return self.axes @ np.diag(self.sigmas ** 2) @ self.axes.T


def _frame(tangent: np.ndarray, up: np.ndarray, twist_deg: float) -> np.ndarray:
    t = tangent / np.linalg.norm(tangent)
    u = up - (up @ t) * t
    u /= np.linalg.norm(u)
    w = np.cross(t, u)
    c, s = np.cos(np.radians(twist_deg)), np.sin(np.radians(twist_deg))
    # tilt the long axis out of the scaffold plane so views are distinguishable
    t2 = c * t + s * u
    u2 = -s * t + c * u
    return np.stack([t2, u2, w], axis=1)


def scaffold_domains() -> list[Blob]:
    """The 15 scaffold domains (deliberately non-symmetric)."""
    blobs = []
    thetas = np.radians(np.linspace(-128.0, 132.0, 15))
    amps = 0.85 + 0.3 * np.abs(np.sin(np.arange(15) * 1.7))
    for j, th in enumerate(thetas):
        radius = 36.0 + 4.0 * np.sin(1.5 * th)          # slightly oval, uneven
        z = 10.0 * np.sin(0.9 * th + 0.4) + 0.05 * np.degrees(th)  # helical rise
        center = np.array([radius * np.cos(th), radius * np.sin(th), z])
        tangent = np.array([-np.sin(th), np.cos(th), 0.12 + 0.002 * np.degrees(th)])
        twist = 28.0 if j % 2 == 0 else -22.0
        axes = _frame(tangent, np.array([0.0, 0.0, 1.0]), twist)
        sigmas = np.array([6.6, 3.6, 4.3]) * (1.0 + 0.04 * np.cos(2.3 * j))
        blobs.append(Blob(center, axes, sigmas, float(amps[j]), "domain", j))
    return blobs


APPENDAGE_CENTERS = {
    # position A: tucked against the inner face near domains 10-12
    "A": np.array([-8.0, 20.0, 16.0]),
    # position B: swung ~22 angstrom towards the open mouth of the horseshoe
    "B": np.array([12.0, 12.0, 4.0]),
}


def appendage(state: str = "A") -> Blob:
    if state not in APPENDAGE_CENTERS:
        raise ValueError(f"state must be 'A' or 'B', got {state!r}")
    axes = _frame(np.array([1.0, 0.4, 0.2]), np.array([0.0, 0.0, 1.0]), 10.0)
    return Blob(APPENDAGE_CENTERS[state], axes, np.array([8.5, 6.0, 5.5]), 1.05,
                "appendage", 15)


def blob_params(state: str | None = "A") -> list[Blob]:
    """Domains (+ appendage unless state is None) — used for drawing in Manim."""
    blobs = scaffold_domains()
    if state is not None:
        blobs.append(appendage(state))
    return blobs


def pseudo_atoms(state: str | None = "A") -> tuple[np.ndarray, np.ndarray]:
    """Deterministic pseudo-atom centres (N, 3) and the blob index they belong to.

    These are the 'toy atomic model' points: sampled inside each blob with a
    fixed seed. They give the density fine (high-frequency) structure.
    """
    rng = np.random.default_rng(SEED)
    pts, owner = [], []
    # the same seeded offsets are drawn for every state; the appendage atoms
    # simply move rigidly with the appendage centre
    for b in scaffold_domains() + [appendage(state or "A")]:
        local = rng.standard_normal((ATOMS_PER_DOMAIN, 3)) * (0.8 * b.sigmas)
        if b.kind == "appendage" and state is None:
            continue
        pts.append(b.center + local @ b.axes.T)
        owner.append(np.full(ATOMS_PER_DOMAIN, b.index))
    return np.concatenate(pts), np.concatenate(owner)


def grid_coords(n: int = BOX_N, voxel: float = VOXEL_A) -> np.ndarray:
    ax = (np.arange(n) - n // 2) * voxel
    return ax


def _eval_blob(b: Blob, X, Y, Z) -> np.ndarray:
    d = np.stack([X - b.center[0], Y - b.center[1], Z - b.center[2]], axis=-1)
    loc = d @ b.axes                                  # coordinates in principal frame
    q = np.sum((loc / b.sigmas) ** 2, axis=-1)
    return b.amplitude * np.exp(-0.5 * q)


def density(points: np.ndarray, state: str | None = "A", atoms: bool = True,
            atom_weight: float = 0.35) -> np.ndarray:
    """Evaluate rho at arbitrary points (..., 3) in angstrom."""
    X, Y, Z = points[..., 0], points[..., 1], points[..., 2]
    out = np.zeros(points.shape[:-1])
    for b in blob_params(state):
        out += _eval_blob(b, X, Y, Z)
    if atoms:
        pos, _ = pseudo_atoms(state)
        for p in pos:
            r2 = (X - p[0]) ** 2 + (Y - p[1]) ** 2 + (Z - p[2]) ** 2
            out += atom_weight * np.exp(-0.5 * r2 / ATOM_SIGMA_A ** 2)
    return out


def _stamp_atoms(vol: np.ndarray, pos: np.ndarray, voxel: float, weight: float):
    n = vol.shape[0]
    rad = int(np.ceil(4 * ATOM_SIGMA_A / voxel))
    offs = np.arange(-rad, rad + 1)
    for p in pos:
        c = np.round(p / voxel).astype(int) + n // 2
        idx = [np.clip(c[a] + offs, 0, n - 1) for a in range(3)]
        gx, gy, gz = [(idx[a] - n // 2) * voxel - p[a] for a in range(3)]
        g = np.exp(-0.5 * (gx[:, None, None] ** 2 + gy[None, :, None] ** 2
                           + gz[None, None, :] ** 2) / ATOM_SIGMA_A ** 2)
        vol[np.ix_(*idx)] += weight * g


def phantom_volume(n: int = BOX_N, voxel: float = VOXEL_A, state: str | None = "A",
                   atoms: bool = True, atom_weight: float = 0.35) -> np.ndarray:
    """Sample the phantom on an n^3 grid (array axes x, y, z)."""
    ax = grid_coords(n, voxel)
    X, Y, Z = np.meshgrid(ax, ax, ax, indexing="ij")
    vol = np.zeros((n, n, n))
    for b in blob_params(state):
        vol += _eval_blob(b, X, Y, Z)
    if atoms:
        pos, _ = pseudo_atoms(state)
        _stamp_atoms(vol, pos, voxel, atom_weight)
    return vol


def ground_truth(n: int = BOX_N, voxel: float = VOXEL_A, **kw) -> np.ndarray:
    """The single ground-truth density used in clips 01-04 and 06 (appendage at A)."""
    return phantom_volume(n, voxel, state="A", **kw)


def conformer(state: str, n: int = BOX_N, voxel: float = VOXEL_A, **kw) -> np.ndarray:
    """The two related densities of clip 05 (appendage at A or at B)."""
    return phantom_volume(n, voxel, state=state, **kw)


def appendage_mask(state: str, n: int = BOX_N, voxel: float = VOXEL_A,
                   nsig: float = 1.5) -> np.ndarray:
    """Boolean mask around the appendage in the given state (for contrast tests)."""
    ax = grid_coords(n, voxel)
    X, Y, Z = np.meshgrid(ax, ax, ax, indexing="ij")
    b = appendage(state)
    d = np.stack([X - b.center[0], Y - b.center[1], Z - b.center[2]], axis=-1) @ b.axes
    return np.sum((d / b.sigmas) ** 2, axis=-1) < nsig ** 2
