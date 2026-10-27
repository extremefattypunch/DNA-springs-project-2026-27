"""The abstract source: five arbitrary, dimensionless points in 3D.

The coordinates have no meaning beyond being strongly asymmetric, so that
rotations, projections, alignments and state changes are visibly distinct.
On screen they are only ever drawn as dots. For numerical images and maps the
same points carry isotropic Gaussian kernels,

    rho(r) = sum_j a_j exp(-|r - x_j|^2 / (2 sigma_j^2)),

which are never displayed as a 3D surface.

State "B" moves only dot 5 (index 4); the other four are shared.
"""
from __future__ import annotations

import numpy as np

SEED = 2027

POINTS = np.array([
    [-1.1, -0.4, 0.2],
    [0.7, -0.9, 0.6],
    [0.2, 0.6, -0.7],
    [-0.6, 0.9, 0.4],
    [1.0, 0.5, 0.9],
])
MOVED = 4                                   # index of the dot that differs in state B
POINT_B5 = np.array([1.15, -0.35, 0.3])     # its state-B position
AMPLITUDES = np.array([1.0, 0.8, 1.2, 0.9, 0.7])
SIGMAS = np.array([0.16, 0.20, 0.14, 0.18, 0.22])

# The asymmetric 2D pattern used by the alignment clip: the (x, y) coordinates
# of dots 1-4 (same notation x_1..x_4, viewed down z).
PATTERN_2D = POINTS[:4, :2].copy()
PATTERN_AMPLITUDES = AMPLITUDES[:4].copy()
PATTERN_SIGMAS = SIGMAS[:4].copy()

# numerical grids (dimensionless units)
FIELD_N = 64            # 3D grid edge
FIELD_EXTENT = 3.0      # grid covers [-3, 3)^3
IMAGE_N = 64            # 2D image edge
IMAGE_EXTENT = 2.0      # images cover [-2, 2)^2


def points(state: str = "A") -> np.ndarray:
    """(5, 3) coordinates for state "A" (default) or "B" (dot 5 moved)."""
    p = POINTS.copy()
    if state == "B":
        p[MOVED] = POINT_B5
    elif state != "A":
        raise ValueError(f"unknown state {state!r}")
    return p


def spacing(n: int, extent: float) -> float:
    return 2 * extent / n


def grid_coords(n: int, extent: float) -> np.ndarray:
    """Centred 1D coordinates: index i -> (i - n//2) * spacing (origin at n//2)."""
    return (np.arange(n) - n // 2) * spacing(n, extent)


def density_at(pts: np.ndarray, state: str = "A", pts_override: np.ndarray | None = None
               ) -> np.ndarray:
    """rho evaluated at points pts[..., 3]."""
    X = points(state) if pts_override is None else np.asarray(pts_override, float)
    out = np.zeros(pts.shape[:-1])
    for x, a, s in zip(X, AMPLITUDES, SIGMAS):
        r2 = np.sum((pts - x) ** 2, axis=-1)
        out += a * np.exp(-0.5 * r2 / s ** 2)
    return out


def field(state: str = "A", n: int = FIELD_N, extent: float = FIELD_EXTENT,
          pts_override: np.ndarray | None = None) -> np.ndarray:
    """rho sampled on an n^3 centred grid, indexed [x, y, z]."""
    ax = grid_coords(n, extent)
    G = np.stack(np.meshgrid(ax, ax, ax, indexing="ij"), -1)
    return density_at(G, state, pts_override)
