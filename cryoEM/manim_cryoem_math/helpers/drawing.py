"""Light-weight 3D Manim geometry for the phantom (ThreeDScene clips).

The phantom is drawn as a handful of translucent ellipsoid Surfaces (one per
Gaussian blob, ~1 sigma iso-surface scaled by `nsig`), never as voxels.
"""
from __future__ import annotations

import numpy as np
from manim import Surface, VGroup

from . import styling as st
from . import synthetic_phantom as sp

A_TO_UNITS = 0.032     # 100 angstrom ~ 3.2 scene units


def ellipsoid(center, axes, radii, color=st.TRUTH, opacity=0.55, resolution=(14, 10),
              scale=A_TO_UNITS) -> Surface:
    c = np.asarray(center) * scale
    Ax = np.asarray(axes) * (np.asarray(radii) * scale)[None, :]

    def f(u, v):
        p = np.array([np.cos(u) * np.sin(v), np.sin(u) * np.sin(v), np.cos(v)])
        return c + Ax @ p

    s = Surface(f, u_range=[0, 2 * np.pi], v_range=[0, np.pi], resolution=resolution,
                fill_opacity=opacity, stroke_width=0.3, stroke_opacity=0.25,
                checkerboard_colors=[color, color])
    s.set_fill(color, opacity=opacity)
    s.set_stroke(color, width=0.3, opacity=0.25)
    return s


def phantom_3d(state: str | None = "A", nsig: float = 1.6, color=st.TRUTH,
               appendage_color=None, opacity=0.55, scale=A_TO_UNITS,
               resolution=(14, 10)) -> VGroup:
    """VGroup of ellipsoids; the last one is the appendage when state is given."""
    g = VGroup()
    for b in sp.blob_params(state):
        col = appendage_color if (b.kind == "appendage" and appendage_color) else color
        g.add(ellipsoid(b.center, b.axes, nsig * b.sigmas, col, opacity, resolution, scale))
    return g


def to_scene(R: np.ndarray) -> np.ndarray:
    """Matrix to apply to Manim points so the 3D drawing matches projection R.

    Manim mobjects live in the object frame; applying R gives the microscope
    frame, whose z axis is the beam. Use mob.apply_matrix(R) about ORIGIN.
    """
    return np.asarray(R)
