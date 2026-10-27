"""Shared look: palette, typography, badges, cards and array -> image conversion.

Semantic colours are reserved across all six clips:
    cyan   = observation / measured image      green  = ground truth / signal
    orange = estimate (pose, reconstruction)    purple = alternate conformer
    red    = residual / noise                   white/grey = text
"""
from __future__ import annotations

import numpy as np
from manim import (DOWN, LEFT, RIGHT, UL, UP, UR, ImageMobject, MathTex, Rectangle,
                   RoundedRectangle, Text, VGroup, config)
from PIL import Image

BG = "#0E1116"
PANEL = "#1A1F27"
OBS = "#4FD1E8"        # cyan
TRUTH = "#5BD67A"      # green
EST = "#F5A142"        # orange
ALT = "#B48CF2"        # purple
RESID = "#F2645A"      # red
TEXT = "#ECEDEF"
MUTED = "#9AA3AE"
BADGE = "#FFD166"

FONT = "Nimbus Sans"
TITLE_SIZE = 44
BODY_SIZE = 32          # >= 32 px-equivalent at 1080p
SMALL_SIZE = 26         # captions/badges only
MATH_SIZE = 44
MARGIN = 0.45           # outer frame margin (scene units)

TEX_TO_COLOR = {}


def T(text: str, size: int = BODY_SIZE, color: str = TEXT, **kw) -> Text:
    return Text(text, font=FONT, font_size=size, color=color, **kw)


def M(*tex: str, size: int = MATH_SIZE, color: str = TEXT, **kw) -> MathTex:
    return MathTex(*tex, font_size=size, color=color, **kw)


def title_card(title: str, question: str | None = None) -> VGroup:
    """Clip title (top-left) and optional guiding question beneath it."""
    t = T(title, TITLE_SIZE, weight="BOLD")
    g = VGroup(t)
    if question:
        q = T(question, BODY_SIZE - 2, MUTED)
        q.next_to(t, DOWN, aligned_edge=LEFT, buff=0.15)
        g.add(q)
    g.to_corner(UL, buff=MARGIN)
    return g


def simulation_badge() -> VGroup:
    """Small 'SIMULATION' tag for the upper-right corner (keep it on whenever
    synthetic images or maps are visible). In ThreeDScene add it with
    add_fixed_in_frame_mobjects."""
    t = Text("SIMULATION", font=FONT, font_size=SMALL_SIZE, color=BADGE, weight="BOLD")
    box = RoundedRectangle(corner_radius=0.08, width=t.width + 0.35, height=t.height + 0.22,
                           stroke_color=BADGE, stroke_width=2, fill_color=BG, fill_opacity=0.85)
    t.move_to(box)
    g = VGroup(box, t)
    g.to_corner(UR, buff=0.3)
    return g


def phantom_disclaimer() -> Text:
    return T("synthetic PP2A-inspired shape; not a measured structure", SMALL_SIZE, MUTED,
             slant="ITALIC")


def caption(text: str, color: str = MUTED, size: int = BODY_SIZE - 4) -> Text:
    """Bottom-centre caption for the current shot."""
    c = T(text, size, color)
    c.to_edge(DOWN, buff=MARGIN)
    return c


def takeaway_card(math: str | None, text: str, xray: str | None = None,
                  width: float = 12.5) -> VGroup:
    """Summary freeze: one mathematical takeaway and an optional X-ray analogy."""
    rows = []
    head = T("Takeaway", BODY_SIZE - 2, BADGE, weight="BOLD")
    rows.append(head)
    if math:
        rows.append(M(math, size=MATH_SIZE))
    rows.append(T(text, BODY_SIZE))
    if xray:
        rows.append(T("X-ray analogy: " + xray, BODY_SIZE - 4, MUTED, slant="ITALIC"))
    body = VGroup(*rows).arrange(DOWN, buff=0.28)
    if body.width > width - 0.6:
        body.scale_to_fit_width(width - 0.6)
    box = RoundedRectangle(corner_radius=0.15, width=max(body.width + 0.8, 6),
                           height=body.height + 0.6, stroke_color=MUTED, stroke_width=1.5,
                           fill_color=PANEL, fill_opacity=0.95)
    body.move_to(box)
    return VGroup(box, body)


def legend(items: list[tuple[str, str]], size: int = SMALL_SIZE) -> VGroup:
    """Colour legend: [(color, label), ...] stacked vertically."""
    rows = VGroup()
    for col, lab in items:
        sw = Rectangle(width=0.32, height=0.2, fill_color=col, fill_opacity=1, stroke_width=0)
        rows.add(VGroup(sw, T(lab, size)).arrange(RIGHT, buff=0.15))
    return rows.arrange(DOWN, aligned_edge=LEFT, buff=0.12)


# ---------------------------------------------------------------- images
def _hex_to_rgb(h: str) -> np.ndarray:
    h = h.lstrip("#")
    return np.array([int(h[i:i + 2], 16) for i in (0, 2, 4)], float)


GRAY = [(0.0, "#000000"), (1.0, "#FFFFFF")]
DIVERGING = [(0.0, "#3B6FD8"), (0.5, "#F2F2F2"), (1.0, "#D8433B")]   # CTF sign
MAGMA = [(0.0, "#000004"), (0.25, "#3B0F70"), (0.5, "#8C2981"), (0.75, "#DE4968"),
         (0.9, "#FE9F6D"), (1.0, "#FCFDBF")]
LOCRES = [(0.0, "#2F7DE1"), (0.5, "#F2F2F2"), (1.0, "#E8513B")]      # good -> poor


def colorize(a01: np.ndarray, stops=GRAY) -> np.ndarray:
    """Map an array in [0, 1] through colour stops -> uint8 RGB."""
    pos = np.array([p for p, _ in stops])
    cols = np.stack([_hex_to_rgb(c) for _, c in stops])
    a = np.clip(a01, 0, 1)
    out = np.stack([np.interp(a, pos, cols[:, ch]) for ch in range(3)], -1)
    return out.astype(np.uint8)


def to_display(arr_xy: np.ndarray) -> np.ndarray:
    """(x, y)-indexed array -> image rows (row 0 = top = largest y)."""
    return np.flipud(np.swapaxes(arr_xy, 0, 1))


def array_to_image_mobject(arr_xy: np.ndarray, height: float = 3.0, vmin=None, vmax=None,
                           stops=GRAY, nearest: bool = False, alpha=None,
                           symmetric: bool = False) -> ImageMobject:
    """Scalar (x, y) array -> ImageMobject. Precompute arrays; never call per frame.

    symmetric=True maps [-m, m] -> [0, 1] (for signed data such as the CTF).
    alpha (same shape, 0..1) makes the background transparent.
    """
    a = np.asarray(arr_xy, float)
    if symmetric:
        m = np.max(np.abs(a)) if vmax is None else vmax
        a01 = 0.5 + 0.5 * a / (m if m > 0 else 1)
    else:
        lo = np.min(a) if vmin is None else vmin
        hi = np.max(a) if vmax is None else vmax
        a01 = (a - lo) / (hi - lo if hi > lo else 1)
    rgb = colorize(to_display(a01), stops)
    if alpha is not None:
        al = (np.clip(to_display(np.asarray(alpha, float)), 0, 1) * 255).astype(np.uint8)
        rgb = np.concatenate([rgb, al[..., None]], -1)
    im = ImageMobject(rgb)
    if nearest:
        from manim import RESAMPLING_ALGORITHMS
        im.set_resampling_algorithm(RESAMPLING_ALGORITHMS["nearest"])
    im.height = height
    return im


def density_image(render: dict, color: str = TRUTH, height: float = 3.0,
                  label_colors: dict | None = None) -> ImageMobject:
    """Shaded isodensity render (from projections.render_density) -> RGBA image.

    label_colors maps blob index -> hex colour (e.g. local resolution); others
    use `color`.
    """
    shade, mask, label = render["shade"], render["mask"], render["label"]
    base = np.broadcast_to(_hex_to_rgb(color), shade.shape + (3,)).copy()
    if label_colors:
        for idx, col in label_colors.items():
            base[label == idx] = _hex_to_rgb(col)
    rgb = np.clip(base * shade[..., None], 0, 255)
    a = mask.astype(float)
    rgba = np.concatenate([rgb, 255 * a[..., None]], -1).astype(np.uint8)
    rgba = np.stack([to_display(rgba[..., c]) for c in range(4)], -1)
    im = ImageMobject(rgba)
    im.height = height
    return im


def framed(img: ImageMobject, color: str = MUTED, buff: float = 0.0) -> VGroup:
    """Thin frame around an image panel."""
    r = Rectangle(width=img.width + 2 * buff, height=img.height + 2 * buff,
                  stroke_color=color, stroke_width=2).move_to(img)
    return VGroup(r)
