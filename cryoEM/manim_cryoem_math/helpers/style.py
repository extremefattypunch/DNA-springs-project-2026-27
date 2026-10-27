"""Shared look: palette, typography, label, cards and array -> image conversion.

Semantic colours are reserved across all six clips:
    cyan   = observed data                  green  = underlying signal / ground truth
    orange = estimates and alignment arrows purple = alternate mathematical state
    red    = noise / residuals              grey   = axes, grids, secondary text
"""
from __future__ import annotations

import numpy as np
from manim import (DOWN, LEFT, RIGHT, UL, UR, ImageMobject, MathTex, Rectangle,
                   RoundedRectangle, Text, VGroup)

BG = "#0E1116"
PANEL = "#1A1F27"
GRID = "#3A424D"
OBS = "#4FD1E8"        # cyan
TRUTH = "#5BD67A"      # green
EST = "#F5A142"        # orange
ALT = "#B48CF2"        # purple
RESID = "#F2645A"      # red
TEXT = "#ECEDEF"
MUTED = "#9AA3AE"
LABEL = "#FFD166"      # the MATHEMATICAL ILLUSTRATION tag and "Takeaway" heading

FONT = "Nimbus Sans"
TITLE_SIZE = 44
BODY_SIZE = 32
SMALL_SIZE = 26         # tag, footnotes and tick labels only
MATH_SIZE = 44
MARGIN = 0.45           # outer frame margin (scene units)
FRAME_W, FRAME_H = 14.222, 8.0
TITLE_MAX_W = 7.9       # title starts at x = -6.66; the tag's left edge is at x = 1.71


def T(text: str, size: int = BODY_SIZE, color: str = TEXT, **kw) -> Text:
    return Text(text, font=FONT, font_size=size, color=color, **kw)


def M(*tex: str, size: int = MATH_SIZE, color: str = TEXT, **kw) -> MathTex:
    return MathTex(*tex, font_size=size, color=color, **kw)


def fit_width(mob, max_w: float):
    """Scale down (never up) so mob.width <= max_w."""
    if mob.width > max_w:
        mob.scale_to_fit_width(max_w)
    return mob


def title_card(title: str, question: str | None = None) -> VGroup:
    """Clip title (top-left) and optional question beneath it. Long text is
    broken by the caller with '\n'; the title is kept clear of the tag."""
    t = fit_width(T(title, TITLE_SIZE, weight="BOLD"), TITLE_MAX_W)
    g = VGroup(t)
    if question:
        q = fit_width(T(question, BODY_SIZE - 2, MUTED), FRAME_W - 2 * MARGIN)
        q.next_to(t, DOWN, aligned_edge=LEFT, buff=0.15)
        g.add(q)
    g.to_corner(UL, buff=MARGIN)
    return g


def illustration_label() -> VGroup:
    """Unobtrusive 'MATHEMATICAL ILLUSTRATION' tag, upper-right, on for the whole
    clip. In ThreeDScene add it with add_fixed_in_frame_mobjects."""
    t = Text("MATHEMATICAL ILLUSTRATION", font=FONT, font_size=SMALL_SIZE - 4, color=LABEL,
             weight="BOLD")
    box = RoundedRectangle(corner_radius=0.08, width=t.width + 0.3, height=t.height + 0.2,
                           stroke_color=LABEL, stroke_width=1.5, stroke_opacity=0.8,
                           fill_color=BG, fill_opacity=0.85)
    t.move_to(box)
    g = VGroup(box, t)
    g.to_corner(UR, buff=MARGIN)
    return g


def simplified_note(text: str = "deliberately simplified equation") -> Text:
    """Small italic note to mark a simplified / illustrative equation."""
    return T(text, SMALL_SIZE, MUTED, slant="ITALIC")


def caption(text: str, color: str = MUTED, size: int = BODY_SIZE - 4) -> Text:
    """Bottom-centre caption for the current shot."""
    c = fit_width(T(text, size, color), FRAME_W - 2 * MARGIN)
    c.to_edge(DOWN, buff=MARGIN)
    return c


def takeaway_card(lines: list[str], math: str | None = None, width: float = 12.4,
                  line_spacing: float = 0.16) -> VGroup:
    """Summary freeze card: heading, optional MathTex row, then text lines.

    Pass text as separate short lines (<= ~60 characters at BODY_SIZE). Raises
    if the body would need shrinking below the readable minimum.
    """
    head = T("Takeaway", BODY_SIZE - 2, LABEL, weight="BOLD")
    rows = [head]
    if math:
        rows.append(M(math, size=MATH_SIZE - 4))
    txt = VGroup(*[T(s, BODY_SIZE - 2) for s in lines]).arrange(DOWN, buff=line_spacing)
    rows.append(txt)
    body = VGroup(*rows).arrange(DOWN, buff=0.25)
    if body.width > width - 0.6:
        scale = (width - 0.6) / body.width
        if scale < 0.9:
            raise ValueError("takeaway_card: a line is too wide; break it into shorter lines")
        body.scale(scale)
    box = RoundedRectangle(corner_radius=0.15, width=max(body.width + 0.7, 6),
                           height=body.height + 0.5, stroke_color=MUTED, stroke_width=1.5,
                           fill_color=PANEL, fill_opacity=0.96)
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
DIVERGING = [(0.0, "#3B6FD8"), (0.5, "#F2F2F2"), (1.0, "#D8433B")]   # signed quantities


def colorize(a01: np.ndarray, stops=GRAY) -> np.ndarray:
    """Map an array in [0, 1] through colour stops -> uint8 RGB."""
    pos = np.array([p for p, _ in stops])
    cols = np.stack([_hex_to_rgb(c) for _, c in stops])
    a = np.clip(a01, 0, 1)
    return np.stack([np.interp(a, pos, cols[:, ch]) for ch in range(3)], -1).astype(np.uint8)


def to_display(arr_xy: np.ndarray) -> np.ndarray:
    """[x, y]-indexed array -> image rows (row 0 = top = largest y)."""
    return np.flipud(np.swapaxes(arr_xy, 0, 1))


def heatmap(arr_xy: np.ndarray, height: float = 2.0, vmin=None, vmax=None, stops=GRAY,
            symmetric: bool = False, smooth: bool = False) -> ImageMobject:
    """Tiny grayscale (or signed) matrix -> ImageMobject with visible pixels.

    Precompute arrays; never call per frame. symmetric=True maps [-m, m] -> [0, 1].
    Pixels are nearest-neighbour by default: these are matrices, not photographs.
    """
    a = np.asarray(arr_xy, float)
    if symmetric:
        m = np.max(np.abs(a)) if vmax is None else vmax
        a01 = 0.5 + 0.5 * a / (m if m > 0 else 1)
    else:
        lo = np.min(a) if vmin is None else vmin
        hi = np.max(a) if vmax is None else vmax
        a01 = (a - lo) / (hi - lo if hi > lo else 1)
    im = ImageMobject(colorize(to_display(a01), stops))
    if not smooth:
        from manim import RESAMPLING_ALGORITHMS
        im.set_resampling_algorithm(RESAMPLING_ALGORITHMS["nearest"])
    im.height = height
    return im


def framed(img: ImageMobject, color: str = GRID, buff: float = 0.0) -> Rectangle:
    """Thin frame around a heatmap (returns the Rectangle; group with Group())."""
    return Rectangle(width=img.width + 2 * buff, height=img.height + 2 * buff,
                     stroke_color=color, stroke_width=2).move_to(img)
