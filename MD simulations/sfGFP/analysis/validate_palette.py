#!/usr/bin/env python3
"""Validate this project's chart palette against the computable colour checks.

Why this exists as a file rather than a one-off
-----------------------------------------------
The rule for a categorical palette is that you compute whether it is
colour-vision-deficiency safe rather than judging it by eye, and the reference
implementation of those checks is JavaScript.  There is no node on this cluster, so this
is the Python twin: same thresholds, same Machado-Oliveira-Fernandes simulation, same
OKLab distance, so the numbers are comparable with anyone else's run rather than being
this project's private opinion.

Checks, and what each is for:

  2.  Lightness band     OKLCH L inside the mode's band, so no series is so dark it
                         reads as ink or so light it disappears into the surface.
  3.  Chroma floor       OKLCH C >= 0.10; below that a hue reads as grey and stops
                         carrying identity.
  4.  CVD separation     OKLab dE x100 between adjacent slots under simulated protanopia
                         and deuteranopia, taking the worse of the two.  8 is the target,
                         6 a floor that is only legal alongside a second encoding.
  4b. Normal-vision floor  the worst adjacent pair under ordinary vision must clear 15.
                         A palette can be CVD-safe and still put two blues next to each
                         other that nobody can tell apart; this catches that.
  5.  Contrast           WCAG ratio of each mark against the chart surface, 3:1 minimum.

Checks 1 and 6 -- that hues are assigned in a fixed order and drawn from the documented
palette -- are structural and cannot be measured from hex values, so they are not here.

    python analysis/validate_palette.py                  # this project's palette
    python analysis/validate_palette.py --mode dark
    python analysis/validate_palette.py --palette "#2a78d6,#eb6834" --pairs all
"""
from __future__ import annotations

import argparse
import math
import re
import sys
from pathlib import Path

BAND = {"light": (0.43, 0.77), "dark": (0.48, 0.67)}     # OKLCH L
CHROMA_FLOOR = 0.10
CVD_TARGET, CVD_FLOOR = 8.0, 6.0                         # OKLab dE x100, adjacent pairs
NORMAL_FLOOR = 15.0
CONTRAST_MIN = 3.0
DEFAULT_SURFACE = {"light": "#fcfcfb", "dark": "#1a1a19"}

# Machado, Oliveira & Fernandes (2009), severity 1.0, in LINEAR RGB.  The thresholds
# above are calibrated to this particular simulation, so substituting another model
# (Vienot 1999, say) would move borderline pairs and require recalibrating them.
MACHADO = {
    "protan": ((0.152286, 1.052583, -0.204868),
               (0.114503, 0.786281, 0.099216),
               (-0.003882, -0.048116, 1.051998)),
    "deutan": ((0.367322, 0.860646, -0.227968),
               (0.280085, 0.672501, 0.047413),
               (-0.011820, 0.042940, 0.968881)),
    "tritan": ((1.255528, -0.076749, -0.178779),
               (-0.078411, 0.930809, 0.147602),
               (0.004733, 0.691367, 0.303900)),
}

# Kept deliberately in step with the JavaScript original: the two runtimes strip
# different edge whitespace, so the shared set is their intersection.
_WS = "[ \t\n\v\f\r   -     　]+"


def _strip(v: str) -> str:
    return re.sub(f"^{_WS}|{_WS}$", "", v)


def _is_hex(v: str) -> bool:
    return bool(re.fullmatch(r"#?[0-9a-fA-F]{6}", v))


def _srgb(h: str):
    h = _strip(h).lstrip("#")
    return [int(h[i:i + 2], 16) / 255 for i in (0, 2, 4)]


def _lin(h: str):
    return [c / 12.92 if c <= 0.04045 else ((c + 0.055) / 1.055) ** 2.4
            for c in _srgb(h)]


def _rel_lum(h: str) -> float:
    r, g, b = _lin(h)
    return 0.2126 * r + 0.7152 * g + 0.0722 * b


def contrast(a: str, b: str) -> float:
    hi, lo = sorted((_rel_lum(a), _rel_lum(b)), reverse=True)
    return (hi + 0.05) / (lo + 0.05)


def _oklab_from_lin(rgb):
    r, g, b = rgb
    l = (0.4122214708 * r + 0.5363325363 * g + 0.0514459929 * b) ** (1 / 3)
    m = (0.2119034982 * r + 0.6806995451 * g + 0.1073969566 * b) ** (1 / 3)
    s = (0.0883024619 * r + 0.2817188376 * g + 0.6299787005 * b) ** (1 / 3)
    return (0.2104542553 * l + 0.7936177850 * m - 0.0040720468 * s,
            1.9779984951 * l - 2.4285922050 * m + 0.4505937099 * s,
            0.0259040371 * l + 0.7827717662 * m - 0.8086757660 * s)


def oklch(h: str):
    L, a, b = _oklab_from_lin(_lin(h))
    return L, math.hypot(a, b)


def _simulate(h: str, kind: str):
    r, g, b = _lin(h)
    M = MACHADO[kind]
    return [min(1.0, max(0.0, M[i][0] * r + M[i][1] * g + M[i][2] * b))
            for i in range(3)]


def delta_e(h1: str, h2: str, kind: str | None = None) -> float:
    a = _oklab_from_lin(_simulate(h1, kind) if kind else _lin(h1))
    b = _oklab_from_lin(_simulate(h2, kind) if kind else _lin(h2))
    return 100 * math.dist(a, b)


def validate(palette, mode="light", surface=None, pairs="adjacent"):
    surface = surface or DEFAULT_SURFACE[mode]
    for c in list(palette) + [surface]:
        if not _is_hex(c):
            sys.exit(f"not a hex colour: {c!r} -- refusing to run, because an unguarded "
                     "parse propagates through every check and the run fails OPEN")
    lo, hi = BAND[mode]
    rows, fails, warns = [], [], []

    for i, c in enumerate(palette):
        L, C = oklch(c)
        cr = contrast(c, surface)
        band_ok = lo <= L <= hi
        chroma_ok = C >= CHROMA_FLOOR
        rows.append({"slot": i, "hex": c, "L": round(L, 3), "C": round(C, 3),
                     "contrast": round(cr, 2),
                     "band": "pass" if band_ok else "FAIL",
                     "chroma": "pass" if chroma_ok else "FAIL",
                     "vs surface": "pass" if cr >= CONTRAST_MIN else "WARN"})
        if not band_ok:
            fails.append(f"slot {i} {c}: OKLCH L {L:.3f} outside the {mode} band "
                         f"[{lo}, {hi}]")
        if not chroma_ok:
            fails.append(f"slot {i} {c}: OKLCH C {C:.3f} below the {CHROMA_FLOOR} "
                         "chroma floor -- it will read as grey")
        if cr < CONTRAST_MIN:
            warns.append(f"slot {i} {c}: contrast {cr:.2f}:1 against the surface is "
                         f"below {CONTRAST_MIN}:1 -- obligates visible labels or a "
                         "table view, and is not dismissable")

    idx = ([(i, i + 1) for i in range(len(palette) - 1)] if pairs == "adjacent"
           else [(i, j) for i in range(len(palette)) for j in range(i + 1, len(palette))])
    pair_rows = []
    for i, j in idx:
        a, b = palette[i], palette[j]
        p, d, tr = (delta_e(a, b, "protan"), delta_e(a, b, "deutan"),
                    delta_e(a, b, "tritan"))
        worst_cvd, normal = min(p, d), delta_e(a, b)
        if worst_cvd >= CVD_TARGET:
            verdict = "pass"
        elif worst_cvd >= CVD_FLOOR:
            verdict = "WARN"
            warns.append(f"{a}/{b}: CVD dE {worst_cvd:.1f} is in the {CVD_FLOOR}-"
                         f"{CVD_TARGET} floor band -- legal only with a second encoding")
        else:
            verdict = "FAIL"
            fails.append(f"{a}/{b}: CVD dE {worst_cvd:.1f} below the {CVD_FLOOR} floor")
        if normal < NORMAL_FLOOR:
            verdict = "FAIL"
            fails.append(f"{a}/{b}: normal-vision dE {normal:.1f} below the "
                         f"{NORMAL_FLOOR} floor -- full-colour readers cannot tell them "
                         "apart, and a second encoding does not excuse this one")
        pair_rows.append({"pair": f"{i}-{j}", "protan": round(p, 1),
                          "deutan": round(d, 1), "tritan": round(tr, 1),
                          "normal": round(normal, 1), "verdict": verdict})
    return rows, pair_rows, fails, warns


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--palette", default=None,
                    help="comma-separated hex; default is this project's SERIES")
    ap.add_argument("--mode", choices=("light", "dark"), default="light")
    ap.add_argument("--surface", default=None)
    ap.add_argument("--pairs", choices=("adjacent", "all"), default="adjacent")
    a = ap.parse_args()

    if a.palette:
        pal = [_strip(x) for x in a.palette.split(",") if _strip(x)]
    else:
        sys.path.insert(0, str(Path(__file__).resolve().parent))
        import vizstyle as VS
        pal = list(VS.SERIES)
        a.surface = a.surface or (VS.SURFACE if a.mode == "light" else None)

    rows, pair_rows, fails, warns = validate(pal, a.mode, a.surface, a.pairs)
    surface = a.surface or DEFAULT_SURFACE[a.mode]
    print(f"palette of {len(pal)} against the {a.mode} surface {surface}, "
          f"{a.pairs} pairs\n")
    print(f"{'slot':>4} {'hex':>9} {'L':>6} {'C':>6} {'contr':>6} "
          f"{'band':>6} {'chroma':>7} {'surface':>8}")
    for r in rows:
        print(f"{r['slot']:>4} {r['hex']:>9} {r['L']:>6.3f} {r['C']:>6.3f} "
              f"{r['contrast']:>6.2f} {r['band']:>6} {r['chroma']:>7} "
              f"{r['vs surface']:>8}")
    print(f"\n{'pair':>6} {'protan':>7} {'deutan':>7} {'tritan':>7} {'normal':>7} "
          f"{'verdict':>8}")
    for r in pair_rows:
        print(f"{r['pair']:>6} {r['protan']:>7.1f} {r['deutan']:>7.1f} "
              f"{r['tritan']:>7.1f} {r['normal']:>7.1f} {r['verdict']:>8}")
    print()
    for w in warns:
        print(f"  WARN  {w}")
    for x in fails:
        print(f"  FAIL  {x}")
    if not fails and not warns:
        print("  all checks pass")
    sys.exit(1 if fails else 0)


if __name__ == "__main__":
    main()
