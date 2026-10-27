#!/usr/bin/env python3
"""Intensity profiles for locating lanes and bands in a gel image.

  profile.py IMAGE cols --y0 Y0 --y1 Y1
      Column-mean profile between rows y0..y1 (find lane centres).
  profile.py IMAGE lane --x X [--w W] [--y0 Y0 --y1 Y1] [--tilt DEG --cy CY]
      Row profile of one lane (mean over a strip of width W centred at X,
      following the gel tilt if given), with detected peaks: y, height above
      the local background, FWHM and a relative intensity (0-1 vs the
      strongest peak found in that lane).
  profile.py IMAGE lanes --xs X1,X2,... [--w W] [...]
      Same for several lanes; intensities are relative to the strongest peak
      across ALL listed lanes, so brightness is comparable between lanes.

Add --invert for gels with dark bands on a light background.
Values are on a 0-255 grey scale of the original JPEG.
"""
import argparse, math
import numpy as np
from PIL import Image

ap = argparse.ArgumentParser()
ap.add_argument("image")
ap.add_argument("mode", choices=["cols", "lane", "lanes"])
ap.add_argument("--x", type=float)
ap.add_argument("--xs")
ap.add_argument("--w", type=int, default=12)
ap.add_argument("--y0", type=int, default=0)
ap.add_argument("--y1", type=int)
ap.add_argument("--tilt", type=float, default=0.0, help="gel tilt in degrees, clockwise positive")
ap.add_argument("--cy", type=float, help="y at which lane x is measured (default image mid-height)")
ap.add_argument("--invert", action="store_true")
ap.add_argument("--min-prom", type=float, default=4.0, help="minimum peak prominence (grey levels)")
a = ap.parse_args()

im = np.asarray(Image.open(a.image).convert("L"), dtype=float)
if a.invert:
    im = 255 - im
H, W = im.shape
y1 = a.y1 if a.y1 is not None else H
cy = a.cy if a.cy is not None else H / 2
t = math.tan(math.radians(a.tilt))


def smooth(v, k=3):
    if k <= 1:
        return v
    ker = np.ones(k) / k
    return np.convolve(v, ker, mode="same")


if a.mode == "cols":
    prof = smooth(im[a.y0:y1, :].mean(axis=0), 3)
    print(f"column means rows {a.y0}..{y1}")
    for x in range(W):
        bar = "#" * int(prof[x] / 4)
        print(f"{x:4d} {prof[x]:6.1f} {bar}")
    raise SystemExit


def lane_profile(xc):
    ys = np.arange(a.y0, y1)
    vals = []
    for y in ys:
        x = xc - (y - cy) * t
        lo, hi = int(round(x - a.w / 2)), int(round(x + a.w / 2))
        lo, hi = max(lo, 0), min(hi, W)
        vals.append(im[y, lo:hi].mean() if hi > lo else 0)
    return ys, smooth(np.array(vals), 3)


def peaks(ys, v):
    # local background: running minimum-ish via wide percentile filter
    n = len(v)
    bg = np.array([np.percentile(v[max(0, i - 25):min(n, i + 26)], 15) for i in range(n)])
    s = v - bg
    out = []
    for i in range(2, n - 2):
        if s[i] >= s[i - 1] and s[i] > s[i + 1] and s[i] >= s[i - 2] and s[i] >= s[i + 2] and s[i] > a.min_prom:
            half = s[i] / 2
            l = i
            while l > 0 and s[l] > half:
                l -= 1
            r = i
            while r < n - 1 and s[r] > half:
                r += 1
            out.append((int(ys[i]), float(s[i]), int(r - l), float(v[i]), float(bg[i])))
    return out


xs = [a.x] if a.mode == "lane" else [float(x) for x in a.xs.split(",")]
allp = []
for xc in xs:
    ys, v = lane_profile(xc)
    allp.append((xc, peaks(ys, v)))
top = max([p[1] for _, ps in allp for p in ps] or [1])
for xc, ps in allp:
    print(f"lane x={xc:.1f} (strip {a.w}px, rows {a.y0}..{y1})")
    print("     y  height  fwhm   raw   bg   rel")
    for (y, h, fw, raw, bg) in ps:
        print(f"  {y:4d}  {h:6.1f}  {fw:4d}  {raw:5.1f} {bg:5.1f}  {h / top:4.2f}")
