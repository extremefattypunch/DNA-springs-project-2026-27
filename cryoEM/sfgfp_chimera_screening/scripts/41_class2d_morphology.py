"""Objective shape descriptors for every 2D class average (feeds the blind class review and the census).

  $DNASPRING_PY -I scripts/41_class2d_morphology.py KEY      (exports from 30_cs_export.py must exist)

For each class average (density = bright, as CryoSPARC writes them):
  mask        pixels > mean + 2 sd of a soft-masked background annulus, largest connected component near centre
  length_A    longest extent of the mask (max Feret diameter, A)
  width_A     mask area / length (mean width, A); DNA ~20 A, sfGFP ~25-40 A
  elong       length / width
  peak_z      brightest compact blob (12 A Gaussian) in sd units: sfGFP stands out strongly
  n_blobs     compact blobs above peak_z/2 (1 = lone dot; 2 = pair/dumbbell)
  hole        1 if the mask encloses a background hole of >= 30 A diameter (closed ring candidate)
  ring_A      if hole: outer diameter of the enclosing mask
Writes results/class2d_morph_KEY.csv. Heuristics only: they describe the averages, they do not label them.
"""
import csv
import glob
import os
import sys

import mrcfile
import numpy as np
from scipy import ndimage as ndi

REPO = os.environ.get("SFGFP_REPO", os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
DATA = os.environ.get("SFGFP_DATA", "/n/netscratch/hekstra_lab/Lab/ian_poon/sfgfp_chimera_screening_20261112")


def feret(mask, apix):
    ys, xs = np.nonzero(mask)
    if len(xs) < 2:
        return 0.0
    pts = np.stack([xs, ys], 1).astype(float)
    if len(pts) > 600:
        pts = pts[np.random.default_rng(0).choice(len(pts), 600, replace=False)]
    d = np.sqrt(((pts[:, None] - pts[None]) ** 2).sum(-1))
    return float(d.max() + 1) * apix


def describe(img, apix):
    n = img.shape[0]
    yy, xx = np.mgrid[:n, :n] - (n - 1) / 2
    r = np.hypot(yy, xx)
    # CryoSPARC flattens everything outside its circular window, so estimate the window radius from ring
    # variance and take robust (median / MAD) background statistics from inside it
    rings = np.array([img[(r >= a) & (r < a + 1)].std() for a in range(n // 2)])
    inside = np.nonzero(rings > 0.1 * rings[: max(2, n // 8)].max())[0]
    rw = (inside.max() + 1) if len(inside) else n // 2
    win = r < rw
    med = np.median(img[win])
    mad = 1.4826 * np.median(np.abs(img[win] - med)) + 1e-9
    z = np.where(win, (img - med) / mad, 0.0)
    lab, nl = ndi.label(z > 2.0)
    if nl == 0:
        return dict(length_A=0, width_A=0, elong=0, peak_z=float(z.max()), n_blobs=0, hole=0, ring_A=0, area_A2=0)
    # component closest to the centre (weighted by size)
    best, score = 0, -1
    for i in range(1, nl + 1):
        m = lab == i
        d = r[m].min()
        s = m.sum() / (1 + d)
        if s > score:
            best, score = i, s
    m = ndi.binary_closing(lab == best, iterations=1)
    area = m.sum() * apix ** 2
    L = feret(m, apix)
    W = area / L if L else 0
    sm = ndi.gaussian_filter(z, 12 / apix / 2.355)
    pk = ndi.maximum_filter(sm, size=max(3, int(20 / apix))) == sm
    vals = sm[pk & m]
    peak = float(vals.max()) if len(vals) else float(sm.max())
    nb = int((vals > peak / 2).sum()) if len(vals) else 0
    filled = ndi.binary_fill_holes(m)
    holes = filled & ~m
    hl, nh = ndi.label(holes)
    hole, ringd = 0, 0.0
    for i in range(1, nh + 1):
        h = hl == i
        if 2 * np.sqrt(h.sum() / np.pi) * apix >= 30:
            hole, ringd = 1, feret(filled, apix)
    return dict(length_A=round(L, 1), width_A=round(W, 1), elong=round(L / W, 2) if W else 0, peak_z=round(peak, 2),
                n_blobs=nb, hole=hole, ring_A=round(ringd, 1), area_A2=round(area, 0))


def main(key):
    exp = os.path.join(DATA, "exports")
    mrcs = sorted(glob.glob(os.path.join(exp, key, "*class_averages*.mrc")))
    with mrcfile.open(mrcs[-1], permissive=True) as mm:
        stack = np.asarray(mm.data, dtype=np.float32)
        apix = float(mm.voxel_size.x)
    p = np.load(os.path.join(exp, f"{key}__particles.npz"))
    counts = np.bincount(p["alignments2D__class"].astype(int), minlength=stack.shape[0])
    rows = []
    for k in range(stack.shape[0]):
        d = describe(stack[k], apix)
        rows.append({"class": k, "particles": int(counts[k]), "fraction": round(counts[k] / counts.sum(), 5), **d})
    rows.sort(key=lambda r: -r["particles"])
    out = os.path.join(REPO, "results", f"class2d_morph_{key}.csv")
    with open(out, "w", newline="") as fh:
        w = csv.DictWriter(fh, list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    print("wrote", out, f"({len(rows)} classes, {apix:.2f} A/px)")


if __name__ == "__main__":
    main(sys.argv[1])
