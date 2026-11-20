"""2D class-average montage (Fig S5B style), sorted by particle count, with per-class counts.

  $DNASPRING_PY -I scripts/40_fig_2d_montage.py KEY [--labels results/class2d_labels_KEY.csv] [--ncol 10] [--top N]

Reads exports written by 30_cs_export.py:  exports/KEY/*class_averages*.mrc  and  exports/KEY__particles.npz
(field alignments2D/class). Writes figures/2d/KEY_montage.png and results/class2d_counts_KEY.csv.
If a labels CSV (class,label) exists, each tile is outlined by population colour.
"""
import argparse
import csv
import glob
import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import mrcfile  # noqa: E402
import numpy as np  # noqa: E402

REPO = os.environ.get("SFGFP_REPO", os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
DATA = os.environ.get("SFGFP_DATA", "/n/netscratch/hekstra_lab/Lab/ian_poon/sfgfp_chimera_screening_20261112")
COLOURS = {"chimera": "#2a9d8f", "lollipop": "#8ab17d", "rod": "#e76f51", "blob": "#f4a261", "junk": "#9e9e9e"}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("key")
    ap.add_argument("--labels")
    ap.add_argument("--ncol", type=int, default=10)
    ap.add_argument("--top", type=int, default=0)
    a = ap.parse_args()
    exp = os.path.join(DATA, "exports")
    mrcs = sorted(glob.glob(os.path.join(exp, a.key, "*class_averages*.mrc")))
    if not mrcs:
        raise SystemExit(f"no class-average stack in {exp}/{a.key}")
    with mrcfile.open(mrcs[-1], permissive=True) as m:  # last iteration's stack
        stack = np.asarray(m.data, dtype=np.float32)
        apix = float(m.voxel_size.x)
    p = np.load(os.path.join(exp, f"{a.key}__particles.npz"))
    cls = p["alignments2D__class"].astype(int)
    K = stack.shape[0]
    counts = np.bincount(cls, minlength=K)
    order = np.argsort(-counts)
    if a.top:
        order = order[: a.top]
    labels = {}
    if a.labels and os.path.exists(a.labels):
        labels = {int(r["class"]): r["label"] for r in csv.DictReader(open(a.labels))}
    os.makedirs(os.path.join(REPO, "results"), exist_ok=True)
    with open(os.path.join(REPO, "results", f"class2d_counts_{a.key}.csv"), "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["class", "particles", "fraction", "label"])
        for k in np.argsort(-counts):
            w.writerow([k, counts[k], round(counts[k] / counts.sum(), 5), labels.get(k, "")])
    n = len(order)
    nrow = int(np.ceil(n / a.ncol))
    fig, axes = plt.subplots(nrow, a.ncol, figsize=(a.ncol * 1.3, nrow * 1.45))
    for ax in np.ravel(axes):
        ax.axis("off")
    for i, k in enumerate(order):
        ax = np.ravel(axes)[i]
        img = stack[k]
        lo, hi = np.percentile(img, [1, 99.5])
        ax.imshow(img, cmap="gray", vmin=lo, vmax=hi)
        ax.set_title(f"{k}: {counts[k]:,}", fontsize=6.5, pad=1.5)
        if k in labels:
            for s in ax.spines.values():
                s.set_visible(True)
                s.set_color(COLOURS.get(labels[k], "k"))
                s.set_linewidth(2.2)
            ax.axis("on")
            ax.set_xticks([])
            ax.set_yticks([])
    # 100 A scale bar on the first tile
    ax0 = np.ravel(axes)[0]
    L = 100.0 / apix
    ax0.plot([4, 4 + L], [stack.shape[1] - 6] * 2, color="w", lw=2)
    ax0.text(4, stack.shape[1] - 10, "100 Å", color="w", fontsize=6)
    fig.suptitle(f"{a.key}: {K} classes, {counts.sum():,} particles, box {stack.shape[1]} px @ {apix:.2f} Å/px "
                 f"(sorted by size)", fontsize=9)
    plt.tight_layout(rect=(0, 0, 1, 0.97))
    out = os.path.join(REPO, "figures", "2d", f"{a.key}_montage.png")
    os.makedirs(os.path.dirname(out), exist_ok=True)
    plt.savefig(out, dpi=150)
    print("wrote", out)


if __name__ == "__main__":
    main()
