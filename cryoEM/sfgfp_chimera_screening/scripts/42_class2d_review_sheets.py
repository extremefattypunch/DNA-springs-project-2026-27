"""Large-tile review sheets of 2D class averages for blind classification (40 classes per page, 100 A bar per tile).

  $DNASPRING_PY -I scripts/42_class2d_review_sheets.py KEY [--min-particles 270]
Writes figures/2d/review/KEY_pNN.png and prints the page -> class-id mapping.
"""
import argparse
import glob
import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import mrcfile  # noqa: E402
import numpy as np  # noqa: E402

REPO = os.environ.get("SFGFP_REPO", os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
DATA = os.environ.get("SFGFP_DATA", "/n/netscratch/hekstra_lab/Lab/ian_poon/sfgfp_chimera_screening_20261112")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("key")
    ap.add_argument("--min-particles", type=int, default=270)
    a = ap.parse_args()
    exp = os.path.join(DATA, "exports")
    with mrcfile.open(sorted(glob.glob(os.path.join(exp, a.key, "*class_averages*.mrc")))[-1], permissive=True) as m:
        stack = np.asarray(m.data, np.float32)
        apix = float(m.voxel_size.x)
    cls = np.load(os.path.join(exp, f"{a.key}__particles.npz"))["alignments2D__class"].astype(int)
    counts = np.bincount(cls, minlength=stack.shape[0])
    order = [k for k in np.argsort(-counts) if counts[k] >= a.min_particles]
    out = os.path.join(REPO, "figures", "2d", "review")
    os.makedirs(out, exist_ok=True)
    n = stack.shape[1]
    for p0 in range(0, len(order), 40):
        page = order[p0:p0 + 40]
        fig, axes = plt.subplots(5, 8, figsize=(16, 10.6))
        for ax in axes.ravel():
            ax.axis("off")
        for ax, k in zip(axes.ravel(), page):
            img = stack[k]
            lo, hi = np.percentile(img, [1, 99.7])
            ax.imshow(img, cmap="gray", vmin=lo, vmax=hi, interpolation="lanczos")
            ax.plot([3, 3 + 100 / apix], [n - 4] * 2, color="#ffd400", lw=2.5)
            ax.set_title(f"class {k}  ({counts[k]:,} ptcls)", fontsize=9)
        fig.suptitle(f"{a.key} review sheet {p0 // 40 + 1}: classes sorted by particle count; box {n * apix:.0f} Å; "
                     f"yellow bar = 100 Å; density is WHITE", fontsize=11)
        plt.tight_layout(rect=(0, 0, 1, 0.97))
        f = os.path.join(out, f"{a.key}_p{p0 // 40 + 1:02d}.png")
        plt.savefig(f, dpi=110)
        plt.close(fig)
        print(f, "classes:", ",".join(map(str, page)))
    print(f"{len(order)} classes with >= {a.min_particles} particles "
          f"({counts[order].sum():,} of {counts.sum():,} particles, {counts[order].sum() / counts.sum():.1%})")


if __name__ == "__main__":
    main()
