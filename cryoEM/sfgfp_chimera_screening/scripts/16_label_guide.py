"""Visual guide for hand-labelling: what each species SHOULD look like in these exact images.

The L6_loop50 atomic model (and its sfGFP / DNA parts, plus a straight 50 bp duplex) is projected at 2.84 A/px,
multiplied by this dataset's CTF (200 kV, Cs 2.7 mm, amp 0.1, defocus of a real pilot micrograph), scaled so the
simulated sfGFP matches the measured contrast of the real sfGFP dots in THAT micrograph, and added onto the real
micrograph (its own noise and ice). The composite then goes through exactly the labeller's two views: 15 A low-pass
and Topaz pretrained denoising. Real crops of automatically detected "bare", "one-arm" and "two-arm closed" dots
(results/G0_dots_pilot.csv.gz) are shown alongside.

This is a viewing aid only: nothing simulated here is used for picking, classification or refinement.
  $DNASPRING_PY -I scripts/16_label_guide.py      -> figures/qc/label_guide.png
"""
import os
import subprocess
import tempfile

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import mrcfile  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

REPO = os.environ.get("SFGFP_REPO", os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
DATA = os.environ.get("SFGFP_DATA", "/n/netscratch/hekstra_lab/Lab/ian_poon/sfgfp_chimera_screening_20261112")
PROJ = "/n/netscratch/hekstra_lab/Lab/ian_poon/cryosparc_projects/CS-cs-sfgfp-chimera-screening/"
TOPAZ = "/n/holylabs/hekstra_lab/Everyone/ianpoon/software/bin/topaz-cs"
APIX, N, KV, CS, AMP, LP = 2.84, 1024, 200.0, 2.7e7, 0.1, 15.0
Z = {"C": 6, "N": 7, "O": 8, "P": 15, "S": 16}
rng = np.random.default_rng(7)


def atoms(path):
    xyz, w = [], []
    for line in open(path):
        if line.startswith(("ATOM", "HETATM")):
            e = (line[76:78].strip() or line[12:16].strip()[0]).upper()[:1]
            if e == "H":
                continue
            xyz.append([float(line[30:38]), float(line[38:46]), float(line[46:54])])
            w.append(Z.get(e, 6))
    return np.array(xyz), np.array(w, float)


def rot(axis, deg):
    a = np.radians(deg)
    c, s = np.cos(a), np.sin(a)
    return {"x": np.array([[1, 0, 0], [0, c, -s], [0, s, c]]),
            "y": np.array([[c, 0, s], [0, 1, 0], [-s, 0, c]]),
            "z": np.array([[c, -s, 0], [s, c, 0], [0, 0, 1]])}[axis]


def ctf_kernel(df_A):
    lam = 12.2643 / np.sqrt(KV * 1e3 * (1 + 0.978466e-6 * KV * 1e3))
    ky = np.fft.fftfreq(N, APIX)[:, None]
    kx = np.fft.rfftfreq(N, APIX)[None, :]
    k2 = ky ** 2 + kx ** 2
    chi = np.pi * lam * df_A * k2 - 0.5 * np.pi * CS * lam ** 3 * k2 ** 2
    return -(np.sqrt(1 - AMP ** 2) * np.sin(chi) + AMP * np.cos(chi)), k2


def lowpass(img, k2):
    return np.fft.irfft2(np.fft.rfft2(img) * np.exp(-k2 * LP ** 2), s=img.shape)


def project(xyz, w, R, centre):
    """Deposit atom weights (projected along z after rotation R) into an NxN image at `centre` (px)."""
    x = (xyz - xyz.mean(0)) @ R.T
    img = np.zeros((N, N))
    ij = np.round(x[:, :2] / APIX).astype(int) + np.asarray(centre)[::-1]
    ok = (ij[:, 0] >= 0) & (ij[:, 0] < N) & (ij[:, 1] >= 0) & (ij[:, 1] < N)
    np.add.at(img, (ij[ok, 1], ij[ok, 0]), w[ok])
    return img


def plane_frame(xyz):
    """Rotation taking the DNA loop's best-fit plane normal onto z (face-on view)."""
    c = xyz - xyz.mean(0)
    _, _, vt = np.linalg.svd(c, full_matrices=False)
    return vt  # rows: in-plane major, in-plane minor, normal


def fcrop(a, n):
    F = np.fft.fftshift(np.fft.fft2(a))
    c = a.shape[0] // 2
    F = F[c - n // 2: c + n // 2, c - n // 2: c + n // 2]
    return np.real(np.fft.ifft2(np.fft.ifftshift(F))).astype(np.float32) * (n / a.shape[0]) ** 2


def depth(img_lp, xs, ys):
    """Centre value minus local 40-60 A annulus median, per point (px coordinates in img)."""
    out = []
    yy, xx = np.mgrid[-25:26, -25:26]
    r = np.hypot(yy, xx) * APIX
    ann = (r > 40) & (r < 60)
    for x, y in zip(xs, ys):
        x, y = int(round(x)), int(round(y))
        if 26 <= x < N - 26 and 26 <= y < N - 26:
            p = img_lp[y - 25:y + 26, x - 25:x + 26]
            out.append(p[25, 25] - np.median(p[ann]))
    return np.array(out)


def main():
    dots = pd.read_csv(os.path.join(REPO, "results", "G0_dots_pilot.csv.gz"))
    # a pilot micrograph near 1.5 um defocus with many dots
    per = dots.groupby("micrograph_path").agg(df=("df_um", "first"), n=("x_px", "size")).reset_index()
    per = per[(per.n > 250)]
    mic = per.iloc[(per.df - 1.5).abs().argsort()].iloc[0]
    raw = fcrop(mrcfile.read(PROJ + mic.micrograph_path).astype(np.float32), N)
    df_A = mic.df * 1e4
    ctf, k2 = ctf_kernel(df_A)
    md = dots[dots.micrograph_path == mic.micrograph_path]
    real_depth = np.median(depth(lowpass(raw, k2), md.x_px / 2, md.y_px / 2))

    gfp = atoms(os.path.join(REPO, "model/L6_loop50/L6_protein.pdb"))
    loop = atoms(os.path.join(REPO, "model/L6_loop50/L6_loop50_solute.pdb"))
    dna_loop = atoms(os.path.join(REPO, "model/L6_loop50/L6_dna.pdb"))
    dna_lin = atoms(os.path.join(REPO, "model/L6_loop50/dna/straight_50bp.pdb"))
    # calibrate: simulated sfGFP (random orientations) depth vs the real dots' median depth
    sims = []
    for _ in range(24):
        R = rot("z", rng.uniform(0, 360)) @ rot("y", rng.uniform(0, 180)) @ rot("x", rng.uniform(0, 360))
        s = lowpass(np.fft.irfft2(np.fft.rfft2(project(*gfp, R, (512, 512))) * ctf, s=(N, N)), k2)
        sims.append(depth(s, [512], [512])[0])
    scale = real_depth / np.median(sims)

    # lollipop: sfGFP with a straight duplex leaving one attachment site radially
    frame = plane_frame(dna_loop[0])
    gc = gfp[0].mean(0)
    d_lin = dna_lin[0] - dna_lin[0].mean(0)
    _, _, vt = np.linalg.svd(d_lin, full_matrices=False)
    axis = vt[0]
    out_dir = (gc - dna_loop[0].mean(0))
    out_dir /= np.linalg.norm(out_dir)
    # rotate duplex axis onto out_dir and place its end 20 A beyond the protein surface
    v = np.cross(axis, out_dir)
    c = float(axis @ out_dir)
    vx = np.array([[0, -v[2], v[1]], [v[2], 0, -v[0]], [-v[1], v[0], 0]])
    Rl = np.eye(3) + vx + vx @ vx / (1 + c)
    lin_xyz = d_lin @ Rl.T + gc + out_dir * (85 + 25)
    lolly = (np.vstack([gfp[0], lin_xyz]), np.concatenate([gfp[1], dna_lin[1]]))

    species = [
        ("Closed loop\nface-on", loop, frame),
        ("Closed loop\ntilted 50°", loop, rot("x", 50) @ frame),
        ("Closed loop\nedge-on", loop, rot("x", 90) @ frame),
        ("Dot + long tail\n(open, 'lollipop')", lolly, rot("x", 20) @ frame),
        ("Free 50 bp DNA\n(~170 Å)", dna_lin, rot("y", 90)),
        ("Free sfGFP\n(side view)", gfp, rot("x", 90)),
    ]
    # insertion sites: a grid, 300 A apart, on the real micrograph
    sites = [(200 + 210 * (i % 3), 260 + 260 * (i // 3)) for i in range(len(species))]
    model_only = np.zeros((N, N))
    clean_tiles = []
    for (name, (xyz, w), R), (cx, cy) in zip(species, sites):
        p = np.fft.irfft2(np.fft.rfft2(project(xyz, w, R, (cy, cx))) * ctf, s=(N, N)) * scale
        model_only += p
        iso = lowpass(np.fft.irfft2(np.fft.rfft2(project(xyz, w, R, (512, 512))) * ctf, s=(N, N)) * scale, k2)
        clean_tiles.append(iso[512 - 64:512 + 64, 512 - 64:512 + 64])
    comp = raw + model_only
    lp_comp = lowpass(comp, k2)
    # Topaz denoise of the composite (same pretrained model as the labeller's 'Denoised' view)
    tmp = tempfile.mkdtemp(dir=os.environ.get("TMPDIR", "/tmp"))
    with mrcfile.new(os.path.join(tmp, "comp.mrc"), overwrite=True) as m:
        m.set_data(comp.astype(np.float32))
    subprocess.run([TOPAZ, "denoise", "--device", "-1", "-o", os.path.join(tmp, "dn"), os.path.join(tmp, "comp.mrc")],
                   check=True, capture_output=True)
    dn_comp = mrcfile.read(os.path.join(tmp, "dn", "comp.mrc")).astype(np.float32)

    def tile(img, cx, cy, h=64):
        return img[cy - h:cy + h, cx - h:cx + h]

    def show(ax, t, ref=None):
        lo, hi = np.percentile(ref if ref is not None else t, [0.5, 99.5])
        ax.imshow(t, cmap="gray", vmin=lo, vmax=hi)
        ax.set_xticks([]); ax.set_yticks([])

    real_rows = []
    for cls, n in (("two_arm_closed", 8), ("one_arm", 4), ("bare", 4)):
        sub = dots[(dots.arm_class == cls) & (dots.df_um.between(1.2, 1.9))].sample(frac=1, random_state=3)
        real_rows.append((cls, sub.head(n)))

    fig = plt.figure(figsize=(16, 15))
    gs = fig.add_gridspec(5, 8, hspace=0.35, wspace=0.08)
    rows = [("Model only\n(no noise)", None), ("In real ice,\n15 Å low-pass view", lp_comp), ("In real ice,\nDenoised view", dn_comp)]
    for r, (rname, img) in enumerate(rows):
        for j, ((name, _, _), (cx, cy)) in enumerate(zip(species, sites)):
            ax = fig.add_subplot(gs[r, j])
            if img is None:
                show(ax, clean_tiles[j], ref=lp_comp)
                ax.set_title(name, fontsize=9)
            else:
                show(ax, tile(img, cx, cy), ref=img)
            if j == 0:
                ax.set_ylabel(rname, fontsize=9)
            if r == 0 and j == 0:
                ax.plot([8, 8 + 100 / APIX], [118, 118], color="#ffd400", lw=3)
                ax.text(8, 110, "100 Å", color="#ffd400", fontsize=8)
    # real crops: auto-detected dot classes (pilot micrographs, 15 A low-pass, labeller scale)
    cache = {}
    k = 0
    for cls, sub in real_rows:
        for _, d in sub.iterrows():
            if d.micrograph_path not in cache:
                im = fcrop(mrcfile.read(PROJ + d.micrograph_path).astype(np.float32), N)
                c2, kk = ctf_kernel(d.df_um * 1e4)
                cache[d.micrograph_path] = lowpass(im, kk)
            img = cache[d.micrograph_path]
            cx, cy = int(d.x_px / 2), int(d.y_px / 2)
            cx, cy = min(max(cx, 64), N - 64), min(max(cy, 64), N - 64)
            ax = fig.add_subplot(gs[3 + k // 8, k % 8])
            show(ax, tile(img, cx, cy), ref=img)
            ax.set_title({"two_arm_closed": "auto: 2 arms, closed?", "one_arm": "auto: 1 arm", "bare": "auto: bare dot"}[cls], fontsize=8)
            k += 1
    fig.suptitle(
        f"What to look for at the labeller's scale (2.84 Å/px; each tile 128 px = 364 Å). Rows 1-3: the L6_loop50 model "
        f"projected with this dataset's CTF ({mic.df:.2f} µm defocus, 200 kV) and contrast-matched to the real sfGFP dots "
        f"of that micrograph, then placed on its real ice.\nRows 4-5: real pilot crops picked by the automatic dot/thread "
        f"detector (15 Å low-pass). Simulations are a viewing aid only.", fontsize=10)
    out = os.path.join(REPO, "figures", "qc", "label_guide.png")
    plt.savefig(out, dpi=90, bbox_inches="tight")
    print("wrote", out, "| scale", round(float(scale), 4), "| real dot depth", round(float(real_depth), 3),
          "| micrograph", os.path.basename(mic.micrograph_path), "df", mic.df)


if __name__ == "__main__":
    main()
