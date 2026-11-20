#!/usr/bin/env python
"""Fig S7 analogue: gallery of the 10 PCA class_3D volumes (P6 J59) under the dot-on-rod NU map (J56).

Stages (run as   $PY -I scripts/50_figS7_class3d.py prep|montage):
  prep     Load J56_003_volume_map.mrc (160^3, 2.84 A/px) and J59_class_0N_00014_volume.mrc (80^3, 5.68 A/px, the same
           454 A box), Fourier-resample everything to 160^3 @ 2.84 A, Gaussian low-pass to 12 A (half amplitude),
           normalise to the background shell (150-178 A from the box centre), and choose ONE contour criterion for all
           maps: the threshold at which all density within 205 A of the box centre encloses 65,000 A^3 (the volume of
           the project's primary contour, scripts/46_dotrod_3d_length.py: sfGFP 33,000 + a 50 bp duplex 32,000 A^3;
           total rather than largest-component volume, because the latter jumps when an arm detaches).
           Class_3D keeps the input (NU) poses, so the class volumes are already in the NU frame; this is checked by
           the real-space correlation of each class with the NU map (both low-passed) at zero shift vs the best
           integer shift (recorded).  One rigid rotation, from the principal axes of the NU map's thresholded main
           component (x = longest, z = thinnest), is applied identically to every map, so every panel has the same
           orientation and scale.  Display maps (level 1.0 = the 65,000 A^3 contour) and the .cxc go to scratch/scripts.
  montage  Crop the ChimeraX renders with one common box per view and assemble figures/final/figS7_class3d.png, plus
           figS7_class3d.json (every plotted number and its source file).
"""
import glob
import json
import os
import sys

import numpy as np
import mrcfile
from scipy import ndimage as ndi

REPO = "/n/home05/ianpoon/DNA-springs-project-2026-27/cryoEM/sfgfp_chimera_screening"
CS = "/n/netscratch/hekstra_lab/Lab/ian_poon/cryosparc_projects/CS-cs-sfgfp-chimera-screening"
SCR = "/n/netscratch/hekstra_lab/Lab/ian_poon/sfgfp_chimera_screening_20261112/figures_scratch/figS7"
CXC = os.path.join(REPO, "scripts", "50_figS7_class3d.cxc")
FIG = os.path.join(REPO, "figures", "final", "figS7_class3d.png")
META = os.path.join(SCR, "prep.json")

APIX, N, LP_RES, PRIMARY = 2.84, 160, 12.0, 65000.0
NU_MAP = f"{CS}/J56/J56_003_volume_map.mrc"
NU_PART = f"{CS}/J56/J56_003_particles.cs"
CLS_MAP = f"{CS}/J59/J59_class_{{:02d}}_00014_volume.mrc"
CLS_PART = f"{CS}/J59/J59_class_{{:02d}}_00014_particles.cs"
REJ_PART = f"{CS}/J59/J59_particles_rejected.cs"
COLOR_NU, COLOR_CLS = "#4a3aa7", "#2a78d6"     # vizstyle SERIES[6], SERIES[0]


# ------------------------------------------------------------------------------------------------ map utilities
def fourier_resample(d, n_out):
    n = d.shape[0]
    F = np.fft.fftshift(np.fft.fftn(d))
    out = np.zeros((n_out,) * 3, complex)
    s = n_out // 2 - n // 2
    out[s:s + n, s:s + n, s:s + n] = F
    return np.real(np.fft.ifftn(np.fft.ifftshift(out))) * (n_out / n) ** 3


def load(path):
    with mrcfile.open(path, permissive=True) as m:
        d = np.asarray(m.data, dtype=np.float64)
        apix = float(m.voxel_size.x)
    if abs(d.shape[0] * apix - N * APIX) > 1.0:
        raise SystemExit(f"box mismatch {path}")
    if d.shape[0] != N:
        d = fourier_resample(d, N)
    return d, apix, d.shape[0]


f = np.fft.fftfreq(N, APIX)
KX, KY, KZ = np.meshgrid(f, f, f, indexing="ij")
K = np.sqrt(KX ** 2 + KY ** 2 + KZ ** 2)
c = (np.arange(N) - N // 2) * APIX
ZZ, YY, XX = np.meshgrid(c, c, c, indexing="ij")
RAD = np.sqrt(XX ** 2 + YY ** 2 + ZZ ** 2)
S26 = np.ones((3, 3, 3), bool)


def lowpass(d):
    return np.real(np.fft.ifftn(np.fft.fftn(d) * np.exp(-np.log(2) * (K * LP_RES) ** 2)))


def normalise(d):
    sh = (RAD > 150) & (RAD < 178)
    return (d - np.median(d[sh])) / d[sh].std()


def largest(m):
    m = m & (RAD < 205)
    lab, n = ndi.label(m, S26)
    if n == 0:
        return m, []
    sz = np.bincount(lab.ravel())[1:]
    big = int(np.argmax(sz)) + 1
    others = sorted([float(s) * APIX ** 3 for i, s in enumerate(sz) if i + 1 != big and s * APIX ** 3 > 2000], reverse=True)
    return lab == big, others


def threshold_for_volume(lp, V):
    """Bisection so that ALL density above threshold within 205 A of the box centre encloses V A^3 (monotone in the
    threshold, unlike the largest-component volume, which jumps when an arm detaches)."""
    inside = RAD < 205
    hi, lo = float(lp.max()), float(np.percentile(lp, 50))
    for _ in range(60):
        t = 0.5 * (hi + lo)
        v = ((lp > t) & inside).sum() * APIX ** 3
        if v > V:
            lo = t
        else:
            hi = t
        if abs(v - V) / V < 0.002:
            break
    m, others = largest(lp > t)
    return t, m, others


def cc(a, b, mask=None):
    if mask is not None:
        a, b = a[mask], b[mask]
    a = a - a.mean()
    b = b - b.mean()
    return float((a * b).sum() / np.sqrt((a * a).sum() * (b * b).sum()))


def best_shift_cc(ref, mov):
    """Cross-correlation over all integer shifts (FFT); returns (cc at 0, best cc, best shift in A)."""
    r = ref - ref.mean()
    m = mov - mov.mean()
    x = np.real(np.fft.ifftn(np.fft.fftn(r) * np.conj(np.fft.fftn(m)))) / np.sqrt((r * r).sum() * (m * m).sum())
    i = np.unravel_index(int(np.argmax(x)), x.shape)
    sh = [((k + N // 2) % N - N // 2) * APIX for k in i]
    return float(x[0, 0, 0]), float(x[i]), sh


# ------------------------------------------------------------------------------------------------ prep
def prep():
    os.makedirs(SCR, exist_ok=True)
    n_nu = len(np.load(NU_PART))
    counts = [len(np.load(CLS_PART.format(i))) for i in range(10)]
    n_rej = len(np.load(REJ_PART))
    maps = [dict(key="nu", label="NU refinement J56 (dot-on-rod)", path=NU_MAP, n=n_nu, n_source=NU_PART)]
    maps += [dict(key=f"c{i}", label=f"class_3D J59 class {i}", path=CLS_MAP.format(i), n=counts[i],
                  n_source=CLS_PART.format(i)) for i in range(10)]
    lps = {}
    for M in maps:
        d, apix, n0 = load(M["path"])
        lp = normalise(lowpass(d))
        t, m, others = threshold_for_volume(lp, PRIMARY)
        M.update(input_apix=apix, input_box=n0, threshold_sigma=round(float(t), 3),
                 main_component_A3=round(float(m.sum() * APIX ** 3)), other_components_A3=[round(o) for o in others],
                 total_above_threshold_A3=round(float((lp > t).sum() * APIX ** 3)))
        lps[M["key"]] = (lp, t, m)
        print(M["key"], M["n"], M["threshold_sigma"], M["main_component_A3"], M["other_components_A3"], flush=True)

    # frame check: class volumes vs NU (both 12 A low-passed)
    nu_lp, nu_t, nu_m = lps["nu"]
    mask = ndi.binary_dilation(nu_m, S26, iterations=4)
    for M in maps[1:]:
        lp = lps[M["key"]][0]
        c0, cbest, sh = best_shift_cc(nu_lp, lp)
        M.update(cc_vs_nu_masked=round(cc(nu_lp, lp, mask), 3), cc_vs_nu_zero_shift=round(c0, 3),
                 cc_vs_nu_best_shift=round(cbest, 3), best_shift_A=[round(s, 2) for s in sh])

    # one rotation for all maps: principal axes of the NU main component
    idx = np.argwhere(nu_m).astype(float)                   # (z, y, x) voxel indices
    xyz = idx[:, ::-1]
    cen = xyz.mean(0)
    w, V = np.linalg.eigh(np.cov((xyz - cen).T))
    order = np.argsort(w)[::-1]                             # longest first
    R = V[:, order].T                                       # rows = new x, y, z axes in old (x, y, z) coords
    if np.linalg.det(R) < 0:
        R[2] *= -1                                          # keep a proper rotation: never mirror the map
    # sign convention: new +x points from the sfGFP-lobe-free end toward ... fixed by skewness so it is reproducible
    proj = (xyz - cen) @ R.T
    for k in range(2):
        if (proj[:, k] ** 3).mean() < 0:
            R[k] *= -1
            R[2] *= -1                                      # flip two axes together: still a proper rotation
            proj = (xyz - cen) @ R.T
    ctr = np.array([N / 2] * 3)                             # CryoSPARC box centre (voxel N/2)

    def rotate(d):
        # output voxel o (x,y,z) samples input at cen + R^T (o - ctr)
        Rzyx = R.T[::-1, ::-1]
        off = cen[::-1] - Rzyx @ ctr[::-1]
        return ndi.affine_transform(d, Rzyx, offset=off, order=3, mode="constant", cval=0.0)

    for M in maps:
        lp, t, _ = lps[M["key"]]
        out = os.path.join(SCR, f"disp_{M['key']}.mrc")
        r = rotate(lp / t).astype(np.float32)
        with mrcfile.new(out, overwrite=True) as mm:
            mm.set_data(r)
            mm.voxel_size = APIX
        M["display_map"] = out
        mk = largest(r > 1.0)[0]
        M["display_main_component_A3_after_rotation"] = round(float(mk.sum() * APIX ** 3))
        xs = np.where(mk.any(axis=(0, 1)))[0]
        M["display_x_extent_A"] = round(float((xs.max() - xs.min() + 1) * APIX), 1)
    meta = dict(script="scripts/50_figS7_class3d.py", lowpass_A=LP_RES, apix=APIX, box=N, contour_A3=PRIMARY,
                contour_rule="threshold at which all voxels above threshold within 205 A of the box centre of the 12 A "
                             "Gaussian-low-passed, background-normalised map enclose 65,000 A^3",
                rotation_rows_newxyz_in_old=R.round(5).tolist(), rotation_centre_xyz_vox=cen.round(3).tolist(),
                pca_eigenvalues_A2=(w[order] * APIX ** 2).round(1).tolist(),
                n_nu_particles=n_nu, n_class3d_assigned=int(sum(counts)), n_class3d_rejected=n_rej,
                rejected_source=REJ_PART, maps=maps)
    json.dump(meta, open(META, "w"), indent=1)
    write_cxc(maps)


def write_cxc(maps):
    L = ["# Fig S7 analogue renders (generated by scripts/50_figS7_class3d.py prep).",
         "# All display maps share one frame and one rotation; level 1.0 = density within 205 A encloses 65,000 A^3.",
         "log clear", "set bgColor white", "lighting soft", "graphics silhouettes true width 1.5",
         ]
    for i, M in enumerate(maps):
        col = COLOR_NU if M["key"] == "nu" else COLOR_CLS
        L.append(f"open {M['display_map']} id {i + 1}")
        L.append(f"volume #{i + 1} level 1.0 style surface color {col} step 1")
    # one camera for every panel: fit the union of all maps once, then never re-fit
    L += ["view orient", "view pad 0.04", f"open {os.path.join(REPO, 'scripts', '50_figS7_ortho.py')}",
          "view pad 0.04", "zoom 0.8"]
    for view, cmd in (("face", None), ("edge", "turn x -90")):
        if cmd:
            L.append(cmd)
        for i, M in enumerate(maps):
            L.append("hide #1-11 models")
            L.append(f"show #{i + 1} models")
            L.append(f"save {SCR}/r_{M['key']}_{view}.png width 1400 height 1400 supersample 3")
    L.append("exit")
    open(CXC, "w").write("\n".join(L) + "\n")
    print("wrote", CXC)


# ------------------------------------------------------------------------------------------------ montage
def montage():
    import textwrap
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import matplotlib.image as mpimg
    sys.path.insert(0, "/n/home05/ianpoon/DNA-springs-project-2026-27/cryoEM/pp2a_day2025_reproduction/scripts")
    import vizstyle
    vizstyle.apply(matplotlib)
    meta = json.load(open(META))
    maps = meta["maps"]
    tot = meta["n_class3d_assigned"]
    # one horizontal crop for both views (the edge view is a 90 deg turn about x, so x stays horizontal at the
    # same scale); one vertical crop per view; identical for every map
    imgs, ys, xs = {}, {"face": [], "edge": []}, []
    for view in ("face", "edge"):
        for M in maps:
            im = mpimg.imread(f"{SCR}/r_{M['key']}_{view}.png")[..., :3]
            imgs[(M["key"], view)] = im
            nz = np.argwhere(im.min(-1) < 0.97)
            ys[view] += [nz[:, 0].min(), nz[:, 0].max()]
            xs += [nz[:, 1].min(), nz[:, 1].max()]
    pad = 14
    X0, X1 = max(min(xs) - pad, 0), max(xs) + pad
    Y = {v: (max(min(ys[v]) - pad, 0), max(ys[v]) + pad) for v in ys}
    crop = lambda k, v: imgs[(k, v)][Y[v][0]:Y[v][1], X0:X1]   # noqa: E731
    wpx = X1 - X0
    af, ae = (Y["face"][1] - Y["face"][0]) / wpx, (Y["edge"][1] - Y["edge"][0]) / wpx

    # manual layout in inches
    W, gap, col_gap, title_h, edge_gap = 2.5, 0.18, 0.08, 0.32, 0.04
    Wn = 2.0 * W
    top_h = title_h + Wn * (af + ae) + edge_gap
    tile_h = title_h + W * (af + ae) + edge_gap
    cap_h, sup_h = 1.25, 0.42
    figw = 5 * W + 4 * col_gap + 0.3
    figh = sup_h + top_h + gap + 2 * tile_h + gap + cap_h
    fig = plt.figure(figsize=(figw, figh), facecolor="white")
    fig.set_layout_engine(None)

    def place(x_in, y_top_in, w_in, h_in):
        return fig.add_axes([x_in / figw, 1 - (y_top_in + h_in) / figh, w_in / figw, h_in / figh])

    def tile(key, x, y, w, title, fs):
        fig.text((x + w / 2) / figw, 1 - (y + title_h * 0.55) / figh, title, ha="center", va="center", fontsize=fs)
        a1 = place(x, y + title_h, w, w * af)
        a1.imshow(crop(key, "face"))
        a1.axis("off")
        a2 = place(x, y + title_h + w * af + edge_gap, w, w * ae)
        a2.imshow(crop(key, "edge"))
        a2.axis("off")
        return a1

    nu = maps[0]
    y = sup_h
    a_nu = tile("nu", (figw - Wn) / 2, y, Wn, f"NU refinement J56 (consensus map), {nu['n']:,} particles", 10.5)
    # 100 A scale bar: render px per A calibrated from the NU face view (x extent of its contoured main component in
    # A, from the display map, vs its non-white pixel width in the ortho render; +-1 voxel = +-2 %)
    nz = np.argwhere(imgs[("nu", "face")].min(-1) < 0.97)
    px_per_A = (nz[:, 1].max() - nz[:, 1].min()) / nu["display_x_extent_A"]
    im0 = crop("nu", "face")
    x0, y0 = 0.02 * im0.shape[1], 0.95 * im0.shape[0]
    a_nu.plot([x0, x0 + 100 * px_per_A], [y0, y0], color=vizstyle.INK, lw=2.5, solid_capstyle="butt")
    a_nu.text(x0 + 50 * px_per_A, y0 - 0.025 * im0.shape[0], "100 Å", ha="center", va="bottom", fontsize=9)
    meta["scalebar_px_per_A"] = round(float(px_per_A), 4)
    y += top_h + gap
    for r in range(2):
        for jc in range(5):
            i = 5 * r + jc
            M = maps[i + 1]
            a_t = tile(M["key"], 0.15 + jc * (W + col_gap), y, W, f"class {i}: {M['n']:,} ({100 * M['n'] / tot:.1f} %)", 9.5)
            if i == 0:   # same crop and render px per A as the NU tile, so the same 100 A in image pixels
                a_t.plot([x0, x0 + 100 * px_per_A], [y0, y0], color=vizstyle.INK, lw=2.0, solid_capstyle="butt")
                a_t.text(x0 + 50 * px_per_A, y0 - 0.03 * im0.shape[0], "100 Å", ha="center", va="bottom", fontsize=8)
        y += tile_h + (gap if r == 0 else 0)
    fig.text(0.15 / figw, 1 - 0.2 / figh, "Fig S7 analogue: 3D classification of the dot-on-rod NU particles "
             "(P6 J59: PCA initialisation, 10 classes, 3 components, 12 Å)", ha="left", va="center", fontsize=11.5)
    zero = all(np.allclose(M["best_shift_A"], 0, atol=APIX + 1e-6) for M in maps[1:])
    frame = ("their correlation with the NU map peaks at zero shift for every class" if zero else
             "see figS7_class3d.json for the class-vs-NU shift check")
    # splitting sentence generated from the per-class components (prep.json = figS7_class3d.json "prep" maps):
    # a separate piece >= 5,000 A^3 is a "large piece", a smaller one a "tip" (components < 2,000 A^3 are not recorded)
    BIG = 5000.0
    big = [(M["key"][1:], [o for o in M["other_components_A3"] if o >= BIG]) for M in maps[1:]]
    big = [(k, o) for k, o in big if o]
    tips = [(M["key"][1:], [o for o in M["other_components_A3"] if o < BIG]) for M in maps[1:]]
    tips = [(k, o) for k, o in tips if o]
    multi = [(M["key"][1:], M["other_components_A3"]) for M in maps[1:] if len(M["other_components_A3"]) > 1]
    big = [(k, o) for k, o in big if k not in [m[0] for m in multi]]   # classes with several pieces are listed once, below
    and_ = lambda xs: ", ".join(xs[:-1]) + " and " + xs[-1] if len(xs) > 1 else xs[0]  # noqa: E731
    split_txt = (f"At this common contour the density splits into separate pieces: classes {and_([k for k, _ in big])} "
                 f"each have a large separate piece ("
                 + "; ".join(f"class {k}: {' + '.join(f'~{round(v, -2):,.0f}' for v in o)} Å³" for k, o in big) + ")")
    if tips:
        split_txt += "; " + "; ".join(f"class {k} has a separate ~{round(o[0], -2):,.0f} Å³ tip" for k, o in tips if k not in
                                      [m[0] for m in multi])
    for k, o in multi:
        nw = {2: "two", 3: "three", 4: "four"}.get(len(o), str(len(o)))
        split_txt += (f"; class {k} has {nw} extra pieces ({' and '.join(f'~{round(v, -2):,.0f}' for v in o)} Å³)")
    split_txt += " (pieces < 2,000 Å³ not counted; all component volumes in figS7_class3d.json)."
    meta["split_sentence"] = split_txt
    cap = (f"Each tile: face view (looking down the thinnest principal axis of the NU map) above the edge view "
           f"(rotated 90° about the horizontal axis). One rotation, one orthographic camera and one crop for every map; "
           f"class_3D keeps the NU poses, so all maps are in the NU frame ({frame}). All maps Gaussian low-passed to "
           f"12 Å and contoured by one rule: density within 205 Å of the box centre encloses 65,000 Å³ (sfGFP + one "
           f"50 bp duplex). {split_txt} Class tiles are drawn at 0.5× the NU panel scale (100 Å bars on the NU panel "
           f"and on class 0). Counts: particles assigned to each class in J59 "
           f"({tot:,} assigned, {meta['n_class3d_rejected']} rejected); % of assigned.")
    fig.text(0.15 / figw, (cap_h - 0.08) / figh, "\n".join(textwrap.wrap(cap, 190)), ha="left", va="top",
             fontsize=8.2, color=vizstyle.INK_2)
    fig.savefig(FIG, dpi=200, facecolor="white")
    plt.close(fig)
    src = dict(figure=os.path.relpath(FIG, REPO), script="scripts/50_figS7_class3d.py",
               renders="scripts/50_figS7_class3d.cxc (ChimeraX, offscreen; ortho camera via scripts/50_figS7_ortho.py)",
               prep=meta,
               plotted_numbers=[dict(panel=M["key"], particles=M["n"], source=M["n_source"],
                                     percent_of_assigned=None if M["key"] == "nu" else round(100 * M["n"] / tot, 2),
                                     map=M["path"], threshold_sigma=M["threshold_sigma"],
                                     main_component_A3=M["main_component_A3"],
                                     other_components_A3=M["other_components_A3"]) for M in maps],
               caption=cap, split_sentence=split_txt,
               scalebar=dict(px_per_A=meta["scalebar_px_per_A"], source="NU display map x extent (prep.json "
                             "display_x_extent_A) vs NU face-render pixel width",
                             drawn_on=["NU face view (100 A)", "class 0 face view (100 A, same render px per A; tiles "
                                       "are 0.5x the NU panel's printed size)"]))
    json.dump(src, open(FIG.replace(".png", ".json"), "w"), indent=1)
    print("wrote", FIG)


if __name__ == "__main__":
    {"prep": prep, "montage": montage}[sys.argv[1]]()
