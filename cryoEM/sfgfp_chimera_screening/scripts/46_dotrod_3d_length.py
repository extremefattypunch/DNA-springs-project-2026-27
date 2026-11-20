"""W4 3D: length and geometry of the DNA attached to sfGFP, measured in 3D (removes 2D projection foreshortening).

Question (report §13a, results/A1_length_vs_hinge_verified.json): in 2D the DNA on sfGFP projects ~75 A past the
sfGFP surface vs ~109 A for free rods. Is the attached duplex genuinely SHORTER (H1, yielded/incomplete loops) or a
FULL-LENGTH duplex tilted out of plane (H2/H3)? A 3D map has no projection foreshortening.

Maps (all Fourier-resampled to 2.84 A/px, 160^3 = 454 A box, then Gaussian low-passed to 12 A for segmentation):
  nu      W4 NU refine J56 (final iteration _003; 54,181 particles = hetero class 0)
  het0/1  W4 hetero refine J55 class 0 (dot-on-rod, 54,181) and class 1 (free-rod-seeded, 42,863)
  c3d_0-9 W4 class_3D PCA J59 classes 0-9 (on the NU particles)
  free    W3 free-duplex control NU refine J39 (final iteration _003; _000 analysed too as 'free_it0')
Per map, for a scan of contour levels (enclosed volume of the largest connected component 40k-120k A^3; primary
65,000 A^3 = sfGFP 33,000 + full 50 bp duplex 32,000, i.e. the generous end: if the DNA is shorter this contour
over-encloses, which can only make the rod look LONGER, so it biases against an H1 call):
  * segmentation: compact lobe = morphological opening of the contour mask with a ball whose radius lies between
    the rod and lobe half-widths (both measured from the distance transform), component holding the density max;
    rod = rest of the mask.
  * skeleton: TEASAR-style medial geodesic paths (Dijkstra on 26-connected voxels, cost penalised near the
    surface) from the lobe seed to every distal tip; each branch >= 25 A is reported.
  * per branch: arc length from lobe centre and from lobe surface, end-to-end distance, kink/bend angles,
    attachment point and angle (vs radial direction and vs lobe long axis), rod diameter (contour: 2x medial EDT;
    contour-free: FWHM of the radial density profile around the axis), second contact with the lobe.
  * topology: number of tunnels (Betti-1 from the cubical Euler characteristic) -> a closed loop shows as >= 1.
  * occupancy profile: density of the 12 A map in an 8 A tube along the main path, extended 60 A past the tip,
    normalised to the proximal rod; the half-occupancy end is a contour-free length (a mixture of short and long
    DNA shows as a gradual decay rather than a step).
  * helical repeat (best maps): cylinder (r 4-10 A) unwrapped along the rod in the unfiltered map and both half
    maps; power at axial period 30-38 A for azimuthal orders n = 1, 2, and half-map phase consistency.
Then ChimeraX renders principal-axis aligned maps (scripts/46_dotrod_render_fit.cxc) and rigid-fits sfGFP
(L6_protein.pdb) into the lobe as validation only (scripts/46_dotrod_fit_sfgfp.cxc; maps are never modified).
The 'poses' stage uses the NU particle poses to test 2D foreshortening and re-projects the map per 2D class.

  python -I scripts/46_dotrod_3d_length.py analyse     # numbers, aligned display maps, QC figure
  python -I scripts/46_dotrod_3d_length.py poses       # pose-based foreshortening + re-projection tests
  $CHIMERAX --nogui --offscreen --silent scripts/46_dotrod_render_fit.cxc
  $CHIMERAX --nogui --offscreen scripts/46_dotrod_fit_sfgfp.cxc > $OUTD/chimerax_fit_log.txt
  python -I scripts/46_dotrod_3d_length.py finalize    # parse fitmap log, montages, final JSON
"""
import json
import os
import re
import sys
import time

import matplotlib
import mrcfile
import numpy as np
from scipy import ndimage as ndi
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import dijkstra

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

REPO = "/n/home05/ianpoon/DNA-springs-project-2026-27/cryoEM/sfgfp_chimera_screening"
CS = "/n/netscratch/hekstra_lab/Lab/ian_poon/cryosparc_projects/CS-cs-sfgfp-chimera-screening"
RELION = "/n/netscratch/hekstra_lab/Lab/ian_poon/sfgfp_chimera_screening_20261112/relion/W4_dotrod"
OUTD = "/n/netscratch/hekstra_lab/Lab/ian_poon/sfgfp_chimera_screening_20261112/dotrod_3d_length"
JSON_OUT = os.path.join(REPO, "results", "dotrod_3d_length.json")
FIG_QC = os.path.join(REPO, "figures", "qc", "dotrod_3d_length.png")
FIG_MAPS = os.path.join(REPO, "figures", "maps")
PDB_SFGFP = os.path.join(REPO, "model", "L6_loop50", "L6_protein.pdb")
PDB_DNA = os.path.join(REPO, "model", "L6_loop50", "dna", "straight_50bp.pdb")

APIX = 2.84
N = 160
LP_RES = 12.0            # segmentation low-pass (A, half amplitude)
CONTOURS = [40000, 50000, 65000, 80000, 100000, 120000]
PRIMARY = 65000
MIN_BRANCH = 25.0        # A of new path beyond already-covered skeleton
DNA_FULL = 170.0         # 50 bp straight B-DNA
SFGFP_VOL = 33000.0

J59_SIZES = [7340, 5528, 2896, 7457, 4881, 5527, 3927, 5599, 4713, 6280]
MAPS = (
    [dict(key="nu", label="NU refine J56 (dot-on-rod)", path=f"{CS}/J56/J56_003_volume_map.mrc",
          halves=(f"{CS}/J56/J56_003_volume_map_half_A.mrc", f"{CS}/J56/J56_003_volume_map_half_B.mrc"),
          sharp=f"{CS}/J56/J56_003_volume_map_sharp.mrc", n=54181, group="nu", helix=True),
     dict(key="het0", label="hetero J55 class 0 (dot-on-rod)", path=f"{CS}/J55/J55_class_00_00102_volume.mrc",
          n=54181, group="hetero", helix=True),
     dict(key="het1", label="hetero J55 class 1 (free-rod seeded)", path=f"{CS}/J55/J55_class_01_00102_volume.mrc",
          n=42863, group="hetero", helix=True),
     dict(key="free", label="free-DNA control J39 (final it.)", path=f"{CS}/J39/J39_003_volume_map.mrc",
          halves=(f"{CS}/J39/J39_003_volume_map_half_A.mrc", f"{CS}/J39/J39_003_volume_map_half_B.mrc"),
          sharp=f"{CS}/J39/J39_003_volume_map_sharp.mrc", n=None, group="free", helix=True),
     dict(key="free_it0", label="free-DNA control J39 _000 (as named in task)",
          path=f"{CS}/J39/J39_000_volume_map.mrc", n=None, group="free", helix=False),
     dict(key="relion_pp", label="RELION postprocess W4 (cs mask, sharpened)",
          path=f"{RELION}/pp_csmask/postprocess.mrc", n=54181, group="nu", helix=True)]
    + [dict(key=f"c3d_{i}", label=f"class_3D J59 class {i}", path=f"{CS}/J59/J59_class_{i:02d}_00014_volume.mrc",
            n=J59_SIZES[i], group="class3d", helix=False) for i in range(10)]
)
RENDER_KEYS = ["nu", "het0", "het1", "free"] + [f"c3d_{i}" for i in range(10)]


# ----------------------------------------------------------------------------------------------- map utilities
def load_map(path):
    with mrcfile.open(path, permissive=True) as m:
        d = np.asarray(m.data, dtype=np.float64)
        apix = float(m.voxel_size.x)
    n = d.shape[0]
    if n != N:
        if abs(apix - APIX) < 0.01 and n < N:
            # same sampling, smaller box (J39: 128 px = 363 A): zero-pad about the box centre
            s = N // 2 - n // 2
            d = np.pad(d - np.median(d[[0, -1]]), [(s, N - n - s)] * 3)
        elif abs(n * apix - N * APIX) > 1.0:
            raise SystemExit(f"box mismatch {path}: {n} x {apix}")
        else:
            d = fourier_resample(d, N)
    return d, apix


def fourier_resample(d, n_out):
    n = d.shape[0]
    F = np.fft.fftshift(np.fft.fftn(d))
    out = np.zeros((n_out,) * 3, complex)
    s = n_out // 2 - n // 2
    out[s:s + n, s:s + n, s:s + n] = F
    return np.real(np.fft.ifftn(np.fft.ifftshift(out))) * (n_out / n) ** 3


_K = None


def kgrid():
    global _K
    if _K is None:
        f = np.fft.fftfreq(N, APIX)
        kx, ky, kz = np.meshgrid(f, f, f, indexing="ij")
        _K = np.sqrt(kx ** 2 + ky ** 2 + kz ** 2)
    return _K


def lowpass(d, res):
    return np.real(np.fft.ifftn(np.fft.fftn(d) * np.exp(-np.log(2) * (kgrid() * res) ** 2)))


_R = None


def radius_grid():
    global _R
    if _R is None:
        c = (np.arange(N) - N // 2) * APIX
        z, y, x = np.meshgrid(c, c, c, indexing="ij")
        _R = np.sqrt(x ** 2 + y ** 2 + z ** 2)
    return _R


def normalise(d):
    """Zero the background median and scale to the background SD (shell 150-178 A from the box centre; inside the 363 A box of J39)."""
    sh = (radius_grid() > 150) & (radius_grid() < 178)
    return (d - np.median(d[sh])) / d[sh].std()


S26 = np.ones((3, 3, 3), bool)


def largest_component(m):
    m = m & (radius_grid() < 205)
    lab, n = ndi.label(m, S26)
    if n == 0:
        return m, 0, []
    sz = np.bincount(lab.ravel())[1:]
    big = int(np.argmax(sz)) + 1
    others = sorted([float(s) * APIX ** 3 for i, s in enumerate(sz) if i + 1 != big and s * APIX ** 3 > 2000],
                    reverse=True)
    return lab == big, n, others


def threshold_for_volume(lp, V):
    """Bisection on the threshold so the largest connected component encloses V A^3."""
    hi, lo = float(lp.max()), float(np.percentile(lp, 50))
    for _ in range(34):
        t = 0.5 * (hi + lo)
        m, _, _ = largest_component(lp > t)
        v = m.sum() * APIX ** 3
        if v > V:
            lo = t
        else:
            hi = t
        if abs(v - V) / V < 0.01:
            break
    m, ncomp, others = largest_component(lp > t)
    return t, m, ncomp, others


def ball(r_vox):
    r = int(np.ceil(r_vox))
    g = np.arange(-r, r + 1)
    z, y, x = np.meshgrid(g, g, g, indexing="ij")
    return x ** 2 + y ** 2 + z ** 2 <= r_vox ** 2 + 1e-6


def euler_tunnels(m):
    """Betti numbers of the union of closed voxel cubes: b0 (26-conn), b2 (cavities), b1 = tunnels."""
    p = np.pad(m, 1).astype(bool)
    # vertices: lattice points touching any of 8 voxels
    v = np.zeros(tuple(s - 1 for s in p.shape), bool)
    for a in (0, 1):
        for b in (0, 1):
            for c in (0, 1):
                v |= p[a:a + p.shape[0] - 1, b:b + p.shape[1] - 1, c:c + p.shape[2] - 1]
    nV = int(v.sum())
    nE = 0
    for ax in range(3):  # edges parallel to axis ax touch 4 voxels (the two other axes, offsets 0/1)
        o = [i for i in range(3) if i != ax]
        e = None
        for a in (0, 1):
            for b in (0, 1):
                sl = [slice(None)] * 3
                sl[ax] = slice(1, p.shape[ax] - 1)
                sl[o[0]] = slice(a, a + p.shape[o[0]] - 1)
                sl[o[1]] = slice(b, b + p.shape[o[1]] - 1)
                t = p[tuple(sl)]
                e = t if e is None else (e | t)
        nE += int(e.sum())
    nF = 0
    for ax in range(3):  # faces normal to ax touch 2 voxels along ax
        sl0 = [slice(1, s - 1) for s in p.shape]
        sl1 = list(sl0)
        sl0[ax] = slice(0, p.shape[ax] - 1)
        sl1[ax] = slice(1, p.shape[ax])
        nF += int((p[tuple(sl0)] | p[tuple(sl1)]).sum())
    nC = int(p.sum())
    chi = nV - nE + nF - nC
    b0 = ndi.label(p, S26)[1]
    comp_lab, ncomp = ndi.label(~p)  # 6-connected complement
    b2 = ncomp - 1
    return dict(b0=int(b0), b1=int(b0 + b2 - chi), b2=int(b2), chi=int(chi))


# ----------------------------------------------------------------------------------------------- skeleton
def build_graph(mask, edt):
    idx = -np.ones(mask.shape, np.int64)
    pts = np.argwhere(mask)
    idx[tuple(pts.T)] = np.arange(len(pts))
    emax = edt[mask].max()
    rows, cols, wl, wp = [], [], [], []
    offs = [(a, b, c) for a in (-1, 0, 1) for b in (-1, 0, 1) for c in (-1, 0, 1) if (a, b, c) > (0, 0, 0)]
    for o in offs:
        q = pts + np.array(o)
        ok = np.all((q >= 0) & (q < N), axis=1)
        p0, q = pts[ok], q[ok]
        j = idx[tuple(q.T)]
        good = j >= 0
        i0 = idx[tuple(p0[good].T)]
        j = j[good]
        ln = APIX * np.sqrt(sum(x * x for x in o))
        em = 0.5 * (edt[tuple(p0[good].T)] + edt[tuple(q[good].T)])
        rows.append(i0)
        cols.append(j)
        wl.append(np.full(len(j), ln))
        wp.append(ln * (1.0 + 10.0 * (1.0 - em / emax) ** 4))
    rows, cols = np.concatenate(rows), np.concatenate(cols)
    n = len(pts)
    Gl = coo_matrix((np.concatenate(wl), (rows, cols)), shape=(n, n)).tocsr()
    Gp = coo_matrix((np.concatenate(wp), (rows, cols)), shape=(n, n)).tocsr()
    return pts, idx, Gl, Gp


def trace(pred, src, dst):
    path = [dst]
    while path[-1] != src:
        p = pred[path[-1]]
        if p < 0:
            return None
        path.append(p)
    return path[::-1]


def teasar(mask, edt, seed_vox, max_iter=80):
    pts, idx, Gl, Gp = build_graph(mask, edt)
    s = int(idx[tuple(seed_vox)])
    dg = dijkstra(Gl, directed=False, indices=s)
    _, pred = dijkstra(Gp, directed=False, indices=s, return_predecessors=True)
    xyz = pts * APIX
    covered = np.zeros(len(pts), bool)
    covered[s] = True
    branches = []
    for _ in range(max_iter):
        cand = np.where(~covered & np.isfinite(dg))[0]
        if len(cand) == 0:
            break
        end = cand[np.argmax(dg[cand])]
        path = trace(pred, s, end)
        if path is None:
            covered[end] = True
            continue
        path = np.array(path)
        incov = covered[path]
        # new part = after the last covered node on the path
        last_cov = np.where(incov)[0].max() if incov.any() else 0
        newpart = path[last_cov:]
        newlen = np.sum(np.linalg.norm(np.diff(xyz[newpart], axis=0), axis=1)) if len(newpart) > 1 else 0.0
        # cover voxels near the path
        rad = np.minimum(edt[tuple(pts[path].T)] * 1.5 + 2.0 * APIX, 20.0)
        for p, r in zip(path[max(0, last_cov - 1):], rad[max(0, last_cov - 1):]):
            covered |= np.sum((xyz - xyz[p]) ** 2, axis=1) <= r * r
        covered[end] = True
        if newlen >= MIN_BRANCH or not branches:
            branches.append(dict(path=path, junction_index=int(last_cov), new_len=float(newlen)))
    return pts, branches, dg


def smooth_path(xyz, sigma_A=5.0):
    if len(xyz) < 4:
        return xyz.copy()
    # resample at ~1 A first so the smoothing is in A
    seg = np.linalg.norm(np.diff(xyz, axis=0), axis=1)
    s = np.concatenate([[0], np.cumsum(seg)])
    ss = np.arange(0, s[-1] + 1e-9, 1.0)
    r = np.stack([np.interp(ss, s, xyz[:, k]) for k in range(3)], 1)
    sm = ndi.gaussian_filter1d(r, sigma_A, axis=0, mode="nearest")
    sm[0], sm[-1] = r[0], r[-1]
    return sm


def arclen(p):
    return np.concatenate([[0], np.cumsum(np.linalg.norm(np.diff(p, axis=0), axis=1))])


def unit(v):
    n = np.linalg.norm(v)
    return v / n if n > 0 else v


def angle(a, b):
    return float(np.degrees(np.arccos(np.clip(np.dot(unit(a), unit(b)), -1, 1))))


def direction(p, s, s0, s1):
    """Mean direction of path p over arc window [s0, s1] (principal axis, oriented with increasing s)."""
    sel = (s >= s0) & (s <= s1)
    q = p[sel]
    if len(q) < 3:
        return None
    u = np.linalg.svd(q - q.mean(0), full_matrices=False)[2][0]
    return u if np.dot(u, q[-1] - q[0]) >= 0 else -u


def kink_profile(p, s, s_start, s_end, w=15.0):
    out = []
    for x in np.arange(s_start + w, s_end - w + 1e-9, 1.0):
        a = direction(p, s, x - w, x)
        b = direction(p, s, x, x + w)
        if a is not None and b is not None:
            out.append((float(x), angle(a, b)))
    return out


def sample_tube(vol, p, s, radii=(0.0, 4.0, 8.0), nang=12):
    """Mean density on discs normal to the path at each path point (vol in zyx, p in xyz A)."""
    T = np.gradient(p, axis=0)
    T /= np.linalg.norm(T, axis=1, keepdims=True) + 1e-12
    vals = []
    for i in range(len(p)):
        t = T[i]
        a = unit(np.cross(t, [1, 0, 0] if abs(t[0]) < 0.9 else [0, 1, 0]))
        b = np.cross(t, a)
        q = [p[i]]
        for r in radii:
            if r > 0:
                for th in np.linspace(0, 2 * np.pi, nang, endpoint=False):
                    q.append(p[i] + r * (np.cos(th) * a + np.sin(th) * b))
        q = np.array(q) / APIX
        vals.append(ndi.map_coordinates(vol, [q[:, 2], q[:, 1], q[:, 0]], order=1, mode="constant").mean())
    return np.array(vals)


def radial_profile(vol, p, s, s0, s1, rmax=20.0):
    sel = np.where((s >= s0) & (s <= s1))[0]
    if len(sel) < 5:
        return None
    T = np.gradient(p, axis=0)
    T /= np.linalg.norm(T, axis=1, keepdims=True) + 1e-12
    rr = np.arange(0, rmax + 0.1, 1.0)
    prof = np.zeros(len(rr))
    for i in sel:
        t = T[i]
        a = unit(np.cross(t, [1, 0, 0] if abs(t[0]) < 0.9 else [0, 1, 0]))
        b = np.cross(t, a)
        for k, r in enumerate(rr):
            th = np.linspace(0, 2 * np.pi, 24, endpoint=False)
            q = (p[i][None] + r * (np.cos(th)[:, None] * a + np.sin(th)[:, None] * b)) / APIX
            prof[k] += ndi.map_coordinates(vol, [q[:, 2], q[:, 1], q[:, 0]], order=1).mean()
    prof /= len(sel)
    base = prof[-4:].mean()
    pk = prof[0] - base
    half = base + 0.5 * pk
    below = np.where(prof < half)[0]
    if len(below) == 0 or pk <= 0:
        return dict(fwhm=None, profile=prof.tolist())
    k = below[0]
    r_half = rr[k - 1] + (prof[k - 1] - half) / (prof[k - 1] - prof[k]) * (rr[k] - rr[k - 1])
    return dict(fwhm=float(2 * r_half), profile=[round(float(x), 4) for x in prof])


# ----------------------------------------------------------------------------------------------- per-contour analysis
def analyse_contour(lp, thr, mask, others):
    edt = ndi.distance_transform_edt(mask) * APIX
    pts_all = np.argwhere(mask)
    seed = pts_all[np.argmax(lp[mask])]
    xyz_seed = seed[::-1] * APIX  # zyx index -> xyz A
    # main path first (for rod radius)
    pts, branches, dg = teasar(mask, edt, seed)
    xyz = pts[:, ::-1] * APIX
    lobe_r = float(edt[tuple(pts[np.sum((xyz - xyz_seed) ** 2, 1) <= 10.0 ** 2].T)].max())
    main = branches[0]["path"]
    dmain = dg[main]
    sel = main[(dmain > lobe_r + 15) & (dmain < dmain.max() - 10)]
    if len(sel) <= 3:  # short main path: use everything beyond the lobe half-width
        sel = main[dmain > lobe_r]
    rod_r = float(np.median(edt[tuple(pts[sel].T)])) if len(sel) > 3 else 10.0
    distinct = bool(np.isfinite(rod_r) and lobe_r / rod_r >= 1.2)
    R_open = float(np.clip(0.5 * (lobe_r + rod_r), rod_r + 1.0, lobe_r - 0.5)) if distinct else rod_r + 1.0
    opened = ndi.binary_opening(mask, structure=ball(R_open / APIX))
    olab, on = ndi.label(opened, S26)
    lobe_id = olab[tuple(seed)]
    if lobe_id == 0 and on > 0:
        # nearest opened component to the seed
        op = np.argwhere(opened)
        lobe_id = olab[tuple(op[np.argmin(np.sum((op - seed) ** 2, 1))])]
    lobe = olab == lobe_id if lobe_id > 0 else (np.sum((np.indices(mask.shape).T - seed).T ** 2, 0) * APIX ** 2
                                                <= lobe_r ** 2) & mask
    # other compact blobs (distal lobes) at the same opening radius
    blobs = []
    for i in range(1, on + 1):
        if i == lobe_id:
            continue
        b = olab == i
        v = b.sum() * APIX ** 3
        if v >= 3000:
            c = np.argwhere(b).mean(0)[::-1] * APIX
            blobs.append(dict(volume_A3=round(float(v)), centre_xyz=[round(float(x), 1) for x in c],
                              dist_from_lobe_centre_A=round(float(np.linalg.norm(c - xyz_seed)), 1),
                              max_edt_A=round(float(edt[b].max()), 1)))
    lv = np.argwhere(lobe)[:, ::-1] * APIX
    lc = lv.mean(0)
    ev, evec = np.linalg.eigh(np.cov((lv - lc).T))
    ev, evec = ev[::-1], evec[:, ::-1]
    proj = (lv - lc) @ evec
    lobe_info = dict(volume_A3=round(float(lobe.sum() * APIX ** 3)),
                     extent_A=[round(float(x), 1) for x in proj.max(0) - proj.min(0)],
                     ellipsoid_axes_A=[round(float(2 * np.sqrt(5 * max(e, 0))), 1) for e in ev],
                     centre_xyz=[round(float(x), 1) for x in lc], seed_xyz=[round(float(x), 1) for x in xyz_seed],
                     lobe_halfwidth_edt_A=round(lobe_r, 1), rod_radius_edt_A=round(rod_r, 1),
                     opening_radius_A=round(R_open, 1), distinct_lobe=distinct)
    lobe_long = evec[:, 0]
    lobe_set = set(map(tuple, np.argwhere(lobe)))
    lobe_d = ndi.distance_transform_edt(~lobe) * APIX  # distance to lobe
    in_lobe = np.array([tuple(p) in lobe_set for p in pts])
    br_out = []
    for bi, b in enumerate(branches):
        path = b["path"]
        P = smooth_path(xyz[path])
        s = arclen(P)
        raw_in = in_lobe[path]
        # lobe exit on the raw path, mapped to arc length on the smoothed path
        out_idx = np.where(~raw_in)[0]
        if len(out_idx) == 0:
            continue
        # last time the path is inside the lobe
        k_exit = (np.where(raw_in)[0].max() + 1) if raw_in.any() else 0
        k_exit = min(k_exit, len(path) - 1)
        raw_s = arclen(xyz[path])
        s_exit = float(raw_s[k_exit] / raw_s[-1] * s[-1]) if raw_s[-1] > 0 else 0.0
        s_junc = float(raw_s[b["junction_index"]] / raw_s[-1] * s[-1]) if raw_s[-1] > 0 else 0.0
        jv = pts[path[b["junction_index"]]]
        starts_at_lobe = bool(lobe_d[tuple(jv)] <= 10.0) or bi == 0
        L_center = float(s[-1])
        L_surface = float(s[-1] - s_exit)
        p_exit = P[np.searchsorted(s, s_exit).clip(0, len(P) - 1)]
        e2e = float(np.linalg.norm(P[-1] - p_exit))
        d_in = direction(P, s, s_exit, s_exit + 25)
        d_far = direction(P, s, s[-1] - 30, s[-1])
        d_near = direction(P, s, s_exit, s_exit + 30)
        att_radial = angle(d_in, p_exit - lc) if d_in is not None else None
        att_long = angle(d_in, lobe_long) if d_in is not None else None
        if att_long is not None:
            att_long = min(att_long, 180 - att_long)
        kinks = kink_profile(P, s, s_exit, s[-1])
        kmax = max(kinks, key=lambda x: x[1]) if kinks else (None, None)
        T = np.gradient(P, axis=0)
        T /= np.linalg.norm(T, axis=1, keepdims=True) + 1e-12
        sel = s >= s_exit
        turn = float(np.degrees(np.sum(np.arccos(np.clip(np.sum(T[sel][1:] * T[sel][:-1], 1), -1, 1)))))
        # medial EDT diameter along rod (excluding 8 A at each end)
        rsel = (raw_s >= raw_s[k_exit] + 8) & (raw_s <= raw_s[-1] - 8)
        diam = float(2 * np.median(edt[tuple(pts[path[rsel]].T)])) if rsel.sum() > 2 else None
        # second contact: distal path points (>40 A past exit) vs lobe
        far = P[s >= s_exit + 40] / APIX
        if len(far):
            dl = ndi.map_coordinates(lobe_d, [far[:, 2], far[:, 1], far[:, 0]], order=1)
            min_back = float(dl.min())
            min_back_s = float(s[s >= s_exit + 40][np.argmin(dl)])
        else:
            min_back, min_back_s = None, None
        br_out.append(dict(
            branch=bi, starts_at_lobe=starts_at_lobe, junction_arc_A=round(s_junc, 1), new_len_A=round(b["new_len"], 1),
            L_from_lobe_centre_A=round(L_center, 1), L_from_lobe_surface_A=round(L_surface, 1),
            end_to_end_from_surface_A=round(e2e, 1), tortuosity=round(L_surface / e2e, 3) if e2e > 0 else None,
            tip_xyz=[round(float(x), 1) for x in P[-1]], attach_xyz=[round(float(x), 1) for x in p_exit],
            attach_angle_vs_radial_deg=None if att_radial is None else round(att_radial, 1),
            attach_angle_vs_lobe_long_axis_deg=None if att_long is None else round(att_long, 1),
            near_vs_far_angle_deg=round(angle(d_near, d_far), 1) if (d_near is not None and d_far is not None
                                                                      and s[-1] - s_exit > 45) else None,
            max_kink_deg=None if kmax[1] is None else round(kmax[1], 1),
            max_kink_at_A_from_surface=None if kmax[0] is None else round(kmax[0] - s_exit, 1),
            total_turning_deg=round(turn, 1), diameter_contour_edt_A=None if diam is None else round(diam, 1),
            min_dist_back_to_lobe_A=None if min_back is None else round(min_back, 1),
            min_dist_back_at_arc_A=None if min_back_s is None else round(min_back_s - s_exit, 1),
            second_contact=bool(min_back is not None and min_back <= 1.5 * APIX),
            _P=P, _s=s, _s_exit=s_exit))
    # tips of distinct branches approaching each other (DNA-DNA closure)
    tips = [b["_P"][-1] for b in br_out]
    tip_pairs = [dict(a=br_out[i]["branch"], b=br_out[j]["branch"], dist_A=round(float(np.linalg.norm(tips[i] - tips[j])), 1))
                 for i in range(len(tips)) for j in range(i + 1, len(tips))]
    # rod contact patches
    rod = mask & ~ndi.binary_dilation(lobe, S26)
    rl, rn = ndi.label(rod, S26)
    lobe_shell = ndi.binary_dilation(lobe, S26, iterations=2) & ~lobe & mask
    multi_patch = []
    for i in range(1, rn + 1):
        comp = rl == i
        if comp.sum() * APIX ** 3 < 1500:
            continue
        touch = ndi.binary_dilation(comp, S26, iterations=2) & lobe_shell
        pl, pn = ndi.label(touch, S26)
        cents = [np.argwhere(pl == k).mean(0) * APIX for k in range(1, pn + 1)]
        far_patches = pn > 1 and max(np.linalg.norm(a - b) for a in cents for b in cents) > 15
        multi_patch.append(dict(rod_component_vol_A3=round(float(comp.sum() * APIX ** 3)), n_patches=int(pn),
                                patches_separated=bool(far_patches)))
    topo = euler_tunnels(mask)
    return dict(threshold_sigma=round(float(thr), 3), volume_A3=round(float(mask.sum() * APIX ** 3)),
                other_components_A3=[round(x) for x in others[:5]], lobe=lobe_info, distal_blobs=blobs,
                branches=br_out, tip_pairs=tip_pairs, rod_components=multi_patch, topology=topo,
                _lobe=lobe, _mask=mask, _edt=edt)


# ----------------------------------------------------------------------------------------------- helical repeat
def rmf_frames(P):
    T = np.gradient(P, axis=0)
    T /= np.linalg.norm(T, axis=1, keepdims=True) + 1e-12
    a = unit(np.cross(T[0], [1, 0, 0] if abs(T[0][0]) < 0.9 else [0, 1, 0]))
    A = [a]
    for i in range(1, len(P)):
        a = A[-1] - np.dot(A[-1], T[i]) * T[i]
        A.append(unit(a))
    A = np.array(A)
    return T, A, np.cross(T, A)


def unwrap_cylinder(vol, P, radii=(4, 6, 8, 10), nphi=32):
    T, A, Bv = rmf_frames(P)
    phi = np.linspace(0, 2 * np.pi, nphi, endpoint=False)
    out = np.zeros((len(radii), len(P), nphi))
    for ir, r in enumerate(radii):
        q = (P[:, None, :] + r * (np.cos(phi)[None, :, None] * A[:, None, :] + np.sin(phi)[None, :, None]
                                  * Bv[:, None, :])) / APIX
        q = q.reshape(-1, 3)
        out[ir] = ndi.map_coordinates(vol, [q[:, 2], q[:, 1], q[:, 0]], order=1).reshape(len(P), nphi)
    return out


def helix_power(cyl, pad=512):
    """Power vs axial frequency for azimuthal orders n=0..3 (cylinder averaged over radii)."""
    c = cyl - cyl.mean(axis=(1, 2), keepdims=True)
    L = c.shape[1]
    w = np.hanning(L)[None, :, None]
    F = np.fft.fft(np.fft.fft(c * w, n=pad, axis=1), axis=2)
    Pw = (np.abs(F) ** 2).mean(0)  # (kz, n)
    kz = np.fft.fftfreq(pad, 1.0)
    pos = kz > 0
    res = {}
    for n in range(4):
        pn = Pw[:, n] + (Pw[:, -n] if n > 0 else 0)
        res[n] = (kz[pos], pn[pos])
    return res


def helix_metrics(vol, P, halves=None):
    cyl = unwrap_cylinder(vol, P)
    pw = helix_power(cyl)
    out = {"segment_length_A": round(float(arclen(P)[-1]), 1)}
    for n in (0, 1, 2):
        kz, p = pw[n]
        band = (kz >= 1 / 38) & (kz <= 1 / 30)
        ctrl = ((kz >= 1 / 70) & (kz < 1 / 45)) | ((kz > 1 / 26) & (kz <= 1 / 16))
        allb = (kz >= 1 / 60) & (kz <= 1 / 16)
        out[f"n{n}_band30_38_over_ctrl"] = round(float(p[band].mean() / p[ctrl].mean()), 3) if ctrl.any() else None
        kk = kz[allb][np.argmax(p[allb])]
        out[f"n{n}_peak_period_A"] = round(float(1 / kk), 1)
    out["_spectra"] = {n: (pw[n][0].tolist(), pw[n][1].tolist()) for n in (0, 1, 2)}
    if halves is not None:
        ca = unwrap_cylinder(halves[0], P)
        cb = unwrap_cylinder(halves[1], P)

        def bandpass(c):
            c = c - c.mean(axis=2, keepdims=True)  # remove n=0
            L = c.shape[1]
            F = np.fft.fft(c, axis=1)
            k = np.abs(np.fft.fftfreq(L, 1.0))
            F[:, (k < 1 / 42) | (k > 1 / 28), :] = 0
            return np.real(np.fft.ifft(F, axis=1))

        a, b = bandpass(ca), bandpass(cb)
        rs = []
        for sh in range(a.shape[2]):
            bb = np.roll(b, sh, axis=2)
            rs.append(float(np.corrcoef(a.ravel(), bb.ravel())[0, 1]))
        out["halfmap_band28_42_corr_at_0"] = round(rs[0], 3)
        out["halfmap_band28_42_corr_rank_of_32_rotations"] = int(1 + np.sum(np.array(rs[1:]) >= rs[0]))
        out["halfmap_band28_42_corr_max_rotation"] = round(max(rs), 3)
    return out


# ----------------------------------------------------------------------------------------------- frames for display
def pca_frame(mask):
    pts = np.argwhere(mask)[:, ::-1] * APIX
    c = pts.mean(0)
    ev, evec = np.linalg.eigh(np.cov((pts - c).T))
    evec = evec[:, ::-1]
    if np.linalg.det(evec) < 0:
        evec[:, 2] *= -1
    return c, evec


def to_frame(vol, c, Rm, order=1):
    """Resample vol (zyx) so that output x,y,z = PC1,PC2,PC3 through centroid c, centred in the box."""
    F = np.eye(3)[::-1]
    Mz = F @ Rm @ F
    cbox = np.array([N // 2] * 3, float)
    cin = c[::-1] / APIX
    return ndi.affine_transform(vol, Mz, offset=cin - Mz @ cbox, order=order, mode="constant")


def xyz_to_frame(p, c, Rm):
    return (p - c) @ Rm + (N // 2) * APIX


def write_mrc(path, d):
    with mrcfile.new(path, overwrite=True) as m:
        m.set_data(d.astype(np.float32))
        m.voxel_size = APIX


def write_cmm(path, P, color=(1.0, 0.55, 0.0), r=2.0, step=4):
    P = P[::step] if len(P) > step else P
    with open(path, "w") as f:
        f.write('<marker_set name="skeleton">\n')
        for i, p in enumerate(P):
            f.write(f'<marker id="{i + 1}" x="{p[0]:.2f}" y="{p[1]:.2f}" z="{p[2]:.2f}" r="{color[0]}" g="{color[1]}" '
                    f'b="{color[2]}" radius="{r}"/>\n')
        for i in range(len(P) - 1):
            f.write(f'<link id1="{i + 1}" id2="{i + 2}" r="{color[0]}" g="{color[1]}" b="{color[2]}" radius="{r * 0.6}"/>\n')
        f.write("</marker_set>\n")


def write_pdb_centred(src, dst, centre_xyz):
    lines = [ln for ln in open(src) if ln.startswith(("ATOM", "HETATM"))]
    xyz = np.array([[float(ln[30:38]), float(ln[38:46]), float(ln[46:54])] for ln in lines])
    sh = np.asarray(centre_xyz) - xyz.mean(0)
    with open(dst, "w") as f:
        for ln, p in zip(lines, xyz + sh):
            f.write(f"{ln[:30]}{p[0]:8.3f}{p[1]:8.3f}{p[2]:8.3f}{ln[54:]}")
        f.write("END\n")


def model_dims(path):
    xyz = np.array([[float(ln[30:38]), float(ln[38:46]), float(ln[46:54])] for ln in open(path)
                    if ln.startswith(("ATOM", "HETATM"))])
    c = xyz - xyz.mean(0)
    p = c @ np.linalg.svd(c, full_matrices=False)[2].T
    return [round(float(x), 1) for x in p.max(0) - p.min(0)]


def save_json(obj):
    def dflt(o):
        if isinstance(o, np.ndarray):
            return o.tolist()
        if isinstance(o, (np.floating, np.integer, np.bool_)):
            return o.item()
        raise TypeError(type(o))
    tmp = JSON_OUT + ".tmp"
    with open(tmp, "w") as fh:
        json.dump(obj, fh, indent=1, default=dflt)
    os.replace(tmp, JSON_OUT)


# ----------------------------------------------------------------------------------------------- driver
def strip(o):
    if isinstance(o, dict):
        return {k: strip(v) for k, v in o.items() if not k.startswith("_")}
    if isinstance(o, list):
        return [strip(v) for v in o]
    if isinstance(o, (np.floating,)):
        return float(o)
    if isinstance(o, (np.integer,)):
        return int(o)
    if isinstance(o, np.bool_):
        return bool(o)
    return o


def main_branch(res):
    brs = [b for b in res["branches"] if b["starts_at_lobe"]]
    return max(brs, key=lambda b: b["L_from_lobe_surface_A"]) if brs else None


def full_path(res):
    """Longest tip-to-tip medial path (for the free rod: the whole duplex)."""
    mask, edt = res["_mask"], res["_edt"]
    mb = main_branch(res)
    tip = np.round(mb["_P"][-1] / APIX).astype(int)[::-1]
    pts = np.argwhere(mask)
    tip = pts[np.argmin(np.sum((pts - tip) ** 2, 1))]
    pts2, br2, _ = teasar(mask, edt, tip, max_iter=1)
    P = smooth_path(pts2[br2[0]["path"]][:, ::-1] * APIX)
    return P


def analyse():
    os.makedirs(OUTD, exist_ok=True)
    t0 = time.time()
    results, store = {}, {}
    only = [k for k in os.environ.get("DOTROD_KEYS", "").split(",") if k]  # debugging subset
    for M in MAPS:
        if only and M["key"] not in only:
            continue
        if not os.path.exists(M["path"]):
            print("missing", M["path"], flush=True)
            continue
        d, apix0 = load_map(M["path"])
        dn = normalise(d)
        lp = normalise(lowpass(d, LP_RES))
        scan = {}
        for V in CONTOURS:
            thr, mask, ncomp, others = threshold_for_volume(lp, V)
            scan[V] = analyse_contour(lp, thr, mask, others)
        prim = scan[PRIMARY]
        mb = main_branch(prim)
        # occupancy profile along the main branch, extended 60 A past the tip
        P, s, s_exit = mb["_P"], mb["_s"], mb["_s_exit"]
        tdir = direction(P, s, s[-1] - 15, s[-1])
        ext = P[-1] + np.arange(1, 61)[:, None] * tdir[None, :]
        Pe = np.vstack([P, ext])
        se = arclen(Pe)
        prof = sample_tube(lp, Pe, se)
        ref_sel = (se >= s_exit + 15) & (se <= s_exit + 45)
        ref = float(np.median(prof[ref_sel])) if ref_sel.any() else float(prof[se >= s_exit].max())
        norm = prof / ref
        s50 = None
        for i in np.where(se >= s_exit + 15)[0]:
            if norm[i] < 0.5 and np.all(norm[i:i + 10] < 0.5):
                s50 = float(se[i])
                break
        s25 = None
        for i in np.where(se >= s_exit + 15)[0]:
            if norm[i] < 0.25 and np.all(norm[i:i + 10] < 0.25):
                s25 = float(se[i])
                break
        occ = dict(ref_density_sigma=round(ref, 3), lobe_surface_arc_A=round(s_exit, 1),
                   half_occupancy_from_lobe_surface_A=None if s50 is None else round(s50 - s_exit, 1),
                   half_occupancy_from_lobe_centre_A=None if s50 is None else round(s50, 1),
                   quarter_occupancy_from_lobe_surface_A=None if s25 is None else round(s25 - s_exit, 1),
                   contour_tip_from_lobe_surface_A=round(float(s[-1] - s_exit), 1),
                   lobe_peak_over_rod_ref=round(float(prof[se <= max(s_exit, 1)].max() / ref), 2))
        # free rod: full tip-to-tip path and both half-occupancy ends
        Pf = full_path(prim)
        sf = arclen(Pf)
        tfa = direction(Pf, sf, 0, 15)
        tfb = direction(Pf, sf, sf[-1] - 15, sf[-1])
        Pfe = np.vstack([Pf[0] - np.arange(60, 0, -1)[:, None] * tfa[None, :], Pf,
                         Pf[-1] + np.arange(1, 61)[:, None] * tfb[None, :]])
        sfe = arclen(Pfe)
        pf = sample_tube(lp, Pfe, sfe)
        mid = (sfe >= 60 + 0.25 * sf[-1]) & (sfe <= 60 + 0.75 * sf[-1])
        reff = float(np.median(pf[mid]))
        above = np.where(pf / reff >= 0.5)[0]
        full_half = float(sfe[above.max()] - sfe[above.min()]) if len(above) else None
        occ_full = dict(tip_to_tip_contour_A=round(float(sf[-1]), 1),
                        tip_to_tip_half_density_A=None if full_half is None else round(full_half, 1),
                        mid_density_sigma=round(reff, 3))
        # diameter: radial FWHM on an 8 A low-pass, rod portion 10 A past exit to 10 A before tip
        lp8 = normalise(lowpass(d, 8.0))
        rp = radial_profile(lp8, P, s, s_exit + 10, s[-1] - 10)
        hel = None
        if M.get("helix"):
            src = d if M["key"] != "relion_pp" else d
            segP = P[(s >= s_exit + 8) & (s <= s[-1] - 8)]
            if M["group"] == "free":
                segP = Pf[(sf >= 10) & (sf <= sf[-1] - 10)]
            if len(segP) >= 50:
                halves = None
                if M.get("halves"):
                    halves = [lowpass(load_map(h)[0], 6.0) for h in M["halves"]]
                hel = helix_metrics(lowpass(src, 6.0), segP, halves)
                if M.get("sharp"):
                    hs = helix_metrics(lowpass(load_map(M["sharp"])[0], 6.0), segP)
                    hel["sharp_map"] = {k: v for k, v in hs.items() if not k.startswith("_")}
        res = dict(label=M["label"], path=M["path"], group=M["group"], n_particles=M["n"],
                   input_apix=round(apix0, 3), primary=prim, scan={str(V): scan[V] for V in CONTOURS},
                   occupancy=occ, occupancy_full_path=occ_full,
                   radial_fwhm_lp8_A=None if rp is None else rp["fwhm"], radial_profile_lp8=None if rp is None else rp,
                   helix=hel)
        results[M["key"]] = res
        store[M["key"]] = dict(d=d, lp=lp, prof=(se - s_exit, norm), proff=(sfe - 60, pf / reff), Pf=Pf)
        nb = len([b for b in prim["branches"] if b["starts_at_lobe"]])
        print(f"{M['key']:10s} thr={prim['threshold_sigma']:.2f} lobe={prim['lobe']['volume_A3']} "
              f"distinct={prim['lobe']['distinct_lobe']} branches={nb} mainL_surf={mb['L_from_lobe_surface_A']} "
              f"L50={occ['half_occupancy_from_lobe_surface_A']} tunnels={prim['topology']['b1']} "
              f"full={occ_full} t={time.time() - t0:.0f}s", flush=True)
    write_display(results, store)
    summary = summarise(results)
    out = dict(script="scripts/46_dotrod_3d_length.py", date=time.strftime("%Y-%m-%d"),
               parameters=dict(apix=APIX, box=N, lowpass_A=LP_RES, contours_A3=CONTOURS, primary_A3=PRIMARY,
                               min_branch_A=MIN_BRANCH, dna_full_A=DNA_FULL,
                               model_dims_A=dict(sfgfp_L6_protein=model_dims(PDB_SFGFP),
                                                 dna_straight_50bp=model_dims(PDB_DNA))),
               summary=summary, maps=strip(results))
    save_json(out)
    qc_figure(results, store, summary)
    print("done", time.time() - t0, flush=True)


def summarise(results):
    S = {}
    for k, r in results.items():
        p = r["primary"]
        mb = main_branch(p)
        lens = {V: (main_branch(sc)["L_from_lobe_surface_A"] if main_branch(sc) else None) for V, sc in r["scan"].items()}
        nbr = {V: len([b for b in sc["branches"] if b["starts_at_lobe"]]) for V, sc in r["scan"].items()}
        S[k] = dict(n_particles=r["n_particles"], lobe_volume_A3=p["lobe"]["volume_A3"],
                    lobe_extent_A=p["lobe"]["extent_A"], distinct_lobe=p["lobe"]["distinct_lobe"],
                    n_branches_from_lobe=nbr[str(PRIMARY)],
                    branch_lengths_from_surface_A=[b["L_from_lobe_surface_A"] for b in p["branches"] if b["starts_at_lobe"]],
                    main_L_from_surface_A=mb["L_from_lobe_surface_A"], main_L_from_centre_A=mb["L_from_lobe_centre_A"],
                    main_end_to_end_A=mb["end_to_end_from_surface_A"],
                    main_L_vs_contour=lens, n_branches_vs_contour=nbr,
                    half_occupancy_from_surface_A=r["occupancy"]["half_occupancy_from_lobe_surface_A"],
                    tip_to_tip_A=r["occupancy_full_path"]["tip_to_tip_contour_A"],
                    tip_to_tip_half_density_A=r["occupancy_full_path"]["tip_to_tip_half_density_A"],
                    max_kink_deg=mb["max_kink_deg"], near_vs_far_deg=mb["near_vs_far_angle_deg"],
                    attach_vs_radial_deg=mb["attach_angle_vs_radial_deg"],
                    diameter_edt_A=mb["diameter_contour_edt_A"], radial_fwhm_A=r["radial_fwhm_lp8_A"],
                    tunnels_vs_contour={V: sc["topology"]["b1"] for V, sc in r["scan"].items()},
                    second_contact_any_contour=any(b["second_contact"] for sc in r["scan"].values() for b in sc["branches"]),
                    rod_multi_patch_any_contour=any(c["patches_separated"] for sc in r["scan"].values()
                                                    for c in sc["rod_components"]),
                    distal_blobs=p["distal_blobs"])
    cls = [k for k in results if k.startswith("c3d_")]
    if cls:
        w = np.array([results[k]["n_particles"] for k in cls], float)
        L = np.array([S[k]["main_L_from_surface_A"] for k in cls])
        H = np.array([S[k]["half_occupancy_from_surface_A"] or np.nan for k in cls])
        S["class3d_weighted"] = dict(
            weighted_mean_L_from_surface_A=round(float(np.sum(w * L) / w.sum()), 1),
            weighted_sd_L_A=round(float(np.sqrt(np.sum(w * (L - np.sum(w * L) / w.sum()) ** 2) / w.sum())), 1),
            min_max_L_A=[round(float(L.min()), 1), round(float(L.max()), 1)],
            weighted_mean_half_occupancy_A=round(float(np.nansum(w * H) / w[~np.isnan(H)].sum()), 1),
            frac_particles_L_ge_free_minus15=None, per_class={k: (int(results[k]["n_particles"]),
                                                                  S[k]["main_L_from_surface_A"],
                                                                  S[k]["half_occupancy_from_surface_A"],
                                                                  S[k]["n_branches_from_lobe"]) for k in cls})
        if "free" in S:
            fl = S["free"]["tip_to_tip_A"]
            S["class3d_weighted"]["frac_particles_L_ge_free_minus15"] = round(float(w[L >= fl - 15].sum() / w.sum()), 3)
    return S


def write_display(results, store):
    """Principal-axis aligned, normalised display maps (contour 1.0 = primary contour) + skeleton markers."""
    frame_nu = None
    for k in RENDER_KEYS:
        if k not in results:
            continue
        r = results[k]
        prim = r["primary"]
        mask, lobe = prim["_mask"], prim["_lobe"]
        if k == "nu" or not k.startswith(("c3d_", "het0")):
            c, Rm = pca_frame(mask)
            if k == "nu":
                frame_nu = (c, Rm)
        else:
            c, Rm = frame_nu
        lp = store[k]["lp"] / prim["threshold_sigma"]
        if r["group"] == "free":  # no compact lobe in the free-duplex control: show it unsegmented
            lobe = np.zeros_like(lobe)
        soft = ndi.gaussian_filter(ndi.binary_dilation(mask, S26, iterations=2).astype(float), 1.0)
        lob = ndi.gaussian_filter(ndi.binary_dilation(lobe, S26, iterations=2).astype(float), 1.0)
        write_mrc(f"{OUTD}/disp_{k}_all.mrc", to_frame(lp * soft, c, Rm))
        write_mrc(f"{OUTD}/disp_{k}_lobe.mrc", to_frame(lp * soft * lob, c, Rm))
        write_mrc(f"{OUTD}/disp_{k}_rod.mrc", to_frame(lp * soft * (1 - lob), c, Rm))
        cols = [(1.0, 0.55, 0.0), (0.85, 0.1, 0.6), (0.1, 0.7, 0.2), (0.5, 0.2, 0.9), (0.9, 0.9, 0.1)]
        with open(f"{OUTD}/disp_{k}_skel.cmm", "w") as f:
            f.write('<marker_set name="skeleton">\n')
            mid = 1
            for bi, b in enumerate([b for b in prim["branches"] if b["starts_at_lobe"] or True][:5]):
                Q = xyz_to_frame(b["_P"], c, Rm)[::4]
                col = cols[bi % len(cols)]
                ids = []
                for p in Q:
                    f.write(f'<marker id="{mid}" x="{p[0]:.2f}" y="{p[1]:.2f}" z="{p[2]:.2f}" r="{col[0]}" '
                            f'g="{col[1]}" b="{col[2]}" radius="1.6"/>\n')
                    ids.append(mid)
                    mid += 1
                for a, bb in zip(ids[:-1], ids[1:]):
                    f.write(f'<link id1="{a}" id2="{bb}" r="{col[0]}" g="{col[1]}" b="{col[2]}" radius="1.0"/>\n')
            f.write("</marker_set>\n")
    # fitmap inputs: original-frame maps, sfGFP placed at the lobe centre and at each distal blob / free terminal lobe
    for k in ("nu", "free", "het1"):
        if k not in results:
            continue
        prim = results[k]["primary"]
        write_mrc(f"{OUTD}/fit_{k}_map.mrc", store[k]["d"])
        write_pdb_centred(PDB_SFGFP, f"{OUTD}/fit_{k}_sfgfp_lobe.pdb", prim["lobe"]["centre_xyz"])
        tips = [b for b in prim["branches"] if b["starts_at_lobe"]]
        mb = max(tips, key=lambda b: b["L_from_lobe_surface_A"])
        # a control placement on the rod: mid-way along the main branch
        P, s = mb["_P"], mb["_s"]
        mid = P[np.searchsorted(s, mb["_s_exit"] + 0.5 * (s[-1] - mb["_s_exit"]))]
        write_pdb_centred(PDB_SFGFP, f"{OUTD}/fit_{k}_sfgfp_rodmid.pdb", mid)
        blob = prim["distal_blobs"][0]["centre_xyz"] if prim["distal_blobs"] else mb["_P"][-1]
        write_pdb_centred(PDB_SFGFP, f"{OUTD}/fit_{k}_sfgfp_distal.pdb", blob)


def qc_figure(results, store, S):
    fig = plt.figure(figsize=(16, 15))
    gs = fig.add_gridspec(4, 4, hspace=0.45, wspace=0.3)
    # row 0: projections with skeleton for nu, het1, free, c3d best
    keys = [k for k in ("nu", "het0", "het1", "free") if k in results]
    for i, k in enumerate(keys):
        ax = fig.add_subplot(gs[0, i])
        prim = results[k]["primary"]
        c, Rm = pca_frame(prim["_mask"])
        v = to_frame(store[k]["lp"], c, Rm)
        img = v.sum(axis=0)  # sum over z' -> (y', x')
        ax.imshow(img, cmap="gray_r", origin="lower", extent=[0, N * APIX, 0, N * APIX])
        for bi, b in enumerate(prim["branches"]):
            Q = xyz_to_frame(b["_P"], c, Rm)
            ax.plot(Q[:, 0], Q[:, 1], "-", lw=1.5, color=f"C{bi + 1}")
        lc = xyz_to_frame(np.array(prim["lobe"]["centre_xyz"])[None], c, Rm)[0]
        ax.plot(lc[0], lc[1], "o", mfc="none", mec="lime", ms=12)
        ax.set_xlim(N * APIX / 2 - 110, N * APIX / 2 + 110)
        ax.set_ylim(N * APIX / 2 - 110, N * APIX / 2 + 110)
        ax.set_title(f"{k}: proj. on PC1-PC2 plane\nmain {S[k]['main_L_from_surface_A']:.0f} A from lobe surface",
                     fontsize=9)
        ax.set_xlabel("A")
    # row 1 left: length vs contour
    ax = fig.add_subplot(gs[1, 0:2])
    for k in [k for k in ("nu", "het0", "het1", "free", "relion_pp") if k in results]:
        L = [S[k]["main_L_vs_contour"][str(V)] for V in CONTOURS]
        ax.plot(np.array(CONTOURS) / 1000, L, "o-", label=k)
    for k in [k for k in results if k.startswith("c3d_")]:
        L = [S[k]["main_L_vs_contour"][str(V)] for V in CONTOURS]
        ax.plot(np.array(CONTOURS) / 1000, L, "-", color="0.7", lw=0.8)
    if "free" in S:
        ax.axhline(S["free"]["tip_to_tip_A"], color="k", ls="--", lw=0.8, label="free rod tip-to-tip (primary)")
    ax.axhline(DNA_FULL, color="r", ls=":", label="170 A (50 bp)")
    ax.axvline(PRIMARY / 1000, color="0.5", ls=":")
    ax.set_xlabel("contour: enclosed volume of main component (x1000 A^3)")
    ax.set_ylabel("main rod arc length from lobe surface (A)")
    ax.legend(fontsize=7, ncol=2)
    ax.set_title("3D rod length vs contour (grey = 10 class_3D maps)", fontsize=9)
    # row 1 right: occupancy profiles
    ax = fig.add_subplot(gs[1, 2:4])
    for k in [k for k in ("nu", "het0", "het1", "relion_pp") if k in results]:
        x, y = store[k]["prof"]
        ax.plot(x, y, label=f"{k} (from lobe surface)")
    if "free" in store:
        x, y = store["free"]["proff"]
        ax.plot(x, y, "k-", lw=1.2, label="free rod (from its first tip)")
    ax.axhline(0.5, color="0.5", ls=":")
    ax.axvline(DNA_FULL, color="r", ls=":")
    ax.set_xlim(-40, 230)
    ax.set_ylim(-0.3, 3)
    ax.set_xlabel("arc length along rod (A)")
    ax.set_ylabel("tube density / proximal rod density")
    ax.legend(fontsize=7)
    ax.set_title("occupancy along rod (12 A map, 8 A tube; extended 60 A past contour tip)", fontsize=9)
    # row 2: class_3D bars
    ax = fig.add_subplot(gs[2, 0:2])
    cls = [k for k in results if k.startswith("c3d_")]
    if cls:
        w = [results[k]["n_particles"] for k in cls]
        L = [S[k]["main_L_from_surface_A"] for k in cls]
        H = [S[k]["half_occupancy_from_surface_A"] or 0 for k in cls]
        x = np.arange(len(cls))
        ax.bar(x - 0.2, L, 0.4, label="contour arc length")
        ax.bar(x + 0.2, H, 0.4, label="half-occupancy length")
        for xi, wi, nb in zip(x, w, [S[k]["n_branches_from_lobe"] for k in cls]):
            ax.text(xi, 5, f"{wi}\n{nb}br", ha="center", fontsize=6, color="w")
        if "free" in S:
            ax.axhline(S["free"]["tip_to_tip_A"], color="k", ls="--", lw=0.8, label="free tip-to-tip")
        ax.axhline(DNA_FULL, color="r", ls=":", label="170 A")
        ax.set_xticks(x)
        ax.set_xticklabels([k.replace("c3d_", "c") for k in cls])
        ax.set_ylabel("A from lobe surface")
        ax.legend(fontsize=7)
        cw = S["class3d_weighted"]
        ax.set_title(f"class_3D (J59): particle-weighted mean {cw['weighted_mean_L_from_surface_A']} A "
                     f"(SD {cw['weighted_sd_L_A']})", fontsize=9)
    # row 2 right: helical spectra
    ax = fig.add_subplot(gs[2, 2:4])
    for k in [k for k in ("nu", "relion_pp", "free", "het1") if k in results and results[k]["helix"]]:
        for n, ls in ((1, "-"), (2, "--")):
            kz, p = results[k]["helix"]["_spectra"][n]
            kz, p = np.array(kz), np.array(p)
            sel = (kz > 1 / 80) & (kz < 1 / 14)
            ax.plot(1 / kz[sel], p[sel] / p[sel].mean(), ls, label=f"{k} n={n}")
    ax.axvspan(30, 38, color="y", alpha=0.3)
    ax.set_xlabel("axial period (A)")
    ax.set_ylabel("power / mean")
    ax.legend(fontsize=6, ncol=2)
    ax.set_title("helical layer-line power along rod (6 A low-pass; yellow = 30-38 A)", fontsize=9)
    # row 3: free control & radial profiles
    ax = fig.add_subplot(gs[3, 0:2])
    for k in [k for k in ("nu", "het1", "free") if k in results and results[k]["radial_profile_lp8"]]:
        pr = results[k]["radial_profile_lp8"]["profile"]
        ax.plot(np.arange(len(pr)), np.array(pr) / pr[0], label=f"{k} FWHM {results[k]['radial_fwhm_lp8_A']:.1f} A")
    ax.set_xlabel("distance from rod axis (A)")
    ax.set_ylabel("density / axis density (8 A low-pass)")
    ax.legend(fontsize=7)
    ax.set_title("rod radial profile (diameter)", fontsize=9)
    ax = fig.add_subplot(gs[3, 2:4])
    ax.axis("off")
    txt = []
    for k in [k for k in ("nu", "het0", "het1", "free", "relion_pp") if k in S]:
        s = S[k]
        txt.append(f"{k}: lobe {s['lobe_volume_A3']} A3 {s['lobe_extent_A']}; branches {s['branch_lengths_from_surface_A']}; "
                   f"L50 {s['half_occupancy_from_surface_A']}; tip-tip {s['tip_to_tip_A']}; tunnels {list(s['tunnels_vs_contour'].values())}")
    ax.text(0, 1, "\n".join(txt), fontsize=7, va="top", family="monospace", wrap=True)
    fig.suptitle("W4 dot-on-rod 3D: attached-DNA length (scripts/46_dotrod_3d_length.py)", fontsize=12)
    fig.savefig(FIG_QC, dpi=110, bbox_inches="tight")


# ----------------------------------------------------------------------------------------------- derived summary
# Reference tip positions (A, map frame shared by J56/J55-class-0/J59/RELION) of the three rod segments seen in the
# NU map at the primary contour: S = straight rod off the lobe; V = second proximal rod ending at the vertex; F = far
# arm that leaves the vertex and ends beyond it. Branches are assigned to the nearest reference tip (<= 35 A).
REF_TIPS = {"S_straight": (211.0, 150.0, 213.0), "V_vertex": (225.0, 181.0, 237.0), "F_far_arm": (226.0, 278.0, 284.0)}


def derived(out):
    D = {}
    for k, r in out["maps"].items():
        if r["group"] == "free":
            continue
        per = {}
        for V, sc in r["scan"].items():
            seg = {}
            for b in sc["branches"]:
                if b["L_from_lobe_surface_A"] < 30:
                    continue
                t = np.array(b["tip_xyz"])
                dists = {n: float(np.linalg.norm(t - np.array(p))) for n, p in REF_TIPS.items()}
                n = min(dists, key=dists.get)
                if dists[n] <= 35 and (n not in seg or b["L_from_lobe_surface_A"] > seg[n]):
                    seg[n] = b["L_from_lobe_surface_A"]
            per[V] = seg
        S = [per[V].get("S_straight") for V in per]
        Sv = [x for x in S if x is not None]
        Fv = [V for V in per if "F_far_arm" in per[V]]
        D[k] = dict(n_particles=r["n_particles"], segments_vs_contour=per,
                    S_straight_from_surface_A_median=round(float(np.median(Sv)), 1) if Sv else None,
                    S_straight_range_A=[min(Sv), max(Sv)] if Sv else None,
                    far_arm_connected_at_contours_A3=[int(V) for V in Fv],
                    far_arm_L_from_surface_A_at_primary=per[str(PRIMARY)].get("F_far_arm"),
                    lowest_contour_far_arm_connected_A3=min([int(V) for V in Fv]) if Fv else None)
    cls = [k for k in D if k.startswith("c3d_")]
    w = np.array([D[k]["n_particles"] for k in cls], float)
    Sm = np.array([D[k]["S_straight_from_surface_A_median"] or np.nan for k in cls])
    conn65 = np.array([D[k]["far_arm_L_from_surface_A_at_primary"] is not None for k in cls])
    ok = ~np.isnan(Sm)
    D["class3d_weighted"] = dict(S_straight_weighted_mean_A=round(float(np.sum(w[ok] * Sm[ok]) / w[ok].sum()), 1),
                                 S_straight_min_max_A=[float(np.nanmin(Sm)), float(np.nanmax(Sm))],
                                 frac_particles_far_arm_connected_at_65k=round(float(w[conn65].sum() / w.sum()), 3))
    # free-duplex control: split the tip-to-tip medial path at its largest kink
    for k in ("free", "free_it0"):
        M = [m for m in MAPS if m["key"] == k][0]
        d, _ = load_map(M["path"])
        lp = normalise(lowpass(d, LP_RES))
        thr, mask, _, others = threshold_for_volume(lp, PRIMARY)
        res = analyse_contour(lp, thr, mask, others)
        P = full_path(res)
        s_ = arclen(P)
        kp = kink_profile(P, s_, 0, s_[-1])
        kx, ka = max(kp, key=lambda x: x[1])
        D[k] = dict(tip_to_tip_A=round(float(s_[-1]), 1), largest_kink_deg=round(ka, 1),
                    kink_at_arc_A=round(kx, 1), segment_lengths_A=sorted([round(kx, 1), round(float(s_[-1] - kx), 1)],
                                                                         reverse=True))
    return D


def finalize():
    out = json.load(open(JSON_OUT))
    fits = {}
    logf = f"{OUTD}/chimerax_fit_log.txt"
    if os.path.exists(logf):
        txt = re.sub(r"<[^>]+>", " ", open(logf).read())
        txt = txt.replace("&nbsp;", " ")
        for m in re.finditer(r"FIT (\S+) start(.*?)ENDFIT", txt, re.S):
            tag, body = m.group(1), m.group(2)
            cors = re.findall(r"correlation\s*=\s*([-0-9.]+)", body)
            cam = re.findall(r"correlation about mean\s*=\s*([-0-9.]+)", body)
            ov = re.findall(r"overlap\s*=\s*([-0-9.]+)", body)
            fits[tag] = dict(correlation=float(cors[-1]) if cors else None,
                             correlation_about_mean=float(cam[-1]) if cam else None,
                             overlap=float(ov[-1]) if ov else None)
    out["sfgfp_rigid_fit_validation"] = fits
    out["derived_segments"] = derived(out)
    # montages
    from PIL import Image
    made = []
    for k in RENDER_KEYS:
        ims = [f"{FIG_MAPS}/dotrod_{k}_v{i}.png" for i in (1, 2, 3)]
        if all(os.path.exists(x) for x in ims):
            I = [Image.open(x).convert("RGB") for x in ims]
            W = sum(i.width for i in I)
            H = max(i.height for i in I)
            canvas = Image.new("RGB", (W, H), "white")
            x = 0
            for i in I:
                canvas.paste(i, (x, 0))
                x += i.width
            canvas.save(f"{FIG_MAPS}/dotrod_{k}_3views.png")
            made.append(f"figures/maps/dotrod_{k}_3views.png")
    cl = [f"{FIG_MAPS}/dotrod_c3d_{i}_v1.png" for i in range(10)]
    if all(os.path.exists(x) for x in cl):
        I = [Image.open(x).convert("RGB") for x in cl]
        w, h = I[0].size
        canvas = Image.new("RGB", (5 * w, 2 * h), "white")
        for i, im in enumerate(I):
            canvas.paste(im, ((i % 5) * w, (i // 5) * h))
        canvas.save(f"{FIG_MAPS}/dotrod_class3d_montage.png")
        made.append("figures/maps/dotrod_class3d_montage.png")
    out["figures"] = made
    save_json(out)
    print(json.dumps(fits, indent=1))


# ----------------------------------------------------------------------------------------------- pose-based tests
def rodrigues(r):
    th = np.linalg.norm(r, axis=1)
    k = r / np.maximum(th, 1e-12)[:, None]
    K = np.zeros((len(r), 3, 3))
    K[:, 0, 1], K[:, 0, 2], K[:, 1, 0] = -k[:, 2], k[:, 1], k[:, 2]
    K[:, 1, 2], K[:, 2, 0], K[:, 2, 1] = -k[:, 0], -k[:, 1], k[:, 0]
    s, c = np.sin(th)[:, None, None], np.cos(th)[:, None, None]
    return np.eye(3)[None] + s * K + (1 - c) * K @ K


def class_radial_profile(img, apix, rmax=230.0, dr=4.0):
    """Annular mass profile of a 2D class average about its brightest (dot) pixel, background-subtracted."""
    sm = ndi.gaussian_filter(img, 1.0)
    n = img.shape[0]
    yy, xx = np.indices(img.shape)
    cen = (yy - n / 2) ** 2 + (xx - n / 2) ** 2 <= (90.0 / apix) ** 2
    iy, ix = np.unravel_index(np.argmax(np.where(cen, sm, -np.inf)), img.shape)
    r = np.hypot(yy - iy, xx - ix) * apix
    rr = np.hypot(yy - n / 2, xx - n / 2) * apix
    bg = np.median(img[(rr > 190) & (rr < 220)])
    bins = np.arange(0, rmax + dr, dr)
    m = np.histogram(r, bins, weights=img - bg)[0]
    return bins, m


def pose_tests():
    """Use the NU particle poses (and the A1 2D poses carried through) to ask how the 3D map projects.

    1. Validate the CryoSPARC 3D pose convention independently of rod tilt: the projected in-plane direction of a
       3D vector must track the particle's 2D in-plane angle (constant offset per 2D class).
    2. Out-of-plane tilt (= 2D foreshortening) of each 3D rod segment per particle.
    3. Re-projection check: annular mass profile about the sfGFP lobe predicted from the 3D map at each particle's
       pose, averaged per 2D class, vs the same profile measured on the real 2D class average.
    4. Composition of the 3D particle sets by A1 2D class group (dot-at-end / haze decoy / free rod / other).
    """
    out = json.load(open(JSON_OUT))
    sel = json.load(open(os.path.join(REPO, "results", "A1_length_vs_hinge.json")))["selection"]
    DOT = [x["cls"] for x in sel if x.get("group") == "DOT"]
    FREE = [x["cls"] for x in sel if x.get("group") == "FREE"]
    HAZE = [105, 12]
    a = np.load(f"{CS}/J56/J56_003_particles.cs")
    p = np.load(f"{CS}/J56/J56_passthrough_particles.cs")
    assert (a["uid"] == p["uid"]).all()
    R = rodrigues(a["alignments3D/pose"].astype(float))
    th2 = p["alignments2D/pose"].astype(float)
    cls = np.asarray(p["alignments2D/class"])
    cmap = dict(zip(p["uid"], cls))

    def comp(uids):
        c = np.array([cmap.get(u, -1) for u in uids])
        return dict(n=int(len(c)), dot=round(float(np.isin(c, DOT).mean()), 3), haze_105_12=round(float(np.isin(c, HAZE).mean()), 3),
                    free=round(float(np.isin(c, FREE).mean()), 3),
                    other=round(float((~np.isin(c, DOT + FREE + HAZE)).mean()), 3))
    composition = {"nu_J56": comp(a["uid"])}
    composition["nu_J56"]["top_2D_classes"] = [(int(k), int(v)) for k, v in
                                               sorted(zip(*np.unique(cls, return_counts=True)), key=lambda x: -x[1])[:8]]
    for i in range(10):
        composition[f"c3d_{i}"] = comp(np.load(f"{CS}/J59/J59_class_{i:02d}_00014_particles.cs")["uid"])
    # hetero classes: map through the J55 passthrough (pooled A1 particles)
    hp = np.load(f"{CS}/J55/J55_passthrough_particles_all_classes.cs")
    hmap = dict(zip(hp["uid"], hp["alignments2D/class"]))
    for i in range(4):
        u = np.load(f"{CS}/J55/J55_class_{i:02d}_00102_particles.cs")["uid"]
        c = np.array([hmap.get(x, -1) for x in u])
        composition[f"het{i}"] = dict(n=int(len(c)), dot=round(float(np.isin(c, DOT).mean()), 3),
                                      haze_105_12=round(float(np.isin(c, HAZE).mean()), 3),
                                      free=round(float(np.isin(c, FREE).mean()), 3),
                                      other=round(float((~np.isin(c, DOT + FREE + HAZE)).mean()), 3))
    # geometry of the NU map at the primary contour (recomputed for the medial paths)
    d, _ = load_map(MAPS[0]["path"])
    lp = normalise(lowpass(d, LP_RES))
    thr, mask, _, others = threshold_for_volume(lp, PRIMARY)
    res = analyse_contour(lp, thr, mask, others)
    lc = np.array(res["lobe"]["centre_xyz"])
    brs = res["branches"]
    longb = max(brs, key=lambda b: b["L_from_lobe_surface_A"])
    straight = [b for b in brs if b is not longb and b["L_from_lobe_surface_A"] > 40]
    straight = max(straight, key=lambda b: b["L_from_lobe_surface_A"]) if straight else None
    P, s_, se = longb["_P"], longb["_s"], longb["_s_exit"]
    kink_s = se + (longb["max_kink_at_A_from_surface"] or 46.0)
    vertex = P[np.searchsorted(s_, kink_s).clip(0, len(P) - 1)]
    vecs = {"lobe_to_long_tip": P[-1] - lc, "lobe_to_vertex": vertex - lc, "vertex_to_long_tip": P[-1] - vertex}
    if straight is not None:
        vecs["lobe_to_straight_tip"] = straight["_P"][-1] - lc
    # 1. convention test
    dotm = np.isin(cls, DOT)
    conv = {}
    for order in ("xyz", "zyx"):
        for tr in (False, True):
            for sgn in (1, -1):
                RR = np.transpose(R, (0, 2, 1)) if tr else R
                vals = {}
                for k, v in vecs.items():
                    vv = v if order == "xyz" else v[::-1]
                    q = RR @ vv
                    phi = np.arctan2(q[:, 1], q[:, 0])
                    rb = []
                    for c in DOT:
                        m = cls == c
                        if m.sum() >= 200:
                            rb.append((m.sum(), abs(np.exp(1j * (phi[m] - sgn * th2[m])).mean())))
                    w = np.array(rb)
                    vals[k] = round(float(np.sum(w[:, 0] * w[:, 1]) / w[:, 0].sum()), 3)
                conv[f"{order}_{'RT' if tr else 'R'}_{'+' if sgn > 0 else '-'}"] = vals
    best = max(conv, key=lambda k: max(conv[k].values()))
    RR = np.transpose(R, (0, 2, 1)) if "_RT_" in best else R
    flip = best.startswith("zyx")
    # 2. tilt / foreshortening per segment
    tilt = {}
    for k, v in vecs.items():
        vv = v[::-1] if flip else v
        q = RR @ vv
        t = np.degrees(np.arcsin(np.clip(np.abs(q[:, 2]) / np.linalg.norm(vv), 0, 1)))
        pl = np.linalg.norm(q[:, :2], axis=1)
        tilt[k] = dict(length_3D_A=round(float(np.linalg.norm(v)), 1),
                       tilt_deg_median_dot=round(float(np.median(t[dotm])), 1),
                       tilt_deg_p10_p90_dot=[round(float(x), 1) for x in np.percentile(t[dotm], [10, 90])],
                       tilt_deg_median_haze=round(float(np.median(t[np.isin(cls, HAZE)])), 1),
                       projected_length_median_dot_A=round(float(np.median(pl[dotm])), 1),
                       projected_length_p10_p90_dot_A=[round(float(x), 1) for x in np.percentile(pl[dotm], [10, 90])],
                       mean_cos_dot=round(float(np.mean(pl[dotm]) / np.linalg.norm(v)), 3))
    # 3. re-projection radial mass profiles about the lobe, per 2D class, vs the real class averages
    big_t, big_m, _, _ = threshold_for_volume(lp, 120000)
    keep = ndi.binary_dilation(big_m, S26, iterations=2)
    X = np.argwhere(keep)[:, ::-1] * APIX
    W = np.clip(lp[keep], 0, None)
    if flip:
        X, lcv = X[:, ::-1], lc[::-1]
    else:
        lcv = lc
    rng = np.random.default_rng(46)
    with mrcfile.open(f"{CS}/J45/J45_030_class_averages.mrc", permissive=True) as m:
        cav = np.asarray(m.data, float)
        capix = float(m.voxel_size.x)
    reproj = {}
    img_store = {}
    n2 = cav.shape[1]
    c2 = n2 / 2
    yy, xx = np.indices((n2, n2))
    for c in DOT[:10] + HAZE:
        idx = np.where(cls == c)[0]
        if len(idx) < 200:
            continue
        idx = rng.choice(idx, min(len(idx), 1500), replace=False)
        # predicted class image: project the masked 12 A NU map at each particle's 3D pose, centre on the projected
        # lobe, rotate by the particle's 2D in-plane angle (validated sign) -> common class frame up to one rotation
        acc = np.zeros((n2, n2))
        for i in idx:
            q = (X - lcv) @ RR[i].T
            ang = th2[i]  # class-frame angle = particle-frame angle + theta2D (+ per-class constant)
            ca, sa = np.cos(ang), np.sin(ang)
            u = ca * q[:, 0] - sa * q[:, 1]
            v = sa * q[:, 0] + ca * q[:, 1]
            acc += np.histogram2d(v / capix + c2, u / capix + c2, bins=n2, range=[[0, n2], [0, n2]], weights=W)[0]
        acc = ndi.gaussian_filter(acc / len(idx), 1.0)
        # real class average, shifted so the dot (brightest pixel near centre) sits at the centre
        real = cav[c].copy()
        sm = ndi.gaussian_filter(real, 1.0)
        rr0 = np.hypot(yy - c2, xx - c2) * capix
        iy, ix = np.unravel_index(np.argmax(np.where(rr0 < 90, sm, -np.inf)), real.shape)
        real = ndi.shift(real, (c2 - iy, c2 - ix), order=1, mode="constant")
        rr = np.hypot(yy - c2, xx - c2) * capix
        disc = rr < 180
        real = real - np.median(real[(rr > 180) & (rr < 220)])
        best_cc, best_fit = -2, None
        for mir in (False, True):
            A = acc[:, ::-1] if mir else acc
            for deg in range(0, 360, 2):
                Bm = ndi.rotate(A, deg, reshape=False, order=1)
                cc = np.corrcoef(Bm[disc], real[disc])[0, 1]
                if cc > best_cc:
                    best_cc, best_fit = cc, (mir, deg, Bm)
        pred = best_fit[2]
        # where does the re-projection put mass that the real average lacks? far-zone masks from the prediction
        far = disc & (rr > 95)
        pm = far & (pred > 0.3 * pred[far].max())
        near = (rr > 25) & (rr < 80)
        pnear = near & (pred > 0.3 * pred[near].max())

        def frac(img, m1, m2):
            return float(img[m1].sum() / img[m2].sum())
        # null: same far mask rotated about the dot by random angles, applied to the real average
        nulls = []
        for deg in range(20, 360, 20):
            mrot = ndi.rotate(pm.astype(float), deg, reshape=False, order=0) > 0.5
            nulls.append(frac(real, mrot & far, pnear))
        reproj[str(c)] = dict(n_particles_nu=int((cls == c).sum()), n_used=int(len(idx)), fit_cc=round(float(best_cc), 3),
                              fit_mirror=bool(best_fit[0]), fit_rot_deg=int(best_fit[1]),
                              predicted_far_arm_over_near_rod=round(frac(pred, pm, pnear), 3),
                              real_far_arm_over_near_rod=round(frac(real, pm, pnear), 3),
                              real_null_rotated_masks_mean_sd=[round(float(np.mean(nulls)), 3), round(float(np.std(nulls)), 3)])
        img_store[c] = (pred, real)
    out["pose_tests"] = dict(
        note="NU J56 particle poses (54,158 particles); 2D poses/classes from the A1 2D (J45) passthrough",
        nu_geometry_vectors={k: [round(float(x), 1) for x in v] for k, v in vecs.items()},
        vertex_xyz=[round(float(x), 1) for x in vertex],
        convention_test_weighted_Rbar=conv, convention_chosen=best, tilt=tilt, reprojection_radial=reproj,
        composition=composition)
    save_json(out)
    # figure: predicted (re-projected 3D) vs real 2D class averages, plus foreshortening bars
    show = [c for c in img_store][:12]
    fig, ax = plt.subplots(3, len(show), figsize=(1.6 * len(show), 6.2))
    lim = (c2 - 200 / capix, c2 + 200 / capix)
    for j, c in enumerate(show):
        pred, real = img_store[c]
        ax[0, j].imshow(pred, cmap="gray", origin="lower")
        ax[1, j].imshow(real, cmap="gray", origin="lower")
        for a_ in (ax[0, j], ax[1, j]):
            a_.set_xlim(*lim)
            a_.set_ylim(*lim)
            a_.set_xticks([])
            a_.set_yticks([])
        r_ = reproj[str(c)]
        ax[0, j].set_title(f"{c} ({r_['n_particles_nu']})\ncc {r_['fit_cc']}", fontsize=7)
        ax[1, j].set_title(f"far arm {r_['predicted_far_arm_over_near_rod']:.2f} vs\n{r_['real_far_arm_over_near_rod']:.2f}",
                           fontsize=7)
        ax[2, j].axis("off")
    ax[0, 0].set_ylabel("re-projected\nNU map", fontsize=8)
    ax[1, 0].set_ylabel("real 2D\nclass avg", fontsize=8)
    axb = fig.add_axes([0.08, 0.04, 0.84, 0.24])
    ks = list(tilt)
    axb.bar(np.arange(len(ks)), [tilt[k]["length_3D_A"] for k in ks], 0.4, label="3D length")
    axb.bar(np.arange(len(ks)) + 0.4, [tilt[k]["projected_length_median_dot_A"] for k in ks], 0.4,
            label="median 2D-projected length (dot-class particles)")
    axb.set_xticks(np.arange(len(ks)) + 0.2)
    axb.set_xticklabels([f"{k.replace('_', ' ')} (tilt {tilt[k]['tilt_deg_median_dot']} deg)" for k in ks], fontsize=7)
    axb.set_ylabel("A")
    axb.legend(fontsize=7)
    fig.suptitle(f"NU map re-projected at each particle's pose ({best}) vs the A1 2D class averages (cols = 2D class)",
                 fontsize=9)
    fig.savefig(FIG_QC.replace(".png", "_poses.png"), dpi=120)
    print("pose tests:", best, {k: v["tilt_deg_median_dot"] for k, v in tilt.items()}, flush=True)


if __name__ == "__main__":
    {"analyse": analyse, "poses": pose_tests, "finalize": finalize}[sys.argv[1] if len(sys.argv) > 1 else "analyse"]()
