#!/usr/bin/env python
"""Fig S5-style panels (Day et al. 2025 Structure, Fig S5) for the two maps of this project: data preparation.

  W3  DNA-control map   : NU J39 (final iteration J39_003), validation J40, orientation J41, RELION relion/W3_duplex
  W4  dot-on-rod map    : NU J56 (final iteration J56_003), validation J57, orientation J58, RELION relion/W4_dotrod

Reads only files on disk (CryoSPARC master is down; .cs files are numpy structured arrays) and writes, per map,
into $SCRATCH/<tag>/:
  fsc.json            every FSC curve (CryoSPARC J40/J57 columns read by header name; RELION postprocess.star
                      columns read by _rln label) and the interpolated 0.5 / 0.143 crossings
  viewdirs.npz        azimuth / elevation of every final particle (convention below)
  classes2d.json/.npz the 32 J27 (W3) or J45 (W4) 2D classes contributing most particles to the final set
  render.json         contour, principal axes, ChimeraX camera, palette used by the .cxc renders
and writes the ChimeraX command files scripts/51_figS5_render_<tag>.cxc, and
results/final/local_resolution_stats.json.

Viewing-direction convention: R = Rodrigues(alignments3D/pose) (scipy Rotation.from_rotvec); viewing direction
d = R @ z_hat; azimuth = atan2(d_y, d_x), elevation = arcsin(d_z). This reproduces CryoSPARC's own "viewing
direction distribution" plot of the same refinement (J39/J56 iteration 3, extracted from the job's gridfs_data and
checked by eye; the transpose convention R^T z_hat does not).

Run:  $DNASPRING_PY -I scripts/50_figS5_prep.py          (everything)
      $DNASPRING_PY -I scripts/50_figS5_prep.py locres   (local_resolution_stats.json only)
"""
import csv
import json
import os

import mrcfile
import numpy as np
from scipy import ndimage
from scipy.spatial.transform import Rotation

REPO = "/n/home05/ianpoon/DNA-springs-project-2026-27/cryoEM/sfgfp_chimera_screening"
DATA = "/n/netscratch/hekstra_lab/Lab/ian_poon/sfgfp_chimera_screening_20261112"
CS = "/n/netscratch/hekstra_lab/Lab/ian_poon/cryosparc_projects/CS-cs-sfgfp-chimera-screening"
SCRATCH = os.path.join(DATA, "figures_scratch")
FIGFINAL = os.path.join(REPO, "figures", "final")
RESFINAL = os.path.join(REPO, "results", "final")
APIX = 2.84
ENCLOSE_A3 = 65000.0          # project-wide primary contour (results/dotrod_3d_length.json parameters.primary_A3)
DUST_A = 15.0                 # render only: hide surface blobs smaller than this (ChimeraX 'surface dust')
# Local-resolution palette shared by both maps (ColorBrewer Spectral-5; red = best). Range set from the
# envelope min/max of both maps (local_resolution_stats.json), rounded outwards to 0.5 A.
PALETTE_HEX = ["#d7191c", "#fdae61", "#ffffbf", "#abdda4", "#2b83ba"]

MAPS = {
    "W3_duplex": dict(
        label="W3 DNA-control map (NU J39)", nu="J39", it="J39_003", val="J40", orient="J41", locres_cs="J42",
        twod_key="W2_2D_rods", twod_job="J27", relion=os.path.join(DATA, "relion", "W3_duplex")),
    "W4_dotrod": dict(
        label="W4 dot-on-rod map (NU J56)", nu="J56", it="J56_003", val="J57", orient="J58", locres_cs=None,
        twod_key="A1_2D", twod_job="J45", relion=os.path.join(DATA, "relion", "W4_dotrod")),
}


def crossing(x, y, thr):
    """Resolution (A) at the first downward crossing of thr, linear interpolation in spatial frequency."""
    x, y = np.asarray(x, float), np.asarray(y, float)
    start = 1 if x[0] == 0 else 0
    for i in range(max(start, 1), len(y)):
        if y[i] < thr <= y[i - 1]:
            xc = x[i - 1] + (thr - y[i - 1]) * (x[i] - x[i - 1]) / (y[i] - y[i - 1])
            return float(1.0 / xc)
    return None


def read_cs_fsc(path):
    """CryoSPARC Validation fsc.txt: line 1 = machine column names, line 2 = human labels, then data."""
    with open(path) as f:
        names = f.readline().strip().split("\t")
        f.readline()
        arr = np.loadtxt(f)
    return {n: arr[:, i] for i, n in enumerate(names)}


def read_star_block(path, block):
    cols, rows, in_block, in_loop = [], [], False, False
    for line in open(path):
        s = line.strip()
        if s.startswith("data_"):
            if in_block:
                break
            in_block, in_loop, cols = (s == block), False, []
            continue
        if not in_block:
            continue
        if s.startswith("loop_"):
            in_loop = True
            continue
        if s.startswith("_rln"):
            cols.append(s.split()[0])
            continue
        if in_loop and s and not s.startswith("#"):
            v = s.split()
            if len(v) == len(cols):
                rows.append([float(t) for t in v])
    arr = np.array(rows)
    return {c: arr[:, i] for i, c in enumerate(cols)}


def read_star_general(path):
    out = {}
    in_gen = False
    for line in open(path):
        s = line.strip()
        if s.startswith("data_"):
            in_gen = s == "data_general"
            continue
        if in_gen and s.startswith("_rln"):
            k, v = s.split()[:2]
            try:
                out[k] = float(v)
            except ValueError:
                out[k] = v
    return out


def events_value(job, key):
    """Pull a 'key: value' line (e.g. 'cFAR: 0.26610') out of a job's events.bson (binary, strings only)."""
    raw = open(os.path.join(CS, job, "events.bson"), "rb").read().decode("latin-1")
    import re
    m = re.findall(re.escape(key) + r":\s*([0-9.]+)", raw)
    if not m:
        return None
    return max(m, key=len)          # the full-precision copy (e.g. 0.26610 rather than 0.27)


def fsc_data(tag, m):
    out = {"curves": {}, "sources": {}}
    p = os.path.join(CS, m["val"], f"{m['val']}_fsc.txt")
    cs = read_cs_fsc(p)
    for col, lab in (("fsc_nomask", "CryoSPARC GSFSC, no mask"),
                     ("fsc_sphericalmask", "CryoSPARC GSFSC, spherical mask"),
                     ("fsc_resmask", "CryoSPARC GSFSC, tight (auto-tightened) mask")):
        x, y = cs["frequency"], cs[col]
        out["curves"][col] = {"label": lab, "x_invA": x.tolist(), "fsc": y.tolist(),
                              "res_0143_A": crossing(x, y, 0.143), "res_05_A": crossing(x, y, 0.5),
                              "source": f"{p} columns 'frequency', '{col}' (read by header name)"}
    for sub, lab in (("pp_csmask", "RELION 3.1, CryoSPARC tight mask"), ("pp_loose", "RELION 3.1, loose mask")):
        sp = os.path.join(m["relion"], sub, "postprocess.star")
        f = read_star_block(sp, "data_fsc")
        g = read_star_general(sp)
        x, y = f["_rlnResolution"], f["_rlnFourierShellCorrelationCorrected"]
        out["curves"][f"relion_{sub}"] = {
            "label": lab + " (phase-randomisation corrected)", "x_invA": x.tolist(), "fsc": y.tolist(),
            "res_0143_A": crossing(x, y, 0.143), "res_05_A": crossing(x, y, 0.5),
            "relion_rlnFinalResolution_A": g.get("_rlnFinalResolution"),
            "relion_bfactor": g.get("_rlnBfactorUsedForSharpening"),
            "relion_randomise_from_A": g.get("_rlnRandomiseFrom"),
            "source": f"{sp} data_fsc _rlnResolution, _rlnFourierShellCorrelationCorrected; data_general _rlnFinalResolution"}
    out["note_relion"] = ("RELION's own _rlnFinalResolution is the last shell above 0.143 (no interpolation); "
                          "res_0143_A is the interpolated crossing used for the plot annotation")
    out["orientation"] = {
        "cFAR": float(events_value(m["orient"], "cFAR")), "SCF_star": float(events_value(m["orient"], "SCF*")),
        "source": f"{CS}/{m['orient']}/events.bson (event-log lines 'cFAR: ...' and 'SCF*: ...')"}
    json.dump(out, open(os.path.join(SCRATCH, tag, "fsc.json"), "w"), indent=1)
    return out


def view_dirs(tag, m):
    p = os.path.join(CS, m["nu"], f"{m['it']}_particles.cs")
    a = np.load(p)
    R = Rotation.from_rotvec(a["alignments3D/pose"]).as_matrix()
    d = R[:, :, 2]                       # R @ z_hat
    az = np.arctan2(d[:, 1], d[:, 0])
    el = np.arcsin(np.clip(d[:, 2], -1, 1))
    np.savez(os.path.join(SCRATCH, tag, "viewdirs.npz"), az=az, el=el)
    return {"n": int(len(a)), "source": p, "n_rejected": int(len(np.load(p.replace("particles.cs", "particles_rejected.cs"))))}


def classes2d(tag, m, ntop=32):
    exp = os.path.join(DATA, "exports")
    z = np.load(os.path.join(exp, f"{m['twod_key']}__particles.npz"))
    final = np.load(os.path.join(CS, m["nu"], f"{m['it']}_particles.cs"))["uid"]
    idx = {u: i for i, u in enumerate(z["uid"])}
    missing = sum(u not in idx for u in final)
    cls_all = z["alignments2D__class"].astype(int)
    cls_fin = cls_all[[idx[u] for u in final if u in idx]]
    stack_p = os.path.join(exp, m["twod_key"], f"{m['twod_job']}_030_class_averages.mrc")
    with mrcfile.open(stack_p, permissive=True) as mm:
        stack = np.asarray(mm.data, np.float32)
        apix2d = float(mm.voxel_size.x)
    K = stack.shape[0]
    n_fin = np.bincount(cls_fin, minlength=K)
    n_all = np.bincount(cls_all, minlength=K)
    order = [int(k) for k in np.argsort(-n_fin, kind="stable")[:ntop] if n_fin[k] > 0]
    # labels
    labels, label_src = {}, ""
    st = json.load(open(os.path.join(REPO, "results", "cs_state.json")))["jobs"]
    if tag == "W3_duplex":
        lp = os.path.join(REPO, "results", "class2d_labels_W2_2D_rods.csv")
        labels = {int(r["class"]): r["consensus_label"] for r in csv.DictReader(open(lp))}
        label_src = lp + " (blind 3-reviewer consensus label)"
        groups = {k: st[k]["classes"] for k in ("W3_duplex_particles", "W3_junk_particles", "W3_ice_particles", "W3_pairs_particles")}
    else:
        groups = {k: st[k]["classes"] for k in ("W4_dot_particles", "W4_free_particles", "W4_haze_particles", "W4_ice_particles")}
        for g, cl in groups.items():
            for c in cl:
                labels[int(c)] = {"W4_dot_particles": "DOT (dot-on-rod)", "W4_free_particles": "FREE (free rod)",
                                  "W4_haze_particles": "haze (decoy)", "W4_ice_particles": "ice (decoy)"}[g]
        label_src = os.path.join(REPO, "results", "cs_state.json") + " jobs W4_{dot,free,haze,ice}_particles.classes (A1 groups)"
    # which hetero-input group did each class feed
    group_of = {int(c): g for g, cl in groups.items() for c in cl}
    tot = int(n_fin.sum())
    by_label = {}
    for k in range(K):
        if n_fin[k]:
            by_label[labels.get(k, "unlabelled")] = by_label.get(labels.get(k, "unlabelled"), 0) + int(n_fin[k])
    info = {"stack": stack_p, "apix_2d": apix2d, "box_px": int(stack.shape[1]), "box_A": round(stack.shape[1] * apix2d, 1),
            "n_final_particles": tot, "n_final_missing_from_2d": int(missing),
            "n_contributing_classes": int((n_fin > 0).sum()),
            "top_classes": [{"class": k, "n_final": int(n_fin[k]), "n_in_2d_class": int(n_all[k]),
                             "label": labels.get(k, "unlabelled"), "hetero_input_group": group_of.get(k)} for k in order],
            "top_fraction_of_final": float(n_fin[order].sum() / tot),
            "final_particles_by_label": by_label,
            "sources": {"final_particles": os.path.join(CS, m["nu"], f"{m['it']}_particles.cs") + " (uid)",
                        "class_of_particle": os.path.join(exp, f"{m['twod_key']}__particles.npz") + " (uid, alignments2D__class)",
                        "labels": label_src}}
    json.dump(info, open(os.path.join(SCRATCH, tag, "classes2d.json"), "w"), indent=1)
    np.savez(os.path.join(SCRATCH, tag, "classes2d.npz"), tiles=stack[order], classes=np.array(order))
    return info


def envelope(vol, apix, enclose):
    n = int(round(enclose / apix ** 3))
    thr = float(np.sort(vol.ravel())[::-1][n])
    env = vol > thr
    lab, nl = ndimage.label(env)
    sizes = ndimage.sum(env, lab, range(1, nl + 1))
    return thr, env, lab, sizes


def render_setup(tag, m, lo, hi):
    """Contour, principal axes and ChimeraX camera for the locres-coloured surface."""
    mp = os.path.join(m["relion"], "pp_loose", "postprocess_masked.mrc")
    with mrcfile.open(mp, permissive=True) as mm:
        vol = np.asarray(mm.data, np.float32)
        apix = float(mm.voxel_size.x)
        org = np.array([float(mm.header.origin.x), float(mm.header.origin.y), float(mm.header.origin.z)])
    thr, env, lab, sizes = envelope(vol, apix, ENCLOSE_A3)
    main = lab == (int(np.argmax(sizes)) + 1)
    kji = np.argwhere(main)
    xyz = kji[:, ::-1] * apix + org      # mrc data are [z, y, x]
    c = xyz.mean(0)
    w, V = np.linalg.eigh(np.cov((xyz - c).T))
    e1, e2 = V[:, 2], V[:, 1]            # longest, middle
    e3 = np.cross(e1, e2)                # thinnest, right-handed
    proj = (xyz - c) @ np.stack([e1, e2, e3], 1)
    ext = proj.max(0) - proj.min(0)
    cen = c + np.stack([e1, e2, e3], 1) @ ((proj.max(0) + proj.min(0)) / 2)
    dist = 1000.0
    pos = cen + e3 * dist
    cam = [e1[0], e2[0], e3[0], pos[0], e1[1], e2[1], e3[1], pos[1], e1[2], e2[2], e3[2], pos[2]]
    field_w = float(np.ceil((ext.max() + 2 * 30.0) / 10.0) * 10.0)   # same scale in all three views, 30 A margin
    step = (hi - lo) / (len(PALETTE_HEX) - 1)
    palette = ":".join(f"{lo + i * step:.2f},{h}" for i, h in enumerate(PALETTE_HEX))
    info = {"map": mp, "locres": os.path.join(m["relion"], "relion_locres.mrc"), "apix": apix,
            "enclose_volume_A3": ENCLOSE_A3, "threshold": thr, "dust_hidden_below_A": DUST_A,
            "threshold_definition": "voxel-count threshold: the top N voxels, N = 65,000 A^3 / apix^3",
            "threshold_sigma": float(thr / vol.std()),
            "components_A3": sorted([float(s * apix ** 3) for s in sizes], reverse=True)[:5],
            "extent_along_axes_A": ext.round(1).tolist(), "axes": [e1.tolist(), e2.tolist(), e3.tolist()],
            "camera_matrix": cam, "field_width_A": field_w, "palette": palette, "palette_range_A": [lo, hi],
            "palette_hex": PALETTE_HEX}
    json.dump(info, open(os.path.join(SCRATCH, tag, "render.json"), "w"), indent=1)
    return info


def write_cxc(tag, r):
    """ChimeraX script (Python form, because 'camera ortho' fails offscreen in this build: the command checks
    opengl_context.stereo, which OffScreenRenderingContext lacks). An orthographic camera with a fixed field
    width gives the same, known A-per-pixel scale in all three views (used for the scale bar)."""
    out = os.path.join(SCRATCH, tag)
    cam = ",".join(f"{v:.6f}" for v in r["camera_matrix"])
    W, H = 1400, 1100
    fw = r["field_width_A"]
    cmds = ["set bgColor #fcfcfb", "lighting soft", "graphics silhouettes true width 1.5",
            f"open {r['map']}", f"open {r['locres']}", "volume #2 hide",
            f"volume #1 level {r['threshold']:.5f} style surface step 1",
            f"surface dust #1 size {DUST_A:.0f}",
            f"color sample #1 map #2 palette {r['palette']}",
            f"view matrix camera {cam}"]
    s = [f"# Fig S5A-style render, {tag}: RELION loose-mask postprocessed map coloured by RELION local resolution.",
         "# Written by scripts/50_figS5_prep.py; run offscreen:  $CHIMERAX --nogui --offscreen --silent <this file>",
         "from chimerax.core.commands import run",
         "from chimerax.graphics import OrthographicCamera",
         "v = session.main_view",
         "def ortho(fw):",
         "    c = OrthographicCamera(fw)",
         "    c.position = v.camera.position",
         "    v.camera.delete()",
         "    v.camera = c",
         "def frame(fw):",
         "    run(session, 'view #1')               # re-centre on the map after every turn",
         "    v.camera.field_width = fw             # then fix the scale",
         "    v.camera.redraw_needed = True"]
    s += [f"run(session, {c!r})" for c in cmds]
    s += [f"ortho({fw})", f"frame({fw})",
          f"run(session, 'save {out}/{tag}_v1.png width {W} height {H} supersample 3')",
          "run(session, 'turn x 90')", f"frame({fw})",
          f"run(session, 'save {out}/{tag}_v2.png width {W} height {H} supersample 3')",
          "run(session, 'turn x -90')", "run(session, 'turn y 90')", f"frame({fw})",
          f"run(session, 'save {out}/{tag}_v3.png width {W} height {H} supersample 3')",
          "import json",
          "from chimerax.map import Volume",
          "vol = [m for m in session.models.list() if isinstance(m, Volume)][0]",
          f"json.dump({{'surface_level': float(vol.minimum_surface_level), 'field_width_A': float(v.camera.field_width), "
          f"'window_px': [{W}, {H}]}}, open('{out}/{tag}_render_log.json', 'w'))",
          "run(session, 'exit')"]
    p = os.path.join(REPO, "scripts", f"51_figS5_render_{tag}.py")
    open(p, "w").write("\n".join(s) + "\n")
    return p


def locres_stats(tag, m):
    lp = os.path.join(m["relion"], "relion_locres.mrc")
    mp = os.path.join(m["relion"], "mask_loose.mrc")
    pp = os.path.join(m["relion"], "pp_loose", "postprocess_masked.mrc")
    lr = mrcfile.read(lp).astype(np.float32)
    mk = mrcfile.read(mp).astype(np.float32)
    vol = mrcfile.read(pp).astype(np.float32)
    thr, env, lab, sizes = envelope(vol, APIX, ENCLOSE_A3)

    def st(v):
        q = np.percentile(v, [0, 25, 50, 75, 100])
        return {"n_voxels": int(v.size), "min": float(q[0]), "p25": float(q[1]), "median": float(q[2]),
                "p75": float(q[3]), "max": float(q[4])}
    inmask = mk >= 0.5
    sampled = lr > 0
    # PRIMARY definition (the one used on every figure and in Table S1): p5 / median / p95 of RELION locres inside the
    # CryoSPARC Validation (FSC) auto-tightened mask, same rule as scripts/50_figS6_tableS1.py (mask value > 0.5,
    # locres finite and > 0)
    vp = os.path.join(CS, m["val"], f"{m['val']}_fsc_mask.mrc")
    vm = mrcfile.read(vp) > 0.5
    lv = lr[vm]
    lv = lv[np.isfinite(lv) & (lv > 0)]
    q = np.percentile(lv, [5, 50, 95])
    out = {
        "map": tag, "locres_map": os.path.realpath(lp), "loose_mask": mp,
        "locres_command": "relion_postprocess --locres --locres_sampling 25 --locres_maskrad 30, loose mask "
                          "(scripts/32_relion_postprocess.sh); RELION 3.1.4",
        "PRIMARY_inside_cryosparc_validation_mask": {
            "n_voxels": int(lv.size), "p5": float(q[0]), "median": float(q[1]), "p95": float(q[2]),
            "mask": vp,
            "definition": "PRIMARY (used on Fig S5 summary boxes and Table S1): 5th / 50th / 95th percentile of RELION "
                          "--locres values in voxels where the CryoSPARC Validation (FSC) auto-tightened mask > 0.5 "
                          "(locres finite and > 0)"},
        "SECONDARY_definitions_note": ("The entries below use other voxel sets and give other ranges; they are kept for "
                                       "reference only and are NOT quoted on the figures or in Table S1. "
                                       "inside_displayed_envelope also sets the Fig S5A colour-bar range."),
        "inside_loose_mask": dict(st(lr[inmask & sampled]),
                                  definition="voxels with loose-mask value >= 0.5 and locres > 0 (0 = not sampled by RELION)",
                                  n_unsampled_zero_excluded=int((inmask & ~sampled).sum()),
                                  n_at_50A_cap=int((inmask & (lr >= 49.99)).sum())),
        "inside_loose_mask_incl_unsampled": dict(st(lr[inmask]), definition="all voxels with loose-mask value >= 0.5"),
        "inside_displayed_envelope": dict(st(lr[env]),
                                          definition=f"voxels above the contour enclosing {ENCLOSE_A3:.0f} A^3 of the "
                                                     f"loose-mask postprocessed map {pp} (threshold {thr:.4f}); this is "
                                                     "what colours the Fig S5A surface"),
        "note": ("The loose mask (15 A low-pass, 98.5th percentile, +6 px, 8 px soft edge) contains much solvent "
                 "(W4: ~70 % of the box), so its percentiles are dominated by solvent; the envelope numbers describe "
                 "the molecule."),
    }
    return out


def main(locres_only=False):
    os.makedirs(FIGFINAL, exist_ok=True)
    os.makedirs(RESFINAL, exist_ok=True)
    stats = {}
    for tag, m in MAPS.items():
        os.makedirs(os.path.join(SCRATCH, tag), exist_ok=True)
        stats[tag] = locres_stats(tag, m)
    lo = np.floor(min(s["inside_displayed_envelope"]["min"] for s in stats.values()) * 2) / 2
    hi = np.ceil(max(s["inside_displayed_envelope"]["max"] for s in stats.values()) * 2) / 2
    stats["palette_range_A"] = [float(lo), float(hi)]
    stats["script"] = os.path.join(REPO, "scripts", "50_figS5_prep.py")
    stats["palette_range_note"] = ("Fig S5A colour-bar range = min/max of inside_displayed_envelope over both maps, "
                                   "rounded outwards to 0.5 A (display only; not a reported resolution range)")
    json.dump(stats, open(os.path.join(RESFINAL, "local_resolution_stats.json"), "w"), indent=1)
    for tag in MAPS:
        print(tag, "PRIMARY", stats[tag]["PRIMARY_inside_cryosparc_validation_mask"])
    if locres_only:
        return
    for tag, m in MAPS.items():
        f = fsc_data(tag, m)
        v = view_dirs(tag, m)
        c = classes2d(tag, m)
        r = render_setup(tag, m, float(lo), float(hi))
        p = write_cxc(tag, r)
        json.dump(v, open(os.path.join(SCRATCH, tag, "viewdirs.json"), "w"), indent=1)
        print(f"== {tag}: n_final {v['n']} (+{v['n_rejected']} rejected)  contour {r['threshold']:.4f} "
              f"({r['threshold_sigma']:.1f} sd) components {r['components_A3'][:3]} extent {r['extent_along_axes_A']} "
              f"field {r['field_width_A']}")
        for k, cv in f["curves"].items():
            print(f"   {k:20s} 0.143: {cv['res_0143_A']}  0.5: {cv['res_05_A']}  relion final: {cv.get('relion_rlnFinalResolution_A')}")
        print("   orientation", f["orientation"])
        print(f"   2D: {c['n_contributing_classes']} classes, top-32 cover {c['top_fraction_of_final']:.3f}; by label {c['final_particles_by_label']}")
        print("   locres", {k: stats[tag][k] for k in ("inside_loose_mask", "inside_displayed_envelope")})
        print("   wrote", p)


if __name__ == "__main__":
    import sys
    main(locres_only="locres" in sys.argv[1:])
