#!/usr/bin/env python
"""Assemble the Fig S5-style figures (Day et al. 2025 Structure, Fig S5) for the two maps.

  figures/final/figS5_W3_duplex.png (+ .json)   W3 DNA-control map, NU J39
  figures/final/figS5_W4_dotrod.png (+ .json)   W4 dot-on-rod map, NU J56

Panels: (A) RELION loose-mask postprocessed map coloured by RELION local resolution, three orthogonal views at one
scale; (B) the 32 2D classes contributing most particles to the final refinement; (C) gold-standard FSC curves
(CryoSPARC Validation J40/J57; RELION 3.1 phase-randomisation-corrected; map-model FSC where a defensible rigid fit
exists); (D) viewing-direction distribution of the final particles; plus a Table S1-style number box.
Inputs are written by scripts/50_figS5_prep.py (+ 51/53/54 for renders, fit, mtriage); every plotted number is
listed with its source file in the JSON next to each figure.
Run:  $DNASPRING_PY -I scripts/55_figS5_assemble.py
"""
import json
import os
import re
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import matplotlib.image as mpimg  # noqa: E402
from matplotlib.colors import LinearSegmentedColormap, LogNorm  # noqa: E402
from matplotlib.patches import Rectangle  # noqa: E402

sys.path.insert(0, "/n/home05/ianpoon/DNA-springs-project-2026-27/cryoEM/pp2a_day2025_reproduction/scripts")
import vizstyle  # noqa: E402
vizstyle.apply(matplotlib)

REPO = "/n/home05/ianpoon/DNA-springs-project-2026-27/cryoEM/sfgfp_chimera_screening"
SCRATCH = "/n/netscratch/hekstra_lab/Lab/ian_poon/sfgfp_chimera_screening_20261112/figures_scratch"
FIG = os.path.join(REPO, "figures", "final")
CS = "/n/netscratch/hekstra_lab/Lab/ian_poon/cryosparc_projects/CS-cs-sfgfp-chimera-screening"
S = vizstyle.SERIES
REL = "/n/netscratch/hekstra_lab/Lab/ian_poon/sfgfp_chimera_screening_20261112/relion"


def title(tag):
    """Figure title; particle composition generated from classes2d.json (final particles by 2D label).
    Attachment of the DNA to sfGFP is NOT shown by either map (report §9), so neither title claims it."""
    info = json.load(open(os.path.join(SCRATCH, tag, "classes2d.json")))
    by, tot = info["final_particles_by_label"], info["n_final_particles"]
    pc = lambda n: f"{100 * n / tot:.0f}%"  # noqa: E731
    if tag == "W3_duplex":
        dec = tot - by["rod_no_dot"]
        return ("W3 rod_no_dot (“free DNA”) control: rod + terminal lobe; NU J39\n"
                f"final particles {pc(by['rod_no_dot'])} rod_no_dot, {pc(dec)} decoy-input "
                f"(noise {by['noise_junk']:,} / dot pairs {by['dot_pair']:,} / ice {by['ice_contaminant']:,})")
    return (f"W4 dot-on-rod map (sfGFP-associated DNA; final particles {pc(by['DOT (dot-on-rod)'])} DOT / "
            f"{pc(by['haze (decoy)'])} haze; NU J56)")


def bfac_loose(tag):
    st = os.path.join(REL, tag, "pp_loose", "postprocess.star")
    m = re.search(r"_rlnBfactorUsedForSharpening\s+([-0-9.eE+]+)", open(st).read())
    return float(m.group(1)), st
LABEL_COLOUR = {"rod_no_dot": S[0], "DOT (dot-on-rod)": S[0], "FREE (free rod)": S[2], "dot_pair": S[4],
                "noise_junk": S[3], "haze (decoy)": S[3], "ice_contaminant": S[6], "ice (decoy)": S[6],
                "dot_long_rod": S[1], "dot_pair_thread": S[5], "unclear": vizstyle.INK_MUTED}
LABEL_TEXT = {"rod_no_dot": "rod, no dot (rod_no_dot target)", "dot_pair": "dot pairs (decoy input)",
              "noise_junk": "noise/junk (decoy input)", "ice_contaminant": "ice (decoy input)"}


def letter(fig, x, y, s):
    fig.text(x, y, s, fontsize=14, fontweight="bold", va="top", ha="left", color=vizstyle.INK)


def crop_box(img, bg, tol=4 / 255):
    d = np.abs(img[..., :3] - np.array(bg)[None, None, :]).max(-1) > tol
    ys, xs = np.where(d)
    return ys.min(), ys.max(), xs.min(), xs.max()


def panel_a(fig, tag, rect, src):
    r = json.load(open(os.path.join(SCRATCH, tag, "render.json")))
    rl = json.load(open(os.path.join(SCRATCH, tag, f"{tag}_render_log.json")))
    a_per_px = rl["field_width_A"] / rl["window_px"][0]
    bg = matplotlib.colors.to_rgb(vizstyle.SURFACE)
    imgs = [mpimg.imread(os.path.join(SCRATCH, tag, f"{tag}_v{i}.png")) for i in (1, 2, 3)]
    boxes = [crop_box(im, bg) for im in imgs]
    h = max(b[1] - b[0] for b in boxes) + 40
    w = max(b[3] - b[2] for b in boxes) + 40
    crops = []
    for im, (y0, y1, x0, x1) in zip(imgs, boxes):
        cy, cx = (y0 + y1) // 2, (x0 + x1) // 2
        pad = np.ones((im.shape[0] + 2 * h, im.shape[1] + 2 * w, im.shape[2]), im.dtype)
        pad[..., :3] = bg
        pad[h:h + im.shape[0], w:w + im.shape[1]] = im
        crops.append(pad[cy + h - h // 2: cy + h + h // 2, cx + w - w // 2: cx + w + w // 2])
    x, y, W, H = rect
    cbw = 0.012
    gap = 0.012
    pw = (W - cbw - 0.06 - 2 * gap) / 3
    names = ["view 1: down the thinnest axis", "view 2: turned 90° about x", "view 3: turned 90° about y (end-on)"]
    for i, (c, nm) in enumerate(zip(crops, names)):
        ax = fig.add_axes([x + i * (pw + gap), y, pw, H])
        ax.imshow(c, interpolation="lanczos")
        ax.set_xticks([]); ax.set_yticks([])
        for s in ax.spines.values():
            s.set_visible(False)
        ax.set_title(nm, fontsize=8.5, color=vizstyle.INK_2, loc="center", pad=3)
        if i == 1:   # scale bar in the empty lower-left corner of view 2
            L = 50.0 / a_per_px
            yb = c.shape[0] - 0.04 * c.shape[0]
            ax.plot([25, 25 + L], [yb, yb], color=vizstyle.INK, lw=2.5, solid_capstyle="butt")
            ax.text(25 + L / 2, yb - 0.035 * c.shape[0], "50 Å", ha="center", va="bottom", fontsize=8, color=vizstyle.INK)
    # colour bar from the exact ChimeraX palette
    stops = [(float(t.split(",")[0]), t.split(",")[1]) for t in r["palette"].split(":")]
    lo, hi = stops[0][0], stops[-1][0]
    cmap = LinearSegmentedColormap.from_list("locres", [((v - lo) / (hi - lo), c) for v, c in stops])
    cax = fig.add_axes([x + W - cbw - 0.045, y + 0.12 * H, cbw, 0.76 * H])
    cb = matplotlib.colorbar.ColorbarBase(cax, cmap=cmap, norm=matplotlib.colors.Normalize(lo, hi), orientation="vertical")
    cb.set_ticks([v for v, _ in stops])
    cb.set_label("RELION local resolution (Å)", fontsize=8.5)
    cb.outline.set_visible(False)
    B, bsrc = bfac_loose(tag)
    why = (", the same enclosed volume as W4 for comparability (not tuned for a single rod)" if tag == "W3_duplex"
           else " (sfGFP 33,000 + one 50 bp duplex 32,000 Å³)")
    cap = (f"Contour: 65,000 Å³ enclosed{why}; level {r['threshold']:.3f} ({r['threshold_sigma']:.1f} SD) on the "
           f"RELION loose-mask postprocessed map, B\u00a0=\u00a0{B:.0f}\u00a0Å²; surface coloured by RELION local resolution "
           f"(colour-bar range = display only; reported range in the summary box).")
    import textwrap
    fig.text(x, y - 0.002, "\n".join(textwrap.wrap(cap, 175)), fontsize=8.2, color=vizstyle.INK_2, va="top", ha="left",
             linespacing=1.35)
    src["A"] = {
        "caption": cap, "sharpening_B_A2": B, "sharpening_B_source": bsrc + " _rlnBfactorUsedForSharpening",
        "map": r["map"], "colouring": r["locres"],
        "contour": {"threshold": r["threshold"], "threshold_sd": r["threshold_sigma"],
                    "definition": f"{r['threshold_definition']} (65,000 A^3 = project primary contour, results/dotrod_3d_length.json "
                                  "parameters.primary_A3)", "surface_level_in_chimerax": rl["surface_level"],
                    "dust_hidden_below_A": r["dust_hidden_below_A"]},
        "palette": r["palette"], "palette_range_A": r["palette_range_A"],
        "palette_range_source": os.path.join(REPO, "results", "final", "local_resolution_stats.json") + " palette_range_A",
        "scale": {"A_per_px": a_per_px, "field_width_A": rl["field_width_A"], "scale_bar_A": 50,
                  "source": os.path.join(SCRATCH, tag, f"{tag}_render_log.json")},
        "renders": [os.path.join(SCRATCH, tag, f"{tag}_v{i}.png") for i in (1, 2, 3)],
        "render_script": os.path.join(REPO, "scripts", f"51_figS5_render_{tag}.py"),
        "views": "camera set on the principal axes of the contoured map (render.json axes); v2 = turn x 90, v3 = turn y 90 from v1; 'view' re-frames after every turn, then the orthographic field width is fixed so all three views share one scale"}


def panel_b(fig, tag, rect, src):
    info = json.load(open(os.path.join(SCRATCH, tag, "classes2d.json")))
    z = np.load(os.path.join(SCRATCH, tag, "classes2d.npz"))
    tiles = z["tiles"]
    x, y, W, H = rect
    nr, nc = 4, 8
    leg_h = 0.085
    tw = W / nc
    th = (H - leg_h) / nr
    side = min(tw * fig.get_figwidth(), th * fig.get_figheight() * 0.82)   # inches, square tile
    for i, t in enumerate(info["top_classes"][: nr * nc]):
        r_, c_ = divmod(i, nc)
        ax = fig.add_axes([x + c_ * tw + (tw - side / fig.get_figwidth()) / 2,
                           y + H - (r_ + 1) * th, side / fig.get_figwidth(), side / fig.get_figheight()])
        img = tiles[i]
        lo, hi = np.percentile(img, [1, 99.5])
        ax.imshow(img, cmap="gray", vmin=lo, vmax=hi, interpolation="lanczos")
        ax.set_xticks([]); ax.set_yticks([])
        col = LABEL_COLOUR.get(t["label"], vizstyle.INK_MUTED)
        for s in ax.spines.values():
            s.set_visible(True); s.set_color(col); s.set_linewidth(2.0)
        ax.set_title(f"{t['class']}: {t['n_final']:,}", fontsize=6.5, pad=1.5, loc="center", color=vizstyle.INK_2)
        if i == 0:
            L = 100.0 / info["apix_2d"]
            n = img.shape[0]
            ax.plot([n * 0.06, n * 0.06 + L], [n * 0.92] * 2, color="w", lw=2, solid_capstyle="butt")
            ax.text(n * 0.06, n * 0.86, "100 Å", color="w", fontsize=6)
    tot = info["n_final_particles"]
    parts = sorted(info["final_particles_by_label"].items(), key=lambda kv: -kv[1])
    for j, (lab, n) in enumerate(parts):
        col = LABEL_COLOUR.get(lab, vizstyle.INK_MUTED)
        lx = x + (j % 2) * W / 2
        ly = y + 0.03 - (j // 2) * 0.02
        fig.patches.append(Rectangle((lx, ly - 0.005), 0.009, 0.009 * fig.get_figwidth() / fig.get_figheight(),
                                     transform=fig.transFigure, facecolor="none", edgecolor=col, lw=2))
        fig.text(lx + 0.014, ly, f"{LABEL_TEXT.get(lab, lab)}: {n:,} final particles ({100 * n / tot:.0f}%)",
                 fontsize=7.6, va="center", color=vizstyle.INK_2)
    stream = "J27 rod-stream round-1" if tag == "W3_duplex" else "J45 A1 window-free"
    fig.text(x, y + leg_h - 0.006,
             f"Top {len(info['top_classes'])} of {info['n_contributing_classes']} contributing {stream} 2D classes "
             f"({100 * info['top_fraction_of_final']:.1f}% of the {tot:,} final particles).\n"
             f"Tile title = class: final particles in it. Box {info['box_A']:.0f} \u00c5 at {info['apix_2d']:.2f} \u00c5/px. "
             f"Outline = 2D label of the class:", fontsize=7.6, color=vizstyle.INK_2, va="top", linespacing=1.4)
    src["B"] = {k: info[k] for k in ("stack", "apix_2d", "box_px", "box_A", "n_final_particles", "n_contributing_classes",
                                     "top_fraction_of_final", "final_particles_by_label", "sources")}
    src["B"]["tiles"] = [{k: t[k] for k in ("class", "n_final", "n_in_2d_class", "label")} for t in info["top_classes"]]
    src["B"]["display"] = "per-tile grey scale 1st-99.5th percentile; 100 A scale bar on the first tile"


def read_two_col(p):
    a = np.loadtxt(p)
    return a[:, 0], a[:, 1]


def panel_c(ax, tag, src, fsc, mm):
    styles = {"fsc_nomask": (vizstyle.INK_MUTED, ":", 1.6), "fsc_sphericalmask": (vizstyle.INK_2, "-.", 1.4),
              "fsc_resmask": (S[0], "-", 2.0), "relion_pp_csmask": (S[1], "--", 1.6), "relion_pp_loose": (S[2], "-", 2.0)}
    src["C"] = {"curves": {}, "note": fsc["note_relion"]}
    for k, (col, ls, lw) in styles.items():
        cv = fsc["curves"][k]
        if k.startswith("relion"):
            res = cv["relion_rlnFinalResolution_A"]
            lab = cv["label"].replace(" (phase-randomisation corrected)", ", phase-rand. corrected") + f"  {res:.1f} Å"
        else:
            res = cv["res_0143_A"]
            lab = cv["label"].replace("CryoSPARC GSFSC, ", "CryoSPARC GSFSC ") + f"  {res:.1f} Å"
        ax.plot(cv["x_invA"], cv["fsc"], color=col, ls=ls, lw=lw, label=lab)
        src["C"]["curves"][k] = {"label": lab, "res_at_0.143_A_annotated": res, "res_0143_interpolated_A": cv["res_0143_A"],
                                 "res_05_interpolated_A": cv["res_05_A"], "source": cv["source"]}
    if mm is not None:
        ax.plot(mm["x"], mm["y"], color=S[6], lw=2.0, ls="-",
                label=f"map vs rigid 32 bp B-DNA rod (mtriage, model mask)  {mm['d05']:.1f} Å at 0.5")
        src["C"]["map_model"] = {k: v for k, v in mm.items() if k not in ("x", "y")}
    for thr in (0.5, 0.143):
        ax.axhline(thr, color=vizstyle.INK_MUTED, lw=0.9, ls="--")
        ax.text(0.1745, thr + 0.012, f"{thr}", fontsize=7.5, color=vizstyle.INK_MUTED, ha="right", va="bottom")
    ax.set_xlim(0, 1 / (2 * 2.84))
    ax.set_ylim(-0.05, 1.03)
    ax.set_xlabel("spatial frequency (1/Å)")
    ax.set_ylabel("Fourier shell correlation")
    sec = ax.secondary_xaxis("top", functions=(lambda f: f, lambda f: f))
    ticks = [40, 20, 12, 9, 7, 5.68]
    sec.set_ticks([1 / t for t in ticks])
    sec.set_xticklabels([f"{t:g}" for t in ticks])
    sec.set_xlabel("resolution (Å)", fontsize=8.5)
    ax.grid(True, axis="y")
    ax.legend(loc="upper left", bbox_to_anchor=(-0.02, -0.2), fontsize=7.2, handlelength=2.6, ncol=1, borderaxespad=0.0)


def panel_d(fig, ax, tag, src, fsc):
    d = np.load(os.path.join(SCRATCH, tag, "viewdirs.npz"))
    vj = json.load(open(os.path.join(SCRATCH, tag, "viewdirs.json")))
    az, el = np.degrees(d["az"]), np.degrees(d["el"])
    H, xe, ye = np.histogram2d(az, el, bins=[72, 36], range=[[-180, 180], [-90, 90]])
    Hm = np.ma.masked_equal(H.T, 0)
    im = ax.pcolormesh(xe, ye, Hm, cmap="viridis", norm=LogNorm(vmin=1, vmax=H.max()), rasterized=True)
    ax.set_facecolor("#f0efe9")
    ax.set_xlim(-180, 180); ax.set_ylim(-90, 90)
    ax.set_xticks(np.arange(-180, 181, 90)); ax.set_yticks(np.arange(-90, 91, 45))
    ax.set_xlabel("azimuth (°)"); ax.set_ylabel("elevation (°)")
    cb = fig.colorbar(im, ax=ax, pad=0.02, fraction=0.05)
    cb.set_label("particles per 5°×5° bin", fontsize=8)
    cb.outline.set_visible(False)
    o = fsc["orientation"]
    ax.text(0.015, 0.97, f"{vj['n']:,} particles   cFAR {o['cFAR']:.2f}   SCF* {o['SCF_star']:.2f}",
            transform=ax.transAxes, fontsize=7.5, va="top", ha="left", color=vizstyle.INK,
            bbox=dict(boxstyle="round,pad=0.25", fc=vizstyle.SURFACE, ec="none", alpha=0.85))
    src["D"] = {"particles": vj["source"], "n_particles": vj["n"], "n_rejected_by_refinement": vj["n_rejected"],
                "convention": ("R = Rodrigues(alignments3D/pose) (scipy Rotation.from_rotvec); viewing direction d = R @ z_hat; "
                               "azimuth = atan2(d_y, d_x), elevation = arcsin(d_z); counts per 5x5 deg bin, not area-normalised; "
                               "white = empty bin. Checked against CryoSPARC's own viewing-direction plot of the same job "
                               "(iteration 3, from the job's gridfs_data): same peaks; the transpose convention does not match."),
                "bins": [72, 36], "max_bin_count": float(H.max()), "empty_bin_fraction": float((H == 0).mean()),
                "cFAR": o["cFAR"], "SCF_star": o["SCF_star"], "cFAR_SCF_source": o["source"]}


def parse_mtriage(log):
    s = open(log).read()
    m = re.search(r"FSC\(map,model map\)=0\.5\s*:\s*([\d.]+)\s+([\d.]+)", s)
    return (float(m.group(1)), float(m.group(2))) if m else (None, None)


def mapmodel(tag):
    """W3 only. Decision rule, fixed before looking at the FSC values:
    include the map-model FSC if (1) the global fitmap search converges on one rod (top-15 fits: centres within 5 A
    of the best, axis lines within 5 deg) and (2) d(FSC=0.5, model mask) changes by <= 1 A between the best fit and
    its register variants (spin 90, spin 180, polarity flip), i.e. the number does not depend on the helical register,
    which the fit does not determine. Otherwise omit it and record why."""
    if tag != "W3_duplex":
        return None
    d = os.path.join(SCRATCH, tag, "fit")
    fj = os.path.join(d, "fit_W3_dna.json")
    if not os.path.exists(fj):
        return None
    fit = json.load(open(fj))
    best_c = np.array(fit["best_centre_A"]); best_a = np.array(fit["best_axis"])
    dc = [float(np.linalg.norm(np.array(t["centre_A"]) - best_c)) for t in fit["search_top"]]
    da = [float(np.degrees(np.arccos(min(1.0, abs(float(np.dot(t["axis"], best_a))))))) for t in fit["search_top"]]
    converged = max(dc) <= 5.0 and max(da) <= 5.0
    d05 = {}
    for v in ("best", "spin090", "spin180", "flip"):
        lg = os.path.join(d, f"mtriage_{v}", "mtriage.log")
        if os.path.exists(lg):
            d05[v] = parse_mtriage(lg)
    if "best" not in d05 or d05["best"][0] is None:
        return None
    masked = [x[0] for x in d05.values() if x[0] is not None]
    spread = max(masked) - min(masked)
    spin = fit["spin_test"]
    dec = {"rule": mapmodel.__doc__.split("Decision rule, fixed before looking at the FSC values:")[1].strip(),
           "search_converged": converged, "max_centre_offset_top15_A": max(dc), "max_axis_angle_top15_deg": max(da),
           "n_opposite_polarity_in_top15": int(sum(np.dot(t["axis"], best_a) < 0 for t in fit["search_top"])),
           "best_correlation": fit["best_correlation"], "flip_correlation_after_local_refit": fit["flip_corr_after_local_refit"],
           "spin_corr_no_refit_range": [min(t["corr_no_refit"] for t in spin), max(t["corr_no_refit"] for t in spin)],
           "d_fsc05_model_mask_by_variant_A": {k: v[0] for k, v in d05.items()},
           "d_fsc05_unmasked_by_variant_A": {k: v[1] for k, v in d05.items()},
           "d05_spread_A": spread, "use_map_model_fsc": bool(converged and spread <= 1.0),
           "sources": {"fit": fj, "variants": os.path.join(d, "register_variants.json"),
                       "mtriage_logs": [os.path.join(d, f"mtriage_{v}", "mtriage.log") for v in d05]}}
    if not dec["use_map_model_fsc"]:
        dec["reason"] = ("Omitted: " + ("the global fit did not converge on one rod placement. " if not converged else "")
                         + (f"d(FSC=0.5) changes by {spread:.1f} A between register variants, so the value depends on the "
                            "helical register (groove phase/polarity), which the 8-9 A map does not determine." if spread > 1.0 else ""))
    json.dump(dec, open(os.path.join(d, "fit_decision.json"), "w"), indent=1)
    if not dec["use_map_model_fsc"]:
        return None
    cur = os.path.join(d, "mtriage_best", "fsc_model.masked.mtriage.log")
    x, y = read_two_col(cur)
    return {"x": x, "y": y, "d05": d05["best"][0], "d05_unmasked": d05["best"][1], "curve_source": cur,
            "log": os.path.join(d, "mtriage_best", "mtriage.log"), "decision": dec,
            "model": os.path.join(d, "dna_32bp_fit.pdb"),
            "model_note": ("straight B-DNA 32 bp (model/L6_loop50/dna/straight_50bp.pdb chain A 10-41 / B -41..-10), rigid "
                           "global fitmap into relion/W3_duplex/raw.mrc at 8.7 A; only the rod placement is determined: "
                           "groove phase and polarity are not (see decision)")}


def numbers_box(fig, rect, tag, src, fsc, mm):
    lr = json.load(open(os.path.join(REPO, "results", "final", "local_resolution_stats.json")))[tag]
    info = json.load(open(os.path.join(SCRATCH, tag, "classes2d.json")))
    c = fsc["curves"]
    pr = lr["PRIMARY_inside_cryosparc_validation_mask"]
    nu = {"W3_duplex": ("J39", "J40", "J41"), "W4_dotrod": ("J56", "J57", "J58")}[tag]
    rows = [
        ("Final particles", f"{json.load(open(os.path.join(SCRATCH, tag, 'viewdirs.json')))['n']:,}"),
        ("Pixel / box", f"2.84 Å / {128 if tag == 'W3_duplex' else 160} px"),
        ("Symmetry", "C1"),
        ("GSFSC 0.143, CryoSPARC", f"{c['fsc_nomask']['res_0143_A']:.1f} / {c['fsc_sphericalmask']['res_0143_A']:.1f} / {c['fsc_resmask']['res_0143_A']:.1f} Å"),
        ("", "no mask / spherical / tight"),
        ("RELION corrected 0.143", f"{c['relion_pp_csmask']['relion_rlnFinalResolution_A']:.1f} / {c['relion_pp_loose']['relion_rlnFinalResolution_A']:.1f} Å"),
        ("", "CryoSPARC mask / loose mask"),
        ("Local resolution range", f"{pr['p5']:.1f}–{pr['p95']:.1f} Å (median {pr['median']:.1f} Å)"),
        ("", "RELION locres, p5–p95 in CryoSPARC validation (tight FSC) mask"),
        ("cFAR / SCF*", f"{fsc['orientation']['cFAR']:.2f} / {fsc['orientation']['SCF_star']:.2f}"),
        ("Map-model FSC 0.5", (f"{mm['d05']:.1f} Å (model mask)" if mm else ("omitted: depends on DNA register" if tag == "W3_duplex" else "omitted: no defensible model"))),
        ("Jobs", f"NU {nu[0]}, validation {nu[1]}, orientation {nu[2]}"),
    ]
    ax = fig.add_axes(rect)
    ax.axis("off")
    ax.set_title("Table S1-style summary", fontsize=10, loc="left", pad=4)
    yy = 0.95
    for k, v in rows:
        ax.text(0.0, yy, k, fontsize=7.8, color=vizstyle.INK_2, va="top", transform=ax.transAxes)
        ax.text(0.47, yy, v, fontsize=7.8 if k else 6.8, color=vizstyle.INK if k else vizstyle.INK_MUTED, va="top",
                transform=ax.transAxes)
        yy -= 0.083 if k else 0.065
    src["summary_box"] = {
        "local_resolution": os.path.join(REPO, "results", "final", "local_resolution_stats.json")
        + f" {tag}.PRIMARY_inside_cryosparc_validation_mask (same definition as Table S1)",
        "local_resolution_values": {k: pr[k] for k in ("p5", "median", "p95", "n_voxels", "mask")},
        "pixel_box": f"{CS}/{nu[0]}/{nu[0]}_003_particles.cs fields blob/psize_A, blob/shape",
        "other_rows": "same sources as panels C and D"}


def make(tag):
    fsc = json.load(open(os.path.join(SCRATCH, tag, "fsc.json")))
    mm = mapmodel(tag)
    src = {"figure": os.path.join(FIG, f"figS5_{tag}.png"), "script": os.path.join(REPO, "scripts", "55_figS5_assemble.py"),
           "prep_script": os.path.join(REPO, "scripts", "50_figS5_prep.py"),
           "style": "Day et al. 2025 Structure Fig S5 layout; house style vizstyle.py (pp2a_day2025_reproduction/scripts)"}
    fig = plt.figure(figsize=(14, 13.0), layout=None)
    fig.patch.set_facecolor(vizstyle.SURFACE)
    ttl = title(tag)
    fig.text(0.03, 0.990, f"Fig S5-style validation \u2014 {ttl}", fontsize=12.5, fontweight="normal", va="top",
             linespacing=1.35)
    src["title"] = ttl
    # A: three views, full width
    letter(fig, 0.02, 0.950, "A")
    panel_a(fig, tag, [0.04, 0.722, 0.94, 0.222], src)
    # B: 2D classes (left) | C: FSC (right)
    letter(fig, 0.02, 0.685, "B")
    panel_b(fig, tag, [0.04, 0.355, 0.56, 0.33], src)
    letter(fig, 0.635, 0.685, "C")
    axc = fig.add_axes([0.68, 0.49, 0.30, 0.155])
    panel_c(axc, tag, src, fsc, mm)
    # D: viewing directions (left) | Table S1-style numbers (right)
    letter(fig, 0.02, 0.325, "D")
    axd = fig.add_axes([0.07, 0.045, 0.40, 0.255])
    panel_d(fig, axd, tag, src, fsc)
    numbers_box(fig, [0.56, 0.035, 0.42, 0.27], tag, src, fsc, mm)
    if mm is None:
        src["C"]["map_model"] = None
        src["C"]["map_model_omitted_because"] = (
            json.load(open(os.path.join(SCRATCH, tag, "fit", "fit_decision.json")))["reason"] if tag == "W3_duplex"
            and os.path.exists(os.path.join(SCRATCH, tag, "fit", "fit_decision.json")) else
            ("No defensible rigid model: the W4 map is an average of heterogeneous, nearly planar sfGFP-DNA assemblies "
             "(report.md section 13b): the DNA outside the sfGFP lobe adds up to ~1.4 duplex-equivalents with an ~85 deg "
             "kink, the lobe is ~10x smaller than rigid sfGFP at 12 A (blurred or partly occupied), and the kink/attachment "
             "geometry is not decidable (cFAR 0.17). Any sfGFP + B-DNA placement would encode an unproven connectivity, so a "
             "map-model FSC would measure the model choice, not the map."))
    out = os.path.join(FIG, f"figS5_{tag}.png")
    fig.savefig(out, dpi=200, facecolor=vizstyle.SURFACE)
    plt.close(fig)
    json.dump(src, open(out.replace(".png", ".json"), "w"), indent=1, default=float)
    print("wrote", out)


if __name__ == "__main__":
    for t in (sys.argv[1:] or ["W3_duplex", "W4_dotrod"]):
        make(t)
