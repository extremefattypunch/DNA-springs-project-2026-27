#!/usr/bin/env python
"""Fig S6 analogue (processing workflow) and Table S1 analogue for the two final maps of the sfGFP screen.

  figures/final/figS6_workflow.png        boxes-and-arrows workflow, movies -> W3 rod_no_dot ("free DNA") control map and
                                          W4 dot-on-rod map, with counts / resolutions read from files
  figures/final/figS6_workflow.json       every number on the figure with the file it came from
  results/final/tableS1.csv / .md         Table S1 analogue, columns W3 (rod_no_dot "free DNA" control) and W4 (dot-on-rod)
  results/final/tableS1_sources.json      per-cell source (file + field + derivation)
  results/final/figS6_tableS1_numbers.json  the full number register both outputs draw on

Conventions follow Day et al. 2025 Structure (Fig S6 workflow; Table S1 rows) and the house style of
../pp2a_day2025_reproduction/scripts/27_em_figures.py + vizstyle.py.  The CryoSPARC master is not contacted:
every CryoSPARC number is read offline from the P6 job directories (job.json output counts, .cs files, FSC txt,
job.log).  Map thumbnails come from scripts/50_figS6_thumbs.cxc (run first, as a CPU job).

Run:  /n/holylabs/hekstra_lab/Everyone/ianpoon/envs/dnaspring/bin/python -I scripts/50_figS6_tableS1.py
"""
import csv
import glob
import json
import os
import re
import sys

import numpy as np
import mrcfile
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch, Polygon  # noqa: E402

REPO = "/n/home05/ianpoon/DNA-springs-project-2026-27/cryoEM/sfgfp_chimera_screening"
DATA = "/n/netscratch/hekstra_lab/Lab/ian_poon/sfgfp_chimera_screening_20261112"
P6 = "/n/netscratch/hekstra_lab/Lab/ian_poon/cryosparc_projects/CS-cs-sfgfp-chimera-screening"
SCR = os.path.join(DATA, "figures_scratch")
FIGDIR = os.path.join(REPO, "figures", "final")
RESDIR = os.path.join(REPO, "results", "final")
RES = os.path.join(REPO, "results")
REL = os.path.join(DATA, "relion")

sys.path.insert(0, "/n/home05/ianpoon/DNA-springs-project-2026-27/cryoEM/pp2a_day2025_reproduction/scripts")
import vizstyle  # noqa: E402
vizstyle.apply(matplotlib)

N = {}  # key -> {"value", "source", "how"}


def put(key, value, source, how=""):
    if isinstance(value, np.generic):
        value = value.item()
    N[key] = {"value": value, "source": source, "how": how}
    return value


def v(key):
    return N[key]["value"]


def jobjson(uid):
    return json.load(open(os.path.join(P6, uid, "job.json")))


def out_count(uid, output):
    """num_items of a CryoSPARC job output, from job.json (spec.outputs.<output>.num_items)."""
    o = jobjson(uid)["spec"]["outputs"][output]
    n = o.get("num_items")
    if n is None:
        n = o["results"][0]["num_items"]
    return int(n)


def put_count(key, uid, output, how=""):
    return put(key, out_count(uid, output), f"{P6}/{uid}/job.json", f"spec.outputs.{output}.num_items" + (f"; {how}" if how else ""))


def read_fsc_txt(path):
    """CryoSPARC Validation (FSC) txt; columns read by header name (first header line)."""
    with open(path) as f:
        names = f.readline().split()
    arr = np.loadtxt(path, skiprows=2)
    return {n: arr[:, i] for i, n in enumerate(names)}


def res_at(x, y, thr=0.143):
    idx = np.where(y < thr)[0]
    idx = idx[idx > 0]
    if len(idx) == 0:
        return float("nan")
    i = idx[0]
    x0, x1, y0, y1 = x[i - 1], x[i], y[i - 1], y[i]
    xc = x0 + (thr - y0) * (x1 - x0) / (y1 - y0) if y1 != y0 else x1
    return float(1.0 / xc)


def star_general(path):
    out = {}
    block = None
    for line in open(path):
        s = line.strip()
        if s.startswith("data_"):
            block = s
            continue
        if block == "data_general" and s.startswith("_rln"):
            k, val = s.split()[:2]
            out[k] = val
    return out


def log_value(path, pattern):
    m = re.search(pattern, open(path, "rb").read().decode("latin-1"))
    return float(m.group(1)) if m else float("nan")


# ----------------------------------------------------------------------------------------------- collect numbers
def collect():
    acq = os.path.join(RES, "g0_acquisition.json")
    g = json.load(open(acq))
    put("movies", g["movies"], acq, "movies")
    put("grid_squares", 20, os.path.join(RES, "raw_inventory_summary.json"), "grid_squares")
    assert json.load(open(os.path.join(RES, "raw_inventory_summary.json")))["grid_squares"] == 20
    put("kv", g["accel_kv"], acq, "accel_kv")
    put("cs_mm", g["cs_mm"], acq, "cs_mm")
    put("apix_phys", g["pixel_size_A"], acq, "pixel_size_A")
    put("dose", g["dose_e_per_A2_used"], acq, "dose_e_per_A2_used (vacuum calibration)")
    put("dose_detected", g["dose_detected_through_sample_e_per_A2"], acq, "dose_detected_through_sample_e_per_A2")
    put("fractions", g["fractions"], acq, "fractions")
    put("exposure_s", g["exposure_s"], acq, "exposure_s")
    put("microscope", g["microscope"], acq, "microscope")
    put("detector", g["detector"], acq, "detector")

    # nominal magnification / calibrated pixel from every EPU Data xml (cached)
    cache = os.path.join(SCR, "epu_xml_mag.json")
    if not os.path.exists(cache):
        mags, kvs, pix = {}, {}, {}
        xmls = glob.glob(os.path.join(DATA, "raw", "Images-Disc1", "GridSquare_*", "Data", "*.xml"))
        for x in xmls:
            s = open(x, errors="replace").read()
            m = re.search(r"<NominalMagnification>(\d+)", s)
            k = re.search(r"<AccelerationVoltage>(\d+)", s)
            p = re.search(r"<pixelSize><x><numericValue>([0-9.Ee+-]+)", s)
            for d, mm in ((mags, m), (kvs, k), (pix, p)):
                key = mm.group(1) if mm else "missing"
                if d is pix and mm:
                    key = f"{float(mm.group(1)) * 1e10:.4f}"
                d[key] = d.get(key, 0) + 1
        json.dump({"n_xml": len(xmls), "NominalMagnification": mags, "AccelerationVoltage_V": kvs,
                   "pixelSize_x_A": pix}, open(cache, "w"), indent=1)
    c = json.load(open(cache))
    mag = max(c["NominalMagnification"], key=c["NominalMagnification"].get)
    put("mag", int(mag), f"{DATA}/raw/Images-Disc1/GridSquare_*/Data/*.xml (tally cached in {cache})",
        f"<NominalMagnification>; {c['NominalMagnification']} over {c['n_xml']} xml")
    pxl = max(c["pixelSize_x_A"], key=c["pixelSize_x_A"].get)
    put("apix_epu", float(pxl), cache, f"EPU <pixelSize><x> calibration; tally {c['pixelSize_x_A']}")
    put("apix_epu_n", int(c["pixelSize_x_A"][pxl]), cache, "number of EPU xml files with that calibration")
    put("apix_epu_nxml", int(c["n_xml"]), cache, "n_xml")
    put("apix_scale", float(pxl) / g["pixel_size_A"], f"{cache} + {acq}", "apix_epu / apix_phys (0.71 used in processing)")

    # curated exposures: defocus range, pixel sizes, dose
    ez = os.path.join(DATA, "exports", "W2_curated__exposures.npz")
    z = np.load(ez)
    df_um = (z["ctf__df1_A"].astype(float) + z["ctf__df2_A"].astype(float)) / 2 / 1e4
    put("n_curated_npz", int(len(df_um)), ez, "rows")
    put("df_p5", float(np.percentile(df_um, 5)), ez, "p5 of (ctf__df1_A+ctf__df2_A)/2, um")
    put("df_p95", float(np.percentile(df_um, 95)), ez, "p95 of (ctf__df1_A+ctf__df2_A)/2, um")
    put("df_median", float(np.median(df_um)), ez, "median of (ctf__df1_A+ctf__df2_A)/2, um")
    put("df_min", float(df_um.min()), ez, "min of mean defocus, um")
    put("df_max", float(df_um.max()), ez, "max of mean defocus, um")
    put("apix_movie_cs", float(np.unique(z["movie_blob__psize_A"])[0]), ez, "movie_blob__psize_A (unique)")
    put("apix_mic", float(np.unique(z["micrograph_blob__psize_A"])[0]), ez, "micrograph_blob__psize_A (unique)")
    put("dose_cs", float(np.unique(z["mscope_params__total_dose_e_per_A2"])[0]), ez, "mscope_params__total_dose_e_per_A2")
    put("ctf_fit_median", float(np.median(z["ctf__ctf_fit_to_A"])), ez, "median ctf__ctf_fit_to_A")
    pz = os.path.join(DATA, "exports", "P0_ctf__exposures.npz")
    put("P0_ctf_fit_median", float(np.median(np.load(pz)["ctf__ctf_fit_to_A"])), pz, "median ctf__ctf_fit_to_A")

    # pilot
    put_count("P0_import", "J2", "imported_movies")
    put_count("P0_mics", "J3", "micrographs")
    put_count("P0_ctf", "J4", "exposures")
    put_count("P0_picks", "J5", "particles")
    put_count("P0_extract", "J6", "particles")
    put_count("P0_2D", "J7", "particles")
    put_count("P0_2D_K", "J7", "class_averages")
    # pre-processing
    put_count("import", "J8", "imported_movies")
    put("motion", out_count("J10", "micrographs") + out_count("J12", "micrographs"),
        f"{P6}/J10/job.json + {P6}/J12/job.json", "sum of spec.outputs.micrographs.num_items")
    put("ctf", out_count("J15", "exposures") + out_count("J16", "exposures"),
        f"{P6}/J15/job.json + {P6}/J16/job.json", "sum of spec.outputs.exposures.num_items")
    put_count("ctf_incomplete", "J16", "exposures_incomplete")
    put_count("curated", "J17", "exposures")
    cs = json.load(open(os.path.join(RES, "cs_state.json")))["jobs"]
    th = cs["W2_curated"]["thresholds"]
    put("cur_thr", f"CTF fit <= {th['ctf_fit_to_A_max']:g} A, astig <= {th['astig_A_max']:g} A, defocus "
        f"{th['defocus_A_min'] / 1e4:g}-{th['defocus_A_max'] / 1e4:g} um, rel. ice <= {th['ice_thickness_rel_max']:g}",
        os.path.join(RES, "cs_state.json"), "jobs.W2_curated.thresholds")
    # streams
    put_count("rod_picks", "J18", "particles")
    put("rod_power", cs["W2_power_rods"]["rows_out"], os.path.join(RES, "cs_state.json"), "jobs.W2_power_rods.rows_out (J22)")
    put_count("rod_extract", "J23", "particles")
    put_count("rod_2D", "J27", "particles")
    put_count("rod_2D_rej", "J27", "particles_rejected")
    put_count("rod_2D_K", "J27", "class_averages")
    put_count("anc_picks", "J19", "particles")
    put("anc_power", cs["W2_power_anchors"]["rows_out"], os.path.join(RES, "cs_state.json"), "jobs.W2_power_anchors.rows_out (J24)")
    put("anc_25", cs["W2_anchors_25pct"]["rows_out"], os.path.join(RES, "cs_state.json"), "jobs.W2_anchors_25pct.rows_out (J25)")
    put("anc_25_mics", cs["W2_anchors_25pct"]["micrographs"], os.path.join(RES, "cs_state.json"), "jobs.W2_anchors_25pct.micrographs")
    put_count("anc_extract", "J26", "particles")
    put_count("anc_2D", "J28", "particles")
    put_count("anc_2D_rej", "J28", "particles_rejected")
    for k, uid in (("rod_2D_window", "J27"), ("anc_2D_window", "J28"), ("A1_2D_window", "J45")):
        put(k, float(jobjson(uid)["spec"]["params"]["class2D_window_inner_A"]), f"{P6}/{uid}/job.json",
            "spec.params.class2D_window_inner_A")

    # blind review: labels for every reviewed class; closed loops
    labs = {}
    for stream in ("rods", "anchors"):
        f = os.path.join(RES, f"class2d_labels_W2_2D_{stream}.csv")
        rows = list(csv.DictReader(open(f)))
        labs[stream] = rows
        put(f"n_reviewed_{stream}", len(rows), f, "rows (classes with >= 270 particles)")
        put(f"reviewed_ptcl_{stream}", int(sum(int(r["particles"]) for r in rows)), f, "sum particles")
        put(f"closed_loop_{stream}", int(sum(r["consensus_label"] == "closed_loop" for r in rows)), f,
            "classes with consensus_label == closed_loop")
    put("n_reviewed", v("n_reviewed_rods") + v("n_reviewed_anchors"), "class2d_labels_W2_2D_{rods,anchors}.csv", "sum of rows")
    put("closed_loop", v("closed_loop_rods") + v("closed_loop_anchors"), "class2d_labels_W2_2D_{rods,anchors}.csv", "sum")
    rnd = [r for r in labs["rods"] if r["consensus_label"] == "rod_no_dot"]
    put("rod_no_dot_classes", len(rnd), os.path.join(RES, "class2d_labels_W2_2D_rods.csv"), "consensus_label == rod_no_dot")
    put("rod_no_dot_ptcl", int(sum(int(r["particles"]) for r in rnd)), os.path.join(RES, "class2d_labels_W2_2D_rods.csv"),
        "sum particles, consensus_label == rod_no_dot")
    rnd_sorted = sorted(rnd, key=lambda r: -int(r["particles"]))
    put("rod_no_dot_top", [int(r["class"]) for r in rnd_sorted[:6]], os.path.join(RES, "class2d_labels_W2_2D_rods.csv"),
        "six largest rod_no_dot classes (shown as averages)")
    br = json.load(open(os.path.join(RES, "class2d_blind_review_W2.json")))
    put("closed_loop_checks", len(br["loop_checks"]), os.path.join(RES, "class2d_blind_review_W2.json"),
        "len(loop_checks): closed-loop votes needing a skeptic")

    # W3
    for k, key in (("W3_duplex", "W3_duplex_particles"), ("W3_junk", "W3_junk_particles"), ("W3_ice", "W3_ice_particles"),
                   ("W3_pairs", "W3_pairs_particles"), ("W3_pool", "W3_pooled_particles")):
        put(k, cs[key]["rows_out"], os.path.join(RES, "cs_state.json"), f"jobs.{key}.rows_out ({cs[key]['uid']})")
    assert v("W3_duplex") == v("rod_no_dot_ptcl")
    for i, nm in enumerate(("duplex", "junk", "ice", "pairs")):  # class order = vols order in 31_duplex_control.py
        put_count(f"W3_het_{nm}", "J38", f"particles_class_{i}", "class order duplex/junk/ice/pairs per scripts/31_duplex_control.py")
    put_count("W3_nu", "J39", "particles")
    # final-particle composition by 2D label (classes2d.json from scripts/50_figS5_prep.py)
    for W, tag in (("W3", "W3_duplex"), ("W4", "W4_dotrod")):
        cj = os.path.join(SCR, tag, "classes2d.json")
        ci = json.load(open(cj))
        put(f"{W}_final_by_label", ci["final_particles_by_label"], cj,
            "final_particles_by_label (final NU particles traced to their 2D class label)")
        put(f"{W}_final_by_label_n", ci["n_final_particles"], cj, "n_final_particles")
    put_count("W3_nu_rej", "J39", "particles_rejected")
    # W4 / A1
    a1 = cs["A1_recentred"]
    put("A1_recentred", a1["rows_out"], os.path.join(RES, "cs_state.json"), "jobs.A1_recentred.rows_out (J43)")
    for lab in ("dot_long_rod", "dot_short_rod", "rod_no_dot"):
        put(f"A1_in_{lab}", a1["by_label"][lab], os.path.join(RES, "cs_state.json"), f"jobs.A1_recentred.by_label.{lab}")
    put_count("A1_extract", "J44", "particles")
    put_count("A1_extract_mics", "J44", "micrographs")
    put_count("A1_2D", "J45", "particles")
    put_count("A1_2D_K", "J45", "class_averages")
    for k, key in (("W4_dot", "W4_dot_particles"), ("W4_free", "W4_free_particles"), ("W4_haze", "W4_haze_particles"),
                   ("W4_ice", "W4_ice_particles"), ("W4_pool", "W4_pooled_particles")):
        put(k, cs[key]["rows_out"], os.path.join(RES, "cs_state.json"), f"jobs.{key}.rows_out ({cs[key]['uid']})")
        put(k + "_classes", len(cs[key].get("classes", [])), os.path.join(RES, "cs_state.json"), f"len(jobs.{key}.classes)")
    sel = json.load(open(os.path.join(RES, "A1_length_vs_hinge.json")))["selection"]
    dots = sorted([s for s in sel if s.get("group") == "DOT"], key=lambda s: -s["n"])
    put("A1_dot_top", [int(s["cls"]) for s in dots[:6]], os.path.join(RES, "A1_length_vs_hinge.json"),
        "selection, group DOT, six largest (shown as averages)")
    assert sum(s["n"] for s in dots) == v("W4_dot")
    for i, nm in enumerate(("dot", "free", "haze", "ice")):  # vols order in 35_dot_rod_3d.py
        put_count(f"W4_het_{nm}", "J55", f"particles_class_{i}", "class order dot/free/haze/ice per scripts/35_dot_rod_3d.py")
    put_count("W4_nu", "J56", "particles")
    put_count("W4_nu_rej", "J56", "particles_rejected")
    put_count("W4_c3d", "J59", "particles_all_classes")
    put_count("W4_c3d_K", "J59", "volumes_all_classes")
    c3 = [out_count("J59", f"particles_class_{i}") for i in range(v("W4_c3d_K"))]
    put("W4_c3d_min", min(c3), f"{P6}/J59/job.json", "min over spec.outputs.particles_class_N.num_items")
    put("W4_c3d_max", max(c3), f"{P6}/J59/job.json", "max over spec.outputs.particles_class_N.num_items")
    p59 = jobjson("J59")["spec"]["params"]
    put("W4_c3d_res", float(p59["class3D_target_res"]), f"{P6}/J59/job.json", "spec.params.class3D_target_res")
    put("W4_c3d_pcs", int(p59["class3D_PCA_num_components"]), f"{P6}/J59/job.json", "spec.params.class3D_PCA_num_components")

    # per-map numbers
    for W, nu, val, od, mtag in (("W3", "J39", "J40", "J41", "W3_duplex"), ("W4", "J56", "J57", "J58", "W4_dotrod")):
        pcs = os.path.join(P6, nu, f"{nu}_003_particles.cs")
        p = np.load(pcs)
        put(f"{W}_final_ptcl", int(len(p)), pcs, "rows")
        put(f"{W}_apix", float(np.unique(p["blob/psize_A"])[0]), pcs, "blob/psize_A (unique)")
        put(f"{W}_box", int(p["blob/shape"][0][0]), pcs, "blob/shape[0]")
        pt = np.load(os.path.join(P6, nu, f"{nu}_passthrough_particles.cs"))
        mu = pt["location/micrograph_uid"][np.isin(pt["uid"], p["uid"])]
        put(f"{W}_mics", int(len(np.unique(mu))), os.path.join(P6, nu, f"{nu}_passthrough_particles.cs"),
            f"unique location/micrograph_uid of the uids in {nu}_003_particles.cs")
        put(f"{W}_sym", jobjson(nu)["spec"]["params"]["refine_symmetry"], f"{P6}/{nu}/job.json", "spec.params.refine_symmetry")
        fsc = read_fsc_txt(os.path.join(P6, val, f"{val}_fsc.txt"))
        for col, tag in (("fsc_nomask", "nomask"), ("fsc_sphericalmask", "sph"), ("fsc_resmask", "tight")):
            put(f"{W}_cs_{tag}", res_at(fsc["frequency"], fsc[col]), f"{P6}/{val}/{val}_fsc.txt",
                f"first 0.143 crossing of column '{col}' vs 'frequency' (linear interpolation)")
        lg = os.path.join(P6, od, "events.bson")
        put(f"{W}_cfar", log_value(lg, r"cFAR:\s*([0-9.]+)"), lg, "event text 'cFAR: x' (job event log)")
        put(f"{W}_scf", log_value(lg, r"SCF\*:\s*([0-9.]+)"), lg, "event text 'SCF*: x' (job event log)")
        for pp, tag in (("pp_csmask", "csmask"), ("pp_loose", "loose")):
            st = os.path.join(REL, mtag, pp, "postprocess.star")
            gen = star_general(st)
            put(f"{W}_rel_{tag}", float(gen["_rlnFinalResolution"]), st, "data_general _rlnFinalResolution")
            put(f"{W}_B_{tag}", float(gen["_rlnBfactorUsedForSharpening"]), st, "data_general _rlnBfactorUsedForSharpening")
            put(f"{W}_guinier_r_{tag}", float(gen["_rlnCorrelationFitGuinierPlot"]), st, "data_general _rlnCorrelationFitGuinierPlot")
        # local resolution range inside the CryoSPARC validation (tight FSC) mask
        mk_path = os.path.join(P6, val, f"{val}_fsc_mask.mrc")
        mk = mrcfile.read(mk_path) > 0.5
        lr_path = os.path.join(REL, mtag, "relion_locres.mrc")
        lr = mrcfile.read(lr_path)[mk]
        lr = lr[np.isfinite(lr) & (lr > 0)]
        put(f"{W}_locres_rel", [float(np.percentile(lr, 5)), float(np.median(lr)), float(np.percentile(lr, 95))],
            f"{lr_path} within {mk_path} > 0.5", "p5 / median / p95 of RELION 3.1 --locres values")
        if W == "W3":
            cl_path = os.path.join(P6, "J42", "J42_map_locres.mrc")
            cl = mrcfile.read(cl_path)[mk]
            cl = cl[np.isfinite(cl) & (cl > 0)]
            put("W3_locres_cs", [float(np.percentile(cl, 5)), float(np.median(cl)), float(np.percentile(cl, 95))],
                f"{cl_path} within {mk_path} > 0.5", "p5 / median / p95 of CryoSPARC Local Resolution (J42)")
    put("curation_rejected", v("ctf") - v("curated"), "derived from keys ctf and curated", "ctf - curated")
    put("power_rej_rod", v("rod_picks") - v("rod_power"), "derived from keys rod_picks and rod_power", "rod_picks - rod_power")
    return N


# ----------------------------------------------------------------------------------------------- Table S1
def fmt(n):
    return f"{n:,}"


def build_table():
    rows, src = [], {}

    def add(section, name, w3, w4, s3, s4):
        rows.append([section, name, w3, w4])
        src[name] = {"W3": s3, "W4": s4}

    def S(*keys, note=""):
        d = {"numbers": {k: N[k] for k in keys}}
        if note:
            d["note"] = note
        return d

    sec = "Data collection and processing"
    mag = f"{fmt(v('mag'))}× (nominal)"
    add(sec, "Magnification", mag, mag, S("mag"), S("mag"))
    kv = f"{v('kv'):g}"
    add(sec, "Voltage (kV)", kv, kv, S("kv"), S("kv"))
    ex = f"{v('dose'):.2f} (vacuum calibration; {v('dose_detected'):.1f} detected through sample); {v('fractions')} fractions, {v('exposure_s'):.2f} s"
    add(sec, "Electron exposure (e–/Å²)", ex, ex, S("dose", "dose_detected", "fractions", "exposure_s", "dose_cs"),
        S("dose", "dose_detected", "fractions", "exposure_s", "dose_cs"))
    df = f"{v('df_p5'):.2f}–{v('df_p95'):.2f} (p5–p95 of {fmt(v('n_curated_npz'))} curated micrographs; median {v('df_median'):.2f})"
    add(sec, "Defocus range (μm, underfocus)", df, df, S("df_p5", "df_p95", "df_median", "n_curated_npz"),
        S("df_p5", "df_p95", "df_median", "n_curated_npz"))
    px3 = (f"{v('apix_phys'):.2f} physical, as used in processing (EPU xml calibration {v('apix_epu'):.4f}) / {v('apix_mic'):.2f} processed "
           f"(motion-corrected, F-crop ½) / {v('W3_apix'):.2f} refined ({v('W3_box')} px box)")
    px4 = (f"{v('apix_phys'):.2f} physical, as used in processing (EPU xml calibration {v('apix_epu'):.4f}) / {v('apix_mic'):.2f} processed "
           f"(motion-corrected, F-crop ½) / {v('W4_apix'):.2f} refined ({v('W4_box')} px box)")
    pxn = f"; if {v('apix_epu'):.4f} is correct, all Å values ×{v('apix_scale'):.4f} (see pixel-size note)"
    add(sec, "Pixel size (Å)", px3 + pxn, px4 + pxn,
        S("apix_phys", "apix_epu", "apix_epu_n", "apix_scale", "apix_movie_cs", "apix_mic", "W3_apix", "W3_box"),
        S("apix_phys", "apix_epu", "apix_epu_n", "apix_scale", "apix_movie_cs", "apix_mic", "W4_apix", "W4_box"))
    add(sec, "Symmetry imposed", v("W3_sym"), v("W4_sym"), S("W3_sym"), S("W4_sym"))
    add(sec, "Micrographs (collected / curated / contributing to map)",
        f"{fmt(v('movies'))} / {fmt(v('curated'))} / {fmt(v('W3_mics'))}",
        f"{fmt(v('movies'))} / {fmt(v('curated'))} / {fmt(v('W4_mics'))}",
        S("movies", "curated", "W3_mics"), S("movies", "curated", "W4_mics"))
    add(sec, "Initial particle images (no.)",
        f"{fmt(v('rod_picks'))} rod-stream blob picks (J18); {fmt(v('rod_extract'))} extracted after power filter",
        f"{fmt(v('rod_picks'))} rod-stream blob picks (J18); {fmt(v('A1_extract'))} re-centred and re-extracted (J44)",
        S("rod_picks", "rod_power", "rod_extract"), S("rod_picks", "A1_recentred", "A1_extract"))
    add(sec, "Final particle images (no.)", fmt(v("W3_final_ptcl")), fmt(v("W4_final_ptcl")),
        S("W3_final_ptcl", "W3_nu"), S("W4_final_ptcl", "W4_nu"))
    add(sec, "Map resolution (Å)",
        f"{v('W3_cs_tight'):.1f} (CryoSPARC GSFSC, tight mask) / {v('W3_rel_loose'):.1f} (RELION 3.1, loose mask, phase-randomised)",
        f"{v('W4_cs_tight'):.1f} (CryoSPARC GSFSC, tight mask) / {v('W4_rel_loose'):.1f} (RELION 3.1, loose mask, phase-randomised)",
        S("W3_cs_tight", "W3_rel_loose", "W3_rel_csmask", "W3_cs_nomask", "W3_cs_sph"),
        S("W4_cs_tight", "W4_rel_loose", "W4_rel_csmask", "W4_cs_nomask", "W4_cs_sph"))
    add(sec, "    FSC threshold", "0.143", "0.143", {"note": "Day et al. 2025 Table S1 convention; gold-standard half maps"},
        {"note": "Day et al. 2025 Table S1 convention; gold-standard half maps"})
    add(sec, "Map resolution, other masks (Å)",
        f"no mask {v('W3_cs_nomask'):.1f} / spherical {v('W3_cs_sph'):.1f} (CryoSPARC); RELION with CryoSPARC mask {v('W3_rel_csmask'):.1f}",
        f"no mask {v('W4_cs_nomask'):.1f} / spherical {v('W4_cs_sph'):.1f} (CryoSPARC); RELION with CryoSPARC mask {v('W4_rel_csmask'):.1f}",
        S("W3_cs_nomask", "W3_cs_sph", "W3_rel_csmask"), S("W4_cs_nomask", "W4_cs_sph", "W4_rel_csmask"))
    l3, c3, l4 = v("W3_locres_rel"), v("W3_locres_cs"), v("W4_locres_rel")
    add(sec, "Map resolution range (Å)",
        f"{l3[0]:.1f}–{l3[2]:.1f} (RELION locres, p5–p95 in tight mask; median {l3[1]:.1f}); "
        f"{c3[0]:.1f}–{c3[2]:.1f} (CryoSPARC J42; median {c3[1]:.1f})",
        f"{l4[0]:.1f}–{l4[2]:.1f} (RELION locres, p5–p95 in tight mask; median {l4[1]:.1f}); no CryoSPARC local-resolution job",
        S("W3_locres_rel", "W3_locres_cs"), S("W4_locres_rel"))
    add(sec, "Orientation: cFAR / SCF*", f"{v('W3_cfar'):.2f} / {v('W3_scf'):.2f}", f"{v('W4_cfar'):.2f} / {v('W4_scf'):.2f}",
        S("W3_cfar", "W3_scf"), S("W4_cfar", "W4_scf"))

    sec = "Refinement"
    nb3 = "not built: ~8–9 Å map of a rod with a terminal lobe, no secondary structure resolved, so an atomic model would be unconstrained"
    nb4 = "not built: ~11 Å anisotropic average (cFAR 0.17) over heterogeneous sfGFP–DNA assemblies; nothing to refine a model against"
    add(sec, "Initial model used (PDB code)", nb3, nb4, {"note": "report.md §10 (map reading)"}, {"note": "report.md §13b (heterogeneity, anisotropy)"})
    add(sec, "Model resolution (Å)", "not built (no model)", "not built (no model)", {"note": "no model"}, {"note": "no model"})
    add(sec, "Map sharpening B factor (Å²)",
        f"{v('W3_B_csmask'):.0f} (RELION auto-B, CryoSPARC mask) / {v('W3_B_loose'):.0f} (loose mask); Guinier fit from 10 Å only, r = {v('W3_guinier_r_csmask'):.2f} / {v('W3_guinier_r_loose'):.2f}",
        f"{v('W4_B_csmask'):.0f} (RELION auto-B, CryoSPARC mask) / {v('W4_B_loose'):.0f} (loose mask); Guinier fit from 10 Å only, r = {v('W4_guinier_r_csmask'):.2f} / {v('W4_guinier_r_loose'):.2f}",
        S("W3_B_csmask", "W3_B_loose", "W3_guinier_r_csmask", "W3_guinier_r_loose"),
        S("W4_B_csmask", "W4_B_loose", "W4_guinier_r_csmask", "W4_guinier_r_loose"))
    for name in ("Model composition: non-hydrogen atoms / residues / ligands", "B factors (Å²)",
                 "R.m.s. deviations: bond lengths (Å) / bond angles (°)", "MolProbity score", "Clashscore",
                 "Poor rotamers (%)", "Ramachandran plot: favored / allowed / disallowed (%)"):
        add("Validation" if name in ("MolProbity score", "Clashscore", "Poor rotamers (%)") or name.startswith("Ramachandran")
            else sec, name, "not built (no model)", "not built (no model)", {"note": "no atomic model; see 'Initial model used'"},
            {"note": "no atomic model; see 'Initial model used'"})
    return rows, src


COL_W3 = "W3 rod_no_dot (\u201cfree DNA\u201d) control, J39"
COL_W4 = "W4 dot-on-rod (J56)"


def pixel_note():
    return (f"Pixel size: processed with {v('apix_phys'):.2f} Å/px (user screenshot); EPU xml calibration "
            f"{v('apix_epu'):.4f} Å/px (all {fmt(v('apix_epu_n'))} movies) — if {v('apix_epu'):.4f} is correct, all Å "
            f"values scale by ×{v('apix_scale'):.4f}; to be confirmed with the facility.")


def write_table(rows, src):
    os.makedirs(RESDIR, exist_ok=True)
    cols = ["section", "row", COL_W3, COL_W4]
    with open(os.path.join(RESDIR, "tableS1.csv"), "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(cols)
        w.writerows(rows)
        w.writerow(["notes", "pixel size", pixel_note(), pixel_note()])
    md = ["# Table S1 analogue: cryo-EM data collection and processing for the two final maps",
          "",
          "sfGFP + 50 bp DNA screening session (12-11-2026), CryoSPARC v5 project P6. Columns follow Day et al. 2025 "
          "*Structure* Table S1. Every cell's source file and field: `tableS1_sources.json`; built by "
          "`scripts/50_figS6_tableS1.py`.",
          "", f"| | {COL_W3} | {COL_W4} |", "|---|---|---|"]
    last = None
    for sec, name, a, b in rows:
        if sec != last:
            md.append(f"| **{sec}** | | |")
            last = sec
        md.append(f"| {name.strip() if not name.startswith('    ') else '&nbsp;&nbsp;' + name.strip()} | {a} | {b} |")
    md += ["",
           "Notes. Tight mask = CryoSPARC Validation (FSC) auto-tightened mask (column `fsc_resmask`); RELION values are "
           "`relion_postprocess` phase-randomisation-corrected FSC (`scripts/32_relion_postprocess.sh`). The W4 loose-mask "
           "value is the honest whole-particle figure (report §13b). Local-resolution ranges are the 5th–95th percentiles of "
           "RELION local-resolution voxels inside the CryoSPARC validation (tight FSC) mask, the same definition as the Fig S5 summary boxes (other voxel sets give other ranges; they are listed, labelled as secondary, in `local_resolution_stats.json`). The W3 column additionally lists CryoSPARC Local Resolution (J42, same tight mask) for reference; no equivalent CryoSPARC job was run for W4. RELION auto-B values are fitted from 10 Å to Nyquist "
           f"({2 * v('W3_apix'):.2f} Å) on "
           "~8–11 Å maps, so they are poorly determined and are listed only because Table S1 asks for the sharpening B "
           "factor.", "", pixel_note()]
    open(os.path.join(RESDIR, "tableS1.md"), "w").write("\n".join(md) + "\n")
    json.dump({"columns": cols[2:], "cells": src, "pixel_size_note": pixel_note(),
               "pixel_size_note_numbers": {k: N[k] for k in ("apix_phys", "apix_epu", "apix_epu_n", "apix_epu_nxml", "apix_scale")},
               "script": os.path.join(REPO, "scripts", "50_figS6_tableS1.py")},
              open(os.path.join(RESDIR, "tableS1_sources.json"), "w"), indent=1, default=str)


# ----------------------------------------------------------------------------------------------- figure
C = {"pre": vizstyle.INK_2, "pilot": vizstyle.INK_MUTED, "anc": vizstyle.SERIES[1], "rod": vizstyle.SERIES[6],
     "W3": vizstyle.SERIES[0], "A1": vizstyle.SERIES[4], "W4": vizstyle.SERIES[2], "gate": vizstyle.SERIES[3],
     "exit": vizstyle.INK_MUTED}
XMAX, YMIN, YMAX = 100.0, 14.0, 120.0
FW, FH = 20.0, 12.04
LH = 1.42          # y units per text line
FS = 7.4
THUMBS = os.path.join(SCR, "renders")  # written by scripts/50_figS6_thumbs.cxc


def tint(hexc, a=0.10):
    c = np.array(matplotlib.colors.to_rgb(hexc))
    s = np.array(matplotlib.colors.to_rgb(vizstyle.SURFACE))
    return tuple(s * (1 - a) + c * a)


def box_h(nlines):
    return 0.9 + (1 + nlines) * LH


class Box:
    """Rounded box: bold title line + body lines.  ytop given, or yc (vertical centre)."""

    def __init__(self, ax, x, ytop, w, title, lines, color, kind="step", fs=FS, yc=None):
        self.x, self.w = x, w
        self.h = box_h(len(lines))
        if yc is not None:
            ytop = yc + self.h / 2
        self.ytop, self.y = ytop, ytop - self.h
        muted = kind in ("exit", "notrun", "input")
        fc = vizstyle.SURFACE if muted else tint(color, 0.22 if kind == "final" else 0.10)
        ax.add_patch(FancyBboxPatch((x, self.y), w, self.h, boxstyle="round,pad=0,rounding_size=0.6",
                                    fc=fc, ec=color, lw=1.7 if kind == "final" else 1.0,
                                    ls=(0, (3, 2)) if muted else "-", zorder=2))
        tc = vizstyle.INK_MUTED if muted else vizstyle.INK
        bc = vizstyle.INK_MUTED if muted else (vizstyle.INK if kind == "final" else vizstyle.INK_2)
        ax.text(x + 0.6, ytop - 0.45, title, fontsize=fs, fontweight="bold", color=tc, va="top", ha="left", zorder=3)
        for i, ln in enumerate(lines):
            ax.text(x + 0.6, ytop - 0.45 - (i + 1) * LH, ln, fontsize=fs, color=bc, va="top", ha="left", zorder=3)

    @property
    def bot(self):
        return (self.x + self.w / 2, self.y)

    @property
    def top(self):
        return (self.x + self.w / 2, self.ytop)

    def left(self, f=0.5):
        return (self.x, self.y + self.h * f)

    def right(self, f=0.5):
        return (self.x + self.w, self.y + self.h * f)


def arrow(ax, p, q, color=vizstyle.INK_2, dashed=False, lw=1.1):
    ax.add_patch(FancyArrowPatch(p, q, arrowstyle="-|>", mutation_scale=9, color=color, lw=lw,
                                 ls=(0, (3, 2)) if dashed else "-", shrinkA=0, shrinkB=1, zorder=1))


def path(ax, pts, color=vizstyle.INK_2, dashed=False, head=True):
    """polyline through pts, arrowhead on the last segment"""
    ls = (0, (3, 2)) if dashed else "-"
    xs, ys = zip(*pts)
    if head:
        ax.plot(xs[:-1], ys[:-1], color=color, lw=1.1, ls=ls, zorder=1, solid_capstyle="butt")
        arrow(ax, pts[-2], pts[-1], color=color, dashed=dashed)
    else:
        ax.plot(xs, ys, color=color, lw=1.1, ls=ls, zorder=1)


def stack(ax, x, ytop, w, items, color, gap=1.8):
    boxes, y = [], ytop
    for it in items:
        b = Box(ax, x, y, w, it[0], it[1], color, it[2] if len(it) > 2 else "step")
        if boxes:
            arrow(ax, boxes[-1].bot, b.top, color=color)
        boxes.append(b)
        y = b.y - gap
    return boxes


def side_box(ax, src, x, w, title, lines, f=0.5, into=False, fs=FS - 0.3):
    """grey dashed box level with src: a side exit (arrow src -> box) or a decoy input (box -> src)"""
    yc = src.y + src.h * f
    b = Box(ax, x, None, w, title, lines, C["exit"], kind="exit", fs=fs, yc=yc)
    near_src = src.right(f) if x > src.x else src.left(f)
    near_box = (x, yc) if x > src.x else (x + w, yc)
    if into:
        arrow(ax, near_box, near_src, color=C["exit"], dashed=True)
    else:
        arrow(ax, near_src, near_box, color=C["exit"], dashed=True)
    return b


def img_axes(fig, x, y, w, h):
    return fig.add_axes([x / XMAX, (y - YMIN) / (YMAX - YMIN), w / XMAX, h / (YMAX - YMIN)])


def class_avg_strip(fig, mrc, classes, x, y, w, h, color):
    st = mrcfile.read(mrc)
    cw = w / len(classes)
    for i, c in enumerate(classes):
        a = img_axes(fig, x + i * cw + 0.1, y, cw - 0.2, h)
        im = st[c]
        lo, hi = np.percentile(im, [1, 99.5])
        a.imshow(im, cmap="gray", vmin=lo, vmax=hi, interpolation="lanczos")
        a.set_xticks([]); a.set_yticks([])
        for s in a.spines.values():
            s.set_visible(True); s.set_color(color); s.set_linewidth(0.8)


def thumb(fig, png, x, y, w, h):
    """auto-cropped ChimeraX render, aspect preserved, centred in the (x, y, w, h) slot"""
    import matplotlib.image as mpimg
    im = mpimg.imread(png)
    nonwhite = np.where(im[..., :3].min(axis=2) < 0.97)
    r0, r1, c0, c1 = nonwhite[0].min(), nonwhite[0].max(), nonwhite[1].min(), nonwhite[1].max()
    pad = 8
    im = im[max(r0 - pad, 0):r1 + pad, max(c0 - pad, 0):c1 + pad, :3].copy()
    bg = im.min(axis=2) > 0.985
    im[bg] = matplotlib.colors.to_rgb(vizstyle.SURFACE)
    a = img_axes(fig, x, y, w, h)
    a.imshow(im, interpolation="lanczos")
    a.set_anchor("C")
    a.axis("off")


def make_figure():
    fig = plt.figure(figsize=(FW, FH), layout="none")
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, XMAX); ax.set_ylim(YMIN, YMAX); ax.axis("off")
    used = set()

    def u(*keys):
        used.update(keys)

    ax.text(1.0, YMAX - 1.0, "Fig S6 analogue: cryo-EM processing workflow, sfGFP + 50 bp DNA screening session "
            "(CryoSPARC v5 project P6)", fontsize=11.5, fontweight="bold", va="top")
    ax.text(1.0, YMAX - 3.1, "Counts are particles unless marked; J numbers are P6 jobs. Solid boxes: path to a final map; "
            "dashed grey: rejected / decoy side exits, decoy inputs, and the planned chimera branch the gate stopped. "
            "Every number: figS6_workflow.json.", fontsize=8, color=vizstyle.INK_2, va="top")

    # columns (x0, width)
    X0, W0 = 1.0, 14.5     # pilot
    XA, WA = 17.5, 14.5    # anchors / W3 side boxes
    XR, WR = 34.0, 16.5    # rods / W3
    XB, WB = 52.5, 15.0    # A1
    XW, WW = 69.5, 15.5    # W4
    XE, WE = 87.0, 12.0    # W4 side exits / not-run branch

    def head(x, y, text, color):
        ax.text(x, y, text, fontsize=8.6, fontweight="bold", color=color, va="bottom", ha="left")

    # ---------------- data collection + pre-processing
    u("movies", "grid_squares", "kv", "mag", "apix_phys", "apix_epu", "dose", "fractions")
    dc = Box(ax, XA, 113.0, 52, "Data collection (12-11-2026)",
             [f"{fmt(v('movies'))} movies in {v('grid_squares')} grid squares; Glacios {v('kv'):g} kV, Falcon 4i behind a 10 eV energy filter",
              f"{fmt(v('mag'))}× nominal, {v('apix_phys'):.2f} Å/px as processed (EPU xml {v('apix_epu'):.4f}, to be confirmed); "
              f"{v('dose'):.2f} e–/Å² in {v('fractions')} fractions"], C["pre"])
    u("import", "motion", "apix_mic", "ctf", "ctf_incomplete", "curated", "cur_thr", "curation_rejected")
    yb = dc.y - 3.0
    m1 = Box(ax, XA, yb, 17, "Import J8 → Patch motion J10+J12",
             [f"{fmt(v('import'))} movies → {fmt(v('motion'))} micrographs", f"F-crop ½ → {v('apix_mic'):.2f} Å/px"], C["pre"])
    m2 = Box(ax, XR + 3.0, yb, 12.5, "Patch CTF J15+J16", [f"{fmt(v('ctf'))} micrographs", f"({v('ctf_incomplete')} incomplete)"], C["pre"])
    m3 = Box(ax, XB, yb, 19.5, "Curation J17 (scripted thresholds)",
             [f"{fmt(v('curated'))} micrographs", "CTF fit ≤ 6 Å, astig ≤ 1000 Å, Δf 0.5–2.8 µm, ice ≤ 1.2"], C["pre"])
    assert v("cur_thr").startswith("CTF fit <= 6 A, astig <= 1000 A, defocus 0.5-2.8 um, rel. ice <= 1.2")
    arrow(ax, (m1.top[0], dc.y), m1.top, color=C["pre"])
    arrow(ax, m1.right(), m2.left(), color=C["pre"])
    arrow(ax, m2.right(), m3.left(), color=C["pre"])
    side_box(ax, m3, XW + 3.5, 12.0, "Rejected", [f"{fmt(v('curation_rejected'))} micrographs"])

    # ---------------- pilot
    u("P0_import", "P0_mics", "P0_ctf_fit_median", "P0_picks", "P0_extract", "P0_2D", "P0_2D_K")
    head(X0, 112.0 - 0.2, "Pilot P0 (412-movie subset)", C["pilot"])
    pil = stack(ax, X0, 111.0, W0, [
        ("Import J2 → motion J3 → CTF J4", [f"{fmt(v('P0_import'))} movies → {fmt(v('P0_mics'))} micrographs",
                                           f"CTF fit median {v('P0_ctf_fit_median'):.1f} Å"]),
        ("Blob pick J5 (80–140 Å)", [f"{fmt(v('P0_picks'))} picks"]),
        ("Extract J6 (2.84 Å/px)", [f"{fmt(v('P0_extract'))}"]),
        (f"2D J7 (K = {v('P0_2D_K')}, 12 Å)", [f"{fmt(v('P0_2D'))} kept; no ring class:", "dark dots, thin threads,",
                                                "chance dot pairs"]),
        ("Decision: redesign picking", ["two streams, sfGFP-sized", "'anchors' + DNA-sized 'rods';",
                                        "local-power filter (pick NCC", "did not separate objects)"], "input"),
    ], C["pilot"], gap=1.6)
    path(ax, [dc.left(0.5), (X0 + W0 + 1.0, dc.left(0.5)[1]), (X0 + W0 + 1.0, pil[0].right(0.5)[1]), pil[0].right(0.5)],
         color=C["pilot"], dashed=True)

    # ---------------- streams
    ystream = m1.y - 6.0
    head(XA, ystream + 0.6, "Anchor stream (sfGFP dots)", C["anc"])
    head(XR, ystream + 0.6, "Rod stream (DNA-sized objects)", C["rod"])
    # curation -> both streams (bus)
    ybus = ystream + 3.6
    dA, dR = XA + WA - 2.5, XR + WR - 2.5   # drop points right of the stream headers
    path(ax, [m3.bot, (m3.bot[0], ybus), (dA, ybus)], color=C["pre"], head=False)
    path(ax, [(dA, ybus), (dA, ystream - 0.1)], color=C["pre"])
    path(ax, [(dR, ybus), (dR, ystream - 0.1)], color=C["pre"])
    # pilot decision informs both stream designs
    path(ax, [pil[-1].right(0.5), (XA - 1.0, pil[-1].right(0.5)[1]), (XA - 1.0, ybus), (dA - 0.3, ybus)],
         color=C["pilot"], dashed=True)
    u("anc_picks", "anc_power", "anc_25", "anc_25_mics", "anc_extract", "anc_2D", "anc_2D_rej", "anc_2D_window")
    anc = stack(ax, XA, ystream - 0.1, WA, [
        ("Blob pick J19 (circle 30–50 Å)", [f"{fmt(v('anc_picks'))} picks"]),
        ("Power filter J24", [f"{fmt(v('anc_power'))} (top 50 % per Δf bin)"]),
        ("Random 25 % of micrographs J25", [f"{fmt(v('anc_25'))} ({fmt(v('anc_25_mics'))} micrographs)"]),
        ("Extract J26 (192 → 96 px)", [f"{fmt(v('anc_extract'))}"]),
        ("2D round 1 J28 (K = 200, 12 Å)", [f"{fmt(v('anc_2D'))} kept, {v('anc_2D_window'):.0f} Å window",
                                            f"({fmt(v('anc_2D_rej'))} rejected)"]),
    ], C["anc"], gap=1.6)
    u("rod_picks", "rod_power", "power_rej_rod", "rod_extract", "rod_2D", "rod_2D_rej", "rod_2D_window")
    rod = stack(ax, XR, ystream - 0.1, WR, [
        ("Blob pick J18 (ellipse 20–180 Å)", [f"{fmt(v('rod_picks'))} picks"]),
        ("Power filter J22", [f"{fmt(v('rod_power'))} (top 50 % per Δf bin)"]),
        ("Extract J23 (256 → 128 px, 2.84 Å/px)", [f"{fmt(v('rod_extract'))}"]),
        ("2D round 1 J27 (K = 200, 12 Å)", [f"{fmt(v('rod_2D'))} kept, {v('rod_2D_window'):.0f} Å window",
                                            f"({fmt(v('rod_2D_rej'))} rejected)"]),
    ], C["rod"], gap=1.6)
    side_box(ax, rod[1], XB, 12.0, "Low local power", [f"{fmt(v('power_rej_rod'))} picks"])

    # ---------------- gate G2
    u("n_reviewed", "closed_loop", "closed_loop_checks")
    gw, gh = 26.0, 9.0
    gcx = XR + WR / 2
    gy = anc[-1].y - 1.0
    gx = gcx - gw / 2
    ax.add_patch(Polygon([(gx, gy), (gcx, gy + gh / 2), (gx + gw, gy), (gcx, gy - gh / 2)], closed=True,
                         fc=tint(C["gate"], 0.18), ec=C["gate"], lw=1.4, zorder=2))
    ax.text(gcx, gy + 1.8, "Gate G2: ≥ 2 ring classes, ≥ 10k particles?", fontsize=FS, fontweight="bold", ha="center", va="center", zorder=3)
    ax.text(gcx, gy, f"blind 3-reviewer consensus, {v('n_reviewed')} classes", fontsize=FS, ha="center", va="center",
            zorder=3, color=vizstyle.INK_2)
    ax.text(gcx, gy - 1.8, f"closed-loop classes: {v('closed_loop')}  →  NO", fontsize=FS, ha="center", va="center",
            zorder=3, fontweight="bold")
    arrow(ax, rod[-1].bot, (gcx, gy + gh / 2), color=C["rod"])
    path(ax, [anc[-1].bot, (anc[-1].bot[0], gy), (gx, gy)], color=C["anc"])
    nr = Box(ax, XE, None, WE, "Planned chimera 3D", ["ab initio → hetero/NU →", "class_3D + Topaz:", "NOT RUN"],
             C["exit"], kind="notrun", yc=gy)
    path(ax, [(gx + gw, gy), (XE, gy)], color=C["exit"], dashed=True)
    ax.text(gx + gw + 1.0, gy + 0.5, "yes (never met)", fontsize=FS - 0.4, color=C["exit"], ha="left", va="bottom")

    # "no" bus -> fallback deliverables
    ynb = gy - gh / 2 - 2.6
    path(ax, [(gcx, gy - gh / 2), (gcx, ynb)], color=C["gate"], head=False)
    drops = (XA + WA - 2.5, XR + WR - 2.5, XB + WB - 2.5)
    path(ax, [(drops[0], ynb), (drops[-1], ynb)], color=C["gate"], head=False)
    ax.text(gcx + 0.6, ynb + 0.5, "NO → fallback (agreed in plan)", fontsize=FS - 0.2, color=C["gate"],
            ha="left", va="bottom", fontweight="bold")
    ytop_br = ynb - 4.2
    for xc in drops:
        arrow(ax, (xc, ynb), (xc, ytop_br - 0.1), color=C["gate"])
    head(XR, ytop_br + 0.5, "W3: rod_no_dot (\u201cfree DNA\u201d) control", C["W3"])
    head(XB, ytop_br + 0.5, "A1: window-free 2D", C["A1"])
    head(XW, ytop_br + 0.5, "W4: dot-on-rod 3D", C["W4"])

    # census
    u("reviewed_ptcl_rods", "reviewed_ptcl_anchors", "n_reviewed_rods", "n_reviewed_anchors")
    Box(ax, XA, ytop_br - 0.1, WA, "Census + 2D gallery", [f"{fmt(v('reviewed_ptcl_rods'))} rod + {fmt(v('reviewed_ptcl_anchors'))}",
                                                          f"anchor particles in {v('n_reviewed_rods')} + {v('n_reviewed_anchors')}",
                                                          "reviewed classes (no 3D)"], C["gate"], kind="final")

    # ---------------- W3
    u("rod_no_dot_classes", "W3_duplex", "W3_junk", "W3_ice", "W3_pairs", "W3_pool", "W3_het_duplex", "W3_het_junk",
      "W3_het_ice", "W3_het_pairs", "W3_nu", "W3_cs_tight", "W3_cs_nomask", "W3_rel_csmask", "W3_rel_loose", "W3_cfar",
      "W3_scf", "W3_sym", "rod_no_dot_top", "W3_final_by_label", "W3_final_by_label_n")
    b3, n3 = v("W3_final_by_label"), v("W3_final_by_label_n")
    pc = lambda n, tot: f"{100 * n / tot:.0f}%"  # noqa: E731
    s3 = Box(ax, XR, ytop_br - 0.1, WR, f"Select J29: rod_no_dot ({v('rod_no_dot_classes')} classes)",
             [f"{fmt(v('W3_duplex'))} particles (from J27)"], C["W3"])
    sh = 4.3
    sy = s3.y - 1.0 - sh
    class_avg_strip(fig, glob.glob(os.path.join(DATA, "exports", "W2_2D_rods", "*class_averages*.mrc"))[0],
                    v("rod_no_dot_top"), XR, sy, WR, sh, C["W3"])
    ax.text(XR, sy - 0.25, "six largest rod_no_dot averages (J27)", fontsize=FS - 1.0, color=vizstyle.INK_MUTED, va="top")
    w3 = stack(ax, XR, sy - 2.3, WR, [
        ("Ab initio K = 1: J34 target", ["+ 3 decoy ab initios J35–J37"]),
        (f"Hetero refine J38 ({fmt(v('W3_pool'))})", [f"class seeded by rod_no_dot ab initio (J34): {fmt(v('W3_het_duplex'))}",
                                                      "composition mixed (report §13b)"]),
        (f"NU refine J39 ({v('W3_sym')})", [f"{fmt(v('W3_nu'))} particles",
                                            f"{pc(b3['rod_no_dot'], n3)} rod_no_dot, {pc(n3 - b3['rod_no_dot'], n3)} decoy-input"]),
        ("Final map W3 (J40 / RELION): rod + terminal lobe", [f"GSFSC 0.143: {v('W3_cs_tight'):.1f} Å tight, {v('W3_cs_nomask'):.1f} Å no mask",
                                         f"RELION: {v('W3_rel_csmask'):.1f} Å CS mask, {v('W3_rel_loose'):.1f} Å loose",
                                         f"cFAR {v('W3_cfar'):.2f}, SCF* {v('W3_scf'):.2f}"], "final"),
    ], C["W3"], gap=1.6)
    # selection -> ab initio, bypassing the strip on the left
    xl = XR - 0.9
    path(ax, [s3.left(0.5), (xl, s3.left(0.5)[1]), (xl, w3[0].left(0.5)[1]), w3[0].left(0.5)], color=C["W3"])
    side_box(ax, w3[0], XA, WA - 1.6, "Decoy inputs J30–J32 (from J27)",
             [f"junk {fmt(v('W3_junk'))} · ice {fmt(v('W3_ice'))} · pairs {fmt(v('W3_pairs'))}"], into=True)
    side_box(ax, w3[1], XA, WA - 1.6, "Other hetero classes (side exit)",
             [f"class seeded by junk ab initio (J35): {fmt(v('W3_het_junk'))}",
              f"class seeded by ice ab initio (J36): {fmt(v('W3_het_ice'))}",
              f"class seeded by pairs ab initio (J37): {fmt(v('W3_het_pairs'))}",
              "composition mixed (report §13b)"], f=0.35, fs=FS - 0.6)
    th = 9.5
    thumb(fig, os.path.join(THUMBS, "W3_v2.png"), XA, w3[-1].y - 0.5, WA - 1.6, th)
    ax.text(XA, w3[-1].y - 1.0, "J39 map, 65,000 Å³ contour (thumbnails not to common scale)", fontsize=FS - 1.0, color=vizstyle.INK_MUTED, va="top")
    arrow(ax, w3[-1].left(0.5), (XA + WA - 1.4, w3[-1].left(0.5)[1]), color=C["W3"])

    # ---------------- A1
    u("A1_recentred", "A1_in_dot_long_rod", "A1_in_dot_short_rod", "A1_in_rod_no_dot", "A1_extract", "A1_extract_mics",
      "A1_2D", "A1_2D_K", "A1_2D_window", "W4_dot", "W4_dot_classes", "W4_free", "W4_free_classes", "W4_haze", "W4_ice",
      "A1_dot_top")
    a1 = stack(ax, XB, ytop_br - 0.1, WB, [
        ("Re-centre J43 (from J27)", [f"{fmt(v('A1_recentred'))} rods: dot_long_rod",
                                      f"{fmt(v('A1_in_dot_long_rod'))}, dot_short_rod {fmt(v('A1_in_dot_short_rod'))},",
                                      f"rod_no_dot {fmt(v('A1_in_rod_no_dot'))}"]),
        ("Extract J44 (454 Å box)", [f"{fmt(v('A1_extract'))} ({fmt(v('A1_extract_mics'))} micrographs)"]),
        (f"2D J45 (K = {v('A1_2D_K')}, {v('A1_2D_window'):.0f} Å window)", [f"{fmt(v('A1_2D'))} kept;",
                                                                            "still no ring or hook class"]),
        ("Groups (A1 analysis)", [f"DOT {fmt(v('W4_dot'))} ({v('W4_dot_classes')} classes)",
                                  f"FREE {fmt(v('W4_free'))} ({v('W4_free_classes')} classes)",
                                  f"haze {fmt(v('W4_haze'))}, ice {fmt(v('W4_ice'))}"]),
    ], C["A1"], gap=1.6)
    sy4 = a1[-1].y - 1.0 - sh
    class_avg_strip(fig, glob.glob(os.path.join(DATA, "exports", "A1_2D", "*class_averages*.mrc"))[0], v("A1_dot_top"),
                    XB, sy4, WB, sh, C["A1"])
    ax.text(XB, sy4 - 0.25, "six largest DOT averages (J45)", fontsize=FS - 1.0, color=vizstyle.INK_MUTED, va="top")

    # ---------------- W4
    u("W4_pool", "W4_het_dot", "W4_het_free", "W4_het_haze", "W4_het_ice", "W4_nu", "W4_cs_tight", "W4_cs_nomask",
      "W4_rel_csmask", "W4_rel_loose", "W4_cfar", "W4_scf", "W4_sym", "W4_c3d", "W4_c3d_K", "W4_c3d_min", "W4_c3d_max",
      "W4_c3d_res", "W4_c3d_pcs", "W4_final_by_label", "W4_final_by_label_n")
    b4, n4 = v("W4_final_by_label"), v("W4_final_by_label_n")
    w4 = stack(ax, XW, ytop_br - 0.1, WW, [
        ("Ab initio K = 1: J51 DOT", ["+ J52 FREE (real competitor)", "+ decoys J53 haze, J54 ice"]),
        (f"Hetero refine J55 ({fmt(v('W4_pool'))})", [f"class seeded by DOT ab initio (J51): {fmt(v('W4_het_dot'))}",
                                                      "composition mixed (report §13b)"]),
        (f"NU refine J56 ({v('W4_sym')})", [f"{fmt(v('W4_nu'))} particles",
                                            f"{pc(b4['DOT (dot-on-rod)'], n4)} DOT, {pc(b4['haze (decoy)'], n4)} haze"]),
        ("Final map W4 (J57 / RELION)", [f"GSFSC 0.143: {v('W4_cs_tight'):.1f} Å tight, {v('W4_cs_nomask'):.1f} Å no mask",
                                         f"RELION: {v('W4_rel_csmask'):.1f} Å CS mask, {v('W4_rel_loose'):.1f} Å loose",
                                         f"cFAR {v('W4_cfar'):.2f}, SCF* {v('W4_scf'):.2f}"], "final"),
        (f"class_3D J59 (PCA {v('W4_c3d_pcs')} comp., K = {v('W4_c3d_K')})", [f"{fmt(v('W4_c3d'))} particles at {v('W4_c3d_res'):.0f} Å,",
                                                                               f"{fmt(v('W4_c3d_min'))}–{fmt(v('W4_c3d_max'))} per class"]),
    ], C["W4"], gap=1.6)
    path(ax, [a1[-1].right(0.5), (XW - 1.2, a1[-1].right(0.5)[1]), (XW - 1.2, w4[0].left(0.5)[1]), w4[0].left(0.5)],
         color=C["A1"])
    ax.text(XW - 1.0, w4[0].left(0.5)[1] - 1.0, "all four\ngroups", fontsize=FS - 1.0, color=C["A1"], ha="left", va="top") \
        if False else None
    side_box(ax, w4[1], XE, WE, "Other hetero classes (exit)",
             [f"class seeded by FREE ab initio (J52):", f"  {fmt(v('W4_het_free'))}",
              f"class seeded by haze ab initio (J53):", f"  {fmt(v('W4_het_haze'))}",
              f"class seeded by ice ab initio (J54):", f"  {fmt(v('W4_het_ice'))}",
              "composition mixed (report §13b)"], fs=FS - 0.6)
    thumb(fig, os.path.join(THUMBS, "W4_v2.png"), XE, w4[3].y - 4.0, WE, 12.0)
    ax.text(XE, w4[3].y - 4.5, "J56 map, 65,000 Å³ contour", fontsize=FS - 1.0, color=vizstyle.INK_MUTED, va="top")
    arrow(ax, w4[3].right(0.5), (XE + 0.6, w4[3].right(0.5)[1]), color=C["W4"])
    return fig, ax, used



def main():
    os.makedirs(FIGDIR, exist_ok=True)
    os.makedirs(RESDIR, exist_ok=True)
    collect()
    rows, src = build_table()
    write_table(rows, src)
    fig, ax, used = make_figure()
    out = os.path.join(FIGDIR, "figS6_workflow.png")
    fig.savefig(out, dpi=200)
    plt.close(fig)
    json.dump({"figure": out, "script": os.path.join(REPO, "scripts", "50_figS6_tableS1.py"),
               "numbers": {k: N[k] for k in sorted(used)}}, open(os.path.join(FIGDIR, "figS6_workflow.json"), "w"),
              indent=1, default=str)
    json.dump(N, open(os.path.join(RESDIR, "figS6_tableS1_numbers.json"), "w"), indent=1, default=str)
    for k in sorted(N):
        print(f"{k:28s} {N[k]['value']}")


if __name__ == "__main__":
    sys.exit(main())
