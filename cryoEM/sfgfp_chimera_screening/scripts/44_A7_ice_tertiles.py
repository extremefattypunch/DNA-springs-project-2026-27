"""A7: ice-thickness split -- an air-water-interface (AWI) damage test.  (CPU only; real coordinates only.)

Hypothesis 4 of results/gel_em_reconciliation.json: closed DNA rings exist in solution but open or deform at the
AWI. If so, thicker ice (a smaller AWI-adsorbed fraction) should show more DNA-attached sfGFP, longer attached DNA or
ring-like class averages. This script splits the curated micrographs into tertiles of CryoSPARC Patch CTF
ctf_stats/ice_thickness_rel (cut points over all 8,036 curated micrographs; the top decile is carried as an extra
sensitivity stratum) and compares, per stratum:

 (a) census of consensus 2D labels, rod stream (W2_2D_rods, all 8,036 micrographs) and anchor stream (W2_2D_anchors,
     the 2,009 micrographs it was picked on): pooled particle fraction per label and per-micrograph count per um^2
     (full 2048^2 x 1.42 A field = 0.0846 um^2);
 (b) detected dot density (dots per um^2 of the detector's valid mask; dots from scripts/17 'detect');
 (c) attached-DNA proxies: (i) 2h of dot_long_rod rods from rods_analysed.npz (class-average geometry, so it is a
     class-mix readout, and window-limited: the 200 A 2D window clips anything longer); (ii) per micrograph, the number
     of ELIGIBLE detected dots (script-17 rule: not inside the middle region of any real rod-set rod) whose centre is
     within 35 A of a rod end point (centre +/- (h + offset) u, offset 0 A primary, 15 A secondary), using the
     rods_analysed geometry, against a shift null (each rod moved 150-300 A in a random direction on the same
     micrograph, valid-mask redraws as in script 17; NREP replicates). The same middle-region rule is applied to the
     shifted (virtual) rod. Reported as obs/null ratio and excess (obs - null) per um^2 and per detected dot.
 (d) ring-like averages: per-stratum reference-free averages of every 2D class, made by applying each particle's own
     2D pose/shift (script-17 empirical convention: rotate by -theta, subtract the shift, particle frame) to the
     phase-flipped (own per-particle CTF) extracted images, Fourier-cropped to the CryoSPARC class-average grid and
     windowed like CryoSPARC (circular window of class2D_window_inner_A, 10 A cosine edge); the hole test is
     describe() of scripts/41_class2d_morphology.py (imported, unmodified: >= 30 A enclosed background hole in the
     z > 2 mask), plus a lenient variant (z > 1.5 mask, hole >= 20 A). Averages of all strata pooled are correlated
     with the CryoSPARC class averages as a convention check.

Uncertainty: 95% bootstrap CIs resampling FOIL HOLES (clusters of ~3 neighbouring exposures that share ice), with
ratio-of-sums statistics; T3 - T1 differences use independent resamples of the two strata. Confounds: Spearman rho
of ice_thickness_rel with defocus, astigmatism, CTF fit resolution, acquisition time, contamination-mask fraction,
dot and particle densities, plus OLS of the main per-micrograph readouts on z-scored ice, CTF fit and defocus.

Run (logs to /n/netscratch/hekstra_lab/Lab/ian_poon/logs/):
  PY=/n/holylabs/hekstra_lab/Everyone/ianpoon/envs/dnaspring/bin/python
  sbatch -p test -c 8 --mem=32G -t 4:00:00 --wrap "$PY -I scripts/44_A7_ice_tertiles.py main"
  sbatch -p test -c 8 --mem=32G -t 4:00:00 --array=0-3 \
      --wrap "$PY -I scripts/44_A7_ice_tertiles.py avg --stream rods --task \$SLURM_ARRAY_TASK_ID --ntasks 4"
  sbatch -p test -c 8 --mem=32G -t 4:00:00 --array=0-3 \
      --wrap "$PY -I scripts/44_A7_ice_tertiles.py avg --stream anchors --task \$SLURM_ARRAY_TASK_ID --ntasks 4"
  sbatch -p test -c 8 --mem=32G -t 4:00:00 --wrap "$PY -I scripts/44_A7_ice_tertiles.py final"
Outputs: results/A7_ice_tertiles.json, figures/qc/A7_ice_tertiles.png; intermediates under WORK.
"""
import argparse
import importlib.util
import json
import os
import re
import time
from multiprocessing import Pool

import mrcfile
import numpy as np
import pandas as pd
from scipy import ndimage as ndi
from scipy.spatial import cKDTree
from scipy.stats import spearmanr

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def _load(name, fn):
    spec = importlib.util.spec_from_file_location(name, os.path.join(REPO, "scripts", fn))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


s17 = _load("s17", "17_rod_end_colocalisation.py")
s41 = _load("s41", "41_class2d_morphology.py")
g0 = s17.g0

PROJ = s17.PROJ
EXPORTS = s17.EXPORTS
CUR_NPZ = s17.CUR_NPZ
RE_WORK = s17.WORK
WORK = "/n/netscratch/hekstra_lab/Lab/ian_poon/sfgfp_chimera_screening_20261112/A7_ice_tertiles"
OUT_JSON = os.path.join(REPO, "results", "A7_ice_tertiles.json")
FIG = os.path.join(REPO, "figures", "qc", "A7_ice_tertiles.png")
STREAMS = {
    "rods": dict(npz=os.path.join(EXPORTS, "W2_2D_rods__particles.npz"),
                 ca=os.path.join(EXPORTS, "W2_2D_rods", "J27_030_class_averages.mrc"),
                 labels=os.path.join(REPO, "results", "class2d_labels_W2_2D_rods.csv"), window_A=200.0),
    "anchors": dict(npz=os.path.join(EXPORTS, "W2_2D_anchors__particles.npz"),
                    ca=os.path.join(EXPORTS, "W2_2D_anchors", "J28_030_class_averages.mrc"),
                    labels=os.path.join(REPO, "results", "class2d_labels_W2_2D_anchors.csv"), window_A=220.0),
}
A_MIC, A_P = 1.42, 2.84
MIC_N = 2048
MIC_AREA_UM2 = MIC_N * MIC_N * A_MIC ** 2 / 1e8
GROUPS = ["T1_thin", "T2_mid", "T3_thick", "top10_thickest"]
NG = len(GROUPS)
NBOOT = 2000
NREP = 40
END_D_A = 35.0
END_OFF = [0.0, 15.0]
NMIN_AVG = 300             # minimum particles for a per-stratum class average to be tested
SEED = 20261114


# ----------------------------------------------------------------------------------------------- shared
def strata():
    """Per curated micrograph: path, ice, tertile index (0..2), top-decile flag, foil-hole cluster, confounds."""
    c = np.load(CUR_NPZ)
    path = c["micrograph_blob__path"].astype(str)
    ice = c["ctf_stats__ice_thickness_rel"].astype(float)
    q = np.percentile(ice, [100 / 3, 200 / 3, 90])
    tert = np.digitize(ice, q[:2])
    top = ice >= q[2]
    hole = np.array([re.search(r"FoilHole_(\d+)", p).group(1) for p in path])
    tt = []
    for p in path:
        m = re.search(r"_(\d{8})_(\d{6})_Fractions", p)
        tt.append(pd.Timestamp(m.group(1) + m.group(2)).value / 3.6e12 if m else np.nan)   # hours
    tt = np.array(tt)
    df = 0.5 * (c["ctf__df1_A"] + c["ctf__df2_A"]).astype(float)
    ast = np.abs(c["ctf__df1_A"] - c["ctf__df2_A"]).astype(float)
    return dict(path=path, ice=ice, cuts=q, tert=tert, top=top, hole=hole, time_h=tt - np.nanmin(tt), defocus_A=df,
                astig_A=ast, ctf_fit_A=c["ctf__ctf_fit_to_A"].astype(float))


def group_masks(S):
    return [S["tert"] == 0, S["tert"] == 1, S["tert"] == 2, S["top"]]


def boot_ratio(num, den, clus, rng, B=NBOOT):
    """Ratio-of-sums statistic sum(num)/sum(den) with a cluster bootstrap. num: (n,) or (n, k). Returns est, lo, hi,
    and the bootstrap draws (B,) or (B, k)."""
    num = np.asarray(num, float)
    one = num.ndim == 1
    if one:
        num = num[:, None]
    u, inv = np.unique(clus, return_inverse=True)
    cn = np.zeros((len(u), num.shape[1]))
    np.add.at(cn, inv, num)
    cd = np.bincount(inv, weights=np.asarray(den, float), minlength=len(u))
    W = rng.multinomial(len(u), np.full(len(u), 1.0 / len(u)), size=B).astype(float)
    bn = W @ cn
    bd = W @ cd
    draws = bn / bd[:, None]
    est = cn.sum(0) / cd.sum()
    lo, hi = np.percentile(draws, [2.5, 97.5], axis=0)
    if one:
        return float(est[0]), float(lo[0]), float(hi[0]), draws[:, 0]
    return est, lo, hi, draws


def ci_dict(est, lo, hi, nd=5):
    return dict(est=round(float(est), nd), lo95=round(float(lo), nd), hi95=round(float(hi), nd))


def diff_ci(d3, d1):
    d = d3 - d1
    return dict(est=None, lo95=float(np.percentile(d, 2.5)), hi95=float(np.percentile(d, 97.5)),
                p_two_sided=float(min(1.0, 2 * min(np.mean(d <= 0), np.mean(d >= 0)))))


def read_labels(fn, ncls):
    lab = pd.read_csv(fn)
    arr = np.full(ncls, "unlabelled", dtype=object)
    arr[lab["class"].values] = lab["consensus_label"].values
    return arr


def load_dots():
    files = sorted(f for f in os.listdir(os.path.join(RE_WORK, "dots")) if f.endswith(".npz"))
    parts = [np.load(os.path.join(RE_WORK, "dots", f)) for f in files]
    nm = int(max(p["mic_index"].max() for p in parts)) + 1
    mpath = np.empty(nm, dtype=object)
    m4p = np.zeros((nm, MIC_N // g0.BIN, MIC_N // g0.BIN // 8), np.uint8)
    nd = np.zeros(nm, int)
    frac_bad = np.full(nm, np.nan)
    xs, ys, zs, ms = [], [], [], []
    for p in parts:
        mi = p["mic_index"]
        mpath[mi] = p["mic_path"]
        m4p[mi] = p["m4p"]
        nd[mi] = p["n"]
        frac_bad[mi] = p["frac_bad"]
        xs.append(p["x"].astype(float))
        ys.append(p["y"].astype(float))
        zs.append(p["z"])
        ms.append(np.repeat(mi, p["n"]))
    return mpath.astype(str), m4p, nd, frac_bad, np.concatenate(xs), np.concatenate(ys), np.concatenate(zs), \
        np.concatenate(ms)


# ----------------------------------------------------------------------------------------------- (c) rod-end dots
G = {}


def end_dots(X, Y, ux, uy, h):
    """Per micrograph: (# eligible dots within END_D_A of a rod end, # rods with such a dot) for each END_OFF."""
    mic = G["mic"]
    rmax = h + (max(END_OFF) + END_D_A) / A_MIC + 1
    out = {}
    lists = G["tree"].query_ball_point(np.c_[X + mic * s17.OFF, Y], r=rmax)
    lens = np.fromiter((len(v) for v in lists), int, len(lists))
    dot = np.fromiter((j for v in lists for j in v), int, int(lens.sum()))
    rod = np.repeat(np.arange(len(X)), lens)
    rx, ry = G["dx"][dot] - X[rod], G["dy"][dot] - Y[rod]
    s = rx * ux[rod] + ry * uy[rod]
    t = -rx * uy[rod] + ry * ux[rod]
    hr = h[rod]
    mid = (np.abs(t) <= s17.MID_T_A / A_MIC) & (np.abs(s) <= hr - s17.MID_INSET_A / A_MIC)
    ok = G["elig"][dot] & ~mid
    nmic = G["nmic"]
    nodot = G["nodot"][rod]
    dotc = ~nodot
    for eo in END_OFF:
        e = hr + eo / A_MIC
        dp = np.hypot(s - e, t) * A_MIC
        dmn = np.hypot(s + e, t) * A_MIC
        # end subsets: all ends; 'light' = rod_no_dot ends + far (-axis) end of dot classes (the ends the 2D
        # classifier did not select on); and the three components separately
        subsets = {"all": (dp <= END_D_A) | (dmn <= END_D_A),
                   "nodot": nodot & ((dp <= END_D_A) | (dmn <= END_D_A)),
                   "far": dotc & (dmn <= END_D_A),
                   "near": dotc & (dp <= END_D_A)}
        subsets["light"] = subsets["nodot"] | subsets["far"]
        for sn, sm in subsets.items():
            hit = ok & sm
            ud = np.unique(dot[hit])
            out[f"dots_off{int(eo)}_{sn}"] = np.bincount(G["dmic"][ud], minlength=nmic)
            ur = np.unique(rod[hit])
            out[f"rods_off{int(eo)}_{sn}"] = np.bincount(mic[ur], minlength=nmic)
    return out


def null_worker(seed):
    rng = np.random.default_rng(seed)
    X, Y, ux, uy, nfail = s17.null_geometry("shift", rng)
    return end_dots(X, Y, ux, uy, G["h"]), nfail


# ----------------------------------------------------------------------------------------------- stage main
def stage_main(a):
    t0 = time.time()
    os.makedirs(WORK, exist_ok=True)
    rng = np.random.default_rng(SEED)
    S = strata()
    lut = {p: i for i, p in enumerate(S["path"])}
    gm = group_masks(S)
    res = dict(created=time.strftime("%Y-%m-%d %H:%M"), n_curated=int(len(S["path"])),
               ice_cuts=dict(t33=float(S["cuts"][0]), t67=float(S["cuts"][1]), p90=float(S["cuts"][2])),
               groups={g: dict(n_mic=int(m.sum()), n_foilholes=int(len(np.unique(S["hole"][m]))),
                               ice_median=float(np.median(S["ice"][m])),
                               ice_range=[float(S["ice"][m].min()), float(S["ice"][m].max())],
                               defocus_um_median=float(np.median(S["defocus_A"][m]) / 1e4),
                               ctf_fit_A_median=float(np.median(S["ctf_fit_A"][m])),
                               astig_A_median=float(np.median(S["astig_A"][m])))
                       for g, m in zip(GROUPS, gm)})
    print("groups", json.dumps(res["groups"], indent=0), flush=True)

    # ---------------- (a) census
    permic = pd.DataFrame(dict(path=S["path"], ice=S["ice"], tert=S["tert"], top=S["top"], hole=S["hole"],
                               defocus_A=S["defocus_A"], astig_A=S["astig_A"], ctf_fit_A=S["ctf_fit_A"],
                               time_h=S["time_h"]))
    census = {}
    for st, cfg in STREAMS.items():
        P = np.load(cfg["npz"])
        cls = P["alignments2D__class"].astype(int)
        labs = read_labels(cfg["labels"], 200)
        names = sorted(set(labs[cls]))
        li = np.array([names.index(x) for x in labs[cls]])
        mi = np.array([lut[p] for p in P["location__micrograph_path"].astype(str)])
        cnt = np.zeros((len(S["path"]), len(names)))
        np.add.at(cnt, (mi, li), 1)
        have = np.zeros(len(S["path"]), bool)
        have[mi] = True
        for j, nm in enumerate(names):
            permic[f"{st}:{nm}"] = np.where(have, cnt[:, j], np.nan)
        permic[f"{st}:total"] = np.where(have, cnt.sum(1), np.nan)
        tot = cnt.sum(1)
        C = dict(labels=names, n_particles=int(len(cls)), n_mic=int(have.sum()), groups={})
        draws_f, draws_d = {}, {}
        for g, m in zip(GROUPS, gm):
            mm = m & have
            ef, lf, hf, df_ = boot_ratio(cnt[mm], tot[mm], S["hole"][mm], rng)
            ed, ld, hd, dd_ = boot_ratio(cnt[mm], np.full(mm.sum(), MIC_AREA_UM2), S["hole"][mm], rng)
            draws_f[g], draws_d[g] = df_, dd_
            et, lt, ht, _ = boot_ratio(tot[mm], np.full(mm.sum(), MIC_AREA_UM2), S["hole"][mm], rng)
            C["groups"][g] = dict(
                n_mic=int(mm.sum()), n_particles=int(tot[mm].sum()),
                particles_per_um2=ci_dict(et, lt, ht, 1),
                fraction={nm: ci_dict(ef[j], lf[j], hf[j]) for j, nm in enumerate(names)},
                per_um2={nm: ci_dict(ed[j], ld[j], hd[j], 2) for j, nm in enumerate(names)})
        C["T3_minus_T1"] = {nm: dict(fraction=diff_ci(draws_f["T3_thick"][:, j], draws_f["T1_thin"][:, j]),
                                     per_um2=diff_ci(draws_d["T3_thick"][:, j], draws_d["T1_thin"][:, j]))
                            for j, nm in enumerate(names)}
        for nm in names:
            for k in ("fraction", "per_um2"):
                C["T3_minus_T1"][nm][k]["est"] = (C["groups"]["T3_thick"][k][nm]["est"]
                                                  - C["groups"]["T1_thin"][k][nm]["est"])
        census[st] = C
        print(f"census {st}: {len(cls)} particles on {have.sum()} micrographs ({time.time() - t0:.0f} s)", flush=True)
    res["census"] = census

    # ---------------- (b) dots
    mpath, m4p, nd, frac_bad, dx, dy, dz, dm = load_dots()
    dlut = np.array([lut[p] for p in mpath])             # dot-file mic index -> curated index
    varea = np.array([np.unpackbits(m4p[i], axis=1).sum() for i in range(len(mpath))]) * 16 * A_MIC ** 2 / 1e8
    n_cur = len(S["path"])
    nd_c = np.full(n_cur, np.nan)
    va_c = np.full(n_cur, np.nan)
    fb_c = np.full(n_cur, np.nan)
    nd_c[dlut], va_c[dlut], fb_c[dlut] = nd, varea, frac_bad
    permic["dots"] = nd_c
    permic["valid_um2"] = va_c
    permic["contam_frac"] = fb_c
    permic["dots_per_um2"] = nd_c / va_c
    hd = np.isfinite(nd_c) & (va_c > 0)
    D = dict(n_mic=int(hd.sum()), n_dots=int(np.nansum(nd_c)), groups={})
    dd = {}
    for g, m in zip(GROUPS, gm):
        mm = m & hd
        e, lo, hi, dr = boot_ratio(nd_c[mm], va_c[mm], S["hole"][mm], rng)
        dd[g] = dr
        D["groups"][g] = dict(n_mic=int(mm.sum()), dots_per_um2_valid=ci_dict(e, lo, hi, 1),
                              valid_area_frac_median=float(np.median(va_c[mm] / MIC_AREA_UM2)))
    D["T3_minus_T1"] = diff_ci(dd["T3_thick"], dd["T1_thin"])
    D["T3_minus_T1"]["est"] = D["groups"]["T3_thick"]["dots_per_um2_valid"]["est"] - \
        D["groups"]["T1_thin"]["dots_per_um2_valid"]["est"]
    res["dots"] = D
    print("dots done", flush=True)

    # ---------------- (c) attached-rod proxies
    R = np.load(os.path.join(RE_WORK, "rods_analysed.npz"))
    X, Y, ux, uy, h, rmic, rlab = (R[k] for k in ("X", "Y", "ux", "uy", "h", "mic", "lab"))
    rlab = rlab.astype(str)
    rcur = dlut[rmic]
    Rres = dict(n_rods=int(len(X)), groups={})
    L_A = 2 * h * A_MIC
    for g, m in zip(GROUPS, gm):
        sel = m[rcur]
        dl = sel & (rlab == "dot_long_rod")
        rn = sel & (rlab == "rod_no_dot")
        Rres["groups"][g] = dict(
            n_rods=int(sel.sum()), n_dot_long_rod=int(dl.sum()),
            dot_long_rod_2h_A=dict(mean=float(L_A[dl].mean()), median=float(np.median(L_A[dl])),
                                   p10=float(np.percentile(L_A[dl], 10)), p90=float(np.percentile(L_A[dl], 90)),
                                   frac_ge_150A=float(np.mean(L_A[dl] >= 150))),
            rod_no_dot_2h_A=dict(mean=float(L_A[rn].mean()), median=float(np.median(L_A[rn]))),
            frac_dot_long_rod_of_R=float(dl.sum() / sel.sum()))
    # boot for dot_long_rod 2h mean
    dl_all = rlab == "dot_long_rod"
    bd = {}
    for g, m in zip(GROUPS, gm):
        sel = m[rcur] & dl_all
        e, lo, hi, dr = boot_ratio(L_A[sel], np.ones(sel.sum()), S["hole"][rcur[sel]], rng)
        bd[g] = dr
        Rres["groups"][g]["dot_long_rod_2h_A"]["mean_ci"] = ci_dict(e, lo, hi, 2)
    Rres["dot_long_rod_2h_T3_minus_T1"] = diff_ci(bd["T3_thick"], bd["T1_thin"])
    Rres["dot_long_rod_2h_T3_minus_T1"]["est"] = Rres["groups"]["T3_thick"]["dot_long_rod_2h_A"]["mean"] - \
        Rres["groups"]["T1_thin"]["dot_long_rod_2h_A"]["mean"]

    # end-dot counts, observed and shift null
    G.update(tree=cKDTree(np.c_[dx + dm * s17.OFF, dy]), dx=dx, dy=dy, dmic=dm, nmic=len(mpath), mic=rmic, h=h,
             X=X, Y=Y, ux=ux, uy=uy, nodot=(rlab == "rod_no_dot"))
    s17.G.update(m4p=m4p, X=X, Y=Y, ux=ux, uy=uy, h=h, mic=rmic)
    lists = G["tree"].query_ball_point(np.c_[X + rmic * s17.OFF, Y], r=h + 1)
    lens = np.fromiter((len(v) for v in lists), int, len(lists))
    dot = np.fromiter((j for v in lists for j in v), int, int(lens.sum()))
    rod = np.repeat(np.arange(len(X)), lens)
    rx, ry = dx[dot] - X[rod], dy[dot] - Y[rod]
    s = rx * ux[rod] + ry * uy[rod]
    t = -rx * uy[rod] + ry * ux[rod]
    onrod = np.zeros(len(dx), bool)
    onrod[dot[(np.abs(t) <= s17.MID_T_A / A_MIC) & (np.abs(s) <= h[rod] - s17.MID_INSET_A / A_MIC)]] = True
    G["elig"] = ~onrod
    obs = end_dots(X, Y, ux, uy, h)
    print(f"observed end dots (off0): {obs['dots_off0_all'].sum()} (light {obs['dots_off0_light'].sum()}); {onrod.sum()} dots on rod middles excluded "
          f"({time.time() - t0:.0f} s)", flush=True)
    procs = int(os.environ.get("SLURM_CPUS_PER_TASK", 4))
    with Pool(procs) as pool:
        nres = pool.map(null_worker, [SEED * 7 + r for r in range(a.nrep)], chunksize=1)
    print(f"{a.nrep} shift-null replicates done ({time.time() - t0:.0f} s)", flush=True)
    nfail = float(np.mean([r[1] for r in nres]))
    nrods_mic = np.bincount(rmic, minlength=len(mpath))
    ndot_elig = np.bincount(dm[~onrod], minlength=len(mpath))
    for k in obs:
        nul = np.mean([r[0][k] for r in nres], axis=0)
        a_o = np.full(n_cur, np.nan)
        a_n = np.full(n_cur, np.nan)
        a_o[dlut], a_n[dlut] = obs[k], nul
        permic[f"end_{k}_obs"] = a_o
        permic[f"end_{k}_null"] = a_n
    tmp = np.full(n_cur, np.nan)
    tmp[dlut] = nrods_mic
    permic["R_rods"] = tmp
    tmp = np.full(n_cur, np.nan)
    tmp[dlut] = ndot_elig
    permic["dots_eligible"] = tmp
    for nm_, msk in (("R_rods_nodot", rlab == "rod_no_dot"), ("R_rods_dotcls", rlab != "rod_no_dot")):
        tmp = np.full(n_cur, np.nan)
        tmp[dlut] = np.bincount(rmic[msk], minlength=len(mpath))
        permic[nm_] = tmp
    rden = {"all": "R_rods", "light": "R_rods", "nodot": "R_rods_nodot", "far": "R_rods_dotcls",
            "near": "R_rods_dotcls"}
    E = dict(definition=f"eligible detected dots within {END_D_A:.0f} A of a rods_analysed rod end point "
                        "(centre +/- (h + offset) u); shift null = each rod moved 150-300 A, NREP replicates",
             nrep=a.nrep, null_unplaced_rods_mean=nfail, n_dots_on_rod_middles_excluded=int(onrod.sum()))
    for k in obs:
        E[k] = dict(groups={})
        dr_ratio, dr_exc = {}, {}
        for g, m in zip(GROUPS, gm):
            mm = m & hd
            o, nu, ar = permic[f"end_{k}_obs"].values[mm], permic[f"end_{k}_null"].values[mm], va_c[mm]
            cl = S["hole"][mm]
            er, lr, hr, drr = boot_ratio(o, nu, cl, rng)
            ee, le, he, dre = boot_ratio(o - nu, ar, cl, rng)
            dr_ratio[g], dr_exc[g] = drr, dre
            den = permic["dots_eligible"].values[mm] if k.startswith("dots") else \
                permic[rden[k.split("_")[-1]]].values[mm]
            ef, lf, hf, _ = boot_ratio(o - nu, den, cl, rng)
            eo_, lo_, ho_, _ = boot_ratio(o, den, cl, rng)
            E[k]["groups"][g] = dict(obs=int(o.sum()), null=float(nu.sum()), ratio=ci_dict(er, lr, hr, 4),
                                     excess_per_um2=ci_dict(ee, le, he, 2),
                                     obs_per_unit=ci_dict(eo_, lo_, ho_, 5),
                                     excess_per_unit=ci_dict(ef, lf, hf, 5),
                                     unit="eligible dot" if k.startswith("dots") else "rod")
        E[k]["ratio_T3_minus_T1"] = diff_ci(dr_ratio["T3_thick"], dr_ratio["T1_thin"])
        E[k]["ratio_T3_minus_T1"]["est"] = E[k]["groups"]["T3_thick"]["ratio"]["est"] - \
            E[k]["groups"]["T1_thin"]["ratio"]["est"]
        E[k]["excess_T3_minus_T1"] = diff_ci(dr_exc["T3_thick"], dr_exc["T1_thin"])
        E[k]["excess_T3_minus_T1"]["est"] = E[k]["groups"]["T3_thick"]["excess_per_um2"]["est"] - \
            E[k]["groups"]["T1_thin"]["excess_per_um2"]["est"]
    Rres["end_dots"] = E
    for g in GROUPS:
        m = gm[GROUPS.index(g)] & hd
        Rres["groups"][g]["R_rods_per_um2_valid"] = float(np.nansum(permic["R_rods"].values[m]) / va_c[m].sum())
    res["attached"] = Rres

    # ---------------- confounds
    permic["rods_stream_per_um2"] = permic["rods:total"] / MIC_AREA_UM2
    permic["anchors_stream_per_um2"] = permic["anchors:total"] / MIC_AREA_UM2
    permic["frac_dot_long_rod_rods"] = permic["rods:dot_long_rod"] / permic["rods:total"]
    permic["end_excess_per_um2"] = (permic["end_dots_off0_light_obs"] - permic["end_dots_off0_light_null"]) / \
        permic["valid_um2"]
    permic["end_excess_per_dot"] = (permic["end_dots_off0_light_obs"] - permic["end_dots_off0_light_null"]) / \
        permic["dots_eligible"]
    permic["frac_rod_no_dot_rods"] = permic["rods:rod_no_dot"] / permic["rods:total"]
    permic["frac_dot_short_rod_anchors"] = permic["anchors:dot_short_rod"] / permic["anchors:total"]
    permic["frac_dot_pair_anchors"] = permic["anchors:dot_pair"] / permic["anchors:total"]
    conf = {}
    for col in ["defocus_A", "astig_A", "ctf_fit_A", "time_h", "contam_frac", "dots_per_um2", "rods_stream_per_um2",
                "anchors_stream_per_um2", "frac_dot_long_rod_rods", "frac_rod_no_dot_rods",
                "frac_dot_short_rod_anchors", "frac_dot_pair_anchors", "end_excess_per_um2", "end_excess_per_dot"]:
        v = permic[col].values
        ok = np.isfinite(v)
        r = spearmanr(permic["ice"].values[ok], v[ok])
        conf[col] = dict(spearman_rho=round(float(r.statistic), 4), p=float(r.pvalue), n=int(ok.sum()),
                         group_medians={g: float(np.nanmedian(v[m])) for g, m in zip(GROUPS, gm)})
    # within-foil-hole ice spread (how much of the ice variation is between holes)
    hh = permic.groupby("hole")["ice"]
    conf["ice_between_hole_var_frac"] = float(hh.transform("mean").var() / permic["ice"].var())
    # OLS: readout ~ ice + ctf_fit + defocus (+ dot density = crowding), z-scored, cluster bootstrap over foil holes
    readouts = ["frac_dot_long_rod_rods", "frac_rod_no_dot_rods", "frac_dot_short_rod_anchors",
                "frac_dot_pair_anchors", "end_excess_per_um2", "end_excess_per_dot", "dots_per_um2",
                "rods_stream_per_um2"]
    ols = {}
    for ycol in readouts:
        for model, xs in (("ice_ctf_defocus", ("ice", "ctf_fit_A", "defocus_A")),
                          ("plus_dot_density", ("ice", "ctf_fit_A", "defocus_A", "dots_per_um2"))):
            if ycol == "dots_per_um2" and model == "plus_dot_density":
                continue
            d = permic[list(dict.fromkeys([*xs, ycol, "hole"]))].dropna()
            d = d[np.isfinite(d[ycol])]
            Z = np.c_[np.ones(len(d)), *[(d[c] - d[c].mean()) / d[c].std() for c in xs]]
            yv = d[ycol].values
            beta = np.linalg.lstsq(Z, yv, rcond=None)[0]
            u, inv = np.unique(d["hole"].values, return_inverse=True)
            idx_by = [np.nonzero(inv == i)[0] for i in range(len(u))]
            bb = []
            for _ in range(500):
                pick = rng.integers(0, len(u), len(u))
                ii = np.concatenate([idx_by[i] for i in pick])
                bb.append(np.linalg.lstsq(Z[ii], yv[ii], rcond=None)[0])
            bb = np.array(bb)
            o = {nm: dict(beta_per_sd=float(beta[j]), lo95=float(np.percentile(bb[:, j], 2.5)),
                          hi95=float(np.percentile(bb[:, j], 97.5)))
                 for j, nm in enumerate(["intercept", *xs])}
            o["n"] = int(len(d))
            o["y_sd"] = float(yv.std())
            ols.setdefault(ycol, {})[model] = o
    # crowding-matched T3 - T1: stratify micrographs by detected-dot density decile, compare the tertiles within
    # each density bin (ratio of sums), average over bins weighted by min(n_T1, n_T3); foil-hole cluster bootstrap
    dens = permic["dots_per_um2"].values
    okd = np.isfinite(dens)
    edges = np.nanpercentile(dens, np.linspace(0, 100, 11))
    dbin = np.clip(np.digitize(dens, edges[1:-1]), 0, 9)
    u, inv = np.unique(S["hole"], return_inverse=True)
    Wc = rng.multinomial(len(u), np.full(len(u), 1.0 / len(u)), size=500).astype(float)
    Wm = Wc[:, inv]                                     # (B, n_mic) bootstrap weights
    matched = {}
    pairs_ = {"frac_dot_long_rod_rods": ("rods:dot_long_rod", "rods:total"),
              "frac_rod_no_dot_rods": ("rods:rod_no_dot", "rods:total"),
              "frac_dot_pair_thread_rods": ("rods:dot_pair_thread", "rods:total"),
              "frac_dot_short_rod_anchors": ("anchors:dot_short_rod", "anchors:total"),
              "frac_dot_pair_anchors": ("anchors:dot_pair", "anchors:total"),
              "end_light_ratio_obs_over_null": ("end_dots_off0_light_obs", "end_dots_off0_light_null"),
              "end_light_excess_per_dot": ("__exc_light", "dots_eligible"),
              "end_light_excess_per_um2": ("__exc_light", "valid_um2")}
    permic["__exc_light"] = permic["end_dots_off0_light_obs"] - permic["end_dots_off0_light_null"]
    for nm, (cn, cd) in pairs_.items():
        num, den = permic[cn].values, permic[cd].values
        ok = okd & np.isfinite(num) & np.isfinite(den)
        est, bw, bins_used, boot = 0.0, 0.0, [], np.zeros(len(Wm))
        for b in range(10):
            m1 = ok & (dbin == b) & (S["tert"] == 0)
            m3 = ok & (dbin == b) & (S["tert"] == 2)
            if m1.sum() < 30 or m3.sum() < 30:
                continue
            w = min(m1.sum(), m3.sum())
            d_ = num[m3].sum() / den[m3].sum() - num[m1].sum() / den[m1].sum()
            bd_ = (Wm[:, m3] @ num[m3]) / (Wm[:, m3] @ den[m3]) - (Wm[:, m1] @ num[m1]) / (Wm[:, m1] @ den[m1])
            est += w * d_
            boot += w * bd_
            bw += w
            bins_used.append(dict(bin=b, dens_range=[float(edges[b]), float(edges[b + 1])], n_T1=int(m1.sum()),
                                  n_T3=int(m3.sum()), T1=float(num[m1].sum() / den[m1].sum()),
                                  T3=float(num[m3].sum() / den[m3].sum())))
        if bw:
            matched[nm] = dict(T3_minus_T1=est / bw, lo95=float(np.percentile(boot / bw, 2.5)),
                               hi95=float(np.percentile(boot / bw, 97.5)), bins=bins_used)
    permic.drop(columns="__exc_light", inplace=True)
    conf["dot_density_matched_T3_minus_T1"] = matched
    conf["ols_zscored"] = ols
    res["confounds"] = conf
    permic.to_csv(os.path.join(WORK, "permic.csv.gz"), index=False)
    with open(os.path.join(WORK, "main.json"), "w") as fh:
        json.dump(res, fh, indent=1)
    print(f"main done ({time.time() - t0:.0f} s)", flush=True)


# ----------------------------------------------------------------------------------------------- stage avg
CONV = None


def avg_worker(args):
    files, st = args
    P = AV["P"]
    box = AV["box"]
    jj = np.arange(box) - box / 2
    YY, XX = np.meshgrid(jj, jj, indexing="ij")
    sums = np.zeros((NG, 200, box, box))
    cnt = np.zeros((NG, 200), int)
    for f in files:
        rows = AV["by_file"][f]
        with mrcfile.mmap(os.path.join(PROJ, f), permissive=True, mode="r") as m:
            stack = np.array(m.data[np.sort(P["blob__idx"][rows])], np.float64)
            order = np.argsort(P["blob__idx"][rows])
            rows = rows[order]
        for img, i in zip(stack, rows):
            gi = AV["gidx"][i]
            if not gi:
                continue
            a = img - img.mean()
            a /= a.std() + 1e-12
            T, _ = g0.ctf_2d(a.shape, P["ctf__df1_A"][i], P["ctf__df2_A"][i], P["ctf__df_angle_rad"][i], apix=A_P)
            pimg = AV["psign"] * np.fft.irfft2(np.fft.rfft2(a) * np.sign(T), s=a.shape)
            px, py = s17.class_to_particle(XX, YY, float(P["alignments2D__pose"][i]),
                                           float(P["alignments2D__shift"][i][0]),
                                           float(P["alignments2D__shift"][i][1]), CONV)
            smp = ndi.map_coordinates(pimg, [box / 2 + py, box / 2 + px], order=1, mode="constant")
            k = int(P["alignments2D__class"][i])
            for g in gi:
                sums[g, k] += smp
                cnt[g, k] += 1
    return sums, cnt


AV = {}


def stage_avg(a):
    global CONV
    t0 = time.time()
    os.makedirs(WORK, exist_ok=True)
    conv = json.load(open(os.path.join(RE_WORK, "conventions.json")))
    CONV = {k: conv["best"][k] for k in ("rot", "shift", "order", "swap")}
    psign = float(conv["particle_contrast_sign_vs_micrograph"])
    S = strata()
    lut = {p: i for i, p in enumerate(S["path"])}
    gm = group_masks(S)
    cfg = STREAMS[a.stream]
    d = np.load(cfg["npz"])
    keys = ["alignments2D__class", "alignments2D__pose", "alignments2D__shift", "blob__path", "blob__idx",
            "location__micrograph_path", "ctf__df1_A", "ctf__df2_A", "ctf__df_angle_rad", "blob__shape"]
    P = {k: d[k] for k in keys}
    mi = np.array([lut[p] for p in P["location__micrograph_path"].astype(str)])
    gidx = [tuple(g for g in range(NG) if gm[g][m]) for m in mi]
    bp = P["blob__path"].astype(str)
    files = np.unique(bp)
    mine = list(files[np.array_split(np.arange(len(files)), a.ntasks)[a.task]])
    order = np.argsort(bp, kind="stable")
    sb = bp[order]
    by_file = {}
    for f in mine:
        lo, hi = np.searchsorted(sb, f, "left"), np.searchsorted(sb, f, "right")
        by_file[f] = order[lo:hi]
    AV.update(P=P, box=int(P["blob__shape"][0][0]), gidx=gidx, by_file=by_file, psign=psign)
    procs = int(os.environ.get("SLURM_CPUS_PER_TASK", 4))
    chunks = [mine[i::procs * 4] for i in range(procs * 4)]
    print(f"{a.stream} task {a.task}/{a.ntasks}: {len(mine)} files, {sum(len(v) for v in by_file.values())} "
          f"particles, {procs} procs", flush=True)
    box = AV["box"]
    sums = np.zeros((NG, 200, box, box))
    cnt = np.zeros((NG, 200), int)
    with Pool(procs) as pool:
        for n, (s_, c_) in enumerate(pool.imap_unordered(avg_worker, [(c, a.stream) for c in chunks if c])):
            sums += s_
            cnt += c_
            print(f"  chunk {n + 1}/{len(chunks)} ({time.time() - t0:.0f} s)", flush=True)
    np.savez(os.path.join(WORK, f"avg_{a.stream}_{a.task:02d}.npz"), sums=sums.astype(np.float32), cnt=cnt,
             groups=np.array(GROUPS))
    print(f"done ({time.time() - t0:.0f} s)", flush=True)


# ----------------------------------------------------------------------------------------------- stage final
def fourier_crop(img, nout):
    n = img.shape[0]
    if n == nout:
        return img
    F = np.fft.fftshift(np.fft.fft2(img))
    c, h = n // 2, nout // 2
    return np.real(np.fft.ifft2(np.fft.ifftshift(F[c - h:c + h, c - h:c + h]))) * (nout / n) ** 2


def window(img, apix, diam_A, edge_A=10.0):
    n = img.shape[0]
    yy, xx = np.mgrid[:n, :n] - (n - 1) / 2
    r = np.hypot(yy, xx) * apix
    r0 = diam_A / 2
    w = np.where(r <= r0, 1.0, np.where(r >= r0 + edge_A, 0.0, 0.5 * (1 + np.cos(np.pi * (r - r0) / edge_A))))
    return (img - np.median(img[r <= r0])) * w


def hole_lenient(img, apix, zt=1.5, dmin=20.0):
    """describe()'s hole logic with a lower mask threshold and a smaller minimum hole (sensitivity variant)."""
    n = img.shape[0]
    yy, xx = np.mgrid[:n, :n] - (n - 1) / 2
    r = np.hypot(yy, xx)
    rings = np.array([img[(r >= a) & (r < a + 1)].std() for a in range(n // 2)])
    inside = np.nonzero(rings > 0.1 * rings[: max(2, n // 8)].max())[0]
    rw = (inside.max() + 1) if len(inside) else n // 2
    win = r < rw
    med = np.median(img[win])
    mad = 1.4826 * np.median(np.abs(img[win] - med)) + 1e-9
    z = np.where(win, (img - med) / mad, 0.0)
    lab, nl = ndi.label(z > zt)
    if nl == 0:
        return 0, 0.0
    best, score = 0, -1
    for i in range(1, nl + 1):
        m = lab == i
        sc = m.sum() / (1 + r[m].min())
        if sc > score:
            best, score = i, sc
    m = ndi.binary_closing(lab == best, iterations=1)
    holes = ndi.binary_fill_holes(m) & ~m
    hl, nh = ndi.label(holes)
    dmax = 0.0
    for i in range(1, nh + 1):
        dmax = max(dmax, 2 * np.sqrt((hl == i).sum() / np.pi) * apix)
    return int(dmax >= dmin), round(float(dmax), 1)


def stage_final(a):
    t0 = time.time()
    res = json.load(open(os.path.join(WORK, "main.json")))
    avgs = {}
    D = {}
    for st, cfg in STREAMS.items():
        fs = sorted(f for f in os.listdir(WORK) if f.startswith(f"avg_{st}_"))
        sums, cnt = None, None
        for f in fs:
            z = np.load(os.path.join(WORK, f))
            sums = z["sums"].astype(np.float64) if sums is None else sums + z["sums"]
            cnt = z["cnt"] if cnt is None else cnt + z["cnt"]
        with mrcfile.open(cfg["ca"], permissive=True) as m:
            ca = np.asarray(m.data, np.float64)
            apix = float(m.voxel_size.x)
        labs = read_labels(cfg["labels"], 200)
        nca = ca.shape[1]
        allsum = sums[:3].sum(0)
        allcnt = cnt[:3].sum(0)
        rows = []
        A = np.zeros((NG + 1, 200, nca, nca))
        yy, xx = np.mgrid[:nca, :nca] - (nca - 1) / 2
        inwin = np.hypot(yy, xx) * apix < cfg["window_A"] / 2
        for k in range(200):
            row = dict(cls=k, label=str(labs[k]), n_all=int(allcnt[k]))
            row["cs"] = s41.describe(ca[k], apix)
            if allcnt[k] > 0:
                im = window(fourier_crop(allsum[k] / allcnt[k], nca), apix, cfg["window_A"])
                A[NG, k] = im
                row["r_vs_cs"] = float(np.corrcoef(im[inwin], ca[k][inwin])[0, 1])
                row["all"] = s41.describe(im, apix)
            for g in range(NG):
                row[f"n_{GROUPS[g]}"] = int(cnt[g, k])
                if cnt[g, k] >= 1:
                    im = window(fourier_crop(sums[g, k] / cnt[g, k], nca), apix, cfg["window_A"])
                    A[g, k] = im
                if cnt[g, k] >= NMIN_AVG:
                    row[GROUPS[g]] = s41.describe(A[g, k], apix)
                    hl, dmax = hole_lenient(A[g, k], apix)
                    row[GROUPS[g]]["hole_lenient"] = hl
                    row[GROUPS[g]]["max_hole_diam_A_lenient"] = dmax
            rows.append(row)
        avgs[st] = dict(A=A, ca=ca, apix=apix, rows=rows)
        tested = [r for r in rows if all(r.get(f"n_{g}", 0) >= NMIN_AVG for g in GROUPS[:3])]
        rv = [r["r_vs_cs"] for r in rows if r["n_all"] >= 1000 and "r_vs_cs" in r]
        summ = dict(n_classes_tested_all3=len(tested),
                    particles_in_tested_classes=int(sum(r["n_all"] for r in tested)),
                    conv_check_median_r_vs_cs=float(np.median(rv)) if rv else None,
                    conv_check_n_classes=len(rv),
                    holes_cs=[r["cls"] for r in rows if r["cs"]["hole"]],
                    holes_all=[r["cls"] for r in rows if r.get("all", {}).get("hole")])
        for g in GROUPS:
            summ[f"holes_{g}"] = [r["cls"] for r in rows if r.get(g, {}).get("hole")]
            summ[f"holes_lenient_{g}"] = [r["cls"] for r in rows if r.get(g, {}).get("hole_lenient")]
            summ[f"n_tested_{g}"] = int(sum(1 for r in rows if g in r))
        th = set(summ["holes_T3_thick"]) | set(summ["holes_top10_thickest"])
        summ["holes_thick_only"] = sorted(th - set(summ["holes_T1_thin"]) - set(summ["holes_T2_mid"]))
        thl = set(summ["holes_lenient_T3_thick"]) | set(summ["holes_lenient_top10_thickest"])
        summ["holes_lenient_thick_only"] = sorted(thl - set(summ["holes_lenient_T1_thin"])
                                                  - set(summ["holes_lenient_T2_mid"]))
        # dot_long_rod mask length by tertile (attached-DNA length proxy from the averages themselves)
        ln = {}
        for lab_ in ("dot_long_rod", "rod_no_dot", "dot_short_rod", "dot_pair_thread"):
            rr = [r for r in tested if r["label"] == lab_]
            if rr:
                w = np.array([r["n_all"] for r in rr], float)
                ln[lab_] = {g: float(np.average([r[g]["length_A"] for r in rr], weights=w)) for g in GROUPS[:3]
                            if all(g in r for r in rr)}
                ln[lab_]["n_classes"] = len(rr)
        summ["mask_length_A_weighted_by_label"] = ln
        D[st] = dict(summary=summ, classes=rows)
        print(st, json.dumps(summ), flush=True)
    res["averages"] = D
    res["averages_method"] = (
        "per-stratum average of each 2D class: own-CTF phase flip, script-17 pose/shift convention (rot -theta, "
        "shift -, particle frame), contrast sign from conventions.json, Fourier-cropped to the CS class-average grid "
        f"and windowed (class2D_window_inner_A, 10 A cosine edge); hole test = scripts/41 describe() on strata with "
        f">= {NMIN_AVG} particles; lenient = z>1.5 mask, hole >= 20 A")
    res["runtime_final_s"] = round(time.time() - t0, 1)
    res = verdict(res)
    with open(OUT_JSON, "w") as fh:
        json.dump(res, fh, indent=1)
    print("wrote", OUT_JSON, flush=True)
    figure(res, avgs)
    print("wrote", FIG, flush=True)


def verdict(res):
    C, E, Av = res["census"], res["attached"]["end_dots"], res["averages"]
    out = {}
    for st in ("rods", "anchors"):
        out[st] = {nm: dict(T1=C[st]["groups"]["T1_thin"]["fraction"][nm]["est"],
                            T3=C[st]["groups"]["T3_thick"]["fraction"][nm]["est"],
                            diff_ci=[C[st]["T3_minus_T1"][nm]["fraction"]["lo95"],
                                     C[st]["T3_minus_T1"][nm]["fraction"]["hi95"]])
                   for nm in C[st]["labels"]}
    res["key_contrasts_T3_vs_T1"] = dict(
        census_fraction=out,
        end_dot_ratio_light=dict(T1=E["dots_off0_light"]["groups"]["T1_thin"]["ratio"],
                                 T3=E["dots_off0_light"]["groups"]["T3_thick"]["ratio"],
                                 top10=E["dots_off0_light"]["groups"]["top10_thickest"]["ratio"],
                                 diff=E["dots_off0_light"]["ratio_T3_minus_T1"]),
        end_dot_excess_per_um2_light=dict(T1=E["dots_off0_light"]["groups"]["T1_thin"]["excess_per_um2"],
                                          T3=E["dots_off0_light"]["groups"]["T3_thick"]["excess_per_um2"],
                                          diff=E["dots_off0_light"]["excess_T3_minus_T1"]),
        end_dot_excess_per_dot_light=dict(T1=E["dots_off0_light"]["groups"]["T1_thin"]["excess_per_unit"],
                                          T3=E["dots_off0_light"]["groups"]["T3_thick"]["excess_per_unit"]),
        dot_density_matched=res["confounds"]["dot_density_matched_T3_minus_T1"],
        holes_thick_only={st: Av[st]["summary"]["holes_thick_only"] for st in Av},
        holes_lenient_thick_only={st: Av[st]["summary"]["holes_lenient_thick_only"] for st in Av})
    return res


def figure(res, avgs):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    permic = pd.read_csv(os.path.join(WORK, "permic.csv.gz"))
    cols = ["#4C72B0", "#55A868", "#C44E52", "#8172B2"]
    fig = plt.figure(figsize=(24, 34), dpi=80)
    sf = fig.subfigures(4, 1, height_ratios=[1.0, 1.05, 1.05, 3.0])
    # row 1: ice distribution + confounds
    ax = sf[0].subplots(1, 4)
    ice = permic["ice"].values
    ax[0].hist(ice, bins=120, range=(0.95, 1.2), color="0.5")
    for k, c in (("t33", "k"), ("t67", "k"), ("p90", "r")):
        ax[0].axvline(res["ice_cuts"][k], color=c, ls="--", lw=1)
    ax[0].set_xlabel("ctf_stats ice_thickness_rel")
    ax[0].set_title(f"{len(ice)} curated micrographs; tertile cuts (k) and top decile (r)", fontsize=10)
    for a_, col, lab in ((ax[1], "ctf_fit_A", "CTF fit resolution (A)"), (ax[2], "defocus_A", "defocus (A)"),
                         (ax[3], "dots_per_um2", "detected dots per um^2 (valid area)")):
        ok = np.isfinite(permic[col].values)
        a_.hexbin(ice[ok], permic[col].values[ok], gridsize=60, extent=(0.95, 1.2, *np.nanpercentile(
            permic[col].values, [0.5, 99.5])), mincnt=1, cmap="viridis", bins="log")
        rho = res["confounds"][col]["spearman_rho"]
        a_.set_title(f"{lab} vs ice: Spearman rho = {rho:+.3f}", fontsize=10)
        a_.set_xlabel("ice_thickness_rel")
    sf[0].suptitle("A7 ice-thickness split: strata and confounds", fontsize=13)

    # row 2: census fractions
    ax = sf[1].subplots(1, 2, gridspec_kw=dict(width_ratios=[1, 1]))
    for a_, st in zip(ax, ("rods", "anchors")):
        C = res["census"][st]
        names = [n for n in C["labels"]]
        w = 0.2
        for gi, g in enumerate(GROUPS):
            e = np.array([C["groups"][g]["fraction"][n]["est"] for n in names])
            lo = np.array([C["groups"][g]["fraction"][n]["lo95"] for n in names])
            hi = np.array([C["groups"][g]["fraction"][n]["hi95"] for n in names])
            x = np.arange(len(names)) + (gi - 1.5) * w
            a_.bar(x, e, w, yerr=[e - lo, hi - e], color=cols[gi], label=f"{g} (n={C['groups'][g]['n_mic']} mic)",
                   capsize=2)
        a_.set_xticks(range(len(names)))
        a_.set_xticklabels(names, rotation=30, ha="right", fontsize=9)
        a_.set_ylabel("fraction of stream particles")
        a_.set_title(f"{st} stream consensus-label census ({C['n_particles']:,} particles, {C['n_mic']} micrographs);"
                     f" 95% foil-hole bootstrap CI", fontsize=10)
        a_.legend(fontsize=8)

    # row 3: dot density, dot_long_rod 2h, end-dot ratio / excess, crowding-matched contrasts
    ax = sf[2].subplots(1, 6, gridspec_kw=dict(width_ratios=[1, 1, 1.3, 1, 1, 1.3]))
    sf[2].subplots_adjust(left=0.04, right=0.99, wspace=0.55)
    D = res["dots"]["groups"]
    A = res["attached"]

    def bars(a_, vals, title, ylab, ref=None, zero=True):
        e = np.array([v["est"] for v in vals])
        lo = np.array([v["lo95"] for v in vals])
        hi = np.array([v["hi95"] for v in vals])
        a_.bar(range(NG), e, yerr=[e - lo, hi - e], color=cols, capsize=3)
        a_.set_xticks(range(NG))
        a_.set_xticklabels([g.replace("_", "\n") for g in GROUPS], fontsize=8)
        if ref is not None:
            a_.axhline(ref, color="k", lw=0.8, ls=":")
        a_.set_title(title, fontsize=9)
        a_.set_ylabel(ylab)
        if not zero:
            rng_ = hi.max() - lo.min()
            a_.set_ylim(lo.min() - 0.6 * rng_, hi.max() + 0.4 * rng_)

    bars(ax[0], [D[g]["dots_per_um2_valid"] for g in GROUPS], "(b) detected dot density", "dots / um^2")
    bars(ax[1], [A["groups"][g]["dot_long_rod_2h_A"]["mean_ci"] for g in GROUPS],
         "(c-i) dot_long_rod 2h (class geometry;\nwindow-limited, 200 A)", "mean 2h (A)", zero=False)
    ED = A["end_dots"]
    w = 0.2
    for si, (sn, lab_) in enumerate((("light", "selection-light ends"), ("nodot", "rod_no_dot ends"),
                                     ("far", "dot-class far end"), ("near", "dot-class near (selected)"))):
        Eg = ED[f"dots_off0_{sn}"]["groups"]
        e = np.array([Eg[g]["ratio"]["est"] for g in GROUPS])
        lo = np.array([Eg[g]["ratio"]["lo95"] for g in GROUPS])
        hi = np.array([Eg[g]["ratio"]["hi95"] for g in GROUPS])
        ax[2].bar(np.arange(NG) + (si - 1.5) * w, e, w, yerr=[e - lo, hi - e], color=cols,
                  alpha=[1, 0.75, 0.5, 0.25][si], edgecolor="k", lw=0.5, capsize=2, label=lab_ if True else None)
    ax[2].axhline(1, color="k", lw=0.8, ls=":")
    ax[2].set_xticks(range(NG))
    ax[2].set_xticklabels([g.replace("_", "\n") for g in GROUPS], fontsize=8)
    ax[2].set_yscale("log")
    ax[2].set_ylabel("observed / shift null")
    ax[2].set_title(f"(c-ii) dots within {END_D_A:.0f} A of a rod end, obs/shift null ({ED['nrep']} reps);\n"
                    "bars per stratum: light | rod_no_dot | far | near(selected)", fontsize=9)
    bars(ax[3], [ED["dots_off0_light"]["groups"][g]["excess_per_unit"] for g in GROUPS],
         "(c-ii) selection-light end-dot excess\nper eligible detected dot", "excess fraction of dots", ref=0.0)
    bars(ax[4], [ED["dots_off0_light"]["groups"][g]["excess_per_um2"] for g in GROUPS],
         "(c-ii) selection-light end-dot excess\nper um^2", "dots / um^2", ref=0.0)
    M = res["confounds"]["dot_density_matched_T3_minus_T1"]
    ks = list(M)
    raw = {"frac_dot_long_rod_rods": res["census"]["rods"]["T3_minus_T1"]["dot_long_rod"]["fraction"],
           "frac_rod_no_dot_rods": res["census"]["rods"]["T3_minus_T1"]["rod_no_dot"]["fraction"],
           "frac_dot_pair_thread_rods": res["census"]["rods"]["T3_minus_T1"]["dot_pair_thread"]["fraction"],
           "frac_dot_short_rod_anchors": res["census"]["anchors"]["T3_minus_T1"]["dot_short_rod"]["fraction"],
           "frac_dot_pair_anchors": res["census"]["anchors"]["T3_minus_T1"]["dot_pair"]["fraction"],
           "end_light_ratio_obs_over_null": ED["dots_off0_light"]["ratio_T3_minus_T1"],
           "end_light_excess_per_dot": None, "end_light_excess_per_um2": ED["dots_off0_light"]["excess_T3_minus_T1"]}
    for i, k in enumerate(ks):
        sc = 1.0 if "ratio" in k or "per_um2" in k else 100.0
        unit = "" if sc == 1.0 else " (x100)"
        if "per_um2" in k:
            sc = 0.01
            unit = " (/100)"
        ax[5].errorbar(M[k]["T3_minus_T1"] * sc, i, xerr=[[(M[k]["T3_minus_T1"] - M[k]["lo95"]) * sc],
                                                          [(M[k]["hi95"] - M[k]["T3_minus_T1"]) * sc]],
                       fmt="o", color="k", capsize=3, label="dot-density matched" if i == 0 else None)
        r_ = raw.get(k)
        if r_ is not None:
            ax[5].errorbar(r_["est"] * sc, i + 0.25, xerr=[[(r_["est"] - r_["lo95"]) * sc], [(r_["hi95"] - r_["est"]) * sc]],
                           fmt="s", color="0.6", capsize=2, label="raw" if i == 0 else None)
        ks[i] = k.replace("frac_", "").replace("_rods", " (rods)").replace("_anchors", " (anch.)") \
            .replace("end_light_", "end ") + unit
    ax[5].axvline(0, color="r", lw=0.8)
    ax[5].set_yticks(range(len(ks)))
    ax[5].set_yticklabels(ks, fontsize=7)
    ax[5].invert_yaxis()
    ax[5].legend(fontsize=8)
    ax[5].set_title("T3 (thick) - T1 (thin), raw vs within dot-density deciles", fontsize=9)

    # row 4: averages montage
    sub = sf[3].subfigures(2, 1)
    for sfi, st in zip(sub, ("rods", "anchors")):
        av = avgs[st]
        rows = av["rows"]
        summ = res["averages"][st]["summary"]
        flag = set(summ["holes_thick_only"]) | set(summ["holes_lenient_thick_only"]) | \
            set(sum([summ[f"holes_{g}"] for g in GROUPS], []))
        keep = [r for r in sorted(rows, key=lambda r: -r["n_all"])
                if r["label"] not in ("ice_contaminant", "noise_junk", "unlabelled")]
        pick, seen = [], {}
        for r in keep:                      # up to 3 classes per label, by size
            if seen.get(r["label"], 0) < 3:
                pick.append(r)
                seen[r["label"]] = seen.get(r["label"], 0) + 1
        pick = (pick + [r for r in rows if r["cls"] in flag and r not in pick])[:14]
        axs = sfi.subplots(NG + 1, len(pick), squeeze=False)
        n = av["ca"].shape[1]
        for j, r in enumerate(pick):
            k = r["cls"]
            ims = [av["ca"][k]] + [av["A"][g, k] for g in range(NG)]
            labs_ = ["CS"] + GROUPS
            for i, (im, lb) in enumerate(zip(ims, labs_)):
                a_ = axs[i, j]
                a_.imshow(im, cmap="gray")
                a_.set_xticks([])
                a_.set_yticks([])
                c = n * 0.18
                a_.set_xlim(c, n - c)
                a_.set_ylim(n - c, c)
                if i == 0:
                    a_.set_title(f"{k} {r['label'].replace('_', ' ')}\nn={r['n_all']}", fontsize=7)
                else:
                    g = GROUPS[i - 1]
                    d = r.get(g)
                    txt = f"n={r[f'n_{g}']}" + ("" if d is None else
                                                f" h{d['hole']}/{d['hole_lenient']} L{d['length_A']:.0f}")
                    a_.text(0.02, 0.02, txt, transform=a_.transAxes, fontsize=6, color="yellow")
                    if d is not None and (d["hole"] or d["hole_lenient"]):
                        for sp in a_.spines.values():
                            sp.set_color("red")
                            sp.set_linewidth(2)
                if j == 0:
                    a_.set_ylabel(labs_[i], fontsize=8)
        sfi.suptitle(f"(d) {st}: CryoSPARC average (top) and per-stratum reference-free averages; "
                     f"h = hole strict/lenient, L = mask length (A); conv. check median r = "
                     f"{summ['conv_check_median_r_vs_cs']:.2f}\nholes: CS {summ['holes_cs']}, all {summ['holes_all']}, "
                     f"T1 {summ['holes_T1_thin']}, T2 {summ['holes_T2_mid']}, T3 {summ['holes_T3_thick']}, top10 "
                     f"{summ['holes_top10_thickest']}; lenient thick-only {summ['holes_lenient_thick_only']}",
                     fontsize=9)
    fig.savefig(FIG)
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("stage", choices=["main", "avg", "final"])
    ap.add_argument("--stream", choices=list(STREAMS), default="rods")
    ap.add_argument("--task", type=int, default=0)
    ap.add_argument("--ntasks", type=int, default=1)
    ap.add_argument("--nrep", type=int, default=NREP)
    a = ap.parse_args()
    {"main": stage_main, "avg": stage_avg, "final": stage_final}[a.stage](a)


if __name__ == "__main__":
    main()
