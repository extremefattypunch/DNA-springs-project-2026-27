#!/usr/bin/env python
"""Plan section 5G (matplotlib part) - cryo-EM plots and the Table S1 comparison.

  figS5_fsc.png           Fig S5C/G analogue: half-map FSC (deposited CryoSPARC curve, RELION 3.1
                          postprocess corrected FSC, Phenix mtriage masked FSC, CryoSPARC re-run) and the
                          map-model FSC (mtriage) per isoform, with the 0.143 / 0.5 thresholds.
  fig3b_com_shift.png     Fig 3B analogue: per-HEAT-repeat CoM distance between Aalpha and Abeta.
  figS5_localres_hist.png local-resolution histograms (RELION relion_locres.mrc, CryoSPARC map) per isoform.
  results/tableS1_comparison.csv  our numbers next to the paper's Table S1 rows.
  results/em_metrics.csv          metrics.csv-schema rows for every number above.
Inputs: $PP2A_DATA/em/{emdb,phenix,relion,cryosparc,geometry}.  Missing inputs are skipped, not fatal.
Run:  $DNASPRING_PY scripts/27_em_figures.py
"""
import csv
import glob
import json
import os
import re
import sys
import xml.etree.ElementTree as ET

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
try:
    import vizstyle  # noqa: E402  (shared palette / chrome from the sfGFP project)
    vizstyle.apply() if hasattr(vizstyle, "apply") else None
except Exception:
    vizstyle = None

REPO = os.environ.get("PP2A_REPO", "/n/home05/ianpoon/DNA-springs-project-2026-27/cryoEM/pp2a_day2025_reproduction")
DATA = os.environ.get("PP2A_DATA", "/n/netscratch/hekstra_lab/Lab/ian_poon/pp2a_day2025")
FIG = os.path.join(REPO, "figures")
RES = os.path.join(REPO, "results")
ISO = {"Aa": {"emd": "48300", "label": "Aα-B56γ-Cα (EMD-48300)", "color": "#c0392b"},
       "Ab": {"emd": "48224", "label": "Aβ-B56γ-Cα (EMD-48224)", "color": "#555555"}}
PAPER_S1 = {  # Table S1 rows we can reproduce
    "Map resolution (FSC 0.143) [A]": {"Aa": 3.4, "Ab": 3.2},
    "Map resolution range [A]": {"Aa": "5.87-3.04", "Ab": "5.22-2.82"},
    "Model resolution (FSC 0.5) [A]": {"Aa": 3.8, "Ab": 3.6},
    "MolProbity score": {"Aa": 1.25, "Ab": 1.15},
    "Clashscore": {"Aa": 4.88, "Ab": 3.59},
    "Poor rotamers [%]": {"Aa": 0, "Ab": 0},
    "Ramachandran favored [%]": {"Aa": 98.08, "Ab": 98.22},
    "Ramachandran allowed [%]": {"Aa": 1.92, "Ab": 1.78},
    "Ramachandran disallowed [%]": {"Aa": 0, "Ab": 0},
    "RMS bond lengths [A]": {"Aa": 0.010, "Ab": 0.009},
    "RMS bond angles [deg]": {"Aa": 1.291, "Ab": 1.241},
}


def read_xml_fsc(path):
    root = ET.parse(path).getroot()
    x = np.array([float(c.find("x").text) for c in root.iter("coordinate")])
    y = np.array([float(c.find("y").text) for c in root.iter("coordinate")])
    return x, y


def read_two_col(path):
    arr = np.loadtxt(path)
    return arr[:, 0], arr[:, 1]


def read_relion_postprocess(path):
    """postprocess.star: data_fsc loop -> resolution (1/A) and corrected FSC; plus FINAL RESOLUTION."""
    x, y, cols, in_loop = [], [], [], False
    with open(path) as f:
        block = None
        for line in f:
            s = line.strip()
            if s.startswith("data_"):
                block = s
                cols, in_loop = [], False
                continue
            if block != "data_fsc":
                continue
            if s.startswith("_rln"):
                cols.append(s.split()[0])
                in_loop = True
                continue
            if in_loop and s and not s.startswith(("#", "loop_")):
                vals = s.split()
                if len(vals) >= len(cols):
                    d = dict(zip(cols, vals))
                    x.append(float(d["_rlnResolution"]))
                    y.append(float(d.get("_rlnFourierShellCorrelationCorrected", d.get("_rlnFourierShellCorrelationMaskedMaps"))))
    return np.array(x), np.array(y)


def res_at(x, y, thr):
    """first crossing below thr (x in 1/A) -> resolution in A"""
    idx = np.where(y < thr)[0]
    if len(idx) == 0 or idx[0] == 0:
        return float("nan")
    i = idx[0]
    x0, x1, y0, y1 = x[i - 1], x[i], y[i - 1], y[i]
    xc = x0 + (thr - y0) * (x1 - x0) / (y1 - y0) if y1 != y0 else x1
    return 1.0 / xc


def parse_mtriage(log):
    out = {}
    s = open(log).read()
    m = re.search(r"FSC\(map,model map\)=0\.5\s*:\s*([\d.]+)\s+([\d.]+)", s)
    if m:
        out["model_fsc05_masked"], out["model_fsc05_unmasked"] = float(m.group(1)), float(m.group(2))
    m = re.search(r"FSC\(half map 1,2\)=0\.143 \(d_fsc\)\s*:\s*([\d.]+)\s+([\d.]+)", s)
    if m:
        out["half_fsc0143_masked"], out["half_fsc0143_unmasked"] = float(m.group(1)), float(m.group(2))
    return out


def parse_molprobity(out_file):
    s = open(out_file).read()
    g = lambda pat: (float(re.search(pat, s).group(1)) if re.search(pat, s) else float("nan"))  # noqa: E731
    return {"rama_outliers": g(r"Ramachandran outliers\s*=\s*([\d.]+)"), "rama_favored": g(r"favored\s*=\s*([\d.]+)"),
            "rotamer_outliers": g(r"Rotamer outliers\s*=\s*([\d.]+)"), "clashscore": g(r"Clashscore\s*=\s*([\d.]+)"),
            "rms_bonds": g(r"RMS\(bonds\)\s*=\s*([\d.]+)"), "rms_angles": g(r"RMS\(angles\)\s*=\s*([\d.]+)"),
            "molprobity": g(r"MolProbity score\s*=\s*([\d.]+)")}


def model_mask(mrc_path, iso, dilate_A=3.0):
    """Voxels within dilate_A of the deposited model's heavy atoms (the paper's range is quoted for the
    molecule, not for the solvent-inclusive RELION mask)."""
    import mdtraj as md
    import mrcfile
    from scipy.ndimage import binary_dilation, generate_binary_structure
    m = mrcfile.open(mrc_path, permissive=True)
    vs = float(m.voxel_size.x)
    origin = np.array([float(m.header.origin.x), float(m.header.origin.y), float(m.header.origin.z)])
    pdb = {"Aa": "9MIP", "Ab": "9MF5"}[iso]
    t = md.load(os.path.join(DATA, "em", "models", pdb + ".pdb"))
    xyz = t.xyz[0, t.topology.select("not element H")] * 10.0
    ijk = np.round((xyz - origin) / vs).astype(int)
    shape = m.data.shape
    ok = np.all((ijk >= 0) & (ijk < np.array(shape)[::-1]), axis=1)
    ijk = ijk[ok]
    mask = np.zeros(shape, dtype=bool)
    mask[ijk[:, 2], ijk[:, 1], ijk[:, 0]] = True
    return binary_dilation(mask, structure=generate_binary_structure(3, 1), iterations=int(round(dilate_A / vs)))


def locres_stats(mrc_path, iso):
    import mrcfile
    m = mrcfile.open(mrc_path, permissive=True).data.astype(np.float32)
    mk = model_mask(mrc_path, iso)
    vals = m[mk & np.isfinite(m) & (m > 0)]
    if len(vals) == 0:
        return None
    return {"min": float(np.percentile(vals, 1)), "max": float(np.percentile(vals, 99)), "median": float(np.median(vals)), "values": vals}


def main():
    os.makedirs(FIG, exist_ok=True)
    os.makedirs(RES, exist_ok=True)
    rows, table = [], {k: {} for k in PAPER_S1}

    def row(system, quantity, method, value, unit, paper="", notes=""):
        rows.append(["em", system, quantity, method, f"{value:.3f}" if isinstance(value, float) else value, unit, paper, notes])

    # ---------------- FSC figure
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.2), sharey=True)
    for ax, (iso, meta) in zip(axes, ISO.items()):
        emd = meta["emd"]
        curves = []
        dep = os.path.join(DATA, "em", "emdb", iso, f"emd_{emd}_fsc.xml")
        if os.path.exists(dep):
            x, y = read_xml_fsc(dep)
            curves.append(("deposited (CryoSPARC GSFSC, corrected)", x, y, "k", "-"))
            r = res_at(x, y, 0.143)
            row(iso, "half-map resolution FSC 0.143", "deposited FSC XML", r, "A", PAPER_S1["Map resolution (FSC 0.143) [A]"][iso], "Table S1")
            table["Map resolution (FSC 0.143) [A]"][iso] = f"{r:.2f} (deposited curve)"
        # RELION: the contour-matched (tight) mask from 22_relion_postprocess.sh, and the loose 15 A-lowpass mask
        # from the sweep (22b) where the resolution has plateaued - the mask dependence is itself a finding
        for pp, label, col, tag in ((os.path.join(DATA, "em", "relion", iso, "mask_sweep", "pp_m5_f0.2_s6.star"),
                                     "RELION 3.1 postprocess, loose mask (15 Å-lowpass, 0.2× thr, corrected)", "#1f77b4", "RELION loose mask"),
                                    (os.path.join(DATA, "em", "relion", iso, "postprocess.star"),
                                     "RELION 3.1 postprocess, tight contour-matched mask (corrected)", "#17becf", "RELION tight mask")):
            if os.path.exists(pp):
                x, y = read_relion_postprocess(pp)
                curves.append((label, x, y, col, "-"))
                r = res_at(x, y, 0.143)
                row(iso, "half-map resolution FSC 0.143", f"RELION 3.1 relion_postprocess ({tag})", r, "A", PAPER_S1["Map resolution (FSC 0.143) [A]"][iso], "Table S1")
                table["Map resolution (FSC 0.143) [A]"][iso] += f"; {r:.2f} ({tag})"
        mt = os.path.join(DATA, "em", "phenix", iso, "fsc_half_maps.masked.mtriage.xml")
        if os.path.exists(mt):
            x, y = read_xml_fsc(mt)
            curves.append(("Phenix mtriage half maps (masked)", x, y, "#2ca02c", "--"))
        mm = os.path.join(DATA, "em", "phenix", iso, "fsc_model.masked.mtriage.log")
        if os.path.exists(mm):
            x, y = read_two_col(mm)
            curves.append(("Phenix mtriage map vs model (masked)", x, y, "#ff7f0e", "-"))
        for cs_txt in sorted(glob.glob(os.path.join(DATA, "em", "cryosparc", f"{iso}_fsc_*fsc*.txt"))):
            try:
                arr = np.loadtxt(cs_txt, skiprows=2)
                # CryoSPARC v5 Validation (FSC) txt: wave_number, frequency (1/A), fsc_nomask, fsc_sphericalmask, fsc_resmask
                x = arr[:, 1]
                y = arr[:, -1]
                curves.append(("CryoSPARC v5 Validation (FSC), auto res-mask", x, y, "#9467bd", ":"))
                r = res_at(x, y, 0.143)
                row(iso, "half-map resolution FSC 0.143", "CryoSPARC v5 Validation (FSC)", r, "A", PAPER_S1["Map resolution (FSC 0.143) [A]"][iso], "Table S1")
                table["Map resolution (FSC 0.143) [A]"][iso] += f"; {r:.2f} (CryoSPARC v5)"
                break
            except Exception as e:  # noqa: BLE001
                print("cryosparc fsc parse failed:", cs_txt, e)
        for label, x, y, col, ls in curves:
            ax.plot(x, y, color=col, ls=ls, lw=1.6, label=label)
        for thr, txt in ((0.5, "0.5"), (0.143, "0.143")):
            ax.axhline(thr, color="grey", lw=0.8, ls=":")
            ax.text(0.005, thr + 0.015, txt, fontsize=8, color="grey")
        ax.set_xlim(0, 0.5)
        ax.set_ylim(-0.05, 1.02)
        ax.set_xlabel("1/resolution (1/Å)")
        ax.set_title(meta["label"], fontsize=10)
        ax.legend(fontsize=7, loc="upper right", frameon=False)
    axes[0].set_ylabel("Fourier shell correlation")
    fig.suptitle("Fig S5C/G analogue: FSC of the deposited half maps and of map vs deposited model", fontsize=11)
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "figS5_fsc.png"), dpi=200)
    plt.close(fig)

    # ---------------- Table S1 numbers from Phenix
    for iso in ISO:
        mt_log = os.path.join(DATA, "em", "phenix", iso, f"mtriage_{iso}.log")
        mp_out = os.path.join(DATA, "em", "phenix", iso, "molprobity.out")
        if os.path.exists(mt_log):
            m = parse_mtriage(mt_log)
            if "model_fsc05_masked" in m:
                table["Model resolution (FSC 0.5) [A]"][iso] = f"{m['model_fsc05_masked']:.2f} (mtriage masked; unmasked {m['model_fsc05_unmasked']:.2f})"
                row(iso, "model resolution FSC 0.5", "phenix.mtriage masked", m["model_fsc05_masked"], "A", PAPER_S1["Model resolution (FSC 0.5) [A]"][iso], "Table S1")
            if "half_fsc0143_masked" in m:
                row(iso, "half-map resolution FSC 0.143", "phenix.mtriage masked", m["half_fsc0143_masked"], "A", PAPER_S1["Map resolution (FSC 0.143) [A]"][iso], "Table S1")
        if os.path.exists(mp_out):
            p = parse_molprobity(mp_out)
            for key, name, unit in (("molprobity", "MolProbity score", ""), ("clashscore", "Clashscore", ""),
                                    ("rotamer_outliers", "Poor rotamers [%]", "%"), ("rama_favored", "Ramachandran favored [%]", "%"),
                                    ("rama_outliers", "Ramachandran disallowed [%]", "%"), ("rms_bonds", "RMS bond lengths [A]", "A"),
                                    ("rms_angles", "RMS bond angles [deg]", "deg")):
                table[name][iso] = f"{p[key]:.4g}"
                row(iso, name, "phenix.molprobity on deposited model", p[key], unit, PAPER_S1[name][iso], "Table S1")
            table["Ramachandran allowed [%]"][iso] = f"{100 - p['rama_favored'] - p['rama_outliers']:.2f}"

    # ---------------- local resolution histograms
    fig, axes = plt.subplots(1, 2, figsize=(10, 3.8))
    any_lr = False
    for ax, (iso, meta) in zip(axes, ISO.items()):
        rel = os.path.join(DATA, "em", "relion", iso, "relion_locres.mrc")
        mask = os.path.join(DATA, "em", "relion", iso, "mask.mrc")
        parts = []
        if os.path.exists(rel):
            st = locres_stats(rel, iso)
            if st:
                any_lr = True
                ax.hist(st["values"], bins=60, range=(2, 8), color=meta["color"], alpha=0.7, label=f"RELION 3.1 locres ({st['min']:.2f}–{st['max']:.2f} Å)")
                parts.append(f"{st['min']:.2f}-{st['max']:.2f} (RELION)")
                row(iso, "local resolution range (1-99 pct within 3 A of model)", "RELION 3.1 --locres", f"{st['min']:.2f}-{st['max']:.2f}", "A", PAPER_S1["Map resolution range [A]"][iso], "Table S1 / Fig S5A,E")
        for cs_lr in sorted(glob.glob(os.path.join(DATA, "em", "cryosparc", f"{iso}_locres_*map_locres*.mrc"))):
            st = locres_stats(cs_lr, iso)
            if st:
                any_lr = True
                ax.hist(st["values"], bins=60, range=(2, 8), histtype="step", color="#9467bd", lw=1.5, label=f"CryoSPARC v5 ({st['min']:.2f}–{st['max']:.2f} Å)")
                parts.append(f"{st['min']:.2f}-{st['max']:.2f} (CryoSPARC)")
                row(iso, "local resolution range (1-99 pct within 3 A of model)", "CryoSPARC v5 Local Resolution", f"{st['min']:.2f}-{st['max']:.2f}", "A", PAPER_S1["Map resolution range [A]"][iso], "Fig S5A,E")
                break
        table["Map resolution range [A]"][iso] = "; ".join(parts)
        ax.set_title(meta["label"], fontsize=10)
        ax.set_xlabel("local resolution (Å)")
        ax.legend(fontsize=7, frameon=False)
    axes[0].set_ylabel("voxels in mask")
    fig.suptitle("Fig S5A/E analogue: local resolution of the deposited half maps", fontsize=11)
    fig.tight_layout()
    if any_lr:
        fig.savefig(os.path.join(FIG, "figS5_localres_hist.png"), dpi=200)
    plt.close(fig)

    # ---------------- Fig 3B: per-HEAT CoM shift
    geo = os.path.join(DATA, "em", "geometry", "geometry.json")
    if os.path.exists(geo):
        g = json.load(open(geo))
        sh = g["com_shift_nominal"]["shift"]
        ks = sorted(int(k) for k in sh)
        fig, ax = plt.subplots(figsize=(6.5, 3.4))
        ax.bar([f"HEAT {k}" for k in ks], [sh[str(k)] for k in ks], color="#333333")
        ax.set_ylabel("Δ centre of mass, Aα vs Aβ (Å)")
        ax.set_title("Fig 3B analogue: HEAT-repeat CoM distance between Aα (9MIP) and Aβ (9MF5)\nafter superposition on HEAT 12–15 Cα", fontsize=10)
        ax.set_ylim(0, 4)
        plt.setp(ax.get_xticklabels(), rotation=60, ha="right", fontsize=8)
        fig.tight_layout()
        fig.savefig(os.path.join(FIG, "fig3b_com_shift.png"), dpi=200)
        plt.close(fig)

    # ---------------- EMReady re-sharpening vs deposited map (report question G)
    for iso in ISO:
        cj = os.path.join(DATA, "em", "emready", iso, "compare.json")
        if os.path.exists(cj):
            c = json.load(open(cj))
            row(iso, "EMReady(raw) vs deposited map real-space CC (all voxels)", "EMReady v1.1 on half-map average", c["cc_all"], "", "", "pre-registered > 0.95")
            row(iso, "EMReady(raw) vs deposited map real-space CC (top 1% of deposited density)", "EMReady v1.1 on half-map average", c["cc_top1pct"], "", "", "question G")
            row(iso, "FSC 0.5 between EMReady(raw) and deposited map", "numpy shells, unmasked", c["fsc05_between_maps_A"], "A", "", "question G")

    # ---------------- write tables
    with open(os.path.join(RES, "tableS1_comparison.csv"), "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["row", "paper Aa", "ours Aa", "paper Ab", "ours Ab"])
        for k, v in PAPER_S1.items():
            w.writerow([k, v["Aa"], table[k].get("Aa", ""), v["Ab"], table[k].get("Ab", "")])
    with open(os.path.join(RES, "em_metrics.csv"), "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["track", "system", "quantity", "selection_or_method", "value", "unit", "paper_value", "notes"])
        w.writerows(rows)
    for k, v in PAPER_S1.items():
        print(f"{k:36s} Aa paper {v['Aa']!s:10s} ours {table[k].get('Aa', '-'):45s} | Ab paper {v['Ab']!s:10s} ours {table[k].get('Ab', '-')}")


if __name__ == "__main__":
    sys.exit(main())
