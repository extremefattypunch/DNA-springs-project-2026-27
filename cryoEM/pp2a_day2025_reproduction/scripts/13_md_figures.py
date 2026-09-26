#!/usr/bin/env python
"""Fig 1A (RMSF per residue) and Fig 1B (RMSD vs time, N-terminus and excluding N-terminus) analogues.

Reads results/md_rmsd_<sys>.csv, results/md_rmsf_gmx_<sys>.csv (GROMACS, preferred) or md_rmsf_<sys>.csv (MDAnalysis).
Paper axis ranges: RMSF 0-25 A over residue 0-600; RMSD 0-30 A over 0-100 ns.  Aalpha red, Abeta grey/black as in the paper.
Also appends summary rows to results/md_metrics.csv.  Usage: $DNASPRING_PY scripts/13_md_figures.py
"""
import csv
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

REPO = os.environ.get("PP2A_REPO", "/n/home05/ianpoon/DNA-springs-project-2026-27/cryoEM/pp2a_day2025_reproduction")
FIG, RES = os.path.join(REPO, "figures"), os.path.join(REPO, "results")
COL = {"Aa": "#d62728", "Ab": "#333333"}
LAB = {"Aa": "Aα", "Ab": "Aβ"}
NT = {"Aa": 13, "Ab": 25}


def read_csv(path):
    with open(path) as f:
        r = csv.DictReader(f)
        rows = list(r)
    return {k: np.array([float(x[k]) if x[k] not in ("", None) else np.nan for x in rows]) for k in rows[0]} if rows else {}


def main():
    rmsd = {s: read_csv(os.path.join(RES, f"md_rmsd_{s}.csv")) for s in COL if os.path.exists(os.path.join(RES, f"md_rmsd_{s}.csv"))}
    rmsf = {}
    for s in COL:
        for cand, tool in ((f"md_rmsf_gmx_{s}.csv", "GROMACS 2023.3 gmx rmsf -fit"), (f"md_rmsf_{s}.csv", "MDAnalysis")):
            p = os.path.join(RES, cand)
            if os.path.exists(p):
                rmsf[s] = (read_csv(p), tool)
                break
    if not rmsd and not rmsf:
        print("no MD results yet")
        return
    rows = []
    # ---- Fig 1A
    if rmsf:
        fig, ax = plt.subplots(figsize=(5.2, 3.4))
        for s, (d, tool) in rmsf.items():
            ax.plot(d["resid"], d["rmsf_A"], color=COL[s], lw=1.2, label=LAB[s])
            nt = d["rmsf_A"][d["resid"] <= NT[s]]
            core = d["rmsf_A"][d["resid"] > NT[s]]
            rows.append(["md", s, "RMSF N-terminus mean", f"Calpha, resid 1-{NT[s]}, {tool}", f"{nt.mean():.2f}", "A", "", "Fig 1A"])
            rows.append(["md", s, "RMSF excluding N-terminus mean", f"Calpha, resid > {NT[s]}, {tool}", f"{core.mean():.2f}", "A", "", "Fig 1A"])
            rows.append(["md", s, "RMSF max (residue)", tool, f"{d['rmsf_A'].max():.2f} (res {int(d['resid'][d['rmsf_A'].argmax()])})", "A", "", "Fig 1A"])
        ax.set_xlim(0, 620)
        ax.set_ylim(0, 25)
        ax.set_xlabel("Residue")
        ax.set_ylabel("RMSF (Å)")
        ax.set_title("Fig 1A analogue: RMSF, Aα vs Aβ (100 ns, one replicate each)", fontsize=10)
        ax.legend(frameon=False)
        fig.tight_layout()
        fig.savefig(os.path.join(FIG, "fig1a_rmsf.png"), dpi=200)
        plt.close(fig)
    # ---- Fig 1B (two panels: N-terminus, excluding N-terminus) + whole protein
    if rmsd:
        fig, axes = plt.subplots(1, 3, figsize=(12, 3.4), sharey=True)
        for ax, key, title in zip(axes, ("rmsd_nterm_A", "rmsd_core_A", "rmsd_all_A"),
                                  ("RMSD – N-terminus (1–13 / 1–25)", "RMSD – excluding N-terminus", "RMSD – whole protein (paper text: 10–15 vs 20–25 Å)")):
            for s, d in rmsd.items():
                ax.plot(d["time_ns"], d[key], color=COL[s], lw=0.9, label=LAB[s])
                half = d[key][len(d[key]) // 2:]
                rows.append(["md", s, f"{title} mean over second half", "backbone, self-fit vs frame 0", f"{np.nanmean(half):.2f}", "A",
                             ("10-15 (Aa) / 20-25 (Ab)" if key == "rmsd_all_A" else ""), "Fig 1B"])
            ax.set_xlim(0, 100)
            ax.set_ylim(0, 30)
            ax.set_xlabel("Time (ns)")
            ax.set_title(title, fontsize=9)
        axes[0].set_ylabel("RMSD (Å)")
        axes[0].legend(frameon=False)
        fig.suptitle("Fig 1B analogue: backbone RMSD vs frame 0 (VMD RMSD Trajectory Tool convention)", fontsize=10)
        fig.tight_layout()
        fig.savefig(os.path.join(FIG, "fig1b_rmsd.png"), dpi=200)
        plt.close(fig)
        # supplementary: N-terminus displacement after fitting on the core
        fig, ax = plt.subplots(figsize=(5.2, 3.4))
        for s, d in rmsd.items():
            ax.plot(d["time_ns"], d["rmsd_nterm_corefit_A"], color=COL[s], lw=0.9, label=LAB[s])
        ax.set_xlabel("Time (ns)")
        ax.set_ylabel("N-terminus RMSD after core fit (Å)")
        ax.set_title("N-terminal tail displacement relative to the solenoid", fontsize=10)
        ax.legend(frameon=False)
        fig.tight_layout()
        fig.savefig(os.path.join(FIG, "fig1b_rmsd_nterm_corefit.png"), dpi=200)
        plt.close(fig)
    with open(os.path.join(RES, "md_metrics.csv"), "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["track", "system", "quantity", "selection_or_method", "value", "unit", "paper_value", "notes"])
        w.writerows(rows)
    for r in rows:
        print(r[1], r[2], r[4], r[5], r[6])


if __name__ == "__main__":
    sys.exit(main())
