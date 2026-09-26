#!/usr/bin/env python
"""Plan section 4 analysis - Fig 1A/1B numbers from the production trajectories.

Per system (Aa: N-terminus 1-13, core 14-589; Ab: N-terminus 1-25, core 26-601):
  * writes a protein-only trajectory  $RUN/protein.xtc  and reference  $RUN/protein_ref.pdb  (frame 0)
    - inputs for GROMACS `gmx rmsf` (11_md_rmsf_gromacs.sh) and for ChimeraX (14, 15);
  * RMSD vs frame 0 (the paper's VMD RMSD Trajectory Tool convention: superpose ON the selection, then RMSD OF the
    same selection, backbone atoms) for N-terminus, core ("excluding N-terminus") and whole protein;
    plus the N-terminus RMSD after fitting on the core (how far the tail moves relative to the body);
  * per-residue Calpha RMSF after fitting the core Calpha to the average structure (MDAnalysis; GROMACS in 11).
Outputs results/md_rmsd_<sys>.csv (time_ns, rmsd_nterm, rmsd_core, rmsd_all, rmsd_nterm_corefit) and
results/md_rmsf_<sys>.csv (resid, rmsf_A).  Usage: $DNASPRING_PY scripts/10_md_rmsd_rmsf.py [Aa Ab]
"""
import csv
import os
import sys

import numpy as np
import MDAnalysis as mda
from MDAnalysis.analysis import align, rms

REPO = os.environ.get("PP2A_REPO", "/n/home05/ianpoon/DNA-springs-project-2026-27/cryoEM/pp2a_day2025_reproduction")
DATA = os.environ.get("PP2A_DATA", "/n/netscratch/hekstra_lab/Lab/ian_poon/pp2a_day2025")
SYSTEMS = {"Aa": {"nterm": (1, 13), "core": (14, 589)}, "Ab": {"nterm": (1, 25), "core": (26, 601)}}
FRAME_PS = 100.0


def main(systems):
    os.makedirs(os.path.join(REPO, "results"), exist_ok=True)
    for name in systems:
        run = os.path.join(DATA, "md", "runs", f"{name}_rep1")
        top, dcd = os.path.join(DATA, "md", "systems", name, "system.prmtop"), os.path.join(run, "traj.dcd")
        if not os.path.exists(dcd):
            print(f"{name}: no trajectory yet at {dcd}")
            continue
        u = mda.Universe(top, dcd, format="DCD")
        prot = u.select_atoms("protein")
        print(f"{name}: {len(u.trajectory)} frames, {prot.n_atoms} protein atoms")
        # protein-only trajectory + reference (make the protein whole and centred first)
        xtc, ref = os.path.join(run, "protein.xtc"), os.path.join(run, "protein_ref.pdb")
        if not os.path.exists(xtc):
            with mda.Writer(xtc, prot.n_atoms) as w:
                for ts in u.trajectory:
                    prot.unwrap(compound="fragments")
                    prot.translate(-prot.center_of_mass())
                    w.write(prot)
            u.trajectory[0]
            prot.unwrap(compound="fragments")
            prot.translate(-prot.center_of_mass())
            prot.write(ref)
        up = mda.Universe(ref, xtc)
        n0, n1 = SYSTEMS[name]["nterm"]
        c0, c1 = SYSTEMS[name]["core"]
        sel_n = f"backbone and resid {n0}:{n1}"
        sel_c = f"backbone and resid {c0}:{c1}"
        sel_all = "backbone"
        # RMSD with self-fit for each selection (VMD RMSD Trajectory Tool default), ref = frame 0
        out = {}
        for key, sel in (("nterm", sel_n), ("core", sel_c), ("all", sel_all)):
            R = rms.RMSD(up, up, select=sel, ref_frame=0).run()
            out[key] = R.results.rmsd[:, 2]
        R = rms.RMSD(up, up, select=sel_c, groupselections=[sel_n], ref_frame=0).run()
        out["nterm_corefit"] = R.results.rmsd[:, 3]
        t = np.arange(len(out["all"])) * FRAME_PS / 1000.0
        with open(os.path.join(REPO, "results", f"md_rmsd_{name}.csv"), "w", newline="") as f:
            w = csv.writer(f)
            w.writerow(["time_ns", "rmsd_nterm_A", "rmsd_core_A", "rmsd_all_A", "rmsd_nterm_corefit_A"])
            for i in range(len(t)):
                w.writerow([f"{t[i]:.2f}", f"{out['nterm'][i]:.3f}", f"{out['core'][i]:.3f}", f"{out['all'][i]:.3f}", f"{out['nterm_corefit'][i]:.3f}"])
        # RMSF: fit core Calpha to the average structure, then per-residue Calpha RMSF over the full run
        avg = align.AverageStructure(up, up, select=f"name CA and resid {c0}:{c1}", ref_frame=0).run()
        ref_u = avg.results.universe
        align.AlignTraj(up, ref_u, select=f"name CA and resid {c0}:{c1}", in_memory=True).run()
        ca = up.select_atoms("name CA")
        F = rms.RMSF(ca).run()
        with open(os.path.join(REPO, "results", f"md_rmsf_{name}.csv"), "w", newline="") as f:
            w = csv.writer(f)
            w.writerow(["resid", "resname", "rmsf_A"])
            for a, v in zip(ca, F.results.rmsf):
                w.writerow([a.resid, a.resname, f"{v:.3f}"])
        nt = F.results.rmsf[(ca.resids >= n0) & (ca.resids <= n1)]
        co = F.results.rmsf[(ca.resids >= c0) & (ca.resids <= c1)]
        last = slice(len(t) // 2, None)
        print(f"{name}: RMSD(all) plateau {np.mean(out['all'][last]):.1f} A, RMSD(N-term, self-fit) {np.mean(out['nterm'][last]):.1f} A, "
              f"RMSD(core) {np.mean(out['core'][last]):.1f} A, RMSD(N-term|core fit) {np.mean(out['nterm_corefit'][last]):.1f} A; "
              f"RMSF N-term mean {nt.mean():.2f} A vs core {co.mean():.2f} A")


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:] or list(SYSTEMS)))
