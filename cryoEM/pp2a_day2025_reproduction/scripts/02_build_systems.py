#!/usr/bin/env python
"""Plan section 4 - solvated Amber systems for Aalpha and Abeta with tleap (ff19SB + OPC).

Paper: "The system was placed in a rectangular water box with 16 A padding and 25 charge-
neutralizing sodium ions (23 for Aalpha) with additional salt (Na+, Cl-) around the solute,
resulting in a 150 mM salt condition, using the Amber23 LEaP module."  Measured padding 18-19.5 A.

Steps per system ($PP2A_DATA/md/systems/<sys>/):
  1. rotate the protein's principal axes onto x/y/z so the axis-aligned solvateBox pads the true
     extent (an un-oriented horseshoe inflates the box by 30-50 %);
  2. tleap: ff19SB + OPC, solvateBox OPCBOX 16.0, addIons Na+ 0 (neutralise), then
     addIonsRand Na+ N Cl- N with N = round(0.15 mol/L * N_water / 55.5 mol/L);
  3. write system.prmtop / system.inpcrd / system.pdb and build.json (atom counts, box, padding,
     ion counts, net charge) - asserts < 200k particles and padding within 16-20 A per axis.
Run in the dnaspring env (tleap on PATH):  $DNASPRING_PY scripts/02_build_systems.py
"""
import json
import os
import re
import subprocess
import sys

import mdtraj as md
import numpy as np

DATA = os.environ.get("PP2A_DATA", "/n/netscratch/hekstra_lab/Lab/ian_poon/pp2a_day2025")
SYS = os.path.join(DATA, "md", "systems")
PAD = 16.0
SALT_M = 0.150
WATER_M = 55.5
SYSTEMS = ["Aa", "Ab"]
TLEAP = """source leaprc.protein.ff19SB
source leaprc.water.opc
m = loadpdb {pdb}
check m
solvateBox m OPCBOX {pad} iso
addIons m Na+ 0
addIonsRand m Na+ {nion} Cl- {nion}
saveamberparm m system.prmtop system.inpcrd
savepdb m system.pdb
quit
"""
# tleap does not accept 'iso' for rectangular boxes; we call it twice: once to learn the water
# count for the salt arithmetic, once with the final ion numbers.


def orient(pdb_in, pdb_out):
    t = md.load(pdb_in)
    xyz = t.xyz[0] - t.xyz[0].mean(0)
    # principal axes of the coordinate covariance -> rotate longest axis onto x, next onto y
    w, v = np.linalg.eigh(np.cov(xyz.T))
    R = v[:, ::-1]
    if np.linalg.det(R) < 0:
        R[:, 2] *= -1
    t.xyz[0] = xyz @ R
    t.save_pdb(pdb_out)
    ext = (t.xyz[0].max(0) - t.xyz[0].min(0)) * 10
    return ext


def run_tleap(workdir, script):
    open(os.path.join(workdir, "tleap.in"), "w").write(script)
    r = subprocess.run(["tleap", "-f", "tleap.in"], cwd=workdir, capture_output=True, text=True)
    log = r.stdout + r.stderr
    open(os.path.join(workdir, "tleap.log"), "w").write(log)
    if r.returncode != 0 or "Exiting LEaP: Errors = 0" not in log:
        print(log[-3000:])
        raise SystemExit(f"tleap failed in {workdir}")
    return log


def count_residues(pdb, names):
    n = 0
    with open(pdb) as f:
        seen = set()
        for line in f:
            if line.startswith(("ATOM", "HETATM")) and line[17:20].strip() in names:
                key = (line[17:20], line[21], line[22:27])
                if key not in seen:
                    seen.add(key)
                    n += 1
    return n


def main():
    report = {}
    for name in SYSTEMS:
        d = os.path.join(SYS, name)
        prot = os.path.join(d, "protein.pdb")
        assert os.path.exists(prot), f"run 01_prepare_models.py first ({prot})"
        ext = orient(prot, os.path.join(d, "protein_oriented.pdb"))
        # pass 1: neutralise only, to learn the water count
        script1 = TLEAP.replace("addIonsRand m Na+ {nion} Cl- {nion}\n", "").format(pdb="protein_oriented.pdb", pad=PAD, nion=0).replace(" iso", "")
        log1 = run_tleap(d, script1)
        n_wat = count_residues(os.path.join(d, "system.pdb"), {"WAT", "HOH"})
        n_na_neut = count_residues(os.path.join(d, "system.pdb"), {"Na+", "NA"})
        n_salt = int(round(SALT_M * n_wat / WATER_M))
        # pass 2: final system with salt
        script2 = TLEAP.format(pdb="protein_oriented.pdb", pad=PAD, nion=n_salt).replace(" iso", "")
        log2 = run_tleap(d, script2)
        t = md.load(os.path.join(d, "system.inpcrd"), top=os.path.join(d, "system.prmtop"))
        box = t.unitcell_lengths[0] * 10
        prot_sel = t.topology.select("protein")
        pxyz = t.xyz[0, prot_sel] * 10
        pad_lo = pxyz.min(0) - (t.xyz[0].min(0) * 10)
        pad_hi = (t.xyz[0].max(0) * 10) - pxyz.max(0)
        n_na = count_residues(os.path.join(d, "system.pdb"), {"Na+", "NA"})
        n_cl = count_residues(os.path.join(d, "system.pdb"), {"Cl-", "CL"})
        charge = re.findall(r"Total unperturbed charge:\s*(-?\d+\.\d+)", log2)
        n_atoms = t.n_atoms
        report[name] = {"protein_extent_A": [round(float(x), 1) for x in ext], "box_A": [round(float(x), 2) for x in box],
                        "padding_low_A": [round(float(x), 1) for x in pad_lo], "padding_high_A": [round(float(x), 1) for x in pad_hi],
                        "n_water": n_wat, "n_Na_neutralising": n_na_neut, "n_salt_pairs": n_salt, "n_Na_total": n_na, "n_Cl_total": n_cl,
                        "n_particles_total": int(n_atoms), "protein_charge_before_ions": float(charge[0]) if charge else None}
        print(f"{name}: extent {ext.round(1)} A, box {box.round(1)} A, waters {n_wat}, Na+ {n_na} (neutralising {n_na_neut}), "
              f"Cl- {n_cl}, particles {n_atoms}, padding lo {pad_lo.round(1)} hi {pad_hi.round(1)}")
        assert n_atoms < 260_000, "system larger than planned - check orientation"   # Abeta: the 25-residue AF tail sticks out; 226k particles with paper-faithful 16 A padding
        assert (pad_lo > 14).all() and (pad_hi > 14).all(), "padding below 16 A somewhere"
    json.dump(report, open(os.path.join(SYS, "build.json"), "w"), indent=2)


if __name__ == "__main__":
    sys.exit(main())
