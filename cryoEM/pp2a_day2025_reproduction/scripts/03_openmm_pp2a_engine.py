#!/usr/bin/env python
"""Plan section 4, Engine B - OpenMM 8.6 port of the paper's NAMD protocol (one full-length A subunit).

Fork of "MD simulations/sfGFP/runs/openmm_engine.py" (same restraint / barostat / checkpoint machinery),
re-parameterised to Day et al. 2025 and collapsed into one resumable command:

    python 03_openmm_pp2a_engine.py --sysdir <dir with system.prmtop/inpcrd> --outdir <run dir> --ns 100

Stages (each writes a marker file, so a re-run after a 12 h Slurm limit continues where it stopped):
  1 minimise   10,000 steps, protein heavy atoms restrained (k = 10 kcal/mol/A^2) to the input
               coordinates, then the restraint reference is moved to the minimised coordinates.
  2 heat       NVT, 30 x 10,000 steps, temperature raised 10 K per stage 0 -> 300 K (NAMD reassignFreq/
               reassignIncr analogue), restraints on.
  3 release    NPT, MonteCarloBarostat(1 bar, 300 K, every 25 steps) as the Langevin-piston substitute;
               k_res x 0.8, 0.6, 0.4, 0.2, 0.0 for 50,000 steps each (paper: "0.2 decrement every 50,000 steps").
  4 produce    NPT, LangevinMiddleIntegrator 300 K, friction 5 ps^-1 (NAMD langevinDamping 5), PME 10 A,
               2 fs with rigid bonds (--hmr switches to 4 fs with solute hydrogen-mass repartitioning),
               DCD + state every 100 ps (1,000 frames / 100 ns), checkpoint every 1 ns, restart-safe.
Deviations from NAMD that the README must list: barostat type (MC vs Langevin piston), cutoff/switching
(Amber 10 A hard cutoff + dispersion correction vs NAMD 12/10 switching), timestep unstated in the paper.
"""
import argparse
import json
import os
import sys
import time
from pathlib import Path

import numpy as np
import openmm as mm
import openmm.app as app
import openmm.unit as unit
import parmed

TEMPERATURE = 300.0 * unit.kelvin
PRESSURE = 1.01325 * unit.bar
FRICTION = 5.0 / unit.picosecond            # paper: langevinDamping 5
CUTOFF = 1.0 * unit.nanometer               # 10 A
RESTRAINT_K = 10.0 * unit.kilocalorie_per_mole / unit.angstrom**2
HEAT_STAGES, HEAT_STEPS, HEAT_INCR = 30, 10_000, 10.0      # 0 -> 300 K
RELEASE_SCALES, RELEASE_STEPS = (0.8, 0.6, 0.4, 0.2, 0.0), 50_000
MIN_STEPS = 10_000


def log(msg):
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def load_structure(sysdir, hmr):
    st = parmed.load_file(str(sysdir / "system.prmtop"), xyz=str(sysdir / "system.inpcrd"))
    if hmr:
        parmed.tools.HMassRepartition(st, "3.024").execute()
        log("HMR applied to solute hydrogens (water untouched)")
    return st


def build_system(st):
    system = st.createSystem(nonbondedMethod=app.PME, nonbondedCutoff=CUTOFF, constraints=app.HBonds,
                             rigidWater=True, removeCMMotion=True)
    for f in system.getForces():
        if isinstance(f, mm.NonbondedForce):
            f.setUseDispersionCorrection(True)
            f.setEwaldErrorTolerance(5e-4)
    return system


def protein_heavy_mask(st):
    solvent = {"WAT", "HOH", "Na+", "Cl-", "NA", "CL", "K+"}
    return [(a.residue.name not in solvent and a.atomic_number > 1) for a in st.atoms]


def add_restraints(system, positions_nm, mask, k=RESTRAINT_K):
    force = mm.CustomExternalForce("k_res*periodicdistance(x, y, z, x0, y0, z0)^2")
    force.addGlobalParameter("k_res", k)
    for p in ("x0", "y0", "z0"):
        force.addPerParticleParameter(p)
    idx = [i for i, keep in enumerate(mask) if keep]
    for i in idx:
        force.addParticle(i, [float(v) for v in positions_nm[i]])
    system.addForce(force)
    log(f"restrained {len(idx)} protein heavy atoms at k = {k}")
    return force, idx


def make_sim(st, system, dt, platform_name):
    integrator = mm.LangevinMiddleIntegrator(TEMPERATURE, FRICTION, dt)
    if platform_name:
        platform = mm.Platform.getPlatformByName(platform_name)
        props = {"Precision": "mixed"} if platform_name == "CUDA" else {}
        sim = app.Simulation(st.topology, system, integrator, platform, props)
    else:
        sim = app.Simulation(st.topology, system, integrator)
    log(f"platform {sim.context.getPlatform().getName()}, dt {dt}")
    return sim, integrator


def save_state(sim, path):
    state = sim.context.getState(getPositions=True, getVelocities=True, getEnergy=True, enforcePeriodicBox=True)
    Path(path).write_text(mm.XmlSerializer.serialize(state))
    return state


def load_state(sim, path):
    sim.context.setState(mm.XmlSerializer.deserialize(Path(path).read_text()))


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--sysdir", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--ns", type=float, default=100.0)
    ap.add_argument("--hmr", action="store_true", help="4 fs with hydrogen mass repartitioning (fallback only)")
    ap.add_argument("--platform", default="CUDA")
    ap.add_argument("--seed", type=int, default=2001)
    ap.add_argument("--report-ps", type=float, default=100.0)
    args = ap.parse_args()
    sysdir, out = Path(args.sysdir), Path(args.outdir)
    out.mkdir(parents=True, exist_ok=True)
    dt_prod = (4.0 if args.hmr else 2.0) * unit.femtosecond
    dt_equil = 2.0 * unit.femtosecond
    log(f"host {os.uname().nodename}  sys {sysdir}  out {out}  ns {args.ns}  dt_prod {dt_prod}")
    st = load_structure(sysdir, args.hmr)
    log(f"{len(st.atoms)} particles, {len(st.residues)} residues, box {st.box}")
    mask = protein_heavy_mask(st)
    pos0 = np.array(st.positions.value_in_unit(unit.nanometer))
    system = build_system(st)
    force, ridx = add_restraints(system, pos0, mask)
    meta_path = out / "progress.json"
    meta = json.loads(meta_path.read_text()) if meta_path.exists() else {}

    # ---------- 1 minimise (restrained) ----------
    if not (out / "01_minimised.xml").exists():
        sim, integ = make_sim(st, system, dt_equil, args.platform)
        sim.context.setPositions(st.positions)
        sim.context.setPeriodicBoxVectors(*st.box_vectors)
        e0 = sim.context.getState(getEnergy=True).getPotentialEnergy()
        sim.minimizeEnergy(maxIterations=MIN_STEPS, tolerance=5.0 * unit.kilojoule_per_mole / unit.nanometer)
        state = save_state(sim, out / "01_minimised.xml")
        log(f"minimised: {e0} -> {state.getPotentialEnergy()}")
        app.PDBFile.writeFile(sim.topology, state.getPositions(), open(out / "01_minimised.pdb", "w"), keepIds=True)
        del sim, integ
    # restraint reference = minimised coordinates (paper restrains heavy atoms during heating)
    min_state = mm.XmlSerializer.deserialize((out / "01_minimised.xml").read_text())
    pmin = np.array(min_state.getPositions().value_in_unit(unit.nanometer))
    for n, i in enumerate(ridx):
        force.setParticleParameters(n, i, [float(v) for v in pmin[i]])

    # ---------- 2 heat (NVT, 0 -> 300 K) ----------
    if not (out / "02_heated.xml").exists():
        sim, integ = make_sim(st, system, dt_equil, args.platform)
        load_state(sim, out / "01_minimised.xml")
        sim.context.setVelocitiesToTemperature(HEAT_INCR * unit.kelvin, args.seed)
        sim.reporters.append(app.StateDataReporter(str(out / "02_heat.log"), 1000, step=True, time=True,
                                                   potentialEnergy=True, temperature=True, speed=True, separator="\t"))
        for s in range(1, HEAT_STAGES + 1):
            T = s * HEAT_INCR * unit.kelvin
            integ.setTemperature(T)
            sim.step(HEAT_STEPS)
        save_state(sim, out / "02_heated.xml")
        log(f"heated to {T} in {HEAT_STAGES} stages of {HEAT_STEPS} steps")
        del sim, integ

    # ---------- 3 release restraints (NPT) ----------
    system.addForce(mm.MonteCarloBarostat(PRESSURE, TEMPERATURE, 25))
    if not (out / "03_released.xml").exists():
        sim, integ = make_sim(st, system, dt_equil, args.platform)
        load_state(sim, out / "02_heated.xml")
        sim.reporters.append(app.StateDataReporter(str(out / "03_release.log"), 1000, step=True, time=True,
                                                   potentialEnergy=True, temperature=True, density=True, volume=True,
                                                   speed=True, separator="\t"))
        k0 = RESTRAINT_K.value_in_unit(unit.kilocalorie_per_mole / unit.angstrom**2)
        for scale in RELEASE_SCALES:
            sim.context.setParameter("k_res", k0 * scale * unit.kilocalorie_per_mole / unit.angstrom**2)
            log(f"release: k_res = {k0 * scale:.1f} kcal/mol/A^2 for {RELEASE_STEPS} steps")
            sim.step(RELEASE_STEPS)
        state = save_state(sim, out / "03_released.xml")
        app.PDBFile.writeFile(sim.topology, state.getPositions(), open(out / "03_released.pdb", "w"), keepIds=True)
        del sim, integ

    # ---------- 4 production ----------
    sim, integ = make_sim(st, system, dt_prod, args.platform)
    sim.context.setParameter("k_res", 0.0)
    total = int(round(args.ns * unit.nanosecond / dt_prod))
    interval = int(round(args.report_ps * unit.picosecond / dt_prod))
    ckpt = out / "state.chk"
    done = int(meta.get("steps_done", 0))
    if ckpt.exists() and done > 0:
        sim.loadCheckpoint(str(ckpt))
        log(f"resumed production at {done} steps ({done * dt_prod.value_in_unit(unit.nanosecond):.2f} ns)")
    else:
        load_state(sim, out / "03_released.xml")
        done = 0
        log("production started from the released state")
    if done >= total:
        log("production already complete")
        return
    append = done > 0
    sim.reporters.append(app.DCDReporter(str(out / "traj.dcd"), interval, append=append))
    sim.reporters.append(app.StateDataReporter(str(out / "04_production.log"), interval, step=True, time=True,
                                               potentialEnergy=True, kineticEnergy=True, temperature=True, volume=True,
                                               density=True, speed=True, remainingTime=True, totalSteps=total,
                                               separator="\t", append=append))
    chunk = interval * 10          # checkpoint every 1 ns
    t0, d0 = time.time(), done
    while done < total:
        n = min(chunk, total - done)
        sim.step(n)
        done += n
        sim.saveCheckpoint(str(ckpt))
        ns_per_day = (done - d0) * dt_prod.value_in_unit(unit.nanosecond) / max(1e-9, (time.time() - t0) / 86400)
        meta.update({"steps_done": done, "steps_total": total, "ns_done": round(done * dt_prod.value_in_unit(unit.nanosecond), 4),
                     "ns_total": args.ns, "dt_fs": dt_prod.value_in_unit(unit.femtosecond), "hmr": args.hmr, "seed": args.seed,
                     "ns_per_day_this_job": round(ns_per_day, 2), "n_particles": len(st.atoms)})
        meta_path.write_text(json.dumps(meta, indent=2))
    log(f"production complete: {done} steps = {args.ns} ns")


if __name__ == "__main__":
    sys.exit(main())
