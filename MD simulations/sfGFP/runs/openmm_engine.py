#!/usr/bin/env python3
"""OpenMM engine for the sfGFP-DNA spring chimera: minimise, equilibrate, produce.

One module, three stages, restart-safe.  Runs from Amber prmtop/inpcrd so that the
exotic residues (the CRO chromophore, the Tet2-Et amino acid, the click adduct and
the DNA-side linker, plus four inter-residue bonds tleap makes explicitly) come in
with a topology that is already fully parameterised.

Design notes
------------
* Hydrogen-mass repartitioning is done with ParmEd rather than
  ``createSystem(hydrogenMass=...)`` because ParmEd leaves *water* alone.  Water is
  already rigid, and inflating its hydrogens changes its dynamics -- and therefore
  its viscosity, which matters here: the observable is a sub-angstrom mechanical
  response, so solvent friction is not a free parameter.
* Production uses LangevinMiddleIntegrator, which has markedly better configurational
  sampling accuracy at long timesteps than the older LangevinIntegrator.
* The force clamp is a genuine constant force, ``E = -f*r`` between two atoms, so
  ``F = -dE/dr = +f`` pushes them apart -- the sign a compressed leaf spring applies.
  A moving harmonic restraint would instead impose a velocity and report a
  force that depends on the pulling rate.
* Restarts resume from a checkpoint and append, so ``gpu_requeue`` preemption is
  survivable.  Step accounting lives in a small JSON sidecar rather than being
  inferred from file sizes.
"""
from __future__ import annotations

import argparse
import json
import os
import sys
import time
from pathlib import Path

import openmm as mm
import openmm.app as app
import openmm.unit as unit
import parmed

# ---------------------------------------------------------------------------
# defaults, all in one place so the SLURM templates stay dumb
# ---------------------------------------------------------------------------
TEMPERATURE = 300.0 * unit.kelvin
PRESSURE = 1.0 * unit.bar
FRICTION = 1.0 / unit.picosecond
CUTOFF = 0.9 * unit.nanometer
HMR_MASS = 3.024 * unit.dalton      # factor-3 repartitioning, solute only
DT_EQUIL = 2.0 * unit.femtosecond
DT_PROD = 4.0 * unit.femtosecond
RESTRAINT_K = 10.0 * unit.kilocalorie_per_mole / unit.angstrom**2


def log(msg: str) -> None:
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


# ---------------------------------------------------------------------------
def load_structure(prmtop: str, inpcrd: str, hmr: bool = True):
    """Read the Amber topology, optionally repartitioning solute hydrogen masses."""
    st = parmed.load_file(prmtop, xyz=inpcrd)
    if hmr:
        # dowater=False: leave the 3-site water model's hydrogens at 1.008 Da.
        parmed.tools.HMassRepartition(
            st, str(HMR_MASS.value_in_unit(unit.dalton))).execute()
        log(f"HMR applied to solute hydrogens ({HMR_MASS}), water left untouched")
    return st


def build_system(st, *, rigid_water=True, dispersion_correction=True):
    system = st.createSystem(
        nonbondedMethod=app.PME,
        nonbondedCutoff=CUTOFF,
        constraints=app.HBonds,
        rigidWater=rigid_water,
        removeCMMotion=True,
    )
    for f in system.getForces():
        if isinstance(f, mm.NonbondedForce):
            f.setUseDispersionCorrection(dispersion_correction)
            f.setEwaldErrorTolerance(5e-4)
    return system


def add_positional_restraints(system, st, selection_mask, k=RESTRAINT_K):
    """Flat harmonic restraints on the selected atoms' crystallographic positions.

    ``periodicdistance`` keeps the restraint correct if an atom is wrapped across a
    periodic boundary mid-run, which a plain (x-x0)^2 would not.
    """
    force = mm.CustomExternalForce(
        "k_res*periodicdistance(x, y, z, x0, y0, z0)^2")
    force.addGlobalParameter("k_res", k)
    for p in ("x0", "y0", "z0"):
        force.addPerParticleParameter(p)
    n = 0
    for i, keep in enumerate(selection_mask):
        if not keep:
            continue
        pos = st.positions[i].value_in_unit(unit.nanometer)
        force.addParticle(i, [pos[0], pos[1], pos[2]])
        n += 1
    system.addForce(force)
    log(f"restrained {n} atoms at k = {k}")
    return force


def solute_heavy_mask(st):
    solvent = {"WAT", "HOH", "Na+", "Cl-", "K+", "MG", "Mg+2"}
    return [(a.residue.name not in solvent and a.atomic_number > 1) for a in st.atoms]


def add_force_clamp(system, atom_i: int, atom_j: int, force_pN: float):
    """Constant force pushing atoms i and j apart, in piconewtons.

    E = -f*r, so F = +f along the i->j axis on each atom, outward.  This is the
    load a compressed DNA leaf spring applies at its attachment points.
    """
    # 1 pN = 1e-12 N; in OpenMM's kJ/mol/nm: 1 pN = 0.6022 kJ/mol/nm
    f_openmm = force_pN * 0.60221408  # kJ/mol/nm
    clamp = mm.CustomBondForce("-f_clamp*r")
    clamp.addGlobalParameter("f_clamp", f_openmm)
    clamp.addBond(int(atom_i), int(atom_j), [])
    clamp.setUsesPeriodicBoundaryConditions(False)
    system.addForce(clamp)
    log(f"force clamp: {force_pN:.2f} pN between atoms {atom_i} and {atom_j} "
        f"({f_openmm:.4f} kJ/mol/nm)")
    return clamp


def make_simulation(st, system, integrator, platform_name="CUDA", precision="mixed"):
    if platform_name:
        platform = mm.Platform.getPlatformByName(platform_name)
        props = {}
        if platform_name == "CUDA":
            props = {"Precision": precision, "DeterministicForces": "false"}
        sim = app.Simulation(st.topology, system, integrator, platform, props)
    else:
        sim = app.Simulation(st.topology, system, integrator)
    log(f"platform: {sim.context.getPlatform().getName()}")
    return sim


# ---------------------------------------------------------------------------
def stage_minimise(st, outdir: Path, args):
    system = build_system(st)
    integrator = mm.LangevinMiddleIntegrator(TEMPERATURE, FRICTION, DT_EQUIL)
    sim = make_simulation(st, system, integrator, args.platform)
    sim.context.setPositions(st.positions)
    if st.box_vectors is not None:
        sim.context.setPeriodicBoxVectors(*st.box_vectors)
    e0 = sim.context.getState(getEnergy=True).getPotentialEnergy()
    log(f"initial potential energy: {e0}")
    sim.minimizeEnergy(maxIterations=args.min_steps,
                       tolerance=5.0 * unit.kilojoule_per_mole / unit.nanometer)
    state = sim.context.getState(getPositions=True, getEnergy=True)
    log(f"minimised potential energy: {state.getPotentialEnergy()}")
    app.PDBFile.writeFile(sim.topology, state.getPositions(),
                          open(outdir / "minimised.pdb", "w"), keepIds=True)
    with open(outdir / "minimised.xml", "w") as fh:
        fh.write(mm.XmlSerializer.serialize(state))
    return state


def stage_equilibrate(st, outdir: Path, args):
    """Restrained NVT heating, then staged release under NPT.

    NVT and NPT use two separate contexts rather than one context with the barostat
    toggled.  A MonteCarloBarostat cannot be "paused" by setting its pressure to
    zero -- at zero external pressure it happily drives the box toward collapse --
    and changing its frequency needs a context reinitialisation anyway.  Two
    contexts is the honest way to say "no barostat, then barostat".
    """
    mask = solute_heavy_mask(st)

    # ---- NVT, restrained, no barostat ----
    system = build_system(st)
    add_positional_restraints(system, st, mask)
    integrator = mm.LangevinMiddleIntegrator(TEMPERATURE, FRICTION, DT_EQUIL)
    sim = make_simulation(st, system, integrator, args.platform)

    prev = outdir.parent / "01_minimise" / "minimised.xml"
    if prev.exists():
        sim.context.setState(mm.XmlSerializer.deserialize(prev.read_text()))
        log(f"resumed from {prev}")
    else:
        sim.context.setPositions(st.positions)
        if st.box_vectors is not None:
            sim.context.setPeriodicBoxVectors(*st.box_vectors)
    sim.context.setVelocitiesToTemperature(TEMPERATURE, args.seed)
    sim.reporters.append(app.StateDataReporter(
        str(outdir / "nvt.log"), 500, step=True, time=True, potentialEnergy=True,
        temperature=True, speed=True, separator="\t"))
    nvt = int(args.nvt_ps * unit.picosecond / DT_EQUIL)
    log(f"NVT with restraints: {args.nvt_ps} ps ({nvt} steps)")
    sim.step(nvt)
    nvt_state = sim.context.getState(getPositions=True, getVelocities=True)
    del sim, integrator

    # ---- NPT, restraint ramped off ----
    # Reuse the same System object and just add the barostat to it.  Rebuilding it
    # costs ~75 s of ParmEd work for a 44k-atom system, twice per run, for nothing:
    # a System can be extended and handed to a second Context as long as the first
    # Context is gone.
    system.addForce(mm.MonteCarloBarostat(PRESSURE, TEMPERATURE, 100))
    integrator = mm.LangevinMiddleIntegrator(TEMPERATURE, FRICTION, DT_EQUIL)
    sim = make_simulation(st, system, integrator, args.platform)
    sim.context.setState(nvt_state)
    sim.reporters.append(app.StateDataReporter(
        str(outdir / "npt.log"), 500, step=True, time=True, potentialEnergy=True,
        temperature=True, density=True, volume=True, speed=True, separator="\t"))

    per_stage = int(args.npt_ps * unit.picosecond / DT_EQUIL / max(1, args.ramp_stages))
    k0 = RESTRAINT_K.value_in_unit(unit.kilocalorie_per_mole / unit.angstrom**2)
    for s in range(args.ramp_stages):
        k = k0 * (1.0 - (s + 1) / args.ramp_stages)
        sim.context.setParameter(
            "k_res", k * unit.kilocalorie_per_mole / unit.angstrom**2)
        log(f"NPT stage {s + 1}/{args.ramp_stages}: k_res = {k:.2f} kcal/mol/A^2, "
            f"{per_stage} steps")
        sim.step(per_stage)

    state = sim.context.getState(getPositions=True, getVelocities=True,
                                 enforcePeriodicBox=True)
    with open(outdir / "equilibrated.xml", "w") as fh:
        fh.write(mm.XmlSerializer.serialize(state))
    app.PDBFile.writeFile(sim.topology, state.getPositions(),
                          open(outdir / "equilibrated.pdb", "w"), keepIds=True)
    log("equilibration complete")
    return state


def stage_produce(st, outdir: Path, args):
    system = build_system(st)
    system.addForce(mm.MonteCarloBarostat(PRESSURE, TEMPERATURE, 100))
    if args.clamp_pN is not None:
        if args.clamp_atoms is None:
            sys.exit("--clamp-pN requires --clamp-atoms i,j")
        i, j = (int(v) for v in args.clamp_atoms.split(","))
        add_force_clamp(system, i, j, args.clamp_pN)

    integrator = mm.LangevinMiddleIntegrator(TEMPERATURE, FRICTION, DT_PROD)
    sim = make_simulation(st, system, integrator, args.platform)

    ckpt = outdir / "state.chk"
    meta_path = outdir / "progress.json"
    done = 0
    if ckpt.exists() and meta_path.exists():
        sim.loadCheckpoint(str(ckpt))
        done = json.loads(meta_path.read_text())["steps_done"]
        log(f"resumed from checkpoint at {done} steps "
            f"({done * DT_PROD.value_in_unit(unit.nanosecond):.2f} ns)")
    else:
        eq = outdir.parent / "02_equilibrate" / "equilibrated.xml"
        if not eq.exists():
            sys.exit(f"no checkpoint and no equilibrated state at {eq}")
        sim.context.setState(mm.XmlSerializer.deserialize(eq.read_text()))
        sim.context.setVelocitiesToTemperature(TEMPERATURE, args.seed)
        log(f"started from {eq} (seed {args.seed})")

    total = int(round(args.ns * unit.nanosecond / DT_PROD))
    interval = int(round(args.report_ps * unit.picosecond / DT_PROD))
    if done >= total:
        log(f"already at {done}/{total} steps; nothing to do")
        return
    append = done > 0
    sim.reporters.append(app.DCDReporter(str(outdir / "traj.dcd"), interval,
                                         append=append))
    sim.reporters.append(app.StateDataReporter(
        str(outdir / "production.log"), interval, step=True, time=True,
        potentialEnergy=True, kineticEnergy=True, temperature=True,
        volume=True, density=True, speed=True, remainingTime=True,
        totalSteps=total, separator="\t", append=append))
    sim.reporters.append(app.CheckpointReporter(str(ckpt), interval * 10))

    log(f"production: {args.ns} ns total, {total} steps at {DT_PROD}, "
        f"writing every {args.report_ps} ps")
    chunk = interval * 10
    while done < total:
        n = min(chunk, total - done)
        sim.step(n)
        done += n
        sim.saveCheckpoint(str(ckpt))
        meta_path.write_text(json.dumps({
            "steps_done": done, "steps_total": total,
            "ns_done": round(done * DT_PROD.value_in_unit(unit.nanosecond), 4),
            "dt_fs": DT_PROD.value_in_unit(unit.femtosecond),
            "clamp_pN": args.clamp_pN, "seed": args.seed}, indent=2))
    log(f"production complete: {done} steps")


# ---------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("stage", choices=["minimise", "equilibrate", "produce"])
    ap.add_argument("--prmtop", required=True)
    ap.add_argument("--inpcrd", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--platform", default="CUDA",
                    help="CUDA, OpenCL, CPU, or empty for automatic")
    ap.add_argument("--no-hmr", action="store_true")
    ap.add_argument("--min-steps", type=int, default=10000)
    ap.add_argument("--nvt-ps", type=float, default=100.0)
    ap.add_argument("--npt-ps", type=float, default=900.0)
    ap.add_argument("--ramp-stages", type=int, default=3)
    ap.add_argument("--ns", type=float, default=50.0)
    ap.add_argument("--report-ps", type=float, default=20.0)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--clamp-pN", type=float, default=None,
                    help="constant outward force between --clamp-atoms, in pN")
    ap.add_argument("--clamp-atoms", default=None, help="two 0-based atom indices, i,j")
    args = ap.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    log(f"host {os.uname().nodename}  stage {args.stage}  out {outdir}")
    st = load_structure(args.prmtop, args.inpcrd, hmr=not args.no_hmr)
    log(f"system: {len(st.atoms)} atoms, {len(st.residues)} residues, "
        f"box {st.box}")

    {"minimise": stage_minimise,
     "equilibrate": stage_equilibrate,
     "produce": stage_produce}[args.stage](st, outdir, args)


if __name__ == "__main__":
    main()
