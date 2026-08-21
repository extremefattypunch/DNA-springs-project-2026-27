#!/usr/bin/env python3
"""Drive tleap to produce a solvated, neutralised prmtop/inpcrd for one system.

Two passes, because the ion count depends on the water count and the water count
is only known after solvation:

  pass 1  solvate, count waters, throw the result away
  pass 2  solvate again and add exactly the ions that give the target molality

Salt is specified as a concentration and converted with n = C * N_water / 55.5,
the standard "molality of the water box" convention.  Neutralising counterions are
added on top, so the requested concentration is the *added salt*, not the total
ionic strength -- stated here because the two differ for a protein at net -6e.

Every build ends with assertions on the resulting topology: integral net charge,
the residues that should exist do, and any inter-residue bond that tleap was asked
to make is present in the prmtop.  A silent failure here would produce a system
that runs fine and means nothing.
"""
from __future__ import annotations

import argparse
import json
import re
import shutil
import subprocess
import sys
from pathlib import Path

WATER_MOLALITY = 55.5      # mol water per litre of water

FORCE_FIELD_LINES = [
    "source leaprc.protein.ff14SB",
    "source leaprc.DNA.OL15",
    "source leaprc.water.tip3p",
    # The GFP chromophore.  Breyfogle et al., J Phys Chem B 2023, 127:5772.
    # Must come after the protein force field: it supplies only the GAFF half of
    # the chromophore description and relies on parm94-family types for the rest.
    "source leaprc.xFPchromophores",
    # Li/Merz 12-6 divalent parameters for Mg2+.  Deliberately the 12-6 set and
    # not 12-6-4: OpenMM has no native C4 term, so a 12-6-4 prmtop would silently
    # lose the polarisation correction on the way in.
    "loadamberparams frcmod.ions234lm_126_tip3p",
    # The custom residues are GAFF2-typed, so gaff2.dat has to be present for the
    # bonds, angles and torsions that cross their boundaries -- the per-residue
    # frcmods only cover what parmchk2 saw inside each capped model.
    "loadamberparams gaff2.dat",
]


def run_tleap(tleap: str, script: str, cwd: Path, tag: str) -> str:
    inp = cwd / f"leap_{tag}.in"
    inp.write_text(script)
    out = cwd / f"leap_{tag}.out"
    p = subprocess.run([tleap, "-f", inp.name], cwd=cwd, capture_output=True, text=True)
    out.write_text(p.stdout + p.stderr)
    text = out.read_text()
    fatal = [l for l in text.splitlines() if "FATAL" in l or "Failed to" in l]
    if fatal:
        print("\n".join(text.splitlines()[-40:]))
        sys.exit(f"tleap failed in pass '{tag}':\n  " + "\n  ".join(fatal))
    return text


def count_waters(text: str) -> int:
    m = re.findall(r"Added (\d+) residues", text)
    if m:
        return int(m[-1])
    m = re.findall(r"WAT\s+(\d+)", text)
    return int(m[-1]) if m else 0


def build(args) -> dict:
    out = Path(args.outdir)
    out.mkdir(parents=True, exist_ok=True)
    solute = Path(args.solute).resolve()
    shutil.copy(solute, out / solute.name)

    extra_params = "\n".join(f"loadamberparams {Path(p).name}" for p in args.frcmod)
    extra_libs = "\n".join(f"loadamberprep {Path(p).name}" for p in args.prep)
    for p in list(args.frcmod) + list(args.prep):
        shutil.copy(p, out / Path(p).name)
    bonds = "\n".join(f"bond mol.{a} mol.{b}" for a, b in
                      (s.split(":") for s in args.bond))

    head = "\n".join(FORCE_FIELD_LINES + [extra_params, extra_libs])
    common = f"""{head}
mol = loadpdb {solute.name}
{bonds}
"""

    # ---- pass 1: how much water does a 12 A shell hold? ----
    t1 = run_tleap(args.tleap, common + f"""
solvateoct mol TIP3PBOX {args.pad}
quit
""", out, "pass1")
    n_wat = count_waters(t1)
    if not n_wat:
        sys.exit("could not determine the water count from pass 1")
    n_salt = round(args.salt_M * n_wat / WATER_MOLALITY)
    n_mg = round(args.mg_M * n_wat / WATER_MOLALITY)
    print(f"  pass 1: {n_wat} waters -> {n_salt} NaCl pairs "
          f"({args.salt_M} M) + {n_mg} MgCl2 ({args.mg_M} M)")

    # ---- pass 2: solvate and ionise for real ----
    ions = [f"addionsrand mol Na+ 0"]           # neutralise first
    if n_mg:
        ions.append(f"addionsrand mol MG {n_mg} Cl- {2 * n_mg}")
    if n_salt:
        ions.append(f"addionsrand mol Na+ {n_salt} Cl- {n_salt}")
    t2 = run_tleap(args.tleap, common + f"""
solvateoct mol TIP3PBOX {args.pad}
{chr(10).join(ions)}
savepdb mol system.pdb
saveamberparm mol system.prmtop system.inpcrd
quit
""", out, "pass2")

    charge = re.findall(r"Total unperturbed charge:\s+(-?[\d.]+)", t2)
    report = {
        "system": args.name,
        "solute_pdb": str(solute),
        "pad_A": args.pad,
        "waters_pass1": n_wat,
        "nacl_pairs": n_salt,
        "mgcl2": n_mg,
        "salt_M_requested": args.salt_M,
        "mg_M_requested": args.mg_M,
        "leap_total_charge": float(charge[-1]) if charge else None,
        "close_contacts_reported": t2.count("Close contact"),
        "extra_bonds": args.bond,
    }
    return validate(out, report, args)


def validate(out: Path, report: dict, args) -> dict:
    import parmed
    p = parmed.load_file(str(out / "system.prmtop"), xyz=str(out / "system.inpcrd"))
    q = sum(a.charge for a in p.atoms)
    report.update({
        "atoms": len(p.atoms), "residues": len(p.residues), "bonds": len(p.bonds),
        "net_charge": round(q, 6),
        "box": [round(v, 3) for v in p.box] if p.box is not None else None,
        "waters": sum(1 for r in p.residues if r.name in ("WAT", "HOH")),
        "ions": {n: sum(1 for r in p.residues if r.name == n)
                 for n in ("Na+", "Cl-", "MG") if any(r.name == n for r in p.residues)},
    })
    problems = []
    if abs(q - round(q)) > 1e-4:
        problems.append(f"net charge {q} is not integral")
    if abs(round(q)) > 1e-6:
        problems.append(f"system is not neutral: net charge {round(q)}")

    # every residue that should carry parameters, does
    for resname in args.expect_residue:
        n = sum(1 for r in p.residues if r.name == resname)
        report.setdefault("expected_residues", {})[resname] = n
        if n == 0:
            problems.append(f"expected residue {resname} is absent")

    # every requested inter-residue bond survived into the topology
    made = set()
    for spec in args.bond:
        a, b = spec.split(":")
        ra, aa = a.split("@"); rb, ab = b.split("@")
        ia, ib = int(ra) - 1, int(rb) - 1
        try:
            at = next(x for x in p.residues[ia].atoms if x.name == aa)
            bt = next(x for x in p.residues[ib].atoms if x.name == ab)
        except (IndexError, StopIteration):
            problems.append(f"bond {spec}: atom not found")
            continue
        if any((bd.atom1 is at and bd.atom2 is bt) or (bd.atom1 is bt and bd.atom2 is at)
               for bd in p.bonds):
            made.add(spec)
        else:
            problems.append(f"bond {spec} is missing from the prmtop")
    report["bonds_verified"] = sorted(made)

    # a missing dihedral type would show up as a zero-parameter term
    report["dihedrals"] = len(p.dihedrals)
    report["angles"] = len(p.angles)
    report["problems"] = problems
    (out / "build_report.json").write_text(json.dumps(report, indent=2))

    print(f"  {report['atoms']} atoms, {report['residues']} residues, "
          f"net charge {report['net_charge']:+.4f}, box {report['box']}")
    print(f"  waters {report['waters']}, ions {report.get('ions')}")
    if problems:
        for x in problems:
            print(f"  PROBLEM: {x}")
    else:
        print("  all assertions passed")
    return report


SOLVENT_NAMES = {"WAT", "HOH", "DOD", "Na+", "Cl-", "K+", "MG", "NA", "CL"}


def solute_residue_order(solute_pdb: Path):
    """(residue number, residue name) per solute residue, in file order.

    Solvent is excluded deliberately.  tleap emits the prmtop as
    solute -> ions -> water, which is not the order the input PDB had (crystal
    waters came straight after the chain), so only the solute prefix has a
    position-for-position correspondence worth recording.
    """
    order, seen = [], set()
    for line in solute_pdb.read_text().splitlines():
        if line[:6] not in ("ATOM  ", "HETATM"):
            continue
        resn = line[17:20].strip()
        if resn in SOLVENT_NAMES:
            continue
        # Key on the chain too.  Both DNL linkers are residue 0, one per strand, so a
        # (number, name) key silently collapses them into one and shifts every
        # subsequent residue in the map by one.
        key = (line[21], int(line[22:26]), resn)
        if key not in seen:
            seen.add(key)
            order.append(key)
    return order


def write_residue_map(out: Path, solute_pdb: Path):
    """Map prmtop residue index -> the crystallographic numbering it came from.

    tleap renumbers residues sequentially and mdtraj further normalises HID/HIE/HIP
    to HIS, so neither the number nor the name in a prmtop can be trusted to identify
    a residue.  Writing the correspondence out at build time means analysis code
    never has to reconstruct it with offset arithmetic -- which is exactly where an
    off-by-one would hide, and this project already has one such trap in the .pse.
    """
    import parmed
    p = parmed.load_file(str(out / "system.prmtop"))
    order = solute_residue_order(solute_pdb)
    # tleap legitimately renames residues to their terminal or tautomer variants:
    # a chain-terminal DG becomes DG3, HIS becomes HID/HIE/HIP.  Accept those; a name
    # change that is not one of them means the map really is misaligned.
    def compatible(prmtop_name, input_name):
        if prmtop_name == input_name:
            return True
        if prmtop_name in (input_name + "3", input_name + "5"):
            return True
        return {prmtop_name, input_name} <= {"HIS", "HID", "HIE", "HIP"}

    rows, mismatched = [], []
    for i, (chain, num, name) in enumerate(order):
        res = p.residues[i]
        if not compatible(res.name, name):
            mismatched.append((i, res.name, name))
        # The chain is recorded because an Amber prmtop has no chain field at all:
        # mdtraj reads the whole system as one chain, so DNA strands cannot be told
        # apart downstream without it.
        rows.append({"index": i, "orig_chain": chain, "orig_resnum": num,
                     "orig_resname": name, "prmtop_resname": res.name})
    if mismatched:
        sys.exit(f"residue map mismatch (first 5): {mismatched[:5]}")
    (out / "residue_map.json").write_text(json.dumps(rows))
    print(f"  residue map: {len(order)} solute residues verified name-for-name")


def write_clamp_atoms(out: Path, solute_pdb: Path, sites: list[int], anchor="CB"):
    """Record the 0-based prmtop indices of the two attachment anchor atoms.

    Resolved by *ordinal position* among the solute residues, because tleap
    renumbers residues sequentially and the crystallographic numbering does not
    survive.  The residue name is checked, so a silent off-by-one cannot pass.
    """
    import parmed
    order = solute_residue_order(solute_pdb)
    p = parmed.load_file(str(out / "system.prmtop"))
    idx = []
    for resi in sites:
        pos = next(i for i, (_, n, _) in enumerate(order) if n == resi)
        want_name = order[pos][2]
        res = p.residues[pos]
        if res.name != want_name:
            sys.exit(f"clamp-atom resolution failed: residue ordinal {pos} is "
                     f"{res.name} in the prmtop but {want_name} in {solute_pdb.name}")
        atom = next((a for a in res.atoms if a.name == anchor), None)
        if atom is None:
            sys.exit(f"{res.name}{resi} has no {anchor} atom")
        idx.append(atom.idx)
        print(f"  clamp anchor: {res.name}{resi} {anchor} -> prmtop atom {atom.idx}")
    (out / "clamp_atoms.txt").write_text(",".join(str(i) for i in idx) + "\n")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--name", required=True)
    ap.add_argument("--solute", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--pad", type=float, default=12.0)
    ap.add_argument("--salt-M", type=float, default=0.150)
    ap.add_argument("--mg-M", type=float, default=0.0)
    ap.add_argument("--frcmod", nargs="*", default=[])
    ap.add_argument("--prep", nargs="*", default=[])
    ap.add_argument("--bond", nargs="*", default=[],
                    help="inter-residue bonds as RESI@ATOM:RESI@ATOM using the "
                         "solute PDB's own residue numbering (1-based ordinal)")
    ap.add_argument("--expect-residue", nargs="*", default=["CRO"])
    ap.add_argument("--clamp-sites", nargs="*", type=int, default=[],
                    help="two residue numbers whose CB atoms the force clamp uses")
    ap.add_argument("--tleap", default=shutil.which("tleap") or "tleap")
    args = ap.parse_args()
    if not shutil.which(args.tleap):
        sys.exit("tleap not found -- source env/activate.sh and use $DNASPRING_ENV/bin")

    print(f"=== building {args.name} ===")
    build(args)
    write_residue_map(Path(args.outdir), Path(args.solute))
    if args.clamp_sites:
        write_clamp_atoms(Path(args.outdir), Path(args.solute), args.clamp_sites)


if __name__ == "__main__":
    main()
