#!/usr/bin/env python3
"""antechamber -> prepgen -> backbone harmonisation -> parmchk2, for TET/TDP/DNL.

Charge model
------------
AM1-BCC for the side chains.  That is the self-consistent choice for GAFF2-typed
atoms: GAFF's parameterisation targets AM1-BCC charges, which are themselves fitted
to reproduce HF/6-31G(d) RESP.  Mixing an independently derived RESP set into GAFF
torsions buys nothing here.

Backbone harmonisation
----------------------
The six backbone atoms are then overwritten with ff14SB's types and charges.  Without
this, the residue's own backbone would carry GAFF types (ns, c3, c, o) and AM1-BCC
charges while both of its peptide neighbours carry ff14SB ones, so the two peptide
bonds into and out of the modified site would be described inconsistently -- and one
of those sites, Asn149, is four residues from His148 and the chromophore.  The
residual charge is spread evenly over the side chain so the residue stays exactly
neutral.  Same idea as the "amide fix" used for the xFP chromophore parameters.

frcmod audit
------------
parmchk2 marks parameters it had to guess with "ATTN, need revision".  For a
1,2,4,5-tetrazine and a 4,5-dihydropyridazine it will guess, so the audit is
reported per residue instead of being buried in a file nobody reads.
"""
from __future__ import annotations

import argparse
import json
import os
import shutil
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from model_builder import FORMAL_CHARGE, build as build_model  # noqa: E402
from residue_defs import RESIDUES, formula, mass, validate  # noqa: E402

# ff14SB backbone: (Amber atom type, charge).  These are the values shared by the
# generic residues; their sum is -0.0984, so a neutral residue's side chain must
# sum to +0.0984.
FF14SB_BACKBONE = {
    "N":  ("N",  -0.4157), "H":  ("H",   0.2719),
    "CA": ("CX",  0.0337), "HA": ("H1",  0.0823),
    "C":  ("C",   0.5973), "O":  ("O",  -0.5679),
}
BACKBONE_SUM = sum(q for _, q in FF14SB_BACKBONE.values())

MC = {
    "TET": {"pre_head": "C", "post_tail": "N", "aa": True},
    "TDP": {"pre_head": "C", "post_tail": "N", "aa": True},
    # DNL is not an amino acid: it runs from the carbamate nitrogen (whose partner is
    # TDP's GAFF carbonyl carbon, type c) to the bridging oxygen that bonds the first
    # nucleotide's phosphorus (Amber DNA type P).
    "DNL": {"pre_head": "c", "post_tail": "P", "aa": False},
    "DNH": {"pre_head": "c", "post_tail": None, "aa": False},
}


def main_chain_path(res: dict) -> list:
    """Interior heavy atoms on the shortest head->tail path through the bonds.

    Derived from the tables rather than listed by hand.  A hand-written list is a
    trap: DNL's carbons were renamed CL1..CL6 -> C1..C6 to stop antechamber reading
    them as chlorine, and the stale list left prepgen unable to find any main-chain
    atom.  It silently fell back to making every one of them a child of the head
    nitrogen, so the bridging oxygen came out bonded to N instead of C6 -- a residue
    with the right atoms, the right charge, and the wrong topology.
    """
    adj = {}
    for a, b, _ in res["bonds"]:
        adj.setdefault(a, set()).add(b)
        adj.setdefault(b, set()).add(a)
    head, tail = res["head"], res["tail"]
    prev, queue = {head: None}, [head]
    while queue:
        x = queue.pop(0)
        if x == tail:
            break
        for y in sorted(adj.get(x, ())):
            if y not in prev:
                prev[y] = x
                queue.append(y)
    if tail not in prev:
        sys.exit(f"{res['name']}: no bonded path from {head} to {tail}")
    path, x = [], tail
    while x is not None:
        path.append(x)
        x = prev[x]
    path.reverse()
    return path[1:-1]          # interior atoms only; prepgen infers head and tail
PREP_FMT = ("{idx:>4}  {name:<4}  {typ:<4}  {topo:<1}  {na:>3} {nb:>3} {nc:>3}  "
            "{r:>8.3f}  {th:>8.3f}  {phi:>9.3f} {q:>9.6f}")


def run(cmd, cwd, logname):
    p = subprocess.run([str(c) for c in cmd], cwd=cwd, capture_output=True, text=True)
    (Path(cwd) / logname).write_text(p.stdout + p.stderr)
    if p.returncode:
        print((p.stdout + p.stderr)[-2500:])
        sys.exit(f"failed: {' '.join(str(c) for c in cmd)}")
    return p.stdout


def parse_prep(path: Path):
    """Return (lines, indices of atom records). Prep atom records have 11 fields."""
    lines = path.read_text().splitlines()
    rows = [i for i, l in enumerate(lines)
            if len(l.split()) == 11 and l.split()[0].isdigit()
            and l.split()[1] != "DUMM"]
    if not rows:
        sys.exit(f"{path}: no atom records found")
    return lines, rows


def harmonise_backbone(prepin: Path) -> dict:
    lines, rows = parse_prep(prepin)

    def put(i, f):
        lines[i] = PREP_FMT.format(idx=int(f[0]), name=f[1], typ=f[2], topo=f[3],
                                   na=int(f[4]), nb=int(f[5]), nc=int(f[6]),
                                   r=float(f[7]), th=float(f[8]), phi=float(f[9]),
                                   q=float(f[10]))

    before, changed, side = {}, [], []
    for i in rows:
        f = lines[i].split()
        if f[1] in FF14SB_BACKBONE:
            typ, q = FF14SB_BACKBONE[f[1]]
            before[f[1]] = {"gaff_type": f[2], "am1bcc_charge": round(float(f[10]), 6)}
            f[2], f[10] = typ, f"{q:.6f}"
            put(i, f)
            changed.append(f[1])
        else:
            side.append(i)

    side_sum = sum(float(lines[i].split()[10]) for i in side)
    target = -BACKBONE_SUM
    delta = (target - side_sum) / len(side)
    for i in side:
        f = lines[i].split()
        f[10] = f"{float(f[10]) + delta:.6f}"
        put(i, f)

    prepin.write_text("\n".join(lines) + "\n")
    total = sum(float(lines[i].split()[10]) for i in rows)
    return {"atoms_reset": changed, "replaced": before,
            "sidechain_atoms": len(side),
            "sidechain_sum_before": round(side_sum, 6),
            "sidechain_target": round(target, 6),
            "per_atom_shift": round(delta, 6),
            "residue_total_charge_after": round(total, 6)}


# GAFF type -> element.  Only the multi-letter halogen/metal types are ambiguous;
# everything else is keyed by its first character.
_GAFF_MULTI = {"cl": "Cl", "br": "Br", "cu": "Cu", "fe": "Fe", "zn": "Zn",
               "mg": "Mg", "na": "Na", "li": "Li", "ca": "C", "cc": "C", "cd": "C",
               "ce": "C", "cf": "C", "cg": "C", "ch": "C", "cp": "C", "cq": "C",
               "cx": "C", "cy": "C", "cu1": "Cu", "cs": "C", "c1": "C", "c2": "C",
               "c3": "C", "c5": "C", "c6": "C", "ni": "Ni", "co": "Co", "ns": "N"}


def element_from_gaff(t: str) -> str:
    """Element implied by a GAFF atom type."""
    tl = t.lower()
    if tl in _GAFF_MULTI:
        return _GAFF_MULTI[tl]
    if tl in ("f", "i"):
        return tl.upper()
    return t[0].upper()


def assert_elements(prepin: Path, res: dict) -> dict:
    """Every atom's assigned type must imply the element the tables declare.

    This is the guard for a whole class of silent failure: antechamber infers
    elements from atom *names*, so a badly chosen name turns a carbon into a
    chlorine and the run proceeds without complaint.  Checking the round trip
    name -> declared element vs. type -> implied element catches it at build time.
    """
    from residue_defs import element_of
    lines, rows = parse_prep(prepin)
    declared = {a: element_of(a) for a in res["heavy"]}
    for hs in res["hydrogens"].values():
        declared.update({h: "H" for h in hs})
    bad = []
    for i in rows:
        f = lines[i].split()
        name, typ = f[1], f[2]
        want = declared.get(name)
        if want is None:
            bad.append((name, typ, "atom not in the residue tables"))
            continue
        got = element_from_gaff(typ)
        if got != want:
            bad.append((name, typ, f"type implies {got}, tables say {want}"))
    if bad:
        for n, ty, why in bad:
            print(f"    ELEMENT MISMATCH {n} (type {ty}): {why}")
        sys.exit(f"{res['name']}: {len(bad)} atom(s) were assigned the wrong element")
    return {"atoms_checked": len(rows), "element_mismatches": 0}


def verify_tree(prepin: Path, res: dict) -> None:
    """Every prep tree parent link must be a real bond in the residue tables."""
    lines, rows = parse_prep(prepin)
    names = {}
    for i in rows:
        f = lines[i].split()
        names[int(f[0])] = f[1]
    # DUMM atoms occupy the first three indices and have no chemical meaning
    for i in rows:
        f = lines[i].split()
        child, parent = f[1], names.get(int(f[4]))
        if parent is None:
            continue
        bonds = {frozenset((a, b)) for a, b, _ in res["bonds"]}
        bonds |= {frozenset((h, par)) for par, hs in res["hydrogens"].items()
                  for h in hs}
        bonded = frozenset((child, parent)) in bonds
        if not bonded and parent not in ("DUMM",):
            sys.exit(f"{res['name']}: prep tree says {child} hangs off {parent}, "
                     f"but that is not a bond in the residue tables")
    print(f"  prep tree: every parent link matches a declared bond")


def audit_frcmod(path: Path) -> dict:
    flagged, section = [], None
    for line in path.read_text().splitlines():
        s = line.strip()
        if s in ("MASS", "BOND", "ANGLE", "DIHE", "IMPROPER", "NONBON"):
            section = s
            continue
        if "ATTN" in line:
            flagged.append({"section": section, "line": " ".join(line.split())})
    by_section = {}
    for f in flagged:
        by_section[f["section"]] = by_section.get(f["section"], 0) + 1
    return {"n_flagged": len(flagged), "by_section": by_section, "flagged": flagged}


def write_mc(name: str, cap_atoms, path: Path):
    res, cfg = RESIDUES[name], MC[name]
    out = [f"HEAD_NAME {res['head']}"]
    if res.get("tail"):
        out.append(f"TAIL_NAME {res['tail']}")
        out += [f"MAIN_CHAIN {a}" for a in main_chain_path(res)]
    out.append(f"PRE_HEAD_TYPE {cfg['pre_head']}")
    if cfg.get("post_tail"):
        out.append(f"POST_TAIL_TYPE {cfg['post_tail']}")
    out += [f"CHARGE {float(res['net_charge']):.4f}"]
    out += [f"OMIT_NAME {a}" for a in cap_atoms]
    path.write_text("\n".join(out) + "\n")


def build_one(name: str, outdir: Path, args) -> dict:
    res = RESIDUES[name]
    if (errs := validate(res)):
        sys.exit(f"{name}: definition invalid: {errs}")
    work = outdir / name
    work.mkdir(parents=True, exist_ok=True)
    print(f"\n=== {name} ===\n  {res['description']}")
    print(f"  target: {formula(res)}, {mass(res):.2f} Da, net charge {res['net_charge']}")

    mol2, caps, fc, smiles = build_model(name, work)
    write_mc(name, caps, work / f"{name}.mc")

    ac = work / f"{name}.ac"
    if ac.exists() and not args.force:
        print(f"  reusing {ac.name} (AM1-BCC already done; --force to redo)")
    else:
        print("  antechamber (AM1-BCC via sqm)...")
        run([args.antechamber, "-i", mol2.name, "-fi", "mol2", "-o", ac.name,
             "-fo", "ac", "-c", "bcc", "-nc", str(fc), "-at", "gaff2",
             "-rn", name, "-pf", "y", "-dr", "no"], work, "antechamber.log")

    prepin = work / f"{name}.prepin"
    run([args.prepgen, "-i", ac.name, "-o", prepin.name, "-f", "prepi",
         "-m", f"{name}.mc", "-rn", name], work, "prepgen.log")
    _, rows = parse_prep(prepin)
    chain = (f"{res['head']} -> {'-'.join(main_chain_path(res))} -> {res['tail']}"
             if res.get("tail") else f"{res['head']} (no tail: terminal residue)")
    print(f"  prepgen: {len(rows)} atoms retained ({len(caps)} cap atoms omitted); "
          f"main chain {chain}")
    verify_tree(prepin, res)
    elem_check = assert_elements(prepin, res)
    print(f"  element check: {elem_check['atoms_checked']} atoms, "
          f"every assigned type matches the declared element")

    report = {"residue": name, "description": res["description"],
              "smiles_from_tables": smiles, "formula": formula(res),
              "mass_Da": round(mass(res), 3), "net_charge": res["net_charge"],
              "model_formal_charge": fc, "cap_atoms": caps,
              "atoms_in_residue": len(rows), "charge_method": "AM1-BCC (GAFF2 types)",
              "element_check": elem_check}

    if MC[name]["aa"]:
        h = harmonise_backbone(prepin)
        report["backbone_harmonisation"] = h
        print(f"  backbone -> ff14SB types/charges; side chain shifted by "
              f"{h['per_atom_shift']:+.6f} e/atom over {h['sidechain_atoms']} atoms; "
              f"residue total {h['residue_total_charge_after']:+.6f}")

    frcmod = work / f"{name}.frcmod"
    run([args.parmchk2, "-i", prepin.name, "-f", "prepi", "-o", frcmod.name,
         "-s", "gaff2", "-a", "Y"], work, "parmchk2.log")
    report["frcmod_audit"] = audit_frcmod(frcmod)
    a = report["frcmod_audit"]
    print(f"  frcmod: {a['n_flagged']} parameters flagged for review {a['by_section']}")
    for f in a["flagged"][:6]:
        print(f"    [{f['section']}] {f['line'][:96]}")
    return report


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", default=str(HERE))
    ap.add_argument("--residues", nargs="+", default=["TET", "TDP", "DNL"])
    ap.add_argument("--force", action="store_true", help="redo the AM1-BCC step")
    for tool in ("antechamber", "prepgen", "parmchk2"):
        ap.add_argument(f"--{tool}", default=shutil.which(tool) or tool)
    args = ap.parse_args()
    for tool in (args.antechamber, args.prepgen, args.parmchk2):
        if not shutil.which(tool):
            sys.exit(f"{tool} not on PATH -- source env/activate.sh and add "
                     "$DNASPRING_ENV/bin")
    if "AMBERHOME" not in os.environ:
        os.environ["AMBERHOME"] = str(Path(shutil.which(args.antechamber)).parents[1])

    reports = [build_one(n, Path(args.outdir), args) for n in args.residues]
    out = Path(args.outdir) / "PARAMS_REPORT.json"
    out.write_text(json.dumps(reports, indent=2))
    print(f"\nwrote {out}")
    print(f"total parameters needing review: "
          f"{sum(r['frcmod_audit']['n_flagged'] for r in reports)}")


if __name__ == "__main__":
    main()
