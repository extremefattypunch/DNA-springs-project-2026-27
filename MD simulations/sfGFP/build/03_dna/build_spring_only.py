#!/usr/bin/env python3
"""Solvate a spring on its own, for measuring the force it exerts directly.

The point of these systems is to answer "is the force from a model or from the force
field?" without a model in the way.  No protein and no linkers: a constant force is
applied between the two atoms the linkers would attach to -- the 5' end and the 3' end
of the *same* strand, which is the pair the corrected topology uses -- and the
equilibrium separation is measured.  Sweeping the force gives the spring's own
force-extension isotherm straight out of the force field.

The 5' phosphates the chimera needs are stripped here: tleap's 5'-terminal template
(DA5/DC5/...) has no phosphate, so leaving them in makes P/OP1/OP2 untyped and the
build dies.  Without them the 5' anchor is O5' instead of P, half an angstrom away,
which is irrelevant to the mechanics of a 90 A duplex.
"""
from __future__ import annotations

import argparse
import json
import subprocess
import sys
from pathlib import Path

# parents[2] is already the sfGFP directory: 03_dna -> build -> sfGFP
SF = Path(__file__).resolve().parents[2]
FIVE_PRIME_PHOSPHATE = {"P", "OP1", "OP2"}


def strip_5p(src: Path, dst: Path) -> dict:
    """Drop the 5'-terminal phosphate of every strand; report the anchor atoms."""
    lines = [l for l in src.read_text().splitlines()
             if l[:6] in ("ATOM  ", "HETATM")]
    first = {}
    last = {}
    for l in lines:
        ch, ri = l[21], int(l[22:26])
        first[ch] = min(first.get(ch, ri), ri)
        last[ch] = max(last.get(ch, ri), ri)
    out, dropped = [], 0
    for l in lines:
        ch, ri, an = l[21], int(l[22:26]), l[12:16].strip()
        if ri == first[ch] and an in FIVE_PRIME_PHOSPHATE:
            dropped += 1
            continue
        out.append(l)
    # keep the TER records so tleap sees two molecules
    body, prev = [], None
    for l in out:
        if prev is not None and l[21] != prev:
            body.append("TER")
        body.append(l)
        prev = l[21]
    dst.write_text("\n".join(body) + "\nTER\nEND\n")
    chains = sorted(first)
    return {"chains": chains, "phosphate_atoms_dropped": dropped,
            "anchor_5prime": {"chain": chains[0], "resi": first[chains[0]],
                              "atom": "O5'"},
            "anchor_3prime": {"chain": chains[0], "resi": last[chains[0]],
                              "atom": "O3'"}}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--spring", required=True, help="build/03_dna/spring_*.pdb")
    ap.add_argument("--name", required=True)
    ap.add_argument("--outdir", default=None)
    a = ap.parse_args()
    sf = SF
    out = Path(a.outdir or sf / "build" / "systems" / a.name)
    out.mkdir(parents=True, exist_ok=True)
    # Not inside `out`: build_system copies the solute into its own output
    # directory, and copying a file onto itself raises SameFileError.
    solute = SF / "build" / "03_dna" / f"solute_{a.name}.pdb"
    info = strip_5p(Path(a.spring), solute)
    print(f"=== {a.name} ===")
    print(f"  stripped {info['phosphate_atoms_dropped']} 5'-phosphate atoms; "
          f"anchors {info['anchor_5prime']} and {info['anchor_3prime']}")

    r = subprocess.run([sys.executable,
                        str(sf / "build/04_assemble/build_system.py"),
                        "--name", a.name, "--solute", str(solute),
                        "--outdir", str(out), "--mg-M", "0.005",
                        "--expect-residue"], capture_output=True, text=True)
    print("  " + "\n  ".join(r.stdout.strip().splitlines()[-4:]))
    if r.returncode:
        print(r.stdout[-2000:], r.stderr[-2000:])
        sys.exit("build_system failed")

    # resolve the two anchor atoms in the finished topology and record them
    import parmed
    p = parmed.load_file(str(out / "system.prmtop"))
    rows = json.loads((out / "residue_map.json").read_text())
    idx = []
    for key in ("anchor_5prime", "anchor_3prime"):
        want = info[key]
        row = next(r for r in rows if r["orig_chain"] == want["chain"]
                   and r["orig_resnum"] == want["resi"])
        res = p.residues[row["index"]]
        atom = next(x for x in res.atoms if x.name == want["atom"])
        idx.append(atom.idx)
        print(f"  {key}: {res.name}{want['resi']} chain {want['chain']} "
              f"{want['atom']} -> prmtop atom {atom.idx}")
    (out / "clamp_atoms.txt").write_text(",".join(str(i) for i in idx) + "\n")
    (out / "spring_only.json").write_text(json.dumps(info, indent=2))


if __name__ == "__main__":
    main()
