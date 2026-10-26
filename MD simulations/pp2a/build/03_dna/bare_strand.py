#!/usr/bin/env python3
"""Solvate the 50 nt strand on its own, to measure the tension it exerts directly.

Same idea as the sfGFP bare-duplex systems (build/03_dna/build_spring_only.py): no
protein, no linkers, a harmonic restraint between the two atoms the linkers would attach
to -- the 5' end and the 3' end of the one strand -- and the mean restraint force read
off at a series of separations.  That is the force field's own force-extension curve,
to set against Marko-Siggia.  The sfGFP work never ran this for its ssDNA, so this is
the first direct test of the "slack, zero tension" statement.

The 5' phosphate is stripped (tleap's 5'-terminal template has none), so the 5' anchor
is O5' rather than P -- half an angstrom, irrelevant at 60-100 A.

The box has to hold the longest window: force_extension.py discards any window whose
separation exceeds 0.85 x (2 x min-image limit), and a truncated octahedron's limit is
about 0.4 x its box edge, so --pad is chosen for the longest x0 rather than the 12 A
the chimera uses.  Check build_report.json["box"] against the ladder before submitting.
"""
from __future__ import annotations

import argparse
import json
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
PP2A = HERE.parents[1]
SFGFP = PP2A.parent / "sfGFP"
sys.path.insert(0, str(SFGFP / "build" / "03_dna"))
from build_spring_only import strip_5p  # noqa: E402


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--spring", default=str(HERE / "spring_50nt.pdb"))
    ap.add_argument("--name", default="D_ss50")
    ap.add_argument("--pad", type=float, default=24.0)
    a = ap.parse_args()
    out = PP2A / "build" / "systems" / a.name
    out.mkdir(parents=True, exist_ok=True)
    solute = HERE / f"solute_{a.name}.pdb"
    info = strip_5p(Path(a.spring), solute)
    print(f"=== {a.name} ===  stripped {info['phosphate_atoms_dropped']} 5'-phosphate "
          f"atoms; anchors {info['anchor_5prime']} / {info['anchor_3prime']}")
    r = subprocess.run([sys.executable, str(SFGFP / "build/04_assemble/build_system.py"),
                        "--name", a.name, "--solute", str(solute), "--outdir", str(out),
                        "--mg-M", "0.005", "--pad", str(a.pad), "--expect-residue"],
                       capture_output=True, text=True)
    print("  " + "\n  ".join(r.stdout.strip().splitlines()[-4:]))
    if r.returncode:
        print(r.stdout[-2000:], r.stderr[-2000:])
        sys.exit("build_system failed")

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
        print(f"  {key}: {res.name}{want['resi']} {want['atom']} -> prmtop atom {atom.idx}")
    (out / "clamp_atoms.txt").write_text(",".join(str(i) for i in idx) + "\n")
    (out / "spring_only.json").write_text(json.dumps(info, indent=2))


if __name__ == "__main__":
    main()
