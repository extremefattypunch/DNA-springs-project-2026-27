#!/usr/bin/env python3
"""Export viewable PDBs of every construct, for opening off the cluster.

What makes these different from the assembly intermediates in build/04_assemble:

* **They have hydrogens and real geometry.** The assembly files are heavy-atom
  scaffolds with linker bonds still ~2 Å long; these come out of the simulation, so
  every bond is at its equilibrium length.
* **Crystallographic numbering is restored.** tleap renumbers residues sequentially,
  so the prmtop calls the attachment sites 131 and 147. These files call them
  Asp133/Asn149 again, and put the protein, each DNA strand and each linker on its own
  chain, so a selection like ``chain A and resi 133`` means what you expect.
* **CONECT records for every non-standard residue.** No viewer has a template for
  TET, TDP, DNL or DNH, so without explicit connectivity PyMOL and ChimeraX guess
  bonds by distance — which quietly mis-draws the fused bicyclononane cage and can
  miss the two carbamate bonds that make the chimera a closed loop.

Two files per system: the equilibrated starting structure and the last production
frame. Solvent is stripped by default; ``--with-solvent`` keeps it.
"""
from __future__ import annotations

import argparse
import json
import shutil
from pathlib import Path

import mdtraj as md
import numpy as np
import parmed
import warnings

warnings.filterwarnings("ignore", module="parmed")

ROOT = Path(__file__).resolve().parent.parent
SOLVENT = {"WAT", "HOH", "Na+", "Cl-", "MG", "K+"}
NONSTANDARD = {"TET", "TDP", "DNL", "DNH", "CRO"}
TITLE = {
    "S0_wt": "sfGFP, wild type (PDB 2B3P prepared)",
    "S1_tet": "sfGFP with two Tet2-Et residues at 133/149, unclicked",
    "S2_clicked": "sfGFP with both tethers clicked, no DNA (unloaded control)",
    "S3_spring27": "sfGFP-DNA chimera, 27 bp spring",
    "S4_spring40": "sfGFP-DNA chimera, 40 bp spring",
    "S5_spring40nick": "sfGFP-DNA chimera, 40 bp spring nicked at the centre",
    "S6_clamp": "sfGFP with two Tet2-Et residues (force-clamp topology)",
}


def load_map(sysdir: Path):
    rows = json.loads((sysdir / "residue_map.json").read_text())
    return {r["index"]: r for r in rows}


def relabel(struct, resmap):
    """Restore crystallographic numbering and per-molecule chain IDs."""
    for res in struct.residues:
        row = resmap.get(res.idx)
        if row is None:                       # solvent, past the mapped prefix
            res.chain = "W" if res.name in ("WAT", "HOH") else "I"
            continue
        res.number = row["orig_resnum"]
        res.chain = row.get("orig_chain") or "A"
    return struct


def conect_lines(struct, serial_of):
    """CONECT records for every bond touching a non-standard residue."""
    partners = {}
    for b in struct.bonds:
        a1, a2 = b.atom1, b.atom2
        if a1.residue.name not in NONSTANDARD and a2.residue.name not in NONSTANDARD:
            continue
        for x, y in ((a1, a2), (a2, a1)):
            if x.idx in serial_of and y.idx in serial_of:
                partners.setdefault(serial_of[x.idx], []).append(serial_of[y.idx])
    out = []
    for s in sorted(partners):
        p = sorted(set(partners[s]))
        for i in range(0, len(p), 4):         # four partners per CONECT line
            out.append("CONECT" + f"{s:5d}" + "".join(f"{q:5d}" for q in p[i:i + 4]))
    return out


def write_pdb(struct, path: Path, title: str, note: str):
    # renumber=False keeps the crystallographic residue numbers we just restored, but
    # it also makes ParmEd write each atom's own `number` -- which a prmtop never sets,
    # so every serial came out as -1.  That breaks CONECT (all bonds collapse onto one
    # serial) and confuses viewers.  Assign the serials explicitly.
    for i, atom in enumerate(struct.atoms):
        atom.number = i + 1
    tmp = path.with_suffix(".tmp.pdb")
    struct.save(str(tmp), overwrite=True, renumber=False)
    body = [l for l in tmp.read_text().splitlines()
            if not l.startswith(("END", "CONECT", "MASTER"))]
    tmp.unlink()
    serial_of = {}
    for line in body:
        if line.startswith(("ATOM", "HETATM")):
            # ParmEd writes atoms in structure order, so the n-th record is atom n
            serial_of[len(serial_of)] = int(line[6:11])
    if len(set(serial_of.values())) != len(serial_of):
        raise RuntimeError(f"{path.name}: atom serials are not unique "
                           f"({len(set(serial_of.values()))} distinct for "
                           f"{len(serial_of)} atoms) -- CONECT would be meaningless")
    header = [f"TITLE     {title}", f"REMARK   1 {note}",
              "REMARK   1 CONECT records cover every bond involving TET/TDP/DNL/DNH/CRO,",
              "REMARK   1 including the carbamate and phosphoester bonds that close the loop."]
    path.write_text("\n".join(header + body + conect_lines(struct, serial_of)
                             + ["END"]) + "\n")
    return len(serial_of)


def export_system(name: str, sysdir: Path, runs: Path, out: Path, keep_solvent: bool):
    top = sysdir / "system.prmtop"
    if not (top.exists() and (sysdir / "residue_map.json").exists()):
        return []
    resmap = load_map(sysdir)
    made = []

    sources = []
    # The equilibrated starting structure comes from the serialised OpenMM State, not
    # from equilibrated.pdb.  Above 99,999 atoms OpenMM writes hybrid-36 atom serials
    # ("A0000"), which ParmEd's PDB reader rejects with thousands of warnings and no
    # coordinates -- and the two 121k-atom chimeras are exactly that size.  The XML
    # state is exact and format-proof.
    eq = sorted(runs.glob(f"{name}*/rep*/02_equilibrate/equilibrated.xml"))
    if eq:
        sources.append(("start", eq[0], "equilibrated starting structure"))
    # the last frame of the longest production run
    trajs = sorted(runs.glob(f"{name}*/rep*/03_production/traj.dcd"),
                   key=lambda p: p.stat().st_size, reverse=True)
    if trajs:
        sources.append(("final", trajs[0], "last frame of the longest production run"))

    for tag, src, note in sources:
        struct = parmed.load_file(str(top))
        if src.suffix == ".dcd":
            # mdtraj's DCD reader cannot seek to a negative index, so ask the file
            # how long it is first rather than counting on -1 to mean "last".
            with md.formats.DCDTrajectoryFile(str(src)) as fh:
                n_frames = len(fh)
            if n_frames == 0:
                continue
            t = md.load_frame(str(src), index=n_frames - 1, top=str(top))
            struct.coordinates = t.xyz[0] * 10.0
            ns = None
            prog = src.parent / "progress.json"
            if prog.exists():
                ns = json.loads(prog.read_text())["ns_done"]
            note = f"{note} ({ns:g} ns)" if ns else note
        elif src.suffix == ".xml":
            import openmm
            import openmm.unit as ommunit
            state = openmm.XmlSerializer.deserialize(src.read_text())
            struct.coordinates = np.array(
                state.getPositions().value_in_unit(ommunit.angstrom))
        else:
            struct.coordinates = md.load(str(src)).xyz[0] * 10.0
        if not keep_solvent:
            struct.strip(f":{','.join(SOLVENT)}")
        relabel(struct, resmap)
        p = out / f"{name}_{tag}.pdb"
        n = write_pdb(struct, p, TITLE.get(name, name),
                      f"{note}; solvent {'kept' if keep_solvent else 'stripped'}")
        made.append((p, n))
        print(f"  {p.name:<32} {n:>6} atoms  ({note})")
    return made


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", default=str(ROOT))
    ap.add_argument("--out", default=None)
    ap.add_argument("--with-solvent", action="store_true")
    a = ap.parse_args()
    root = Path(a.root)
    out = Path(a.out or root / "pdb_exports")
    out.mkdir(parents=True, exist_ok=True)
    runs = root / "data" / "runs"

    print("=== simulated constructs ===")
    made = []
    for sysdir in sorted((root / "build" / "systems").glob("*")):
        if sysdir.name == "S6_clamp":         # identical topology to S1_tet
            continue
        made += export_system(sysdir.name, sysdir, runs, out, a.with_solvent)

    print("\n=== isolated DNA springs (as built, before assembly) ===")
    for p in sorted((root / "build" / "03_dna").glob("spring_*.pdb")):
        shutil.copy(p, out / f"dna_{p.name}")
        n = sum(1 for l in p.read_text().splitlines() if l.startswith("ATOM"))
        print(f"  dna_{p.name:<28} {n:>6} atoms")

    print("\n=== the custom residues, as parameterised ===")
    for d in sorted((root / "build" / "02_params").glob("*/")):
        m = d / f"{d.name}_model.mol2"
        if m.exists():
            shutil.copy(m, out / f"residue_{d.name}_model.mol2")
            print(f"  residue_{d.name}_model.mol2")

    (out / "README.md").write_text(README)
    print(f"\nwrote {len(list(out.glob('*')))} files to {out}")
    tar = root / "pdb_exports.tar.gz"
    shutil.make_archive(str(tar).replace(".tar.gz", ""), "gztar",
                        root_dir=out.parent, base_dir=out.name)
    print(f"and {tar} ({tar.stat().st_size / 1e6:.1f} MB) for a single download")


README = """# Viewable structures

Every construct in the study, as coordinates you can open anywhere. Written by
`build/export_pdbs.py`.

## What is here

| file | what it is |
|---|---|
| `S0_wt_*.pdb` | wild-type sfGFP |
| `S1_tet_*.pdb` | two Tet2-Et residues at 133/149, unclicked — the fluorometry sample |
| `S2_clicked_*.pdb` | both tethers clicked, no DNA — the unloaded control |
| `S3_spring27_*.pdb` | the chimera with a 27 bp spring (~4.9 pN) |
| `S4_spring40_*.pdb` | the chimera with a 40 bp spring (~3.0 pN) |
| `S5_spring40nick_*.pdb` | the same, nicked at the centre (~2.4 pN) |
| `dna_spring_*.pdb` | the isolated duplexes, bent, before assembly |
| `residue_*_model.mol2` | the four custom residues as parameterised, with bond orders |

`_start` is the equilibrated structure at the beginning of production; `_final` is the
last frame of the longest run. Solvent is stripped — the solvated systems are 44k–121k
atoms and live in `build/systems/*/system.prmtop` with the trajectories on netscratch.

## Things worth knowing before you look

Residues use **2B3P numbering**, so the attachment sites are `Asp133` and `Asn149`
(D134/N150 in the construct's own numbering, which runs one higher). Chains: `A` the
protein, `C`/`D` the two DNA strands, `E`/`F` the capped tethers in `S2_clicked`.
The linker residues are `TDP` (the clicked Tet2-Et/sTCO adduct, part of the protein
chain) and `DNL` (the amino-C6 arm, leading each DNA strand); `TET` is the unclicked
tetrazine and `DNH` the tether capped as a free alcohol. `CRO` is the chromophore.

CONECT records are included for every bond involving those residues. Keep them — no
viewer has a template for the fused bicyclononane cage, and distance-based bond
guessing mis-draws it.

## Opening them

PyMOL:

    pymol S3_spring27_final.pdb
    # then, to see the construct the way the report shows it:
    hide everything
    show cartoon, polymer.protein
    show cartoon, polymer.nucleic
    set cartoon_ring_mode, 3
    show sticks, resn CRO+TDP+DNL
    color grey80, polymer.protein
    color skyblue, polymer.nucleic
    color orange, resn TDP+DNL
    color limegreen, resn CRO
    orient

ChimeraX:

    open S3_spring27_final.pdb
    cartoon; nucleotides stubs
    show :CRO,TDP,DNL atoms; style :CRO,TDP,DNL stick
    color /A grey; color /C,D cornflowerblue; color :TDP,DNL orange; color :CRO green

To measure the spring the way the analysis does — the distance between the two 5′
phosphates, which is what sets the force:

    # PyMOL
    distance span, chain C and resi 1 and name P, chain D and resi 1 and name P

And the deformation coordinate:

    distance sites, chain A and resi 133 and name CB, chain A and resi 149 and name CB
"""

if __name__ == "__main__":
    main()
