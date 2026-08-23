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
  miss the bonds that close the covalent loop.  Note that only some of these
  constructs ARE closed loops: L3/L4/L5, C1 and U2/U3 are, while S3/S4/S5 and C2/U4
  are two-arm or two-strand constructs by design.  Each file's REMARK says which.

Two files per system: the equilibrated starting structure and the last production
frame. Solvent is stripped by default; ``--with-solvent`` keeps it.
"""
from __future__ import annotations

import argparse
import json
import re
import shutil
from pathlib import Path

import mdtraj as md
import numpy as np
import parmed
import warnings

warnings.filterwarnings("ignore", module="parmed")

ROOT = Path(__file__).resolve().parent.parent
SOLVENT = {"WAT", "HOH", "Na+", "Cl-", "MG", "K+"}
# Every residue no viewer has a template for.  Missing one is not a crash, it is a
# silently mis-drawn molecule: without explicit connectivity PyMOL and ChimeraX guess
# bonds by distance, which mangles the fused polycyclic cages in particular.
NONSTANDARD = {"TET", "TDP", "DNL", "DNP", "DNH", "CRO", "CMS", "CDB"}
TITLE = {
    "S0_wt": "sfGFP, wild type (PDB 2B3P prepared)",
    "S1_tet": "sfGFP with two Tet2-Et residues at 133/149, unclicked",
    "S2_clicked": "sfGFP with both tethers clicked, no DNA (unloaded control)",
    "S3_spring27": "sfGFP-DNA chimera, 27 bp spring",
    "S4_spring40": "sfGFP-DNA chimera, 40 bp spring",
    "S5_spring40nick": "sfGFP-DNA chimera, 40 bp spring nicked at the centre",
    "S6_clamp": "sfGFP with two Tet2-Et residues (force-clamp topology)",
    # The corrected builds: ONE strand modified at both ends, so the chimera really is a
    # covalent circle.  S3-S5 above put a linker on each of two strands and are kept as
    # the two-arm construct they actually are.
    "L3_loop27": "sfGFP-DNA chimera, 27 bp spring, CLOSED covalent loop",
    "L4_loop40": "sfGFP-DNA chimera, 40 bp spring, closed covalent loop",
    "L5_loop40nick": "sfGFP-DNA chimera, 40 bp closed loop, complement nicked at centre",
    "C1_ss50": "sfGFP-DNA chimera, 50 nt single-stranded spring, closed covalent loop",
    "C2_ssarms": "sfGFP with two 25 nt single strands, non-circular control",
    "D_27bp": "27 bp duplex alone, for the force-extension measurement",
    "D_40bp": "40 bp duplex alone, for the force-extension measurement",
    "D_40bp_nick": "40 bp duplex alone, nicked, for the force-extension measurement",
    "U0_apo": "ubiquitin F45W/T9C/K63C, apo (PDB 1UBQ prepared)",
    "U1_clicked": "ubiquitin with both DBCO crosslinkers clicked, no DNA",
    "U2_ss50": "ubiquitin-DNA chimera, 50 nt ssDNA spring, closed covalent loop",
    "U3_ds50": "ubiquitin-DNA chimera, 50 bp dsDNA spring, closed covalent loop",
    "U4_ssarms": "ubiquitin with two 25 nt single strands, non-circular control",
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


# A run directory is the system name plus, optionally, the physical condition the run
# was done at -- md_chain.sbatch builds TAG="$SYS" + "_<f>pN" / "_x<x0>A" / "_k<k>".  The
# suffixes have to be matched precisely rather than with a trailing "*", because
# "D_40bp*" also matches "D_40bp_nick_x50A": pairing that trajectory with D_40bp's
# topology is a 77397-against-77407 atom mismatch, which is how this was found.
RUN_SUFFIX = re.compile(r"^(_[0-9.]+pN)?(_x[0-9.]+A)?(_k[0-9.]+)?$")


def runs_for(runs: Path, name: str, tail: str):
    """Run files belonging to exactly this system, not to one whose name extends it."""
    out = []
    for d in sorted(runs.glob(f"{name}*")):
        if not d.is_dir() or not RUN_SUFFIX.match(d.name[len(name):]):
            continue
        out += sorted(d.glob(tail))
    return out


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
    eq = runs_for(runs, name, "rep*/02_equilibrate/equilibrated.xml")
    if eq:
        sources.append(("start", eq[0], "equilibrated starting structure"))
    else:
        # A system that has not run yet still has a built structure worth exporting --
        # the ubiquitin set was queued when this was written.
        inp = sysdir / "system.inpcrd"
        if inp.exists():
            sources.append(("built", inp, "as built, before equilibration"))
    trajs = sorted(runs_for(runs, name, "rep*/03_production/traj.dcd"),
                   key=lambda q: q.stat().st_size, reverse=True)
    if trajs:
        sources.append(("final", trajs[0], "last frame of the longest production run"))

    for tag, src, note in sources:
        struct = parmed.load_file(str(top))
        n_expect = len(struct.atoms)
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
        elif src.suffix == ".inpcrd":
            rst = parmed.load_file(str(src))
            if len(rst.coordinates[0] if rst.coordinates.ndim == 3
                   else rst.coordinates) != n_expect:
                print(f"  SKIP {name} {tag}: inpcrd atom count mismatch")
                continue
            struct.coordinates = rst.coordinates
        elif src.suffix == ".xml":
            import openmm
            import openmm.unit as ommunit
            state = openmm.XmlSerializer.deserialize(src.read_text())
            xyz = np.array(state.getPositions().value_in_unit(ommunit.angstrom))
            # Belt and braces alongside runs_for(): if a source ever gets paired with the
            # wrong topology again, say so instead of raising a reshape error thirty
            # lines away from the cause.
            if len(xyz) != n_expect:
                print(f"  SKIP {name} {tag}: {src} has {len(xyz)} atoms, "
                      f"{top.name} expects {n_expect}")
                continue
            struct.coordinates = xyz
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
    show sticks, resn CRO+TET+TDP+DNL+DNP+DNH+CMS+CDB
    color grey80, polymer.protein
    color skyblue, polymer.nucleic
    color orange, resn TET+TDP+DNL+DNP+DNH+CMS+CDB
    color limegreen, resn CRO
    orient

ChimeraX:

    open S3_spring27_final.pdb
    cartoon; nucleotides stubs
    show :CRO,TET,TDP,DNL,DNP,DNH,CMS,CDB atoms; style :CRO,TET,TDP,DNL,DNP,DNH,CMS,CDB stick
    color /A grey; color /C,D,E,F cornflowerblue; color :TET,TDP,DNL,DNP,DNH,CMS,CDB orange; color :CRO green

To measure the spring the way the analysis does — the distance between the two 5′
phosphates, which is what sets the force:

    # PyMOL
    distance span, chain C and resi 1 and name P, chain D and resi 1 and name P

And the deformation coordinate:

    distance sites, chain A and resi 133 and name CB, chain A and resi 149 and name CB
"""

if __name__ == "__main__":
    main()
