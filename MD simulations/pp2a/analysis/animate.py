#!/usr/bin/env python3
"""Render trajectory animations with PyMOL, then encode with ffmpeg.

The sfGFP renderer (../../sfGFP/analysis/animate.py) with its three views rewritten for
a HEAT solenoid carrying a single strand: the PML scaffold, frame loop, ffmpeg settings
and still rendering are imported from it unchanged.

  overview   the horseshoe, the strand across its mouth, the two tethers
  sites      close-up of the two adducts with the Cb-Cb distance drawn
  strain     cartoon putty painted by per-residue C-alpha RMSF (real values: the sfGFP
             version computed the RMSF and never wrote it into the B-factor column)

Water and ions are stripped and frames are superposed on the core C-alpha (14-589) so
the animation shows the protein moving under the strand rather than the box tumbling.
Residue numbers are restored from residue_map.json (UniProt numbering for the protein;
the DNA chain is offset to 301-350 with its linkers at 300 and 351).
"""
from __future__ import annotations

import argparse
import shutil
import subprocess
import sys
from pathlib import Path

import mdtraj as md
import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from common import SFGFP_ANALYSIS, SITES, CORE, load_resmap, load as load_dedup  # noqa: E402
sys.path.insert(0, str(SFGFP_ANALYSIS))
import animate as sf  # noqa: E402

SA, SB = SITES


def prepare(traj_path, top_path, out: Path, n_frames: int, resmap=None, stride_first=1):
    t = load_dedup(traj_path, top_path, stride=stride_first)
    top = t.topology
    if resmap and Path(resmap).exists():
        sf.renumber(top, resmap)
    else:
        for r in top.residues:          # prmtop: 0-based -> UniProt
            r.resSeq = r.index + 1
    rmap = load_resmap(top, resmap)
    keep = [a.index for a in top.atoms
            if a.residue.name not in ("WAT", "HOH", "Na+", "Cl-", "MG", "K+")]
    t = t.atom_slice(keep)
    top = t.topology
    ca = {}
    for i, (num, name, chain) in rmap.items():
        if chain == "A" and name not in ("DA", "DT", "DG", "DC"):
            pass
    core = [a.index for a in top.atoms if a.name == "CA" and a.residue.chain.index == 0
            and CORE[0] <= a.residue.resSeq <= CORE[1] and a.residue.resSeq < 300]
    core = np.array(core)
    t.superpose(t, 0, atom_indices=core)
    step = max(1, t.n_frames // n_frames)
    t = t[::step]
    ca_all = np.array([a.index for a in top.atoms if a.name == "CA" and a.residue.resSeq < 300])
    rmsf = md.rmsf(t, t, 0, atom_indices=ca_all) * 10.0
    per_res = {top.atom(i).residue.index: float(v) for i, v in zip(ca_all, rmsf)}
    # the PDB B-factor field is (-10, 100): a floppy strand or tail can exceed 100 A of RMSF
    bfac = np.clip(np.array([per_res.get(a.residue.index, 0.0) for a in top.atoms]), 0.0, 99.9)
    t[0].save_pdb(str(out / "solute.pdb"), bfactors=bfac)
    t.save_dcd(str(out / "solute.dcd"))
    return t, rmsf


BODIES = {
    "overview": f"""
show cartoon, polymer.protein
color grey80, polymer.protein
spectrum resi, grey90 grey60, polymer.protein and name CA, minimum=1, maximum=589
show cartoon, polymer.nucleic
set cartoon_nucleic_acid_mode, 4
set cartoon_ring_mode, 3
color skyblue, polymer.nucleic
show sticks, resn TDP+TET+DNL+DNP+DNH and not hydro
color orange, resn TDP+TET+DNL+DNP+DNH and elem C
show spheres, (resn TDP+TET) and name CB
set sphere_scale, 0.45
color firebrick, (resn TDP+TET) and name CB
pseudoatom lab_dna, pos=[0,0,0]
python
from pymol import cmd
if cmd.count_atoms("polymer.nucleic"):
    cmd.pseudoatom("lab_dna", selection="polymer.nucleic", label="50 nt ssDNA spring")
    cmd.set("label_position", (0, 12, 10), "lab_dna")
cmd.pseudoatom("lab_a", selection="resi {SA} and name CB", label="D6X (Asp7)")
cmd.pseudoatom("lab_b", selection="resi {SB} and name CB", label="N552X (Asn553)")
cmd.set("label_position", (0, 10, 12), "lab_a")
cmd.set("label_position", (0, -12, 12), "lab_b")
python end
hide everything, lab_dna or lab_a or lab_b
show labels, lab_dna or lab_a or lab_b
orient polymer
zoom polymer, 6
turn x, -20
""",
    "sites": f"""
show cartoon, polymer.protein
set cartoon_transparency, 0.75
color grey80, polymer.protein
show cartoon, polymer.nucleic
color skyblue, polymer.nucleic
set cartoon_transparency, 0.5, polymer.nucleic
show sticks, (resn TDP+TET+DNL+DNP+DNH or (byres (polymer.nucleic within 4 of resn DNL+DNP))) and not hydro
color orange, resn TDP+TET and elem C
color gold, resn DNL+DNP+DNH and elem C
show sticks, (resi {SA}+{SB}) and not hydro
show spheres, (resn TDP+TET) and name CB
set sphere_scale, 0.4
color firebrick, (resn TDP+TET) and name CB
distance dCB, resi {SA} and name CB and polymer.protein, resi {SB} and name CB and polymer.protein
hide labels, dCB
python
from pymol import cmd
cmd.pseudoatom("lab_a", selection="resi {SA} and name CB and polymer.protein", label="Tet2-Et + sTCO at Asp7 (D6X)")
cmd.pseudoatom("lab_b", selection="resi {SB} and name CB and polymer.protein", label="Tet2-Et + sTCO at Asn553 (N552X)")
cmd.set("label_position", (0, 6, 10), "lab_a")
cmd.set("label_position", (0, -6, 10), "lab_b")
python end
hide everything, lab_a or lab_b
show labels, lab_a or lab_b
orient resn TDP+TET+DNL+DNP+DNH or (resi {SA}+{SB} and polymer.protein)
zoom resn TDP+TET+DNL+DNP+DNH or (resi {SA}+{SB} and polymer.protein), 8
""",
    "strain": f"""
show cartoon, polymer.protein
cartoon putty, polymer.protein
set cartoon_putty_scale_min, 0.6
set cartoon_putty_scale_max, 2.5
set cartoon_putty_transform, 0
spectrum b, white_yellow_orange_red, polymer.protein and name CA, minimum=0.5, maximum=6
show cartoon, polymer.nucleic
color grey70, polymer.nucleic
show spheres, (resn TDP+TET) and name CB
set sphere_scale, 0.45
color black, (resn TDP+TET) and name CB
python
from pymol import cmd
cmd.pseudoatom("lab_a", selection="resi {SA} and name CB and polymer.protein", label="Asp7 attachment")
cmd.pseudoatom("lab_b", selection="resi {SB} and name CB and polymer.protein", label="Asn553 attachment")
cmd.set("label_position", (0, 10, 12), "lab_a")
cmd.set("label_position", (0, -12, 12), "lab_b")
python end
hide everything, lab_a or lab_b
show labels, lab_a or lab_b
orient polymer.protein
zoom polymer.protein, 4
turn x, -20
""",
}
sf.BODIES = BODIES


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--traj", required=True)
    ap.add_argument("--top", required=True)
    ap.add_argument("--name", required=True)
    ap.add_argument("--out", default=str(HERE.parent / "figures"))
    ap.add_argument("--views", nargs="+", default=list(BODIES))
    ap.add_argument("--frames", type=int, default=80)
    ap.add_argument("--fps", type=int, default=12)
    ap.add_argument("--width", type=int, default=900)
    ap.add_argument("--height", type=int, default=700)
    ap.add_argument("--still-width", type=int, default=1500)
    ap.add_argument("--still-height", type=int, default=1150)
    ap.add_argument("--resmap", default=None)
    ap.add_argument("--stride-first", type=int, default=5,
                    help="pre-stride when loading (a 100 ns, 20 ps trajectory is 5000 frames)")
    ap.add_argument("--keep-work", action="store_true")
    a = ap.parse_args()
    pymol = shutil.which("pymol") or sys.exit("pymol not found")
    ffmpeg = shutil.which("ffmpeg") or sys.exit("ffmpeg not found")
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    work = out / f"_work_{a.name}"
    if work.exists():
        shutil.rmtree(work)
    work.mkdir(parents=True)
    print(f"=== {a.name} ===")
    t, rmsf = prepare(a.traj, a.top, work, a.frames, a.resmap, a.stride_first)
    print(f"  {t.n_frames} frames, {t.n_atoms} solute atoms, RMSF {rmsf.mean():.2f} A mean")
    made = []
    for v in a.views:
        made.append(sf.render(v, work, out, a.name, a.width, a.height, pymol, ffmpeg, a.fps))
        sf.render_still(v, work, out, a.name, a.still_width, a.still_height, pymol,
                        max(1, t.n_frames // 2))
    if not a.keep_work:
        shutil.rmtree(work)
    print("  wrote: " + ", ".join(m for m in made if m))


if __name__ == "__main__":
    main()
