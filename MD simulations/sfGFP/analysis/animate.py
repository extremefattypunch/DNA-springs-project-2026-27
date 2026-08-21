#!/usr/bin/env python3
"""Render trajectory animations with PyMOL, then encode with ffmpeg.

Water and ions are stripped and the frames are superposed on the barrel core first,
so the animation shows the spring working against the protein rather than the whole
box tumbling.  PyMOL cannot read an Amber prmtop, so a solute-only PDB plus DCD pair
is written for it.

Three views, each answering a different question:
  overview     what the chimera looks like and how the duplex moves
  chromophore  what the chromophore and its hydrogen-bond partners do
  strain       where the protein moves, painted by per-residue RMSF
"""
from __future__ import annotations

import argparse
import shutil
import subprocess
import sys
from pathlib import Path

import mdtraj as md
import numpy as np


def prepare(traj_path, top_path, out: Path, n_frames: int, sites):
    t = md.load(str(traj_path), top=str(top_path))
    # Selected by iterating residues, not with a selection string: mdtraj's parser
    # chokes on the "+" in the Amber ion names Na+ and Cl-.
    solvent = {"WAT", "HOH", "Na+", "Cl-", "MG", "K+"}
    keep = np.array([a.index for a in t.topology.atoms
                     if a.residue.name not in solvent])
    t = t.atom_slice(keep)
    ca = t.topology.select("name CA")
    if len(ca):
        t.superpose(t, 0, atom_indices=ca)
    step = max(1, t.n_frames // n_frames)
    t = t[::step]
    t[0].save_pdb(str(out / "solute.pdb"))
    t.save_dcd(str(out / "solute.dcd"))
    rmsf = None
    if len(ca):
        rmsf = md.rmsf(t, t, 0, atom_indices=ca) * 10.0
    return t, rmsf


# Labels are set up once here rather than per view: a molecular figure that does not
# say which residue is which is decoration, not evidence.  White-backed outlined text
# survives being drawn over both the pale cartoon and the dark sticks, and
# label_position pushes the text toward the camera so side chains do not occlude it.
PML = """
load {pdb}, mol
load_traj {dcd}, mol
hide everything
bg_color white
set ray_opaque_background, 1
set antialias, 2
set ambient_occlusion_mode, 1
set ray_shadow, 0
set cartoon_transparency, 0.0
set movie_fps, 20
set label_size, {label_size}
set label_font_id, 7
set label_color, black
set label_outline_color, white
set label_bg_color, white
set label_bg_transparency, 0.15
set label_bg_outline, 1
set label_position, (0, 0, 3)
set label_distance_digits, 2
set dash_color, grey20
set dash_width, 3.0
set dash_gap, 0.28
set dash_length, 0.30
set dash_radius, 0.035
{body}
set all_states, off
python
from pymol import cmd
n = cmd.count_states("mol")
for i in range(1, n + 1):
    cmd.frame(i)
    cmd.png("{prefix}_%04d.png" % i, width={w}, height={h}, dpi=150, ray=0)
python end
"""

BODIES = {
    "overview": """
show cartoon, polymer.protein
color grey80, polymer.protein
show cartoon, polymer.nucleic
set cartoon_nucleic_acid_mode, 4
set cartoon_ring_mode, 3
set cartoon_ring_finder, 1
color skyblue, polymer.nucleic
show sticks, resn CRO
color limegreen, resn CRO
show sticks, resn TDP+DNL
color orange, resn TDP+DNL
show spheres, resn TDP and name CB
set sphere_scale, 0.38, resn TDP and name CB
color firebrick, resn TDP and name CB
python
from pymol import cmd
# Pseudoatoms at each component's centre of mass, so the parts are named on the
# picture instead of only in the caption.
parts = [("lbl_dna",  "polymer.nucleic",                     "27 bp DNA spring"),
         ("lbl_prot", "polymer.protein and not resn TDP",    "sfGFP beta-barrel"),
         ("lbl_cro",  "resn CRO",                            "chromophore"),
         ("lbl_teth", "resn TDP+DNL",                        "Tet2-Et + sTCO tether")]
for name, sele, text in parts:
    if cmd.count_atoms(sele):
        cmd.pseudoatom(name, pos=list(cmd.centerofmass(sele)), label=text)
cmd.set("label_color", "grey20", "lbl_prot")
python end
orient polymer
zoom all, 6
turn x, -15
""",
    "chromophore": """
show cartoon, polymer.protein
set cartoon_transparency, 0.82
color grey85, polymer.protein
show sticks, resn CRO
color limegreen, resn CRO
util.cnc("resn CRO")
color limegreen, resn CRO and elem C
show sticks, resi 148+203+205+96+94+222 and (sidechain or name CA)
color salmon, resi 148+203+205+96+94+222 and elem C
util.cnc("resi 148+203+205+96+94+222 and not elem C")
# the three bonds that hold the chromophore planar, plus the pair behind it
distance hb_his148,  resn CRO and name OH, resi 148 and name ND1
distance hb_thr203,  resn CRO and name OH, resi 203 and name OG1
distance hb_arg96,   resn CRO and name O2, resi 96  and name NH2
distance hb_gln94,   resn CRO and name O2, resi 94  and name NE2
distance hb_ser205,  resi 205 and name OG, resi 222 and name OE1
color grey20, hb_*
label resi 148 and name NE2, "His148"
label resi 203 and name CG2, "Thr203"
label resi 96  and name NH1, "Arg96"
label resi 94  and name OE1, "Gln94"
label resi 205 and name OG,  "Ser205"
label resi 222 and name OE2, "Glu222"
label resn CRO and name CZ,  "chromophore (CRO)"
orient resn CRO or resi 148+203+96
zoom resn CRO or resi 148+203+96, 3.5
turn x, -10
""",
    "strain": """
show cartoon, polymer.protein
show cartoon, polymer.nucleic
color grey70, polymer.nucleic
spectrum b, white_yellow_orange_red, polymer.protein and name CA
cartoon putty
set cartoon_putty_scale_min, 0.6
set cartoon_putty_scale_max, 2.5
set cartoon_putty_transform, 0
show sticks, resn CRO
color limegreen, resn CRO
show spheres, resn TDP and name CB
set sphere_scale, 0.4, resn TDP and name CB
color black, resn TDP and name CB
python
from pymol import cmd
for resi, text in ((133, "Asp133 attachment"), (149, "Asn149 attachment")):
    sele = f"resi {resi} and name CB"
    if cmd.count_atoms(sele):
        cmd.label(sele, f'"{text}"')
if cmd.count_atoms("resn CRO"):
    cmd.label("resn CRO and name CZ", '"chromophore"')
python end
orient polymer.protein
zoom all, 5
""",
}


STILL_TAIL = """
set all_states, off
frame {frame}
ray {w}, {h}
png {prefix}_still.png, dpi=200
"""


def render_still(view, work: Path, out: Path, name, w, h, pymol, frame):
    """One ray-traced, labelled frame.  A moving GIF is good for showing motion and
    bad for reading labels; the report needs both."""
    pml = work / f"{view}_still.pml"
    pml.write_text(PML.format(pdb="solute.pdb", dcd="solute.dcd", body=BODIES[view],
                              prefix="_x", w=w, h=h, label_size=20)
                   .split("set all_states, off")[0]
                   + STILL_TAIL.format(frame=frame, w=w, h=h, prefix=view))
    r = subprocess.run([pymol, "-cq", pml.name], cwd=work, capture_output=True,
                       text=True)
    src = work / f"{view}_still.png"
    if not src.exists():
        print(f"  {view} still: no output")
        print((r.stdout + r.stderr)[-1200:])
        return None
    dst = out / f"{name}_{view}_still.png"
    shutil.copy(src, dst)
    print(f"  {view} still: {dst.name} ({dst.stat().st_size // 1024} kB, {w}x{h})")
    return dst.name


def render(view, work: Path, out: Path, name, w, h, pymol, ffmpeg, fps):
    prefix = work / f"{view}"
    pml = work / f"{view}.pml"
    pml.write_text(PML.format(pdb="solute.pdb", dcd="solute.dcd",
                              body=BODIES[view], prefix=prefix.name, w=w, h=h,
                              label_size=15))
    r = subprocess.run([pymol, "-cq", pml.name], cwd=work, capture_output=True,
                       text=True)
    frames = sorted(work.glob(f"{view}_*.png"))
    if not frames:
        print(f"  {view}: PyMOL produced no frames")
        print((r.stdout + r.stderr)[-1500:])
        return None
    mp4 = out / f"{name}_{view}.mp4"
    gif = out / f"{name}_{view}.gif"
    subprocess.run([ffmpeg, "-y", "-loglevel", "error", "-framerate", str(fps),
                    "-i", str(work / f"{view}_%04d.png"),
                    "-vf", "pad=ceil(iw/2)*2:ceil(ih/2)*2", "-pix_fmt", "yuv420p",
                    "-crf", "23", str(mp4)], check=True)
    subprocess.run([ffmpeg, "-y", "-loglevel", "error", "-framerate", str(fps),
                    "-i", str(work / f"{view}_%04d.png"),
                    "-vf", "fps=12,scale=560:-1:flags=lanczos,split[a][b];"
                           "[a]palettegen[p];[b][p]paletteuse",
                    str(gif)], check=True)
    print(f"  {view}: {len(frames)} frames -> {mp4.name} "
          f"({mp4.stat().st_size // 1024} kB), {gif.name} "
          f"({gif.stat().st_size // 1024} kB)")
    return mp4.name


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--traj", required=True)
    ap.add_argument("--top", required=True)
    ap.add_argument("--name", required=True)
    ap.add_argument("--out", default="figures")
    ap.add_argument("--views", nargs="+", default=list(BODIES))
    ap.add_argument("--frames", type=int, default=60)
    ap.add_argument("--fps", type=int, default=12)
    ap.add_argument("--width", type=int, default=900)
    ap.add_argument("--height", type=int, default=700)
    ap.add_argument("--still-width", type=int, default=1500)
    ap.add_argument("--still-height", type=int, default=1150)
    ap.add_argument("--sites", nargs=2, type=int, default=[133, 149])
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
    t, rmsf = prepare(a.traj, a.top, work, a.frames, a.sites)
    print(f"  {t.n_frames} frames, {t.n_atoms} solute atoms")
    made = []
    for v in a.views:
        made.append(render(v, work, out, a.name, a.width, a.height, pymol, ffmpeg,
                           a.fps))
        render_still(v, work, out, a.name, a.still_width, a.still_height, pymol,
                     max(1, t.n_frames // 2))
    if not a.keep_work:
        shutil.rmtree(work)
    print("  wrote: " + ", ".join(m for m in made if m))


if __name__ == "__main__":
    main()
