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
set sphere_scale, 0.35, resn TDP and name CB
orient
zoom all, 4
turn x, -15
""",
    "chromophore": """
show cartoon, polymer.protein
set cartoon_transparency, 0.75
color grey90, polymer.protein
show sticks, resn CRO
color limegreen, resn CRO
show sticks, byres (polymer.protein within 4.5 of resn CRO) and sidechain
color salmon, byres (polymer.protein within 4.5 of resn CRO) and sidechain
show sticks, resi 148+203+205+96+94+222 and sidechain
distance hb1, resn CRO and name OH, resi 148 and name ND1
distance hb2, resn CRO and name OH, resi 203 and name OG1
distance hb3, resn CRO and name O2, resi 96 and name NH2
color black, hb*
set dash_gap, 0.35
set dash_width, 2.0
set label_size, 0
orient resn CRO
zoom resn CRO, 6
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
orient polymer.protein
zoom all, 4
""",
}


def render(view, work: Path, out: Path, name, w, h, pymol, ffmpeg, fps):
    prefix = work / f"{view}"
    pml = work / f"{view}.pml"
    pml.write_text(PML.format(pdb="solute.pdb", dcd="solute.dcd",
                              body=BODIES[view], prefix=prefix.name, w=w, h=h))
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
    made = [render(v, work, out, a.name, a.width, a.height, pymol, ffmpeg, a.fps)
            for v in a.views]
    if not a.keep_work:
        shutil.rmtree(work)
    print("  wrote: " + ", ".join(m for m in made if m))


if __name__ == "__main__":
    main()
