#!/usr/bin/env python
"""Videos S1/S2 analogues: smoothed 100 ns trajectory movies with the N-terminal residues in blue.

Paper: "Smoothed trajectory movies were recorded using VMD ... smoothing window of 10".  Here: frames superposed on
the core Calpha of frame 0, coordinates smoothed with a 10-frame running mean (MDAnalysis), written as a
multi-model PDB, and recorded with ChimeraX offscreen (`coordset` playback -> PNG frames -> ffmpeg H.264, 30 fps).
1,000 frames (100 ps each) -> 33 s at 30 fps.  Usage: $DNASPRING_PY scripts/15_md_movies.py [Aa Ab]
"""
import os
import shlex
import subprocess
import sys

import numpy as np
import MDAnalysis as mda
from MDAnalysis.analysis import align
from scipy.ndimage import uniform_filter1d

REPO = os.environ.get("PP2A_REPO", "/n/home05/ianpoon/DNA-springs-project-2026-27/cryoEM/pp2a_day2025_reproduction")
DATA = os.environ.get("PP2A_DATA", "/n/netscratch/hekstra_lab/Lab/ian_poon/pp2a_day2025")
CHIMERAX = os.environ.get("CHIMERAX")
FFMPEG = os.path.join(os.path.dirname(os.environ.get("DNASPRING_PY", sys.executable)), "ffmpeg")
FIG = os.path.join(REPO, "figures")
SYS = {"Aa": dict(nt=13, core=(14, 589), color="#c0392b", video="videoS1_aa"), "Ab": dict(nt=25, core=(26, 601), color="#7f7f7f", video="videoS2_ab")}
SMOOTH = 10
FPS = 30
W, H = 1280, 960


def main(systems):
    for name in systems:
        m = SYS[name]
        run = os.path.join(DATA, "md", "runs", f"{name}_rep1")
        xtc, ref = os.path.join(run, "protein.xtc"), os.path.join(run, "protein_ref.pdb")
        if not os.path.exists(xtc):
            print(f"{name}: protein.xtc missing (run 10 first)")
            continue
        u = mda.Universe(ref, xtc)
        c0, c1 = m["core"]
        align.AlignTraj(u, mda.Universe(ref), select=f"name CA and resid {c0}:{c1}", in_memory=True).run()
        xyz = np.array([u.atoms.positions.copy() for _ in u.trajectory])          # (nframes, natoms, 3)
        xyz = uniform_filter1d(xyz, size=SMOOTH, axis=0, mode="nearest")           # VMD 'smoothing window 10'
        out = os.path.join(run, "movie")
        os.makedirs(out, exist_ok=True)
        smoothed = os.path.join(out, "smoothed.dcd")
        with mda.Writer(smoothed, u.atoms.n_atoms) as w:
            for i in range(len(xyz)):
                u.trajectory[i]
                u.atoms.positions = xyz[i]
                w.write(u.atoms)
        frames_dir = os.path.join(out, "frames")
        os.makedirs(frames_dir, exist_ok=True)
        n = len(xyz)
        cx = f"open {ref} id 1\nopen {smoothed} structureModel #1\n"
        cx += f"hide atoms\ncartoon\ncolor #1 {m['color']}\ncolor #1:1-{m['nt']} blue\ncartoon style #1:1-{m['nt']} width 1.4 thick 1.4\n"
        cx += "set bgColor white\nlighting soft\ngraphics silhouettes true\nview orient\nview\n"
        cx += f"2dlabels create title text 'PP2A {'Aα' if name == 'Aa' else 'Aβ'} — 100 ns, N-terminus 1-{m['nt']} blue (10-frame smoothing)' xpos 0.03 ypos 0.94 size 16 color black\n"
        cx += "2dlabels create tlab text 't = 0 ns' xpos 0.03 ypos 0.05 size 18 color black\n"
        for i in range(n):
            cx += f"coordset #1 {i+1}\n2dlabels change tlab text 't = {i*0.1:5.1f} ns'\nsave {frames_dir}/f{i:04d}.png width {W} height {H}\n"
        cx += "exit\n"
        cxc = os.path.join(out, "movie.cxc")
        open(cxc, "w").write(cx)
        if not CHIMERAX:
            print(f"{name}: wrote {cxc} (CHIMERAX not set)")
            continue
        r = subprocess.run(shlex.split(CHIMERAX) + ["--nogui", "--offscreen", cxc], capture_output=True, text=True, timeout=7200)
        open(os.path.join(out, "movie.log"), "w").write(r.stdout + r.stderr)
        if r.returncode != 0:
            print(f"{name}: ChimeraX failed\n{r.stdout[-800:]}")
            continue
        mp4 = os.path.join(FIG, m["video"] + ".mp4")
        subprocess.run([FFMPEG, "-y", "-loglevel", "error", "-framerate", str(FPS), "-i", os.path.join(frames_dir, "f%04d.png"),
                        "-c:v", "libx264", "-pix_fmt", "yuv420p", "-crf", "20", mp4], check=True)
        print(f"{name}: {mp4} ({n} frames, {n / FPS:.0f} s)")


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:] or list(SYS)))
