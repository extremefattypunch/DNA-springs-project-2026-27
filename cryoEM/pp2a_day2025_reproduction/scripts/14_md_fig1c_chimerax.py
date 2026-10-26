#!/usr/bin/env python
"""Fig 1C analogue: the N-terminus at 0, 20, 40, 60, 80, 100 ns, coloured red -> blue, over the core from frame 0.

Frames are extracted from protein.xtc (10_md_rmsd_rmsf.py) with MDAnalysis after superposing every frame on the
core Calpha of frame 0, written as six PDBs, and rendered with ChimeraX offscreen: core cartoon in the paper's
colour (Aalpha red, Abeta grey) from frame 0, N-terminal residues (1-13 / 1-25) of each time point as a cartoon
tube coloured by time.  Usage: $DNASPRING_PY scripts/14_md_fig1c_chimerax.py [Aa Ab]
"""
import os
import shlex
import subprocess
import sys

import MDAnalysis as mda
from MDAnalysis.analysis import align

REPO = os.environ.get("PP2A_REPO", "/n/home05/ianpoon/DNA-springs-project-2026-27/cryoEM/pp2a_day2025_reproduction")
DATA = os.environ.get("PP2A_DATA", "/n/netscratch/hekstra_lab/Lab/ian_poon/pp2a_day2025")
CHIMERAX = os.environ.get("CHIMERAX")
FIG = os.path.join(REPO, "figures")
SYS = {"Aa": dict(nt=13, core=(14, 589), color="#c0392b"), "Ab": dict(nt=25, core=(26, 601), color="#7f7f7f")}
TIMES_NS = [0, 20, 40, 60, 80, 100]
FRAME_PS = 100.0
# red (0 ns) -> blue (100 ns), as the paper's colour bar
RAMP = ["#ff0000", "#cc0033", "#990066", "#660099", "#3300cc", "#0000ff"]


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
        out = os.path.join(run, "fig1c")
        os.makedirs(out, exist_ok=True)
        pdbs = []
        for t in TIMES_NS:
            i = min(int(round(t * 1000 / FRAME_PS)), len(u.trajectory) - 1)
            u.trajectory[i]
            p = os.path.join(out, f"frame_{t:03d}ns.pdb")
            u.select_atoms("protein").write(p)
            pdbs.append(p)
        cx = ""
        for k, p in enumerate(pdbs):
            cx += f"open {p} id {k+1}\n"
        cx += "hide atoms\ncartoon\n"
        cx += f"color #1 {m['color']}\nhide #2-{len(pdbs)}/A:{m['nt']+1}-1000 cartoon\nhide #2-{len(pdbs)} atoms\n"
        for k, col in enumerate(RAMP[:len(pdbs)]):
            cx += f"color #{k+1}:1-{m['nt']} {col}\ncartoon style #{k+1}:1-{m['nt']} width 1.2 thick 1.2\n"
        cx += "set bgColor white\nlighting soft\ngraphics silhouettes true\nview orient\nview\n"
        cx += f"2dlabels text 'N-terminus 1-{m['nt']}: 0 ns (red) -> 100 ns (blue), 20 ns steps' xpos 0.03 ypos 0.94 size 18 color black\n"
        cx += f"save {FIG}/fig1c_nterm_overlay_{name.lower()}.png width 1600 height 1200 supersample 2\n"
        cx += f"view #1:1-{m['nt']+30}\nzoom 0.8\nsave {FIG}/fig1c_nterm_overlay_{name.lower()}_zoom.png width 1600 height 1200 supersample 2\nexit\n"
        cxc = os.path.join(out, "fig1c.cxc")
        open(cxc, "w").write(cx)
        if CHIMERAX:
            r = subprocess.run(shlex.split(CHIMERAX) + ["--nogui", "--offscreen", cxc], capture_output=True, text=True, timeout=1800)
            open(os.path.join(out, "fig1c.log"), "w").write(r.stdout + r.stderr)
            print(f"{name}: rendered (rc {r.returncode})" + ("" if r.returncode == 0 else r.stdout[-800:]))
        else:
            print(f"{name}: wrote {cxc} (CHIMERAX not set)")


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:] or list(SYS)))
