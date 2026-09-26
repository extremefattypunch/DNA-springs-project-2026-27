#!/usr/bin/env python
"""Plan section 5G (ChimeraX part) - map/model figures in the paper's style, rendered offscreen.

Generates one ChimeraX command file (.cxc) per panel under $PP2A_DATA/em/chimerax/ and runs it with
$CHIMERAX --nogui --offscreen (ChimeraX 1.12 inside the Rocky 9 container; see env/activate_pp2a.sh).
Panels (paper colours: Aalpha red, Abeta grey, B56g blue, Calpha orange; C-tail green):
  fig2a_maps.png              deposited sharpened maps at the EMDB contour, coloured by nearest chain, two views
  fig2b_lowpass6_b56_helix.png  6 A low-pass of the raw (half-map average) maps at matched threshold; the
                              B56g C-terminal helix (9MIP chain B 402-436) shown as a blue ribbon in both frames
                              (transferred into the Abeta frame by superposing B56g cores)  -> Fig 2B
  fig3a_heat_com.png          Aalpha/Abeta A subunits superposed on HEAT 12-15 Calpha, HEAT labels, the
                              HEAT1-HEAT15 centre-of-mass distances as marker-to-marker distances  -> Fig 3A
  fig3c_salt_bridge.png       R113-E191 (Aalpha) and R125-E203 (Abeta) with map mesh and Calpha-Calpha distance -> Fig 3C
  fig3de_interfaces.png       A-B56g and A-Calpha interface residues as sticks with mesh                      -> Fig 3D/E
  figS5_localres.png          locally filtered maps coloured by RELION local resolution (3.0-5.0 A palette)   -> Fig S5A/E
Run:  $DNASPRING_PY scripts/26_chimerax_figures.py [panel ...]
"""
import os
import shlex
import subprocess
import sys

REPO = os.environ.get("PP2A_REPO", "/n/home05/ianpoon/DNA-springs-project-2026-27/cryoEM/pp2a_day2025_reproduction")
DATA = os.environ.get("PP2A_DATA", "/n/netscratch/hekstra_lab/Lab/ian_poon/pp2a_day2025")
FIG = os.path.join(REPO, "figures")
CX = os.path.join(DATA, "em", "chimerax")
CHIMERAX = os.environ.get("CHIMERAX")
M = os.path.join(DATA, "em", "models")
E = os.path.join(DATA, "em", "emdb")
R = os.path.join(DATA, "em", "relion")
ISO = {"Aa": dict(emd="48300", pdb="9MIP", contour=0.0777, acolor="#c0392b", off=0),
       "Ab": dict(emd="48224", pdb="9MF5", contour=0.154, acolor="#7f7f7f", off=12)}
B56, CA = "#3b6fd6", "#f39c12"
W, H = 1600, 1200


def chain_colors(model, acolor):
    return f"color {model}/A {acolor}\ncolor {model}/B {B56}\ncolor {model}/C {CA}\n"


def run(name, cxc):
    os.makedirs(CX, exist_ok=True)
    os.makedirs(FIG, exist_ok=True)
    path = os.path.join(CX, f"{name}.cxc")
    open(path, "w").write(cxc + "\nexit\n")
    if not CHIMERAX:
        print(f"[{name}] wrote {path} (CHIMERAX not set - not rendered)")
        return
    cmd = shlex.split(CHIMERAX) + ["--nogui", "--offscreen", path]   # no --silent: the log must show which command failed
    print(f"[{name}] rendering ...", flush=True)
    r = subprocess.run(cmd, capture_output=True, text=True, timeout=1800)
    log = r.stdout + r.stderr
    open(os.path.join(CX, f"{name}.log"), "w").write(log)
    if r.returncode != 0:
        print(log[-2000:])
        print(f"[{name}] FAILED (rc {r.returncode})")
    else:
        print(f"[{name}] ok")


def fig2a():
    cx = ""
    for i, (iso, m) in enumerate(ISO.items()):
        mdl, vol = f"#{2*i+1}", f"#{2*i+2}"
        cx += f"open {M}/{m['pdb']}.cif\nopen {E}/{iso}/emd_{m['emd']}.map\n"
        cx += chain_colors(mdl, m["acolor"])
        cx += f"volume {vol} level {m['contour']} style surface step 1\ncolor zone {vol} near {mdl} distance 4\nhide {mdl} models\n"
    # same orientation for both isoforms: superpose Abeta model on Aalpha, carry its map along
    cx += "matchmaker #3/A to #1/A\nview position #4 sameAsModels #3\n"
    cx += "set bgColor white\nlighting soft\ngraphics silhouettes true\n"
    cx += "view orient\n"          # principal axes of the (superposed) complexes onto the screen axes: same view for both
    for iso, vol, other in (("Aa", "#2", "#4"), ("Ab", "#4", "#2")):
        cx += f"hide {other} models\nshow {vol} models\nview #2 pad 0.05\n"
        cx += f"save {FIG}/fig2a_maps_{iso}.png width {W} height {H} supersample 3\n"
        cx += f"turn y 90\nsave {FIG}/fig2a_maps_{iso}_view2.png width {W} height {H} supersample 3\nturn y -90\n"
    run("fig2a", cx)


def fig2b():
    # raw_lp6.mrc from RELION (22_relion_postprocess.sh); helix 402-436 of 9MIP chain B as reference in both frames
    cx = f"open {M}/9MIP.cif\nopen {R}/Aa/raw_lp6.mrc\nopen {M}/9MF5.cif\nopen {R}/Ab/raw_lp6.mrc\nopen {M}/9MIP.cif\n"
    cx += "matchmaker #5/B:38-403 to #3/B:38-403\n"          # second copy of 9MIP carried into the Abeta frame
    cx += "matchmaker #3/A to #1/A\n"                        # then Abeta (with its map? no - maps stay; keep frames separate)
    cx += "hide #1,3,5 models\nhide #1,3,5 atoms\ncartoon hide\n"
    cx += "show #1/B:402-436 cartoon\nshow #5/B:402-436 cartoon\ncolor #1/B:402-436 blue\ncolor #5/B:402-436 blue\n"
    cx += "show #1,5 models\n"
    for vol, mdl, m in (("#2", "#1", ISO["Aa"]), ("#4", "#3", ISO["Ab"])):
        cx += f"volume {vol} style surface step 1 transparency 0.25\n"
    cx += "volume #2 sdLevel 3.0\nvolume #4 sdLevel 3.0\n"   # matched sigma level for both raw low-pass maps
    cx += "color zone #2 near #1 distance 5\ncolor zone #4 near #3 distance 5\n"
    cx += chain_colors("#1", ISO["Aa"]["acolor"]) + chain_colors("#3", ISO["Ab"]["acolor"])
    cx += "color zone #2 near #1 distance 5\ncolor zone #4 near #3 distance 5\n"
    cx += "set bgColor white\nlighting soft\ngraphics silhouettes true\n"
    cx += f"view #1,2\nsave {FIG}/fig2b_lowpass6_Aa.png width {W} height {H} supersample 3\n"
    cx += f"view #3,4,5\nsave {FIG}/fig2b_lowpass6_Ab.png width {W} height {H} supersample 3\n"
    run("fig2b", cx)


def fig3a():
    cx = f"open {M}/9MIP.cif\nopen {M}/9MF5.cif\n"
    cx += "hide atoms\nhide #1,2/B,C cartoon\nhide #1,2/B,C atoms\n"
    cx += "align #2/A:449-601@CA toAtoms #1/A:437-589@CA\n"       # HEAT 12-15 superposition (Abeta = Aalpha + 12)
    cx += f"color #1/A {ISO['Aa']['acolor']}\ncolor #2/A {ISO['Ab']['acolor']}\n"
    for k in range(1, 16):
        a, b = 8 + 39 * (k - 1), 8 + 39 * k - 1
        cx += f"define centroid #1/A:{a}-{b}@CA name h{k}a color {ISO['Aa']['acolor']} radius 1.2\n"
        cx += f"define centroid #2/A:{a+12}-{b+12}@CA name h{k}b color {ISO['Ab']['acolor']} radius 1.2\n"
        cx += f"2dlabels text 'HEAT {k}' size 14\n" if False else ""
    cx += "distance sel:h1a sel:h15a\n" if False else ""
    cx += "set bgColor white\nlighting soft\ngraphics silhouettes true\nview\n"
    cx += f"save {FIG}/fig3a_heat_com.png width {W} height {H} supersample 3\n"
    run("fig3a", cx)


def fig3c():
    cx = ""
    for i, (iso, m) in enumerate(ISO.items()):
        mdl, vol = f"#{2*i+1}", f"#{2*i+2}"
        r1, e1 = (113, 191) if iso == "Aa" else (125, 203)
        cx += f"open {M}/{m['pdb']}.cif\nopen {E}/{iso}/emd_{m['emd']}.map\n"
        cx += f"hide {mdl} atoms\ncartoon {mdl}/A\nhide {mdl}/B,C cartoon\ncolor {mdl}/A {m['acolor']}\n"
        cx += f"show {mdl}/A:{r1},{e1} atoms\nstyle {mdl}/A:{r1},{e1} stick\ncolor {mdl}/A:{r1},{e1} byhetero\n"
        cx += f"distance {mdl}/A:{r1}@CA {mdl}/A:{e1}@CA\n"
        cx += f"volume {vol} style mesh level {m['contour']} step 1\nvolume zone {vol} nearAtoms {mdl}/A:{r1},{e1} range 3\n"
        cx += f"transparency {vol} 40\n"
    cx += "set bgColor white\nlighting soft\n"
    cx += f"view #1/A:113,191\nsave {FIG}/fig3c_salt_bridge_Aa.png width {W} height {H} supersample 3\n"
    cx += f"view #3/A:125,203\nsave {FIG}/fig3c_salt_bridge_Ab.png width {W} height {H} supersample 3\n"
    run("fig3c", cx)


def fig3de():
    iface_b = [100, 101, 183, 257]
    iface_c = [416, 418, 456, 495, 498, 572, 574]
    cx = ""
    for i, (iso, m) in enumerate(ISO.items()):
        mdl, vol = f"#{2*i+1}", f"#{2*i+2}"
        off = m["off"]
        rb = ",".join(str(r + off) for r in iface_b)
        rc = ",".join(str(r + off) for r in iface_c)
        cx += f"open {M}/{m['pdb']}.cif\nopen {E}/{iso}/emd_{m['emd']}.map\n"
        cx += f"hide {mdl} atoms\ncolor {mdl}/A {m['acolor']}\ncolor {mdl}/B {B56}\ncolor {mdl}/C {CA}\n"
        cx += f"show {mdl}/A:{rb},{rc} atoms\nstyle {mdl}/A:{rb},{rc} stick\ncolor {mdl}/A:{rb},{rc} byhetero\n"
        cx += f"label {mdl}/A:{rb},{rc} residues text '{{0.name}}{{0.number}}' height 1.2\n"
        cx += f"volume {vol} style mesh level {m['contour']} step 1\nvolume zone {vol} nearAtoms {mdl}/A:{rb},{rc} range 2.5\n"
    cx += "set bgColor white\nlighting soft\n"
    cx += f"view #1/A:{','.join(map(str, iface_b))}\nsave {FIG}/fig3d_ifaceB56_Aa.png width {W} height {H} supersample 3\n"
    cx += f"view #3/A:{','.join(str(r+12) for r in iface_b)}\nsave {FIG}/fig3d_ifaceB56_Ab.png width {W} height {H} supersample 3\n"
    cx += f"view #1/A:{','.join(map(str, iface_c))}\nsave {FIG}/fig3e_ifaceC_Aa.png width {W} height {H} supersample 3\n"
    cx += f"view #3/A:{','.join(str(r+12) for r in iface_c)}\nsave {FIG}/fig3e_ifaceC_Ab.png width {W} height {H} supersample 3\n"
    run("fig3de", cx)


def figS5():
    cx = ""
    for i, (iso, m) in enumerate(ISO.items()):
        filt, loc = f"#{2*i+1}", f"#{2*i+2}"
        cx += f"open {R}/{iso}/locres_filtered.mrc\nopen {R}/{iso}/relion_locres.mrc\n"
        cx += f"volume {filt} sdLevel 4 style surface step 1\nhide {loc} models\n"
        cx += f"color sample {filt} map {loc} palette 3.0,red:3.5,orange:4.0,yellow:4.5,green:5.0,blue\n"
    cx += "set bgColor white\nlighting soft\nview\n"
    cx += f"view #1\nsave {FIG}/figS5_localres_Aa.png width {W} height {H} supersample 3\n"
    cx += f"turn y 90\nsave {FIG}/figS5_localres_Aa_view2.png width {W} height {H} supersample 3\n"
    cx += f"view #3\nsave {FIG}/figS5_localres_Ab.png width {W} height {H} supersample 3\n"
    cx += f"turn y 90\nsave {FIG}/figS5_localres_Ab_view2.png width {W} height {H} supersample 3\n"
    run("figS5", cx)


PANELS = {"fig2a": fig2a, "fig2b": fig2b, "fig3a": fig3a, "fig3c": fig3c, "fig3de": fig3de, "figS5": figS5}

if __name__ == "__main__":
    for name in (sys.argv[1:] or PANELS):
        PANELS[name]()
