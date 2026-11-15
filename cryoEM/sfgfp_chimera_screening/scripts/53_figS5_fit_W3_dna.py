# ChimeraX script (Python): rigid fit of a straight 32 bp B-DNA (model/L6_loop50/dna/straight_50bp.pdb trimmed to
# chain A 10-41 / chain B -41..-10 = 108 A, the ~110 A W3 rod body of report.md section 10) into the W3 DNA-control
# map, for the Fig S5C map-model FSC. Run offscreen:  $CHIMERAX --nogui --offscreen --silent scripts/53_figS5_fit_W3_dna.py
#
#  1. global search (fitmap search 300, simulated map at 8.7 A = RELION loose-mask resolution) into the RELION
#     half-map average raw.mrc, contour enclosing 65,000 A^3 (project primary contour) -> top fits logged
#  2. defensibility tests on the best fit:
#       spin   : rotate the duplex about its own helix axis in 15 deg steps; correlation without refitting, and
#                after local rigid refinement from each start (does the register return to one basin?)
#       slide  : translate along the helix axis in 3.4 A steps (one bp) without refitting
#       flip   : reverse the axis direction (polarity) and refine locally
#  3. save the best fit in the map frame (relModel #1) for phenix.mtriage
import json
import numpy as np
from chimerax.core.commands import run

S = "/n/netscratch/hekstra_lab/Lab/ian_poon/sfgfp_chimera_screening_20261112/figures_scratch/W3_duplex/fit"
MAP = "/n/netscratch/hekstra_lab/Lab/ian_poon/sfgfp_chimera_screening_20261112/relion/W3_duplex/raw.mrc"
RES = 8.7
out = {"map": MAP, "model_in": f"{S}/dna_32bp.pdb", "resolution_sim_A": RES}

vol = run(session, f"open {MAP}")[0]
run(session, "volume #1 encloseVolume 65000")
out["map_level"] = float(vol.minimum_surface_level)
dna = run(session, f"open {S}/dna_32bp.pdb")[0]

xyz = dna.atoms.coords
P = dna.atoms.filter(dna.atoms.names == "P").coords
c0 = xyz.mean(0)
# helix axis in model coordinates: principal axis of the phosphates
w, V = np.linalg.eigh(np.cov((P - P.mean(0)).T))
ax0 = V[:, 2]


def corr(move=False):
    opt = "" if move else " shift false rotate false"
    fl = run(session, f"fitmap #2 inMap #1 resolution {RES} metric correlation envelope false{opt}")
    return float(fl[0].correlation()) if fl else float("nan")


def state():
    p = dna.position
    return {"centre": p * c0, "axis": p.transform_vector(ax0)}


fl = run(session, f"fitmap #2 inMap #1 resolution {RES} metric correlation envelope false search 300 placement sr "
                  f"clusterAngle 6 clusterShift 3 logFits {S}/fitmap_search.txt")
tops = []
for f in fl[:15]:
    f.place_models(session)
    st = state()
    tops.append({"correlation": float(f.correlation()), "hits": int(f.hits()),
                 "centre_A": [round(float(v), 2) for v in st["centre"]], "axis": [round(float(v), 4) for v in st["axis"]]})
out["search_top"] = tops
fl[0].place_models(session)
best_pos = dna.position
best = corr()
out["best_correlation"] = best
st = state()
out["best_centre_A"] = [float(v) for v in st["centre"]]
out["best_axis"] = [float(v) for v in st["axis"]]
run(session, f"save {S}/dna_32bp_fit.pdb models #2 relModel #1")

from chimerax.geometry import rotation, translation
# spin about own helix axis
spin = []
for ang in range(0, 360, 15):
    dna.position = best_pos
    c, a = state()["centre"], state()["axis"]
    dna.position = rotation(a, ang, center=c) * best_pos
    c_fixed = corr(False)
    c_ref = corr(True)
    a2 = state()["axis"]
    # register after refinement: angle between refined and best frames about the axis
    R = (dna.position * best_pos.inverse()).rotation_angle()
    spin.append({"spin_deg": ang, "corr_no_refit": c_fixed, "corr_after_local_refit": c_ref,
                 "residual_rotation_deg": float(R)})
out["spin_test"] = spin
# slide along the axis
slide = []
for k in range(-6, 7):
    dna.position = best_pos
    a = state()["axis"]
    dna.position = translation(a * 3.375 * k) * best_pos
    slide.append({"shift_A": round(3.375 * k, 3), "corr_no_refit": corr(False)})
out["slide_test"] = slide
# polarity flip (rotate 180 about an axis perpendicular to the helix axis, through the centre)
dna.position = best_pos
c, a = state()["centre"], state()["axis"]
perp = np.cross(a, [0, 0, 1.0])
if np.linalg.norm(perp) < 1e-3:
    perp = np.cross(a, [0, 1.0, 0])
perp /= np.linalg.norm(perp)
dna.position = rotation(perp, 180, center=c) * best_pos
out["flip_corr_after_local_refit"] = corr(True)
dna.position = best_pos
out["best_correlation_check"] = corr(False)
json.dump(out, open(f"{S}/fit_W3_dna.json", "w"), indent=1)
run(session, "exit")
