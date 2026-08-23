#!/usr/bin/env python3
"""Build a type II kink into a duplex, so we can ask whether it survives.

Why
---
Doeke asked whether "it does not kink" reflects the force field, the timescale, or an
error in Zocchi's model of when kinking occurs.  The compression windows address the
second and third: they squeeze the duplex past the curvature at which all-atom MD is
reported to nucleate kinks, and if none appears there, spontaneous nucleation is not the
explanation.  But that still leaves the first objection open -- perhaps this force field
simply cannot hold a kink at all, in which case never seeing one says nothing.

This closes it off from the other side.  Put a kink in by hand, release it, and watch:

  * if it is metastable over tens of nanoseconds, the force field CAN support a kink, and
    the reason we never see one spontaneously is the barrier and the timescale -- answers
    (b) and (c);
  * if it heals in a nanosecond or two, that is a force-field result in its own right and
    worth reporting as one -- answer (a).

Either way the test is informative, which is what makes it worth the GPU time.  Threshold
stated in advance, before any run: the seed is judged to have SURVIVED if, averaged over
the last two thirds of the trajectory, the seeded step still carries a roll of more than
30 deg AND the central base pair's N1-N3 distance is still beyond 4 A.  Anything less is
healing.

What a type II kink is
----------------------
From Lankas, Lavery & Maddocks (Structure 2006) and Curuksu et al. (NAR 2009), the type II
kink is a three-base-pair motif:

  * the two flanking pairs keep their Watson-Crick hydrogen bonds;
  * the CENTRAL pair is broken, and its bases stack on their 5' neighbours instead of on
    each other;
  * the junction step carries a very large roll -- Mitchell & Harris report -76 to -92 deg
    -- giving a local bend of 80-95 deg, and the broken pair a propeller of order 120 deg.

So the seed here is two operations, both rigid-body and both built from the same Rodrigues
rotation the rest of the pipeline uses:

  1. ROLL.  Rotate everything downstream of the chosen step about the base-pair LONG axis
     through that step's centre.  Roll is rotation about the long axis by convention
     (Olson/Dickerson), and it is the mode that bends the helix into the grooves, which is
     what a kink does.  Tilt and twist do not produce a kink.
  2. PROPELLER.  Rotate one base of the central pair about its own glycosidic bond, which
     breaks the Watson-Crick pairing and swings it towards its neighbour to stack.

Neither operation tries to place the broken bases perfectly; minimisation and a restrained
equilibration do that.  The point of the seed is to start inside the kinked basin, not to
guess its exact bottom.
"""
from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from build_dna import PRIME3, read_pdb, renumber, residues, write_pdb  # noqa: E402

# The base atoms, i.e. everything on the far side of the glycosidic bond.  Rotating these
# and not the sugar or phosphate is what makes the propeller a propeller rather than a
# backbone distortion.
BASE_ATOMS = {
    "N1", "C2", "N3", "C4", "C5", "C6", "N6", "N7", "C8", "N9", "O2", "O4", "O6",
    "N2", "N4", "C7", "C5M", "H1", "H2", "H3", "H5", "H6", "H8", "H41", "H42",
    "H21", "H22", "H61", "H62", "H71", "H72", "H73", "H51", "H52", "H53",
}
PURINES = {"DA", "DG"}


def rodrigues(axis, angle):
    a = np.asarray(axis, float)
    a = a / np.linalg.norm(a)
    K = np.array([[0, -a[2], a[1]], [a[2], 0, -a[0]], [-a[1], a[0], 0]])
    return np.eye(3) + math.sin(angle) * K + (1 - math.cos(angle)) * (K @ K)


def bp_frame(res, a, b, i, n_bp):
    """(centre, long axis, local helical axis) for base pair i, from C1' and the rise."""
    ra, rb = res[(a, i)], res[(b, n_bp + 1 - i)]
    c1a, c1b = ra["C1'"]["xyz"], rb["C1'"]["xyz"]
    centre = 0.5 * (c1a + c1b)
    long_axis = c1b - c1a
    long_axis = long_axis / np.linalg.norm(long_axis)
    # the local helical direction, from the neighbouring pairs' centres
    lo = max(1, i - 1)
    hi = min(n_bp, i + 1)
    def ctr(k):
        return 0.5 * (res[(a, k)]["C1'"]["xyz"] + res[(b, n_bp + 1 - k)]["C1'"]["xyz"])
    z = ctr(hi) - ctr(lo)
    z = z / np.linalg.norm(z)
    # make the long axis perpendicular to the helical axis, so a "roll" is pure roll and
    # does not smuggle in a twist component
    long_axis = long_axis - (long_axis @ z) * z
    long_axis /= np.linalg.norm(long_axis)
    return centre, long_axis, z


def seed(atoms, step: int, roll_deg: float, propeller_deg: float):
    """Roll everything past `step`, then break the pair AT `step` with a propeller."""
    res = residues(atoms)
    chains = sorted({c for c, _ in res})
    if len(chains) != 2:
        sys.exit(f"need a duplex, found chains {chains}")
    a, b = chains
    n_bp = max(n for _, n in res)
    if not 2 <= step <= n_bp - 2:
        sys.exit(f"step {step} needs at least two base pairs on each side of it "
                 f"(duplex is {n_bp} bp)")

    centre, long_axis, z = bp_frame(res, a, b, step, n_bp)
    report = {"step": step, "n_bp": n_bp, "roll_deg": roll_deg,
              "propeller_deg": propeller_deg}

    # ---- 1. the roll, about the axis that runs THROUGH both backbone junctions ----
    # The roll axis of a base-pair STEP is not the line through the base-pair centre.  It
    # runs across the duplex between the two backbones, and the two O3'-P junctions being
    # rolled about lie on it -- which is the whole point: atoms on the rotation axis do not
    # move, so the backbone stays closed while the helix bends.
    #
    # Rolling about the base-pair centre instead tore both junctions open, to 7.3 A at the
    # middle of the duplex and 35 A near its end, and no rigid translation can repair that:
    # the two junctions are ~18 A apart across the helix, so a roll moves them in different
    # directions and a single shift can only fix their average.  Putting the axis through
    # them removes the problem rather than patching it.
    up_a, dn_a = res[(a, step)].get(PRIME3), res[(a, step + 1)].get("P")
    up_b, dn_b = res[(b, n_bp - step)].get(PRIME3), res[(b, n_bp + 1 - step)].get("P")
    if all(x is not None for x in (up_a, dn_a, up_b, dn_b)):
        j_a = 0.5 * (up_a["xyz"] + dn_a["xyz"])       # strand A's junction midpoint
        j_b = 0.5 * (up_b["xyz"] + dn_b["xyz"])       # strand B's
        pivot, axis = j_a, (j_b - j_a)
        report["roll_axis"] = "through both O3'-P junctions"
    else:
        pivot, axis = centre, long_axis
        report["roll_axis"] = "base-pair long axis (junction atoms not found)"
    R = rodrigues(axis, math.radians(roll_deg))
    moved = 0
    for i in range(step + 1, n_bp + 1):
        for key in ((a, i), (b, n_bp + 1 - i)):
            for at in res[key].values():
                at["xyz"] = pivot + R @ (at["xyz"] - pivot)
                moved += 1
    report["atoms_rolled"] = moved

    # ---- 2. the propeller, breaking the central pair ----
    # Applied to the pyrimidine of the pair, arbitrarily but consistently: it is the
    # smaller base, so swinging it disturbs less of the surrounding stack, and the
    # literature motif has the broken bases stacking on their 5' neighbours either way.
    ra, rb = res[(a, step)], res[(b, n_bp + 1 - step)]
    pick = rb if next(iter(rb.values()))["resn"].rstrip("35") not in PURINES else ra
    nm = next(iter(pick.values()))["resn"].rstrip("35")
    glyco = "N9" if nm in PURINES else "N1"
    if "C1'" not in pick or glyco not in pick:
        sys.exit(f"cannot find the glycosidic bond of the central pair ({nm})")
    p0, p1 = pick["C1'"]["xyz"], pick[glyco]["xyz"]
    Rp = rodrigues(p1 - p0, math.radians(propeller_deg))
    spun = 0
    for name, at in pick.items():
        if name in BASE_ATOMS:
            at["xyz"] = p1 + Rp @ (at["xyz"] - p1)
            spun += 1
    report["base_rotated"] = {"residue": nm, "atoms": spun, "about": f"C1'-{glyco}"}

    # ---- did the seed do what we asked? ----
    res = residues(atoms)
    def wc(i):
        r1, r2 = res[(a, i)], res[(b, n_bp + 1 - i)]
        n1 = "N1" if next(iter(r1.values()))["resn"].rstrip("35") in PURINES else "N3"
        n2 = "N1" if next(iter(r2.values()))["resn"].rstrip("35") in PURINES else "N3"
        if n1 not in r1 or n2 not in r2:
            return None
        return float(np.linalg.norm(r1[n1]["xyz"] - r2[n2]["xyz"]))
    report["N1_N3_A"] = {str(i): (round(wc(i), 2) if wc(i) else None)
                         for i in (step - 1, step, step + 1)}
    # overall bend: angle between the helical directions well before and well after
    def axis_at(k):
        _, _, zz = bp_frame(res, a, b, k, n_bp)
        return zz
    before, after = axis_at(max(2, step - 3)), axis_at(min(n_bp - 1, step + 3))
    report["local_bend_deg"] = round(
        float(np.degrees(np.arccos(np.clip(before @ after, -1, 1)))), 1)
    # the backbone had better still be a backbone
    dev = [float(np.linalg.norm(res[(ch, i)][PRIME3]["xyz"] - res[(ch, i + 1)]["P"]["xyz"]))
           for ch in (a, b) for i in range(1, n_bp)
           if PRIME3 in res[(ch, i)] and "P" in res[(ch, i + 1)]]
    report["O3prime_P"] = [round(float(np.mean(dev)), 3), round(float(np.max(dev)), 3)]
    return report


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--duplex", required=True, help="a built duplex PDB")
    ap.add_argument("--step", type=int, required=True,
                    help="base pair to break; the roll goes at the step after it")
    ap.add_argument("--roll", type=float, default=-85.0,
                    help="roll in degrees; Mitchell & Harris report -76 to -92 for a kink")
    ap.add_argument("--propeller", type=float, default=120.0)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    atoms = renumber(read_pdb(Path(a.duplex)))
    rep = seed(atoms, a.step, a.roll, a.propeller)
    write_pdb(atoms, Path(a.out))
    Path(a.out).with_suffix(".json").write_text(json.dumps(rep, indent=2))
    print(f"seeded a kink at base pair {rep['step']} of {rep['n_bp']}")
    print(f"  rolled {rep['atoms_rolled']} atoms by {a.roll:+.0f} deg about an axis "
          f"{rep.get('roll_axis','')}")
    print(f"  spun {rep['base_rotated']['atoms']} base atoms of "
          f"{rep['base_rotated']['residue']} by {a.propeller:+.0f} deg about "
          f"{rep['base_rotated']['about']}")
    print(f"  N1-N3 at pairs {list(rep['N1_N3_A'])}: "
          f"{list(rep['N1_N3_A'].values())} A  (paired is ~2.8, broken is >4)")
    print(f"  local bend across the motif: {rep['local_bend_deg']} deg")
    print(f"  backbone O3'-P: mean {rep['O3prime_P'][0]}, max {rep['O3prime_P'][1]} A "
          f"(ideal 1.61; minimisation takes the strain out)")
    print(f"  wrote {a.out}")


if __name__ == "__main__":
    main()
