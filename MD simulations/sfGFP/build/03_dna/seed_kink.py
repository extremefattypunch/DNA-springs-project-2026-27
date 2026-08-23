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
# Any two atoms in different residues closer than this are overlapping, not contacting:
# a real stacking contact is 3.3-3.6 A and a hydrogen bond 1.8-2.0, so 1.5 flags only
# genuine interpenetration.  A handful can be relaxed away; sixty cannot.
CLASH_A = 1.5
MAX_CLASHES = 4
WC_BROKEN_A = 4.0     # the central pair counts as broken beyond this
WC_PAIRED_A = 3.2     # the flanking pairs must still be within this


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


def count_clashes(res):
    """Inter-residue atom pairs closer than CLASH_A."""
    from scipy.spatial import cKDTree
    items = [(key, nm, at["xyz"]) for key, r in res.items() for nm, at in r.items()]
    xyz = np.array([x for _, _, x in items])
    who = [k for k, _, _ in items]
    tree = cKDTree(xyz)
    return sum(1 for i, j in tree.query_pairs(CLASH_A, output_type="ndarray")
               if who[i] != who[j])


def wc_distance(res, a, b, i, n_bp):
    r1, r2 = res.get((a, i)), res.get((b, n_bp + 1 - i))
    if r1 is None or r2 is None:
        return None
    n1 = "N1" if next(iter(r1.values()))["resn"].rstrip("35") in PURINES else "N3"
    n2 = "N1" if next(iter(r2.values()))["resn"].rstrip("35") in PURINES else "N3"
    if n1 not in r1 or n2 not in r2:
        return None
    return float(np.linalg.norm(r1[n1]["xyz"] - r2[n2]["xyz"]))


def motif(res, a, b, step, n_bp):
    """Is this still a type II kink?  central pair broken, flanking pairs intact."""
    c = wc_distance(res, a, b, step, n_bp)
    f1 = wc_distance(res, a, b, step - 1, n_bp)
    f2 = wc_distance(res, a, b, step + 1, n_bp)
    flank = max(x for x in (f1, f2) if x is not None) if (f1 or f2) else None
    return {"central": c, "flank": flank}


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
    # The bend is split across the TWO steps flanking the central pair, half at each,
    # rather than dumped on one.  That is the motif as described -- three base pairs with
    # the broken one acting as the hinge and sitting in the wedge -- and it is also what
    # makes the seed buildable: putting the whole 85 deg on a single step closes that
    # groove completely and drove the two flanking bases on one strand to 0.4 A of each
    # other.  Two half-rolls of 42.5 deg leave room for the broken pair between them.
    def junction_axis(k):
        """Pivot and direction for the step between base pairs k and k+1."""
        ua, da = res[(a, k)].get(PRIME3), res[(a, k + 1)].get("P")
        ub, db = res[(b, n_bp - k)].get(PRIME3), res[(b, n_bp + 1 - k)].get("P")
        if not all(x is not None for x in (ua, da, ub, db)):
            return None, None
        ja = 0.5 * (ua["xyz"] + da["xyz"])
        jb = 0.5 * (ub["xyz"] + db["xyz"])
        return ja, (jb - ja)

    moved = 0
    for k, first_moving in ((step - 1, step), (step, step + 1)):
        pivot, axis = junction_axis(k)
        if pivot is None:
            pivot, axis = centre, long_axis
        R = rodrigues(axis, math.radians(roll_deg / 2.0))
        for i in range(first_moving, n_bp + 1):
            for key in ((a, i), (b, n_bp + 1 - i)):
                for at in res[key].values():
                    at["xyz"] = pivot + R @ (at["xyz"] - pivot)
                    moved += 1
    report["roll_axis"] = ("half at each of the two steps flanking the broken pair, "
                           "about the axis through both O3'-P junctions")
    report["atoms_rolled"] = moved

    # ---- 1c. rise: open the wedge the roll closes ----
    # A large roll at a step necessarily drives the base edges together on the closing
    # side -- at 3.4 A rise and ~4 A of base reaching out from the axis, 42 deg closes
    # more than the whole gap.  Real kinks accommodate this by sliding and rising as well
    # as rolling ("kink-and-slide"), and no amount of rotating the bases about their own
    # glycosidic bonds substitutes for it: the colliding atoms are anchored to backbones
    # whose separation the rotation does not change.
    #
    # So translate the two halves apart along the local helical axis until the clashes
    # clear.  This DOES stretch the two O3'-P bonds at the junction, and that is the right
    # trade: a stretched covalent bond is something minimisation pulls straight back in,
    # whereas two rings threaded through each other is something it relaxes INTO and
    # cannot escape -- which is exactly what produced a guanine ring at 2.66 A and NaN
    # coordinates six minutes into production.
    rise_used = 0.0
    for rise in np.arange(0.0, 12.01, 0.25):
        if rise > 0:
            for i in range(step, n_bp + 1):
                for key in ((a, i), (b, n_bp + 1 - i)):
                    for at in res[key].values():
                        at["xyz"] = at["xyz"] + 0.125 * z
            for i in range(step + 1, n_bp + 1):
                for key in ((a, i), (b, n_bp + 1 - i)):
                    for at in res[key].values():
                        at["xyz"] = at["xyz"] + 0.125 * z
        rise_used = rise
        if count_clashes(res) <= MAX_CLASHES:
            break
    report["rise_A"] = round(float(rise_used), 2)

    # ---- 2. break the central pair: BOTH bases swing out of the stack ----
    # The motif is not a wedge inserted into an intact helix.  The central pair leaves the
    # stack -- "its bases stack on their 5' neighbours" -- and the two flanking pairs
    # close up on each other at an angle.  So both bases rotate about their own glycosidic
    # bonds, not just one; leaving the purine in place is what kept colliding with the
    # flanking base on that strand.
    #
    # The two angles are chosen by a small search rather than asserted, because which way
    # each base should swing depends on the local sequence and on the sign of the roll,
    # and there is no reason to guess.  The search is constrained so it cannot return
    # something that is not a kink: the central pair must be broken, both flanking pairs
    # must still be paired, and the bend must survive.  Among the candidates that qualify
    # it takes the one with fewest clashes.
    def rotate_base(residue, angle):
        nm_ = next(iter(residue.values()))["resn"].rstrip("35")
        gly = "N9" if nm_ in PURINES else "N1"
        if "C1'" not in residue or gly not in residue:
            return 0
        q0, q1 = residue["C1'"]["xyz"], residue[gly]["xyz"]
        Rq = rodrigues(q1 - q0, math.radians(angle))
        n = 0
        for name, at in residue.items():
            if name in BASE_ATOMS:
                at["xyz"] = q1 + Rq @ (at["xyz"] - q1)
                n += 1
        return n

    ra, rb = res[(a, step)], res[(b, n_bp + 1 - step)]
    saved = ({k: v["xyz"].copy() for k, v in ra.items()},
             {k: v["xyz"].copy() for k, v in rb.items()})
    best = None
    for ang_a in (0.0, propeller_deg / 2, -propeller_deg / 2, propeller_deg,
                  -propeller_deg):
        for ang_b in (0.0, propeller_deg / 2, -propeller_deg / 2, propeller_deg,
                      -propeller_deg):
            for k, v in ra.items():
                v["xyz"] = saved[0][k].copy()
            for k, v in rb.items():
                v["xyz"] = saved[1][k].copy()
            na = rotate_base(ra, ang_a)
            nb = rotate_base(rb, ang_b)
            m = motif(res, a, b, step, n_bp)
            if m["central"] is None or m["central"] <= WC_BROKEN_A:
                continue
            if m["flank"] is None or m["flank"] > WC_PAIRED_A:
                continue
            nclash = count_clashes(res)
            if best is None or nclash < best[0]:
                best = (nclash, ang_a, ang_b, na, nb)
    if best is None:
        sys.exit("no combination of base rotations both breaks the central pair and "
                 "leaves the flanking pairs intact -- check the step index")
    for k, v in ra.items():
        v["xyz"] = saved[0][k].copy()
    for k, v in rb.items():
        v["xyz"] = saved[1][k].copy()
    rotate_base(ra, best[1]); rotate_base(rb, best[2])
    report["base_rotations_deg"] = {"strand_A": best[1], "strand_B": best[2]}
    report["search"] = f"least-clashing of 25 combinations that still form the motif"

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

    # ---- refuse to write a structure that is already wrecked ----
    # Seeding into an ALREADY BENT duplex folds the arc onto itself: an 85 deg roll on top
    # of a 122 deg arc drove the downstream half into the upstream half, leaving 60 atom
    # pairs closer than 1.5 A with the worst at 0.58 A.  Minimisation cannot undo that --
    # it relaxes into whatever tangle it is handed, and here it left a guanine ring bond
    # stretched to 2.66 A, which then produced NaN coordinates six minutes into production
    # on all six runs.  Cheaper to check here.
    from scipy.spatial import cKDTree
    xyz = np.array([a["xyz"] for a in atoms])
    who = [(a["chain"], a["resi"]) for a in atoms]
    tree = cKDTree(xyz)
    clashes = [(float(np.linalg.norm(xyz[i] - xyz[j])), i, j)
               for i, j in tree.query_pairs(CLASH_A, output_type="ndarray")
               if who[i] != who[j]]
    clashes.sort()
    report["clashes_under_%.1fA" % CLASH_A] = len(clashes)
    if clashes:
        i, j = clashes[0][1], clashes[0][2]
        report["worst_clash"] = {
            "distance_A": round(clashes[0][0], 2),
            "between": f"{who[i][0]}{who[i][1]}.{atoms[i]['name']}"
                       f" and {who[j][0]}{who[j][1]}.{atoms[j]['name']}"}
    if len(clashes) > MAX_CLASHES:
        sys.exit(
            f"{len(clashes)} inter-residue atom pairs are closer than {CLASH_A} A "
            f"(worst {clashes[0][0]:.2f} A, {report['worst_clash']['between']}).\n"
            "That is a wrecked structure, not a seed.  Two usual causes: seeding into a "
            "duplex that is ALREADY bent, where an 85 deg roll on top of an existing arc "
            "folds the helix into itself (use build/03_dna/straight_*.pdb and impose the "
            "compression afterwards with the harmonic restraint); or a step whose local "
            "sequence simply will not accommodate the motif, in which case try an "
            "adjacent one -- the position along the duplex is what matters, not the "
            "particular step.")
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
    br = rep.get("base_rotations_deg", {})
    print(f"  both central bases swung out of the stack: strand A "
          f"{br.get('strand_A', 0):+.0f} deg, strand B {br.get('strand_B', 0):+.0f} deg "
          f"({rep.get('search','')})")
    print(f"  rise opened at the kink: {rep.get('rise_A', 0):.2f} A; "
          f"residual clashes under 1.5 A: {rep.get('clashes_under_1.5A', 0)}")
    print(f"  N1-N3 at pairs {list(rep['N1_N3_A'])}: "
          f"{list(rep['N1_N3_A'].values())} A  (paired is ~2.8, broken is >4)")
    print(f"  local bend across the motif: {rep['local_bend_deg']} deg")
    print(f"  backbone O3'-P: mean {rep['O3prime_P'][0]}, max {rep['O3prime_P'][1]} A "
          f"(ideal 1.61; minimisation takes the strain out)")
    print(f"  wrote {a.out}")


if __name__ == "__main__":
    main()
