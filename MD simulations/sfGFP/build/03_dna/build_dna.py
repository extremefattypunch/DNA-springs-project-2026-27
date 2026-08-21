#!/usr/bin/env python3
"""Build the dsDNA spring: an ideal B-form duplex, then bend it to a target span.

Straight duplex
---------------
PyMOL's ``fnab`` builds ideal B-form DNA with the correct 5'-phosphates.  Three
things about its output have to be normalised before Amber will take it:

* records are HETATM, not ATOM;
* the phosphate oxygens use the old PDB v2 names O1P/O2P, where Amber wants OP1/OP2;
* the complementary strand is numbered **-N..-1**, which silently breaks any
  residue-number filter written for a 1..N chain.

Both strands run 5'->3' with ascending number, and the two 5' ends sit at opposite
ends of the duplex, each carrying a phosphate whose fourth valence is free.  That is
exactly the Zocchi attachment geometry -- one 5' end to each protein site -- and it
means the linker's bridging oxygen can take the structural role of the preceding
residue's O3' with no phosphate surgery.

Bending
-------
The duplex is bent by treating each base pair as a rigid slab and mapping the
straight helical axis onto a circular arc of radius R, chosen so the end-to-end
distance matches the span the protein plus linkers impose.  Moving whole base pairs
rigidly means intra-pair geometry -- the hydrogen bonds and the stacking within a
pair -- is preserved exactly; only the inter-pair geometry changes, which is what
bending physically *is*.  The backbone takes up the strain, and the resulting
O3'-P deviations are reported rather than hidden.

For the spans this project needs, R comes out near 37-40 A, which is essentially
nucleosomal curvature (~42 A) -- strongly bent but thoroughly precedented, not a
structure held together by hope.  Whether the duplex then relieves that strain by
kinking is left to the simulation to answer, since that softening transition is the
thing Zocchi's model is really about.
"""
from __future__ import annotations

import argparse
import json
import math
import subprocess
import sys
from pathlib import Path

import numpy as np

# Zocchi's 2021 review oligos, verified base-by-base against their own splint and
# cDNA.  The 2013 JACS arm-B string does not complement its own stated 60mer and
# must not be used.  Springs are cut from the centre of this 60mer.
ZOCCHI_60MER = "CAGCTGCTTGGATGGTACCGTGGACTCCTGCCAGCAACTCACGGTCTAGGCTCCACACTC"
PRIME3 = "O3'"


def centre_slice(seq: str, n: int) -> str:
    if n > len(seq):
        sys.exit(f"cannot cut {n} bp from a {len(seq)}-mer")
    start = (len(seq) - n) // 2
    return seq[start:start + n]


def run_fnab(seq: str, out: Path, pymol: str) -> Path:
    pml = out.with_suffix(".pml")
    pml.write_text(f"fnab {seq}, name=spring, mode=DNA, form=B, dbl_helix=1\n"
                   f"save {out.name}, spring\n")
    subprocess.run([pymol, "-cq", pml.name], cwd=out.parent, check=True,
                   capture_output=True, text=True)
    if not out.exists():
        sys.exit("fnab produced no output")
    return out


def read_pdb(path: Path):
    atoms = []
    for l in path.read_text().splitlines():
        if l[:6] not in ("ATOM  ", "HETATM"):
            continue
        elem = (l[76:78].strip() or l[12:16].strip()[0]).upper()
        if elem == "H":
            continue
        name = l[12:16].strip()
        name = {"O1P": "OP1", "O2P": "OP2"}.get(name, name)
        atoms.append({"name": name, "resn": l[17:20].strip(), "chain": l[21],
                      "resi": int(l[22:26]), "elem": elem,
                      "xyz": np.array([float(l[30:38]), float(l[38:46]),
                                       float(l[46:54])])})
    return atoms


def renumber(atoms):
    """Make both strands run 1..N ascending in the 5'->3' direction."""
    for ch in sorted({a["chain"] for a in atoms}):
        nums = sorted({a["resi"] for a in atoms if a["chain"] == ch})
        shift = 1 - nums[0]
        if shift:
            for a in atoms:
                if a["chain"] == ch:
                    a["resi"] += shift
    return atoms


def residues(atoms):
    out = {}
    for a in atoms:
        out.setdefault((a["chain"], a["resi"]), {})[a["name"]] = a
    return out


def verify(atoms, n_bp) -> dict:
    res = residues(atoms)
    chains = sorted({c for c, _ in res})
    if len(chains) != 2:
        sys.exit(f"expected 2 chains, found {chains}")
    a, b = chains
    rep = {"chains": chains, "n_bp": n_bp}

    for ch in chains:
        nums = sorted(n for c, n in res if c == ch)
        if nums != list(range(1, n_bp + 1)):
            sys.exit(f"chain {ch} is not numbered 1..{n_bp}: {nums[:3]}..{nums[-3:]}")
        links = sum(np.linalg.norm(res[(ch, i)][PRIME3]["xyz"]
                                   - res[(ch, i + 1)]["P"]["xyz"]) < 2.0
                    for i in nums[:-1])
        if links != n_bp - 1:
            sys.exit(f"chain {ch}: only {links}/{n_bp - 1} O3'->P links ascending; "
                     "the numbering does not run 5'->3'")
        if "P" not in res[(ch, 1)]:
            sys.exit(f"chain {ch} residue 1 has no 5' phosphate to attach a linker to")

    # antiparallel: A:i pairs with B:n+1-i
    d = [np.linalg.norm(res[(a, i)]["C1'"]["xyz"] - res[(b, n_bp + 1 - i)]["C1'"]["xyz"])
         for i in range(1, n_bp + 1)]
    rep["C1'-C1' mean_A"] = round(float(np.mean(d)), 3)
    if not 9.5 < np.mean(d) < 11.5:
        sys.exit(f"base pairing looks wrong: mean C1'-C1' {np.mean(d):.2f} A "
                 "(Watson-Crick B-form is ~10.5)")

    ctr = np.array([(res[(a, i)]["C1'"]["xyz"]
                     + res[(b, n_bp + 1 - i)]["C1'"]["xyz"]) / 2
                    for i in range(1, n_bp + 1)])
    c0 = ctr.mean(axis=0)
    axis = np.linalg.svd(ctr - c0)[2][0]
    if (ctr[-1] - ctr[0]) @ axis < 0:
        axis = -axis
    z = (ctr - c0) @ axis
    # The 3D centre-to-centre spacing is *not* the helical rise: the centres spiral
    # around the axis, so it is sqrt(rise^2 + lateral^2).  Report the axial rise.
    spacing = np.linalg.norm(np.diff(ctr, axis=0), axis=1)
    rep["bp_centre_spacing_A"] = [round(float(spacing.mean()), 3),
                                  round(float(spacing.std()), 3)]
    rise = np.diff(z)
    rep["axial_rise_A"] = [round(float(rise.mean()), 3), round(float(rise.std()), 3)]
    rep["straightness_rms_A"] = round(
        float(np.linalg.norm((ctr - c0) - np.outer(z, axis), axis=1).mean()), 3)

    v = np.array([res[(a, i)]["C1'"]["xyz"] - res[(b, n_bp + 1 - i)]["C1'"]["xyz"]
                  for i in range(1, n_bp + 1)])
    vp = v - np.outer(v @ axis, axis)
    vp /= np.linalg.norm(vp, axis=1)[:, None]
    tw = np.degrees(np.arccos(np.clip((vp[:-1] * vp[1:]).sum(axis=1), -1, 1)))
    rep["twist_deg"] = [round(float(tw.mean()), 2), round(float(tw.std()), 2)]

    # both 5' ends must be at opposite ends of the duplex
    za = (res[(a, 1)]["P"]["xyz"] - c0) @ axis
    zb = (res[(b, 1)]["P"]["xyz"] - c0) @ axis
    rep["five_prime_z_A"] = [round(float(za), 2), round(float(zb), 2)]
    if za * zb > 0:
        sys.exit("both 5' ends are at the same end of the duplex")
    rep["five_prime_separation_A"] = round(
        float(np.linalg.norm(res[(a, 1)]["P"]["xyz"] - res[(b, 1)]["P"]["xyz"])), 2)
    return rep, ctr, axis, c0


def radius_for_chord(L: float, chord: float) -> float:
    """Solve 2R sin(L/2R) = chord for R by bisection."""
    if chord >= L:
        sys.exit(f"target span {chord:.2f} A is not shorter than the contour "
                 f"{L:.2f} A: the duplex would be stretched, not bent")
    lo, hi = L / (2 * math.pi), 1e6          # lo: a full circle; hi: nearly straight
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        c = 2 * mid * math.sin(L / (2 * mid))
        if c < chord:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def bend(atoms, ctr, axis, c0, n_bp, chord: float, bend_azimuth_deg: float = 0.0):
    """Map the straight helical axis onto a circular arc of the required radius."""
    res = residues(atoms)
    a, b = sorted({c for c, _ in res})
    z = (ctr - c0) @ axis
    L = float(z[-1] - z[0])
    R = radius_for_chord(L, chord)
    theta = L / R

    # local frame: e3 along the axis, e1 the bend direction, e2 = e3 x e1
    tmp = np.array([1.0, 0.0, 0.0])
    if abs(tmp @ axis) > 0.9:
        tmp = np.array([0.0, 1.0, 0.0])
    e1 = tmp - (tmp @ axis) * axis
    e1 /= np.linalg.norm(e1)
    e2 = np.cross(axis, e1)
    phi0 = math.radians(bend_azimuth_deg)
    e1, e2 = math.cos(phi0) * e1 + math.sin(phi0) * e2, \
             -math.sin(phi0) * e1 + math.cos(phi0) * e2
    M = np.stack([e1, e2, axis])          # lab -> local (rows are the basis)

    # Each base pair rotates about the point where its own plane crosses the *helical
    # axis*, not about its centre.  Base-pair centres spiral around the axis at ~1.9 A,
    # so rotating about the centre and then dropping that centre onto the arc silently
    # replaces every pair's own off-axis offset with the first pair's -- a lateral
    # shift of up to 3.8 A that tears the backbone open even at zero curvature.  With
    # the axis point as the pivot the transform is the exact identity as R -> infinity,
    # which is the property to test against.
    z0 = float(z[0])
    for i in range(1, n_bp + 1):
        s = float(z[i - 1] - z0)
        phi = s / R
        axis_point = np.array([0.0, 0.0, z0 + s])         # on the straight axis
        arc_point = np.array([R * (1 - math.cos(phi)), 0.0, z0 + R * math.sin(phi)])
        c, sn = math.cos(phi), math.sin(phi)
        Ry = np.array([[c, 0.0, sn], [0.0, 1.0, 0.0], [-sn, 0.0, c]])
        for key in ((a, i), (b, n_bp + 1 - i)):
            for at in res[key].values():
                loc = M @ (at["xyz"] - c0)                # into the local frame
                at["xyz"] = c0 + M.T @ (arc_point + Ry @ (loc - axis_point))

    # report the strain the backbone had to absorb
    res = residues(atoms)
    dev = []
    for ch in (a, b):
        for i in range(1, n_bp):
            dev.append(np.linalg.norm(res[(ch, i)][PRIME3]["xyz"]
                                      - res[(ch, i + 1)]["P"]["xyz"]))
    dev = np.array(dev)
    new_span = float(np.linalg.norm(res[(a, 1)]["P"]["xyz"] - res[(b, 1)]["P"]["xyz"]))
    return {"contour_A": round(L, 2), "radius_A": round(R, 2),
            "arc_deg": round(math.degrees(theta), 1),
            "target_chord_A": round(chord, 2),
            "five_prime_separation_after_A": round(new_span, 2),
            "O3prime_P_after": [round(float(dev.mean()), 3), round(float(dev.std()), 3),
                                round(float(dev.max()), 3)],
            "bend_azimuth_deg": bend_azimuth_deg}


def write_pdb(atoms, path: Path, nick_after: int | None = None):
    """Write with Amber-friendly naming.  A nick is expressed as 3'/5' termini."""
    res = residues(atoms)
    chains = sorted({c for c, _ in res})
    n_bp = max(n for _, n in res)
    lines, serial = [], 1
    for ch in chains:
        for i in range(1, n_bp + 1):
            rname = res[(ch, i)][next(iter(res[(ch, i)]))]["resn"]
            if nick_after is not None and ch == chains[0]:
                if i == nick_after:
                    rname += "3"          # 3'-OH terminus before the nick
                elif i == nick_after + 1:
                    rname += "5"          # 5'-OH terminus after it
            for at in res[(ch, i)].values():
                if nick_after is not None and ch == chains[0] and i == nick_after + 1 \
                        and at["name"] in ("P", "OP1", "OP2"):
                    continue              # the nicked strand loses that phosphate
                nm = at["name"]
                fld = f" {nm:<3}" if len(nm) < 4 else nm
                lines.append(
                    f"ATOM  {serial:5d} {fld}{rname:>4}{ch:>2}{i:4d}    "
                    f"{at['xyz'][0]:8.3f}{at['xyz'][1]:8.3f}{at['xyz'][2]:8.3f}"
                    f"{1.00:6.2f}{0.00:6.2f}          {at['elem']:>2}")
                serial += 1
            if nick_after is not None and ch == chains[0] and i == nick_after:
                lines.append(f"TER   {serial:5d}")
                serial += 1
        lines.append(f"TER   {serial:5d}")
        serial += 1
    path.write_text("\n".join(lines) + "\nEND\n")


def main():
    ap = argparse.ArgumentParser()
    here = Path(__file__).resolve().parent
    ap.add_argument("--n-bp", type=int, required=True)
    ap.add_argument("--span-A", type=float, required=True,
                    help="required 5'P-to-5'P distance in angstrom")
    ap.add_argument("--sequence", default=None,
                    help="default: the centre n-bp of Zocchi's verified 60mer")
    ap.add_argument("--nick-after", type=int, default=None,
                    help="nick the first strand after this residue (centre = n/2)")
    ap.add_argument("--azimuth", type=float, default=0.0)
    ap.add_argument("--outdir", default=str(here))
    ap.add_argument("--pymol", default=None)
    args = ap.parse_args()

    import shutil
    pymol = args.pymol or shutil.which("pymol")
    if not pymol:
        sys.exit("pymol not found; source env/activate.sh")

    out = Path(args.outdir)
    out.mkdir(parents=True, exist_ok=True)
    seq = args.sequence or centre_slice(ZOCCHI_60MER, args.n_bp)
    tag = f"{args.n_bp}bp" + ("_nick" if args.nick_after else "")
    print(f"=== {tag}: {seq} ===")

    raw = run_fnab(seq, out / f"straight_{tag}.pdb", pymol)
    atoms = renumber(read_pdb(raw))
    rep, ctr, axis, c0 = verify(atoms, args.n_bp)
    print(f"  straight: axial rise {rep['axial_rise_A'][0]} +/- "
          f"{rep['axial_rise_A'][1]} A, "
          f"twist {rep['twist_deg'][0]} +/- {rep['twist_deg'][1]} deg, "
          f"off-axis rms {rep['straightness_rms_A']} A")
    print(f"  5' ends at z = {rep['five_prime_z_A']} A, separation "
          f"{rep['five_prime_separation_A']} A")

    brep = bend(atoms, ctr, axis, c0, args.n_bp, args.span_A, args.azimuth)
    print(f"  bent to R = {brep['radius_A']} A ({brep['arc_deg']} deg of arc); "
          f"5'-5' now {brep['five_prime_separation_after_A']} A "
          f"(target {brep['target_chord_A']})")
    print(f"  backbone strain: O3'-P {brep['O3prime_P_after'][0]} +/- "
          f"{brep['O3prime_P_after'][1]} A (max {brep['O3prime_P_after'][2]}); "
          f"ideal is 1.61")

    pdb = out / f"spring_{tag}.pdb"
    write_pdb(atoms, pdb, nick_after=args.nick_after)
    (out / f"spring_{tag}.json").write_text(json.dumps(
        {"sequence": seq, "straight": rep, "bend": brep,
         "nick_after": args.nick_after}, indent=2))
    print(f"  wrote {pdb.name}")


if __name__ == "__main__":
    main()
