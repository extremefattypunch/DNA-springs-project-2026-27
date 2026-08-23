#!/usr/bin/env python3
"""Distances that are longer than the periodic box's minimum image can express.

The bug this exists to prevent
------------------------------
``mdtraj.compute_distances`` applies the minimum-image convention by default
(``periodic=True``).  That is the right thing for a hydrogen bond and the wrong thing for
the end-to-end distance of a DNA spring, because minimum image is only *guaranteed*
correct out to half the smallest perpendicular width of the cell.  Our solvent boxes are
truncated octahedra whose limits are 38-48 A, and the spring's ends are 60-130 A apart.

Measured on real trajectories, the damage was:

    system          periodic=True     raw        disagreement
    S3_spring27       53.80 A        70.73 A       up to 30.7 A
    S4_spring40       62.13          67.27         up to 13.7
    L3_loop27         45.37          71.36         up to 41.9

and the give-away is the spread: the raw distances have a standard deviation of 2.7-4.2 A
while the minimum-imaged ones have 6.7-8.8, because the reported value flips between
images from frame to frame.  Every force this project quotes from Zocchi's f(x) was
evaluated at one of those wrong distances.

Why the raw coordinates are the right answer here
-------------------------------------------------
OpenMM's reporters write positions with ``enforcePeriodicBox``, which translates whole
*molecules* so that each sits in the box.  It never splits a molecule.  So two atoms in
the same covalently connected molecule are always contiguous in the frame as written, and
their raw separation is the true one, however large.  Two atoms in *different* molecules
are a different matter: one may have been translated by a lattice vector relative to the
other, and then neither convention is right beyond the minimum-image limit.

So the rule this module enforces:

* same molecule, any distance          -> raw coordinates (periodic=False)
* different molecules, short distance  -> minimum image (periodic=True)
* different molecules, long distance   -> refuse; there is no answer to be had from a
                                          frame written this way

and it *checks* the reasoning rather than trusting it: whenever a raw distance comes out
below the minimum-image limit, the two conventions have to agree, and if they do not then
the contiguity assumption is wrong and we want to hear about it immediately.
"""
from __future__ import annotations

import numpy as np


def min_image_limit(traj) -> float:
    """Half the smallest perpendicular width of the unit cell, in angstrom.

    For a general triclinic cell the perpendicular width along one axis is the volume
    divided by the area of the opposite face, which is what this computes; for a cube it
    reduces to the side length.  Minimum image is only guaranteed below half of the
    smallest of the three.
    """
    if traj.unitcell_vectors is None:
        return float("inf")
    uc = np.asarray(traj.unitcell_vectors[0], float) * 10.0      # nm -> A
    a, b, c = uc
    widths = []
    for u, v, w in ((a, b, c), (b, c, a), (c, a, b)):
        n = np.cross(v, w)
        widths.append(abs(float(u @ n)) / float(np.linalg.norm(n)))
    return 0.5 * min(widths)


def molecule_ids(top, cut=()):
    """atom index -> molecule id, molecules being connected components of the bonds.

    `cut` is a set of frozenset({i, j}) bonds to ignore.  Cutting the two C-alpha--C-beta
    bonds at the modified residues turns this into the covalent-cycle test: afterwards the
    two attachment atoms are in one component only if a path runs out along one tether,
    through the DNA and back along the other.  Without the cut every construct looks
    connected, because both tethers hang off the protein -- which is exactly the mistake
    that let the original two-arm build pass as a closed loop.

    Built here rather than taken from ``top.find_molecules()`` so the result is a plain
    array that can be indexed cheaply, and so the same definition is used everywhere.
    """
    n = top.n_atoms
    parent = list(range(n))

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    for b in top.bonds:
        if frozenset((b[0].index, b[1].index)) in cut:
            continue
        i, j = find(b[0].index), find(b[1].index)
        if i != j:
            parent[i] = j
    return np.array([find(i) for i in range(n)])


def long_distance(traj, i: int, j: int, label: str = "", mols=None):
    """Distance between atoms i and j in angstrom, per-frame, done correctly.

    Raises if the pair spans two molecules and is too far apart for the minimum image to
    mean anything -- a silent wrong number there would be worse than a stack trace.
    """
    import mdtraj as md

    pair = np.array([[int(i), int(j)]])
    raw = md.compute_distances(traj, pair, periodic=False)[:, 0] * 10.0
    per = md.compute_distances(traj, pair, periodic=True)[:, 0] * 10.0
    limit = min_image_limit(traj)

    if mols is None:
        mols = molecule_ids(traj.topology)
    same = bool(mols[int(i)] == mols[int(j)])

    if same:
        # The molecule is contiguous as written, so the raw separation is the true one.
        # Where it happens to be short, minimum image must agree; if it does not, the
        # contiguity assumption has failed and every long distance is suspect.
        short = raw < 0.98 * limit
        if short.any():
            bad = float(np.abs(raw[short] - per[short]).max())
            if bad > 1e-3:
                raise AssertionError(
                    f"{label or f'atoms {i},{j}'}: the two atoms are in one molecule and "
                    f"{short.sum()} frames put them inside the minimum-image limit "
                    f"({limit:.1f} A), yet the conventions disagree by up to {bad:.3f} A. "
                    "The frames are not molecule-contiguous, so neither convention can "
                    "be trusted -- image the molecules before analysing.")
        return raw, {"convention": "raw (one molecule, contiguous as written)",
                     "min_image_limit_A": round(limit, 2),
                     "frames_beyond_limit": int((raw >= 0.98 * limit).sum()),
                     "max_disagreement_A": round(float(np.abs(raw - per).max()), 3)}

    if not (per >= 0.98 * limit).any():
        return per, {"convention": "minimum image (separate molecules, within the limit)",
                     "min_image_limit_A": round(limit, 2),
                     "frames_beyond_limit": 0,
                     "max_disagreement_A": round(float(np.abs(raw - per).max()), 3)}

    # Different molecules, and too far apart for the minimum image on the ATOM PAIR to
    # mean anything -- taking it anyway is what turned a nicked duplex held at 98 A into a
    # reported 22.8 A and a force of -1400 pN.
    #
    # There is still a right answer when the two molecules are parts of one assembly, as
    # they are for a nicked strand whose halves are base-paired to a common complement.
    # Choose the periodic image of the second MOLECULE that sits closest to the first --
    # a decision made on the molecules' centres of geometry, which are only a fraction of
    # the box apart and so unambiguous -- and then measure the atoms inside that choice.
    # Doing it atom-by-atom is what fails: the two atoms are at opposite ends of the
    # assembly and their own nearest image is the wrong one.
    mi, mj = int(mols[int(i)]), int(mols[int(j)])
    ai = np.where(mols == mi)[0]
    aj = np.where(mols == mj)[0]
    xyz = traj.xyz * 10.0                                  # nm -> A
    ci = xyz[:, ai, :].mean(axis=1)
    cj = xyz[:, aj, :].mean(axis=1)
    H = np.asarray(traj.unitcell_vectors, float) * 10.0    # (frames, 3, 3), rows a/b/c
    shifts = np.array([[a, b, c] for a in (-1, 0, 1) for b in (-1, 0, 1)
                       for c in (-1, 0, 1)], float)
    out = np.empty(len(xyz))
    used = np.zeros((len(xyz), 3))
    for fr in range(len(xyz)):
        cand = cj[fr] + shifts @ H[fr]                     # every image of B's centre
        best = int(np.argmin(np.linalg.norm(cand - ci[fr], axis=1)))
        used[fr] = shifts[best]
        out[fr] = float(np.linalg.norm(
            (xyz[fr, int(j)] + shifts[best] @ H[fr]) - xyz[fr, int(i)]))
    moved = int((np.abs(used).sum(axis=1) > 0).sum())
    return out, {"convention": "whole-molecule image (separate molecules, assembled by "
                               "their centres of geometry)",
                 "min_image_limit_A": round(limit, 2),
                 "frames_beyond_limit": int((out >= 0.98 * limit).sum()),
                 "frames_needing_a_lattice_shift": moved,
                 "max_disagreement_A": round(float(np.abs(out - per).max()), 3)}


def audit(traj, pairs, mols=None):
    """Report, for a list of (label, i, j), which convention each needs and why."""
    out = []
    if mols is None:
        mols = molecule_ids(traj.topology)
    for label, i, j in pairs:
        try:
            d, info = long_distance(traj, i, j, label, mols)
            out.append({"label": label, "mean_A": round(float(d.mean()), 3), **info})
        except AssertionError as e:
            out.append({"label": label, "error": str(e)})
    return out
