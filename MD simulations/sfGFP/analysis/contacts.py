#!/usr/bin/env python3
"""How much of the protein the DNA actually touches.

Why this is the measurement that matters
----------------------------------------
The report treated the DNA spring as a device that applies a force and otherwise stays
out of the way.  Doeke's objection to the force argument is correct on its own terms --
at equilibrium a couple applied through two tethers averages away, and only the
component of force along the coordinate you are watching can displace it -- so if the
chimeras differ from a matched axial pull, the difference has to come from something
that is not a force.  The obvious candidate is that the duplex is a large charged object
tethered a couple of nanometres from the barrel, and it leans on it.

That is also the conclusion of the ubiquitin preprint from this lab: a control bearing
two separate short strands reproduced most of the circular chimera's NMR spectrum, which
points at local protein-DNA interaction rather than spring action.  So this module
measures contact directly -- how many protein atoms the DNA is within reach of, which
residues, and how persistent it is -- for every system that has DNA.

Contacts are counted with a neighbour search rather than a full pair matrix: the naive
version builds a ~1800 x 1700 distance array per frame and was killed by the OOM
reaper on a full-length trajectory.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import mdtraj as md
import numpy as np
import pandas as pd

CONTACT_A = 4.5
# A contact occupancy is a mean over frames, and a mean does not need every frame.  The
# neighbour search runs once per frame from Python, so the cost grows with trajectory
# length -- 145 s per replicate over 146 of them is six hours to answer a question that
# 120 evenly spaced frames answer to better than a percent.  Capping the frame count makes
# the cost per replicate constant.
MAX_FRAMES = 120
NUCLEOTIDES = ("DA", "DT", "DG", "DC")
LINKERS = ("DNL", "DNP", "DNH")
NON_SOLUTE = ("WAT", "HOH", "Na+", "Cl-", "MG", "K+")


def load_map(path):
    return {r["index"]: (r["orig_resnum"], r["orig_resname"], r.get("orig_chain"))
            for r in json.loads(Path(path).read_text())}


def split_atoms(t, resmap, include_linkers=False):
    """(DNA heavy atoms, protein heavy atoms, protein residue index per atom).

    The linkers are excluded from the DNA side by default.  They are covalently bonded
    to the protein, so they are always "in contact" with the residues they are attached
    to and counting them would swamp the signal being looked for, which is the duplex
    itself reaching the barrel.
    """
    dna, prot, owner = [], [], {}
    for i, (num, name, ch) in resmap.items():
        base = name.rstrip("35")
        heavy = [a.index for a in t.top.residue(i).atoms if a.element.symbol != "H"]
        if base in NUCLEOTIDES or (include_linkers and name in LINKERS):
            dna += heavy
        elif name in NON_SOLUTE or name in LINKERS:
            continue
        else:
            prot += heavy
            for a in heavy:
                owner[a] = (num, name)
    return np.array(dna), np.array(prot), owner


def analyse(traj, top_path, resmap_path, out_prefix, stride=1, skip=0):
    t = md.load(str(traj), top=str(top_path), stride=stride)
    if skip and len(t) > skip:
        t = t[skip:]
    n_all = len(t)
    if n_all > MAX_FRAMES:
        t = t[:: int(np.ceil(n_all / MAX_FRAMES))]
    resmap = load_map(resmap_path)
    dna, prot, owner = split_atoms(t, resmap)
    if len(dna) == 0 or len(prot) == 0:
        return {"has_dna": False}

    per_frame, tally = [], {}
    for f in range(len(t)):
        hit = md.compute_neighbors(t[f], CONTACT_A / 10.0, dna,
                                   haystack_indices=prot)[0]
        per_frame.append(len(hit))
        for a in hit:
            num, name = owner[int(a)]
            tally[(num, name)] = tally.get((num, name), 0) + 1
    per_frame = np.array(per_frame, float)
    n = len(t)
    ranked = sorted(tally.items(), key=lambda kv: -kv[1])
    out = {
        "has_dna": True, "frames": n, "frames_available": n_all,
        "cutoff_A": CONTACT_A,
        "atoms_in_contact_mean": round(float(per_frame.mean()), 2),
        "atoms_in_contact_sd": round(float(per_frame.std()), 2),
        "frames_with_contact": round(float((per_frame > 0).mean()), 4),
        "residues_ever_touched": len(tally),
        # A residue in contact for most of the trajectory is a real interaction site,
        # not a passing collision, so those are what get named.
        "residues_over_half_the_time": [
            f"{nm}{num}" for (num, nm), c in ranked if c / n > 0.5],
        "top_residues": [{"residue": f"{nm}{num}", "occupancy": round(c / n, 3)}
                         for (num, nm), c in ranked[:12]],
    }
    if out_prefix:
        pd.DataFrame({"contact_atoms": per_frame}).to_csv(
            f"{out_prefix}_contacts.csv", index=False)
        Path(f"{out_prefix}_contacts.json").write_text(json.dumps(out, indent=2))
    print(f"  DNA-protein contact: {out['atoms_in_contact_mean']:.1f} "
          f"+/- {out['atoms_in_contact_sd']:.1f} protein heavy atoms within "
          f"{CONTACT_A} A, in {100 * out['frames_with_contact']:.0f}% of frames; "
          f"{len(out['residues_over_half_the_time'])} residues touched over half "
          f"the time")
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--traj", required=True)
    ap.add_argument("--top", required=True)
    ap.add_argument("--resmap", required=True)
    ap.add_argument("--stride", type=int, default=5)
    ap.add_argument("--skip", type=int, default=0)
    ap.add_argument("--out-prefix", default=None)
    a = ap.parse_args()
    print(json.dumps(analyse(a.traj, a.top, a.resmap, a.out_prefix, a.stride, a.skip),
                     indent=2))


if __name__ == "__main__":
    main()
