#!/usr/bin/env python3
"""Backbone conformation at the two modified sites -- what Tet2-Et does locally.

Why this module exists
----------------------
Doeke noticed that putting the two Tet2-Et residues in *shortens* the distance between
the attachment sites, and suggested a local secondary-structure change.  The observation
is real and the question is the right one to ask, so this measures the three things that
could account for it, at both sites, in every system:

  1. the backbone dihedrals phi and psi, which is where a local conformational change
     would show up first;
  2. DSSP secondary structure over a window around each site, which is where his
     specific hypothesis would show up;
  3. contacts between the unclicked tetrazine ring and Tyr200, because one of the two
     sites has an aromatic neighbour and the other does not.

The answer that comes out is sharper than a secondary-structure change and is the reason
this is worth a module rather than a footnote.  Site 133 is an aspartate in the wild
type, and aspartate is one of the few residues that commonly sits in the *left-handed*
(positive-phi) alpha region -- a conformation that a phenylalanine-shaped side chain,
which is what Tet2-Et is, tolerates poorly.  Site 149 is ordinary beta.  So the two
sites are not equivalent targets, and that has a practical consequence: it predicts
which positions in a protein can take this amino acid without local distortion.

Circular statistics throughout.  phi at site 133 sits near +60 and psi near +20, well
away from the +/-180 branch cut, but the same code has to read site 149's psi at +137
and the analysis should not silently depend on which side of the cut a residue falls.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import mdtraj as md
import numpy as np
import pandas as pd

# The tetrazine (TET) and its clicked dihydropyridazine (TDP) both carry the ring that
# could stack on an aromatic neighbour; these are the ring heavy atoms in each.
# Names taken from build/02_params/residue_defs.py, not guessed: TET's tetrazine is
# CT3-NT2-NT1-CT6-NT5-NT4 and clicking it opens that ring to TDP's dihydropyridazine
# CP3-NP2-NP1-CP6-CQ5-CQ4, having expelled N2.  Getting these wrong is silent -- the
# atom filter simply matches nothing and the contact count never appears.
RING_ATOMS = {
    "TET": ("CT3", "NT2", "NT1", "CT6", "NT5", "NT4"),
    "TDP": ("CP3", "NP2", "NP1", "CP6", "CQ5", "CQ4"),
}
AROMATIC_NEIGHBOUR = ("TYR", 200)
STACK_CUTOFF_A = 4.5     # heavy-atom contact distance for a stacking contact
DSSP_WINDOW = (128, 155)


def load_map(path):
    """prmtop residue index -> (original residue number, name, original chain).

    The chain is carried through because the residue *number* is not unique: a 50-mer
    of DNA is numbered 1-50 in its own chain and the protein is numbered from 1 in its
    own, so in the ubiquitin chimeras 52 of 130 residue numbers name two different
    residues.  Dropping the chain here is what made every lookup of ubiquitin's site 9
    return DNA nucleotide 9 instead.
    """
    rows = json.loads(Path(path).read_text())
    return {r["index"]: (r["orig_resnum"], r["orig_resname"], r.get("orig_chain"))
            for r in rows}


def drop_equilibration(t, skip):
    if skip and len(t) > skip:
        return t[skip:]
    return t


def circular_stats(deg):
    """Circular mean and sd in degrees, and the mean resultant length R.

    R runs from 0 (uniform) to 1 (a delta function), so it doubles as a check that the
    'mean' being reported means anything: a bimodal residue flipping between two rotamers
    has a small R and a mean that sits in neither well.
    """
    r = np.radians(np.asarray(deg, float))
    c, s = np.cos(r).mean(), np.sin(r).mean()
    R = float(np.hypot(c, s))
    mean = float(np.degrees(np.arctan2(s, c)))
    sd = float(np.degrees(np.sqrt(-2.0 * np.log(max(R, 1e-12)))))
    return mean, sd, R


# Amber's nucleic residue names, including the 3'/5' terminal variants.  Listed rather
# than sniffed so a residue this project has not met yet fails loudly.
_NUCLEIC = {f"D{b}{s}" for b in "ATGC" for s in ("", "3", "5", "N")}


def _by_number(resmap, prefer_protein=True):
    """residue number -> prmtop index, resolving the number collision explicitly.

    Ubiquitin is modified at residues 9 and 63 and carries 50 nucleotides, so number 9
    names both a cysteine and a guanine.  Whichever the dict built last silently won,
    and it was the DNA: site 63 resolved (no nucleotide 63) while site 9 did not, so
    the three DNA-bearing ubiquitin systems reported no backbone angles at all while
    every sfGFP system was fine -- sfGFP's sites are 133 and 149, above any nucleotide
    number, so the collision could not reach them.

    Every caller in this module wants a protein residue, so protein wins by default and
    the ambiguity is recorded rather than resolved by dict-insertion order.
    """
    out, ambiguous = {}, {}
    for idx, entry in resmap.items():
        num, name = entry[0], entry[1]
        is_dna = name in _NUCLEIC
        if num in out:
            ambiguous.setdefault(num, []).append(name)
            # keep whichever side the caller asked for; on a tie keep the first
            keep_new = (not is_dna) if prefer_protein else is_dna
            prev_name = resmap[out[num]][1]
            prev_is_dna = prev_name in _NUCLEIC
            if keep_new and (prev_is_dna if prefer_protein else not prev_is_dna):
                out[num] = idx
            continue
        out[num] = idx
    return out


def backbone_dihedrals(t, top, resmap, resnum):
    """phi = C(i-1)-N(i)-CA(i)-C(i);  psi = N(i)-CA(i)-C(i)-N(i+1), in degrees.

    Computed from named atoms in explicitly identified residues rather than with
    mdtraj's compute_phi/compute_psi, whose output is an array over *all* residues that
    then has to be indexed -- and indexing it correctly requires knowing how mdtraj
    numbered a topology that contains DNA chains and a custom residue.  Naming the four
    atoms cannot drift.
    """
    idx_of = _by_number(resmap)
    if resnum not in idx_of:
        raise KeyError(f"residue {resnum} is not in the residue map")
    i = idx_of[resnum]
    prev, nxt = idx_of.get(resnum - 1), idx_of.get(resnum + 1)

    def atom(ridx, name):
        if ridx is None:
            return None
        return next((a.index for a in top.residue(ridx).atoms if a.name == name), None)

    quads, keys = [], []
    phi = [atom(prev, "C"), atom(i, "N"), atom(i, "CA"), atom(i, "C")]
    psi = [atom(i, "N"), atom(i, "CA"), atom(i, "C"), atom(nxt, "N")]
    for nm, q in (("phi", phi), ("psi", psi)):
        if all(v is not None for v in q):
            quads.append(q)
            keys.append(nm)
    if not quads:
        raise KeyError(f"residue {resnum}: no complete phi/psi quadruple")
    ang = np.degrees(md.compute_dihedrals(t, np.array(quads)))
    return dict(zip(keys, ang.T)), top.residue(i).name


def stacking(t, top, resmap, resnum):
    """Heavy-atom contacts between the site's ring and Tyr200's ring, per frame."""
    idx_of = _by_number(resmap)
    i = idx_of.get(resnum)
    if i is None:
        return None
    rname = top.residue(i).name
    names = RING_ATOMS.get(rname)
    if names is None:
        return None
    ring = [a.index for a in top.residue(i).atoms if a.name in names]
    tyr = None
    for idx, (num, nm, *_) in resmap.items():
        if num == AROMATIC_NEIGHBOUR[1] and nm == AROMATIC_NEIGHBOUR[0]:
            tyr = [a.index for a in top.residue(idx).atoms
                   if a.name in ("CG", "CD1", "CD2", "CE1", "CE2", "CZ", "OH")]
    if not ring or not tyr:
        return None
    pairs = np.array([[a, b] for a in ring for b in tyr])
    d = md.compute_distances(t, pairs) * 10.0
    return (d < STACK_CUTOFF_A).sum(axis=1).astype(float)


def secondary_structure(t, top, resmap):
    """Per-residue DSSP over the window, simplified to H / E / C."""
    idx_of = _by_number(resmap)
    want = [(n, idx_of[n]) for n in range(DSSP_WINDOW[0], DSSP_WINDOW[1] + 1)
            if n in idx_of]
    if not want:
        return {}
    # DSSP needs protein only; a topology carrying DNA and water makes mdtraj's
    # implementation either fail or return codes for residues that have no backbone.
    prot = t.atom_slice(t.top.select("protein"))
    keep = {top.residue(i).resSeq for _, i in want}
    codes = md.compute_dssp(prot, simplified=True)
    out = {}
    for col, r in enumerate(prot.top.residues):
        if r.resSeq not in keep:
            continue
        num = next((n for n, i in want if top.residue(i).resSeq == r.resSeq), None)
        if num is None:
            continue
        c = codes[:, col]
        out[str(num)] = {"H": round(float((c == "H").mean()), 3),
                         "E": round(float((c == "E").mean()), 3),
                         "C": round(float((c == "C").mean()), 3)}
    return out


def analyse(traj, top_path, resmap_path, sites, out_prefix, stride=1, skip=0):
    t = md.load(str(traj), top=str(top_path), stride=stride)
    t = drop_equilibration(t, skip)
    top = t.topology
    resmap = load_map(resmap_path)

    rows, summary = {}, {}
    for s in sites:
        ang, rname = backbone_dihedrals(t, top, resmap, s)
        summary[f"resname_{s}"] = rname
        for nm, v in ang.items():
            rows[f"{nm}_{s}_deg"] = v
            mean, sd, R = circular_stats(v)
            summary[f"{nm}_{s}"] = {"mean": round(mean, 2), "sd": round(sd, 2),
                                    "R": round(R, 4)}
        st = stacking(t, top, resmap, s)
        if st is not None:
            rows[f"stack_Tyr200_{s}"] = st
            summary[f"stack_Tyr200_{s}"] = {
                "contacts_per_frame": round(float(st.mean()), 3),
                "frames_with_contact": round(float((st > 0).mean()), 3)}
        print(f"  site {s} ({rname}): phi {summary[f'phi_{s}']['mean']:+.1f} "
              f"+/- {summary[f'phi_{s}']['sd']:.1f}, "
              f"psi {summary[f'psi_{s}']['mean']:+.1f} "
              f"+/- {summary[f'psi_{s}']['sd']:.1f}"
              + (f", Tyr200 contacts {summary[f'stack_Tyr200_{s}']['contacts_per_frame']:.2f}/frame"
                 if f"stack_Tyr200_{s}" in summary else ""))

    try:
        summary["dssp"] = secondary_structure(t, top, resmap)
    except Exception as e:                                        # noqa: BLE001
        print(f"  DSSP failed: {type(e).__name__}: {e}")
        summary["dssp"] = {}

    if out_prefix:
        pd.DataFrame(rows).to_csv(f"{out_prefix}_sites.csv", index=False)
        Path(f"{out_prefix}_sites.json").write_text(json.dumps(summary, indent=2))
    return summary


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--traj", required=True)
    ap.add_argument("--top", required=True)
    ap.add_argument("--resmap", required=True)
    ap.add_argument("--sites", nargs=2, type=int, default=[133, 149])
    ap.add_argument("--stride", type=int, default=2)
    ap.add_argument("--skip", type=int, default=0)
    ap.add_argument("--out-prefix", default=None)
    a = ap.parse_args()
    print(json.dumps(analyse(a.traj, a.top, a.resmap, a.sites, a.out_prefix,
                             a.stride, a.skip), indent=2))


if __name__ == "__main__":
    main()
