#!/usr/bin/env python3
"""What the modification does at the two sites: backbone phi/psi and ring neighbours.

Reuses the sfGFP machinery (analysis/sites.py: named-atom dihedrals, circular statistics,
the protein-first resolution of a residue number that a nucleotide shares) and replaces
its two sfGFP constants -- Tyr200 as the one aromatic neighbour, a DSSP window around
residues 128-155 -- with a search: every aromatic ring within STACK_CUTOFF of the Tet2-Et
/ dihydropyridazine ring is counted, per frame, for each site.  The sfGFP lesson was
that ring stacking on a neighbour changes between the unclicked and clicked states and
masquerades as a mechanical effect, so it is measured rather than assumed absent.

Native phi/psi of Asp7 (-115, -123) and Asn553 (+118, +153) in the AlphaFold model are
in attachment_sites.json; neither is the left-handed-alpha position that strained sfGFP's
Asp133.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import mdtraj as md
import numpy as np
import pandas as pd

from common import SITES, load, load_resmap, summ  # noqa: F401
import sites as sf_sites                             # sfGFP analysis module
from sites import RING_ATOMS, STACK_CUTOFF_A, backbone_dihedrals, circular_stats

AROMATIC = {"PHE": ("CG", "CD1", "CD2", "CE1", "CE2", "CZ"),
            "TYR": ("CG", "CD1", "CD2", "CE1", "CE2", "CZ"),
            "TRP": ("CD2", "CE2", "CE3", "CZ2", "CZ3", "CH2"),
            "HIS": ("CG", "ND1", "CD2", "CE1", "NE2"), "HID": ("CG", "ND1", "CD2", "CE1", "NE2"),
            "HIE": ("CG", "ND1", "CD2", "CE1", "NE2"), "HIP": ("CG", "ND1", "CD2", "CE1", "NE2")}


def ring_neighbours(t, top, resmap, resnum):
    """Per-frame count of aromatic rings within STACK_CUTOFF of the site's own ring."""
    sf_resmap = {i: v for i, v in resmap.items()}
    idx_of = sf_sites._by_number(sf_resmap)
    ridx = idx_of.get(resnum)
    if ridx is None:
        return None, {}
    res = top.residue(ridx)
    names = RING_ATOMS.get(res.name)
    if not names:
        return None, {}
    ring = [a.index for a in res.atoms if a.name in names]
    partners = {}
    for i, (num, nm, ch) in resmap.items():
        if ch != "A" or i == ridx or nm not in AROMATIC:
            continue
        r = top.residue(i)
        atoms = [a.index for a in r.atoms if a.name in AROMATIC[nm]]
        if not atoms:
            continue
        c = t.xyz[:, atoms].mean(1)
        rc = t.xyz[:, ring].mean(1)
        d = np.linalg.norm(c - rc, axis=1) * 10.0
        if d.min() < 8.0:
            partners[f"{nm}{num}"] = d
    counts = np.zeros(t.n_frames)
    detail = {}
    for k, d in partners.items():
        hit = d < STACK_CUTOFF_A + 1.0        # ring-centroid criterion, a little looser than atom contact
        counts += hit
        detail[k] = {"min_A": round(float(d.min()), 2), "mean_A": round(float(d.mean()), 2),
                     "fraction_within_cutoff": round(float(hit.mean()), 3)}
    return counts, detail


def analyse(traj, top_path, resmap_path, out_prefix, stride=1, skip=0, sites=SITES):
    t = load(traj, top_path, stride, skip)
    top = t.topology
    resmap = load_resmap(top, resmap_path)
    rows, summary = {}, {}
    for s in sites:
        ang, rname = backbone_dihedrals(t, top, resmap, s)
        summary[f"resname_{s}"] = rname
        for nm, v in ang.items():
            rows[f"{nm}_{s}_deg"] = v
            mean, sd, R = circular_stats(v)
            summary[f"{nm}_{s}"] = {"mean": round(mean, 2), "sd": round(sd, 2), "R": round(R, 4)}
        counts, detail = ring_neighbours(t, top, resmap, s)
        if counts is not None:
            rows[f"ring_neighbours_{s}"] = counts
            summary[f"ring_neighbours_{s}"] = {"contacts_per_frame": round(float(counts.mean()), 3),
                                               "frames_with_contact": round(float((counts > 0).mean()), 3),
                                               "partners": detail}
        print(f"  site {s} ({rname}): phi {summary[f'phi_{s}']['mean']:+.1f} +/- {summary[f'phi_{s}']['sd']:.1f}, "
              f"psi {summary[f'psi_{s}']['mean']:+.1f} +/- {summary[f'psi_{s}']['sd']:.1f}"
              + (f", aromatic ring neighbours {counts.mean():.2f}/frame" if counts is not None else ""))
    pd.DataFrame(rows).to_csv(f"{out_prefix}_sites.csv", index=False)
    Path(f"{out_prefix}_sites.json").write_text(json.dumps(summary, indent=2))
    return summary


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--traj", required=True)
    ap.add_argument("--top", required=True)
    ap.add_argument("--resmap", default=None)
    ap.add_argument("--out-prefix", required=True)
    ap.add_argument("--stride", type=int, default=1)
    ap.add_argument("--skip", type=int, default=0)
    a = ap.parse_args()
    analyse(a.traj, a.top, a.resmap, a.out_prefix, a.stride, a.skip)
