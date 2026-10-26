#!/usr/bin/env python3
"""Shared pieces for the PP2A analyses: loading, residue maps, the HEAT grid, block SEM.

Residue numbering throughout is UniProt P30153 (Met1 = 1).  The two attachment sites are
Asp7 and Asn553 (= construct D6X/N552X; see build/01_protein/prepare_pp2a.py).

Two kinds of run are analysed with the same code:
  * the spring systems built here (P1_tet, P2_clicked, P3_ss50), which carry a
    build/systems/<sys>/residue_map.json like every sfGFP system;
  * the cryo-EM reproduction's apo run (Aa_rep1, ff19SB/OPC, 100 ps frames), which has
    no residue map -- its prmtop is plain protein 1-589 and the map is the identity.
`load_resmap` hides that difference; everything else keys on (chain, orig_resnum).
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import mdtraj as md
import numpy as np

HERE = Path(__file__).resolve().parent
SFGFP_ANALYSIS = HERE.parents[1] / "sfGFP" / "analysis"
if str(SFGFP_ANALYSIS) not in sys.path:
    sys.path.insert(0, str(SFGFP_ANALYSIS))

SITES = (7, 553)
LITERATURE_PAIR = (30, 578)         # Yilmaz/Bahar 2025 "N29-F577", Met1-less -> Asn30/Phe578
SALT_BRIDGE = ((113, "ARG"), (191, "GLU"))   # Day 2025 Fig 3C, HEAT3-HEAT5 backbone distance
HEAT_START, HEAT_LEN, N_HEAT, N_RES = 8, 39, 15, 589
N_TAIL = (1, 13)                    # Day 2025: the N-terminus analysed separately
CORE = (14, 589)
REGIONS = {"N (HEAT 1-3)": (1, 3), "central (HEAT 4-12)": (4, 12), "C (HEAT 13-15)": (13, 15)}
PROTEIN_SKIP = {"WAT", "HOH", "Na+", "Cl-", "MG", "K+", "DA", "DT", "DG", "DC",
                "DA5", "DT5", "DG5", "DC5", "DA3", "DT3", "DG3", "DC3", "DNL", "DNP", "DNH"}


def heat_ranges():
    """Nominal 39-residue grid (cryoEM/.../25_geometry.py): HEAT k = 8+39(k-1) .. 8+39k-1."""
    return {k: (HEAT_START + HEAT_LEN * (k - 1), min(HEAT_START + HEAT_LEN * k - 1, N_RES))
            for k in range(1, N_HEAT + 1)}


def heat_of(resnum):
    if resnum < HEAT_START:
        return 0
    return min(N_HEAT, (resnum - HEAT_START) // HEAT_LEN + 1)


def load_resmap(top, path):
    """prmtop residue index -> (orig_resnum, orig_resname, orig_chain).

    Falls back to the identity map (chain A, resSeq, name) when there is no
    residue_map.json, which is the apo run's case.
    """
    if path and Path(path).exists():
        rows = json.loads(Path(path).read_text())
        return {r["index"]: (r["orig_resnum"], r["orig_resname"], r.get("orig_chain", "A"))
                for r in rows}
    # mdtraj numbers prmtop residues from 0 (resSeq == index), so the identity map is
    # index + 1: Met1 is residue 0 in the topology.
    return {r.index: (r.index + 1, r.name, "A") for r in top.residues
            if r.name not in PROTEIN_SKIP}


def protein_residues(top, resmap):
    """[(orig_resnum, prmtop residue)] for the protein, in sequence order."""
    out = [(num, top.residue(i)) for i, (num, name, chain) in resmap.items()
           if chain == "A" and name not in PROTEIN_SKIP]
    out.sort(key=lambda p: p[0])
    return out


def ca_by_resnum(top, resmap):
    return {num: next(a.index for a in res.atoms if a.name == "CA")
            for num, res in protein_residues(top, resmap)
            if any(a.name == "CA" for a in res.atoms)}


def atom_by_resnum(top, resmap, resnum, name):
    for num, res in protein_residues(top, resmap):
        if num == resnum:
            a = next((x for x in res.atoms if x.name == name), None)
            return None if a is None else a.index
    return None


def load(traj, top_path, stride=1, skip=0):
    """Load, stride, and drop ``skip`` strided frames of settling (never below 10 frames)."""
    t = md.load(str(traj), top=str(top_path), stride=stride)
    if skip and t.n_frames - skip >= 10:
        t = t[skip:]
    return t


def block_sem(x, n_blocks=10):
    """SEM from block averages -- the honest error of a single replicate's mean."""
    x = np.asarray(x, dtype=float)
    if len(x) < 2 * n_blocks:
        return float(x.std() / np.sqrt(len(x)))
    blocks = np.array_split(x, n_blocks)
    means = np.array([b.mean() for b in blocks])
    return float(means.std(ddof=1) / np.sqrt(n_blocks))


def angle_deg(u, v):
    u = u / np.linalg.norm(u, axis=-1, keepdims=True)
    v = v / np.linalg.norm(v, axis=-1, keepdims=True)
    return np.degrees(np.arccos(np.clip((u * v).sum(-1), -1.0, 1.0)))


def dihedral_deg(p0, p1, p2, p3):
    b0, b1, b2 = p1 - p0, p2 - p1, p3 - p2
    b1n = b1 / np.linalg.norm(b1, axis=-1, keepdims=True)
    v = b0 - (b0 * b1n).sum(-1, keepdims=True) * b1n
    w = b2 - (b2 * b1n).sum(-1, keepdims=True) * b1n
    x = (v * w).sum(-1)
    y = (np.cross(b1n, v) * w).sum(-1)
    return np.degrees(np.arctan2(y, x))


def summ(x):
    x = np.asarray(x, dtype=float)
    return {"mean": round(float(x.mean()), 3), "sd": round(float(x.std()), 3),
            "sem_block": round(block_sem(x), 4), "min": round(float(x.min()), 3),
            "max": round(float(x.max()), 3)}
