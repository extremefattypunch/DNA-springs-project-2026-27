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


SCRATCH = Path(__import__("os").environ.get("PP2A_MD_SCRATCH",
                                             "/n/netscratch/hekstra_lab/Lab/ian_poon/pp2a-md"))
# Fixed system -> colour slot (vizstyle.SERIES) and display label, so a figure that drops a
# system never repaints the survivors.  P0 is the cryo-EM reproduction's apo run.
SYSTEMS = {
    "P0_apo": {"label": "apo Aα (ff19SB/OPC)", "slot": 0, "marker": "o"},
    "P1_tet": {"label": "2× Tet2-Et, unclicked", "slot": 1, "marker": "s"},
    "P2_clicked": {"label": "clicked, no DNA", "slot": 5, "marker": "D"},
    "P3_ss50": {"label": "50 nt ssDNA, closed loop", "slot": 2, "marker": "^"},
    "P4_ds50": {"label": "50 bp dsDNA, closed loop", "slot": 3, "marker": "v"},
    "P5_ss80": {"label": "80 nt ssDNA, closed loop", "slot": 4, "marker": "<"},
    "P6_ds80": {"label": "80 bp dsDNA, closed loop", "slot": 6, "marker": ">"},
}
CONTROL = "P2_clicked"
# the spring chimeras: the same strand unhybridised (an entropic coil that can only pull)
# and hybridised (Zocchi's leaf spring, which pushes), at two lengths
CHIMERAS = {"P3_ss50": {"n": 50, "topology": "ss"}, "P4_ds50": {"n": 50, "topology": "ds"},
            "P5_ss80": {"n": 80, "topology": "ss"}, "P6_ds80": {"n": 80, "topology": "ds"}}
BARE = {"D_ss50": {"n": 50, "topology": "ss"}, "D_ds50": {"n": 50, "topology": "ds"},
        "D_ss80": {"n": 80, "topology": "ss"}, "D_ds80": {"n": 80, "topology": "ds"}}


def analysis_prefix(system, replicate="rep1"):
    """Where a run's per-frame analysis CSVs live (``<prefix>_rmsd.csv`` etc.)."""
    if system == "P0_apo":
        return SCRATCH / "apo_analysis" / "analysis"
    return SCRATCH / "runs" / system / replicate / "03_production" / "analysis"


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


def frame_mask(traj, log=None):
    """Which DCD frames are genuine: both engines APPEND to traj.dcd on a checkpoint
    restart without truncating, so the frames simulated after the last checkpoint are
    written twice (the apo run has ~17 such duplicates; a preempted gpu_requeue run has
    up to ten per restart).  The StateDataReporter log is written in lockstep, so a frame
    is kept only when its step is below every step after it -- the re-simulated branch
    that the rest of the run continued from, not the abandoned one.  Returns the kept frame
    indices, or every frame (with a warning) if the log and the DCD disagree in length."""
    traj = Path(traj)
    n = len(md.formats.DCDTrajectoryFile(str(traj)))
    if log is None:
        for name in ("production.log", "04_production.log"):
            if (traj.parent / name).exists():
                log = traj.parent / name
                break
    if log is None or not Path(log).exists():
        return np.arange(n)
    steps = np.array([int(float(l.split("\t")[0])) for l in Path(log).read_text().splitlines()
                      if l and not l.startswith("#")])
    if len(steps) != n:
        print(f"  WARNING: {len(steps)} log rows vs {n} DCD frames in {traj.parent.name}; keeping every frame")
        return np.arange(n)
    # keep the LAST occurrence of each step: the branch the run actually continued from
    suffix_min = np.minimum.accumulate(np.concatenate([[np.iinfo(np.int64).max], steps[::-1][:-1]]))[::-1]
    keep = np.flatnonzero(steps < suffix_min)
    if len(keep) < n:
        print(f"  dropped {n - len(keep)} duplicated restart frames of {n}")
    return keep


def load(traj, top_path, stride=1, skip=0, report_ps=None, log=None):
    """Load the genuine frames (restart duplicates removed), stride, and drop ``skip``
    strided frames of settling (never below 10 frames).  t.time is set to real
    picoseconds when report_ps is given -- a DCD carries no clock and mdtraj's time
    would otherwise be the raw frame index."""
    keep_all = frame_mask(traj, log)
    keep = keep_all[::stride]
    chunks, offset, want = [], 0, set(keep.tolist())
    for ch in md.iterload(str(traj), top=str(top_path), chunk=250):
        idx = [i for i in range(ch.n_frames) if offset + i in want]
        if idx:
            chunks.append(ch[idx])
        offset += ch.n_frames
    t = chunks[0] if len(chunks) == 1 else md.join(chunks, check_topology=False)
    # genuine frames are consecutive in step, so time is the rank of the kept frame
    t.time = np.arange(0, len(keep_all), stride, dtype=float)[:t.n_frames] * (float(report_ps) if report_ps else 1.0)
    if skip and t.n_frames - skip >= 10:
        t = t[skip:]
    return t


class dedup_mdtraj_load:
    """Context manager: route every mdtraj.load(traj, top=, stride=) made by the sfGFP
    modules (mechanics, contacts, force_extension) through the loader above, so they too
    skip the duplicated restart frames and carry real times."""

    def __init__(self, report_ps=None, log=None):
        self.report_ps, self.log = report_ps, log

    def __enter__(self):
        self.orig = md.load

        def patched(filename_or_filenames, **kw):
            if isinstance(filename_or_filenames, (str, Path)) and str(filename_or_filenames).endswith(".dcd") \
                    and set(kw) <= {"top", "stride"}:
                return load(filename_or_filenames, kw["top"], kw.get("stride", 1), 0, self.report_ps, self.log)
            return self.orig(filename_or_filenames, **kw)
        md.load = patched
        return self

    def __exit__(self, *exc):
        md.load = self.orig


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
    b0, b1, b2 = p0 - p1, p2 - p1, p3 - p2      # b0 points back, per the standard formula
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
