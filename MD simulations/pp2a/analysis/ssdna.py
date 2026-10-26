#!/usr/bin/env python3
"""What the single strand does: size, shape, where it bends, and whether it is stretched.

The sfGFP analysis gates every DNA shape measurement on there being two strands, so a
single strand got nothing beyond its end-to-end distance.  For a 50-mer that is the
spring, that is not enough.  Measured here, per frame:

  * end-to-end of the strand (the linker-bonded atoms, found from the bond graph exactly
    as mechanics.py does) and its radius of gyration;
  * the Marko-Siggia tension at that end-to-end distance, both as the repository's
    SsSpring returns it (clipped to zero below the chain's rms size) and unclipped, so
    the "slack / matched / stretched" verdict is not an artefact of the clip;
  * a backbone trace through the C4' atoms, its local bend angle per nucleotide
    (angle between successive C4'->C4' vectors) and the count of sharp turns -- the
    single-strand analogue of the duplex kink;
  * the fraction of adjacent bases still stacked (ring-centroid distance < 4.5 A), since
    a stacked run behaves as a stiff segment and an unstacked one as the flexible chain
    the 1.5 nm persistence length assumes;
  * the fraction of nucleotides in contact with the protein (heavy atoms < 4.5 A), the
    Zocchi section 4.5.3 caveat made measurable.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import mdtraj as md
import numpy as np
import pandas as pd

from common import angle_deg, load, load_resmap, protein_residues, summ, SFGFP_ANALYSIS  # noqa: F401
import geometry                                  # sfGFP analysis, on sys.path via common
from spring_model import SsSpring, KT, XI_SS

DNA_NAMES = ("DA", "DT", "DG", "DC")
LINKERS = {"DNL", "DNP"}
RING = {"DA": ("N9", "C8", "N7", "C5", "C6", "N1", "C2", "N3", "C4"),
        "DG": ("N9", "C8", "N7", "C5", "C6", "N1", "C2", "N3", "C4"),
        "DT": ("N1", "C2", "N3", "C4", "C5", "C6"),
        "DC": ("N1", "C2", "N3", "C4", "C5", "C6")}
SHARP_TURN_DEG = 90.0
STACK_A = 4.5
CONTACT_A = 4.5
STRAIGHT_HELIX = Path(__file__).resolve().parents[1] / "build" / "03_dna" / "straight_50nt.pdb"


def axis_window(n):
    """One helical turn, as the duplex analysis smooths (mechanics.py): min(11, max(3, n//3)), odd."""
    w = int(min(11, max(3, n // 3)))
    return w - (w % 2 == 0)


def axis_curvature(xyz, w):
    """Bend per step of a C4' axis smoothed over w nucleotides: (frames, n-w-1) degrees."""
    n = xyz.shape[1]
    sm = np.stack([xyz[:, i:i + w].mean(1) for i in range(n - w + 1)], axis=1)
    sv = sm[:, 1:] - sm[:, :-1]
    return angle_deg(sv[:, :-1], sv[:, 1:])
# backbone torsions of nucleotide i (IUPAC): atom, residue offset
TORSIONS = {
    "alpha": (("O3'", -1), ("P", 0), ("O5'", 0), ("C5'", 0)),
    "beta": (("P", 0), ("O5'", 0), ("C5'", 0), ("C4'", 0)),
    "gamma": (("O5'", 0), ("C5'", 0), ("C4'", 0), ("C3'", 0)),
    "delta": (("C5'", 0), ("C4'", 0), ("C3'", 0), ("O3'", 0)),
    "epsilon": (("C4'", 0), ("C3'", 0), ("O3'", 0), ("P", 1)),
    "zeta": (("C3'", 0), ("O3'", 0), ("P", 1), ("O5'", 1)),
}
CHI = {"DA": ("O4'", "C1'", "N9", "C4"), "DG": ("O4'", "C1'", "N9", "C4"),
       "DT": ("O4'", "C1'", "N1", "C2"), "DC": ("O4'", "C1'", "N1", "C2")}


def circ_mean_deg(a, axis=0):
    """Circular mean in degrees, ignoring NaN (undefined torsions at the strand's ends)."""
    r = np.radians(a)
    return np.degrees(np.arctan2(np.nanmean(np.sin(r), axis), np.nanmean(np.cos(r), axis)))


def torsions(t, nts):
    """Per-frame backbone torsions alpha..zeta and chi, (n_frames, n_nt) each, NaN where an atom
    is missing (alpha of a 5'-terminal nucleotide without a phosphate, epsilon/zeta of the 3' end)."""
    atom = [{a.name: a.index for a in r.atoms} for r in nts]
    out = {}
    for name, spec in list(TORSIONS.items()) + [("chi", None)]:
        quads, cols = [], []
        for i, r in enumerate(nts):
            if name == "chi":
                names = CHI[r.name.rstrip("35")]
                q = [atom[i].get(nm) for nm in names]
            else:
                q = []
                for nm, off in spec:
                    j = i + off
                    q.append(atom[j].get(nm) if 0 <= j < len(nts) else None)
            if all(v is not None for v in q):
                quads.append(q)
                cols.append(i)
        vals = np.full((t.n_frames, len(nts)), np.nan)
        if quads:
            vals[:, cols] = np.degrees(md.compute_dihedrals(t, np.array(quads), periodic=False))
        out[name] = vals
    return out


def strand(top, resmap):
    nts = [(num, top.residue(i)) for i, (num, name, chain) in resmap.items()
           if name.rstrip("35") in DNA_NAMES]
    nts.sort(key=lambda p: p[0])
    return [r for _, r in nts]


def end_atoms(top, resmap, nts):
    """The two DNA atoms bonded to a linker residue; else O5'/P of nt 1 and O3' of nt N."""
    linker_res = {i for i, (num, name, chain) in resmap.items() if name in LINKERS}
    ends = []
    for a, b in top.bonds:
        for x, y in ((a, b), (b, a)):
            if x.residue.index in linker_res and y.residue.name.rstrip("35") in DNA_NAMES:
                ends.append(y.index)
    if len(ends) == 2:
        return sorted(ends), "linker-bonded atoms (5'P of nt 1, 3'O3' of nt N)"
    first = next((a.index for a in nts[0].atoms if a.name in ("P", "O5'")), None)
    last = next(a.index for a in nts[-1].atoms if a.name == "O3'")
    return [first, last], "O5'/P of nt 1 and O3' of nt N (no linkers)"


def analyse(traj, top_path, resmap_path, out_prefix, stride=1, skip=0, report_ps=None):
    t = load(traj, top_path, stride, skip, report_ps)
    top = t.topology
    resmap = load_resmap(top, resmap_path)
    nts = strand(top, resmap)
    if not nts:
        return {"has_dna": False}
    n = len(nts)
    rows = {"time_ps": t.time if t.time is not None else np.arange(t.n_frames)}
    summary = {"has_dna": True, "n_nt": n}

    (i, j), how = end_atoms(top, resmap, nts)
    mols = geometry.molecule_ids(top)
    x, pbc = geometry.long_distance(t, i, j, "ssDNA end-to-end", mols)
    rows["end_to_end_A"] = x
    summary["end_to_end"] = {"atoms": how, "A": summ(x), "pbc": pbc}

    dna_idx = np.array([a.index for r in nts for a in r.atoms if a.element.symbol != "H"])
    rg = md.compute_rg(t.atom_slice(dna_idx)) * 10.0
    rows["rg_A"] = rg
    summary["rg_A"] = summ(rg)

    sp = SsSpring(n)
    clipped = np.array([sp.tension(v / 10.0) for v in x])
    s = np.clip(x / 10.0 / sp.L0, 0, 0.999)
    raw = (KT / XI_SS) * (1.0 / (4.0 * (1.0 - s) ** 2) - 0.25 + s)
    rows["ms_tension_pN"] = clipped
    rows["ms_tension_unclipped_pN"] = raw
    rows["ideal_chain_pN"] = 3.0 * KT * (x / 10.0) / (2.0 * XI_SS * sp.L0)
    summary["marko_siggia"] = {
        "model": "Marko-Siggia WLC, xi 1.5 nm, 0.59 nm/nt; MODELLED, not measured",
        "contour_nm": round(sp.L0, 2), "natural_rms_end_to_end_nm": round(sp.rms_end_to_end, 3),
        "state_at_mean_x": sp.state(float(x.mean()) / 10.0),
        "tension_pN_clipped": summ(clipped), "tension_pN_unclipped": summ(raw),
        "ideal_chain_pN": summ(rows["ideal_chain_pN"]),
        "fraction_of_frames_stretched": round(float((x / 10.0 > sp.rms_end_to_end).mean()), 4)}

    # --- does the DNA kink: bend along the strand, torsions, and where ---
    # Two bend measures.  (i) the raw angle between successive C4'->C4' vectors, whose
    # sharp (> 90 deg) turns are the single-strand analogue of a duplex kink; (ii) the
    # curvature of a smoothed axis (C4' positions averaged over one helical turn before
    # differencing, the duplex analysis's construction in mechanics.py), in degrees per
    # nucleotide.  Analogous to the duplex number but not numerically comparable with it:
    # a WLC of persistence length xi bends thermally by ~sqrt(2a/(xi w)) per step, so the
    # same construction on an ideal straight B-form strand (the baseline recorded below)
    # already gives a few degrees per nucleotide from the helix's own spiral.
    seq = [r.name.rstrip("35")[1] for r in nts]
    c4 = np.array([next(a.index for a in r.atoms if a.name == "C4'") for r in nts])
    xyz = t.xyz[:, c4] * 10.0
    v = xyz[:, 1:] - xyz[:, :-1]
    bend = angle_deg(v[:, :-1], v[:, 1:])            # frames x (n-2), at nt 2..n-1
    w = axis_window(n)
    curv = axis_curvature(xyz, w)                     # frames x (n-w-1), at nt 1+w//2+1 ...
    curv_nt = np.arange(2 + w // 2, 2 + w // 2 + curv.shape[1])
    baseline = None
    if STRAIGHT_HELIX.exists():                        # the same measure on an ideal B-form strand
        try:
            ideal = md.load(str(STRAIGHT_HELIX))
            ch = next(c for c in ideal.topology.chains)
            c4i = [next(a.index for a in r.atoms if a.name == "C4'") for r in ch.residues][:n]
            baseline = float(axis_curvature(ideal.xyz[:, c4i] * 10.0, w).mean())
        except Exception:                              # noqa: BLE001
            baseline = None
    tor = torsions(t, nts)
    ez = (tor["epsilon"] - tor["zeta"] + 180.0) % 360.0 - 180.0     # BI ~ -90, BII ~ +90
    bii = ez > 0
    sharp = bend > SHARP_TURN_DEG
    rows["trace_bend_mean_deg"] = bend.mean(1)
    rows["trace_bend_max_deg"] = bend.max(1)
    rows["axis_curvature_mean_deg_per_nt"] = curv.mean(1)
    rows["axis_curvature_max_deg_per_nt"] = curv.max(1)
    rows["sharp_turns"] = sharp.sum(1)
    rows["bii_fraction"] = np.nanmean(bii[:, ~np.isnan(ez).all(0)], axis=1)
    prof = bend.mean(0)
    steps = [f"{seq[i - 1]}p{seq[i]}" for i in range(1, n)]          # dinucleotide step i-1 -> i
    hot = int(np.argmax(sharp.mean(0)))                               # index into nt 2..n-1
    summary["kink"] = {
        "definition": (f"sharp turn = angle between successive C4'->C4' vectors > {SHARP_TURN_DEG:.0f} deg; "
                       f"curvature = bend of a {w}-nt smoothed C4' axis per nucleotide step"),
        "trace_bend_per_nt_mean_deg": summ(bend.mean(1)),
        "axis_curvature_deg_per_nt": summ(curv.mean(1)),
        "axis_curvature_baseline_straight_helix_deg_per_nt": None if baseline is None else round(baseline, 2),
        "axis_window_nt": int(w),
        "axis_curvature_max_deg_per_nt": summ(curv.max(1)),
        "sharp_turns_per_frame": summ(rows["sharp_turns"]),
        "frames_with_any_sharp_turn": round(float((rows["sharp_turns"] > 0).mean()), 4),
        "hottest_nt": int(hot + 2),
        "hottest_nt_sharp_turn_fraction": round(float(sharp.mean(0)[hot]), 4),
        "hottest_step": f"{steps[hot]}-{steps[hot + 1]}",
        "sharp_turns_at_unstacked_steps": None,     # filled below once stacking is known
        "bii_fraction": summ(rows["bii_fraction"]),
        "torsion_means_deg": {k: {"mean": round(float(circ_mean_deg(v[:, ~np.isnan(v).all(0)].ravel())), 1)}
                              for k, v in tor.items()},
        "profile_mean_deg": [round(float(p), 2) for p in prof]}
    summary["trace_bend"] = {"hottest_nt": summary["kink"]["hottest_nt"],
                             "sharp_turns_per_frame": summary["kink"]["sharp_turns_per_frame"]}
    for k, v in tor.items():
        rows[f"{k}_mean_deg"] = circ_mean_deg(v, axis=1)
    prof_df = pd.DataFrame({"nt": np.arange(1, n + 1), "base": seq,
                            "step_in": ["-"] + steps})
    prof_df["bend_deg_mean"] = np.r_[np.nan, prof, np.nan]
    prof_df["bend_deg_sd"] = np.r_[np.nan, bend.std(0), np.nan]
    prof_df["sharp_turn_fraction"] = np.r_[np.nan, sharp.mean(0), np.nan]
    cv = np.full(n, np.nan); cv[curv_nt - 1] = curv.mean(0)
    prof_df["axis_curvature_deg_per_nt"] = cv
    prof_df["bii_fraction"] = np.where(np.isnan(ez).all(0), np.nan, bii.mean(0))
    for k, v in tor.items():
        prof_df[f"{k}_deg"] = [None if np.isnan(v[:, i]).all() else round(float(circ_mean_deg(v[:, i])), 1)
                               for i in range(n)]

    # --- base stacking ---
    cents = []
    for r in nts:
        names = RING[r.name.rstrip("35")]
        idx = [a.index for a in r.atoms if a.name in names]
        cents.append(t.xyz[:, idx].mean(1) * 10.0)
    cents = np.stack(cents, axis=1)                   # frames x n x 3
    dd = np.linalg.norm(cents[:, 1:] - cents[:, :-1], axis=2)
    stacked = dd < STACK_A
    rows["stacked_fraction"] = stacked.mean(1)
    prof_df["unstacked_fraction_step_in"] = np.r_[np.nan, 1.0 - stacked.mean(0)]
    # do the sharp turns happen where the bases have come apart?  sharp[:, i] is the turn at
    # nucleotide i+2 (1-based); its incoming step is stacked[:, i] and its outgoing stacked[:, i+1]
    both = sharp & (~stacked[:, :-1] | ~stacked[:, 1:])
    summary["kink"]["sharp_turns_at_unstacked_steps"] = round(float(both.sum() / max(sharp.sum(), 1)), 4)
    prof_df.round(4).to_csv(f"{out_prefix}_ssdna_profile.csv", index=False)
    summary["stacking"] = {"criterion": f"adjacent ring centroids < {STACK_A} A",
                           "fraction": summ(rows["stacked_fraction"]),
                           "per_step_fraction": [round(float(f), 3) for f in stacked.mean(0)]}

    # --- contact with the protein ---
    # the residue map's protein (keeps TDP, drops the DNL/DNP linkers that mdtraj calls protein)
    prot_idx = np.array([a.index for _, res in protein_residues(top, resmap)
                         for a in res.atoms if a.element.symbol != "H"], dtype=int)
    if prot_idx.size:
        owner = np.full(top.n_atoms, -1, dtype=int)
        for k, r in enumerate(nts):
            for a in r.atoms:
                if a.element.symbol != "H":
                    owner[a.index] = k
        dna_heavy = np.where(owner >= 0)[0]
        per_nt = np.zeros((t.n_frames, n))
        for fi in range(t.n_frames):
            hit = md.compute_neighbors(t[fi], CONTACT_A / 10.0, prot_idx, haystack_indices=dna_heavy,
                                       periodic=False)[0]
            if len(hit):
                per_nt[fi, np.unique(owner[hit])] = 1.0
        rows["nt_in_contact"] = per_nt.sum(1)
        summary["protein_contact"] = {
            "criterion": f"any nucleotide heavy atom < {CONTACT_A} A from a protein heavy atom (Tet/TDP included, DNA linkers excluded)",
            "nt_in_contact_per_frame": summ(rows["nt_in_contact"]),
            "fraction_of_strand": round(float(per_nt.mean()), 4),
            "per_nt_fraction": [round(float(f), 3) for f in per_nt.mean(0)]}

    print(f"  ssDNA {n} nt: end-to-end {x.mean():.1f} +/- {x.std():.1f} A ({summary['marko_siggia']['state_at_mean_x']}; "
          f"MS tension {clipped.mean():.2f} pN clipped, {raw.mean():.2f} unclipped), Rg {rg.mean():.1f} A, "
          f"sharp turns {rows['sharp_turns'].mean():.1f}/frame (hottest nt {summary['kink']['hottest_nt']}, {summary['kink']['hottest_step']}), "
          f"curvature {curv.mean():.1f} deg/nt, BII {rows['bii_fraction'].mean():.2f}, stacked {rows['stacked_fraction'].mean():.2f}"
          + (f", {rows['nt_in_contact'].mean():.1f} nt on the protein" if 'nt_in_contact' in rows else ""))
    pd.DataFrame(rows).to_csv(f"{out_prefix}_ssdna.csv", index=False)
    Path(f"{out_prefix}_ssdna.json").write_text(json.dumps(summary, indent=2))
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
