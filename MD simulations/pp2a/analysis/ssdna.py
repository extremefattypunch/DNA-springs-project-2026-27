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

from common import angle_deg, load, load_resmap, summ, SFGFP_ANALYSIS  # noqa: F401
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


def analyse(traj, top_path, resmap_path, out_prefix, stride=1, skip=0):
    t = load(traj, top_path, stride, skip)
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

    # --- backbone trace and sharp turns ---
    c4 = np.array([next(a.index for a in r.atoms if a.name == "C4'") for r in nts])
    xyz = t.xyz[:, c4] * 10.0
    v = xyz[:, 1:] - xyz[:, :-1]
    bend = angle_deg(v[:, :-1], v[:, 1:])            # frames x (n-2), at nt 2..n-1
    rows["trace_bend_mean_deg"] = bend.mean(1)
    rows["trace_bend_max_deg"] = bend.max(1)
    rows["sharp_turns"] = (bend > SHARP_TURN_DEG).sum(1)
    prof = bend.mean(0)
    summary["trace_bend"] = {
        "definition": f"angle between successive C4'->C4' vectors; sharp turn > {SHARP_TURN_DEG:.0f} deg",
        "per_nt_mean_deg": summ(bend.mean(1)),
        "sharp_turns_per_frame": summ(rows["sharp_turns"]),
        "hottest_nt": int(np.argmax(prof) + 2),
        "profile_mean_deg": [round(float(p), 2) for p in prof]}
    pd.DataFrame({"nt": np.arange(2, n), "bend_deg_mean": prof.round(3),
                  "bend_deg_sd": bend.std(0).round(3),
                  "sharp_turn_fraction": (bend > SHARP_TURN_DEG).mean(0).round(4)}
                 ).to_csv(f"{out_prefix}_ssdna_profile.csv", index=False)

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
    summary["stacking"] = {"criterion": f"adjacent ring centroids < {STACK_A} A",
                           "fraction": summ(rows["stacked_fraction"]),
                           "per_step_fraction": [round(float(f), 3) for f in stacked.mean(0)]}

    # --- contact with the protein ---
    prot_idx = top.select("protein and not element H")
    if prot_idx.size:
        per_nt = np.zeros((t.n_frames, n))
        for k, r in enumerate(nts):
            heavy = [a.index for a in r.atoms if a.element.symbol != "H"]
            pairs = np.array([(h, p) for h in heavy for p in prot_idx[::3]])  # every 3rd protein heavy atom: cheap, adequate at 4.5 A
            d = md.compute_distances(t, pairs, periodic=False) * 10.0
            per_nt[:, k] = (d.min(1) < CONTACT_A)
        rows["nt_in_contact"] = per_nt.sum(1)
        summary["protein_contact"] = {
            "criterion": f"any nucleotide heavy atom < {CONTACT_A} A from a protein heavy atom",
            "nt_in_contact_per_frame": summ(rows["nt_in_contact"]),
            "fraction_of_strand": round(float(per_nt.mean()), 4),
            "per_nt_fraction": [round(float(f), 3) for f in per_nt.mean(0)]}

    print(f"  ssDNA {n} nt: end-to-end {x.mean():.1f} +/- {x.std():.1f} A ({summary['marko_siggia']['state_at_mean_x']}; "
          f"MS tension {clipped.mean():.2f} pN clipped, {raw.mean():.2f} unclipped), Rg {rg.mean():.1f} A, "
          f"sharp turns {rows['sharp_turns'].mean():.1f}/frame, stacked {rows['stacked_fraction'].mean():.2f}"
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
