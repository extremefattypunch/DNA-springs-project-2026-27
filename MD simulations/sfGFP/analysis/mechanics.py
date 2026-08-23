#!/usr/bin/env python3
"""Mechanics: what the spring actually does, and how hard it pulls.

Three measurements, and one inference that ties them to Zocchi's model.

1. **Protein extension** -- the Cbeta-Cbeta distance between the two attachment
   residues.  This is the deformation coordinate: Zocchi's "enzyme spring" is
   x = (this distance) minus its unstressed value, and his calibration puts the
   protein's linear stiffness near 100 kT/nm^2 with a softening transition around
   3 A.  Predicted displacements are 0.1-3 A, so this is measured against matched
   controls with error bars, never as a single-trajectory difference.

2. **Spring extension** -- the 5'P-to-5'P distance across the duplex.  This is
   exactly the variable x that Zocchi's f(x) takes, so the force the spring delivers
   can be evaluated from the simulation's own geometry rather than assumed.

3. **Where the duplex gives** -- per-base-pair bend angle from local helical-axis
   vectors, plus Watson-Crick hydrogen-bond distances.  A kink shows up as a bend
   spike; melting shows up as base pairs opening.  Whether the duplex relieves its
   strain this way is the softening transition the whole model rests on.

The inference: the force-clamp ladder gives extension as a function of *known*
applied force on the same topology.  Reading a chimera's measured extension against
that curve yields the tension its spring is really delivering, which can then be
compared with f(x) evaluated at the measured x.  Two independent routes to the same
number is the point -- one structural, one from the analytic model.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import mdtraj as md
import numpy as np
import pandas as pd

import sys
sys.path.insert(0, str(Path(__file__).resolve().parent))
import geometry  # noqa: E402
from spring_model import KT, TAU_C_CONTINUOUS, TAU_C_NICKED, Spring  # noqa: E402

WC = {("DA", "DT"): [("N1", "N3"), ("N6", "O4")],
      ("DT", "DA"): [("N3", "N1"), ("O4", "N6")],
      ("DG", "DC"): [("N1", "N3"), ("N2", "O2"), ("O6", "N4")],
      ("DC", "DG"): [("N3", "N1"), ("O2", "N2"), ("N4", "O6")]}


def load_map(path):
    rows = json.loads(Path(path).read_text())
    return {r["index"]: (r["orig_resnum"], r["orig_resname"],
                         r.get("orig_chain", "?")) for r in rows}


def anchor_atoms(top, resmap, sites):
    out = []
    for s in sites:
        idx = next(i for i, v in resmap.items() if v[0] == s)
        res = top.residue(idx)
        a = next((x for x in res.atoms if x.name == "CB"), None)
        if a is None:
            a = next(x for x in res.atoms if x.name == "CA")
        out.append(a.index)
    return out


def dna_residues(top, resmap):
    """[(chain_id, [residues 5'->3'])] per DNA strand, from the build-time chain map.

    An Amber prmtop stores no chain field, so mdtraj presents the entire system as a
    single chain and top.residue(i).chain is useless for telling the two strands
    apart.  The chain recorded by build_system.py at assembly time is the only
    reliable source.
    """
    strands = {}
    for i, (num, name, chain) in resmap.items():
        if name.rstrip("35") in ("DA", "DT", "DG", "DC"):
            strands.setdefault(chain, []).append((num, top.residue(i)))
    for k in strands:
        strands[k].sort(key=lambda p: p[0])
        strands[k] = [r for _, r in strands[k]]
    return sorted(strands.items())


def drop_equilibration(t, skip):
    """Drop the first ``skip`` frames (already strided) as settling time.

    The box is still contracting and the solvent-exposed side chains are still
    finding their positions for the first few ns of production; averaging that into
    an equilibrium observable biases it and inflates its spread.  Never drops so much
    that fewer than 10 frames remain.
    """
    if skip and t.n_frames - skip >= 10:
        return t[skip:]
    return t


def analyse(traj, top_path, resmap_path, sites, out_prefix, stride=1,
            clamp_pN=None, nicked=False, skip=0):
    t = md.load(str(traj), top=str(top_path), stride=stride)
    t = drop_equilibration(t, skip)
    top = t.topology
    resmap = load_map(resmap_path)
    rows = {"time_ps": t.time if t.time is not None else np.arange(t.n_frames)}
    summary = {"n_frames": int(t.n_frames), "clamp_pN": clamp_pN}

    # --- 1. protein extension ---
    i, j = anchor_atoms(top, resmap, sites)
    # Both of the long-range distances in this module go through geometry.long_distance,
    # which picks the periodic convention from the topology and then checks its own
    # reasoning.  mdtraj minimum-images by default, and that is correct for a hydrogen
    # bond and wrong for anything longer than half the cell's smallest width -- 38-48 A
    # in these truncated octahedra.  This pair is ~32 A so both conventions agree, but it
    # is routed through the same code so that the assertion covers it too: if a future
    # box shrinks or a future site pair is farther apart, this fails loudly.
    mols = geometry.molecule_ids(top)
    ext, ext_pbc = geometry.long_distance(t, i, j, "anchor CB-CB", mols)
    summary["pbc"] = {"anchor_CB_CB": ext_pbc}
    rows["anchor_CB_CB_A"] = ext
    summary["anchor_CB_CB_A"] = {"mean": round(float(ext.mean()), 3),
                                 "sd": round(float(ext.std()), 3),
                                 "sem": round(float(ext.std() / np.sqrt(len(ext))), 4)}
    print(f"  anchor Cb-Cb: {ext.mean():.3f} +/- {ext.std():.3f} A")

    # --- 2. spring extension and the force it implies ---
    strands = dna_residues(top, resmap)
    if len(strands) == 2:
        # The spring's end-to-end distance is measured between the two nucleotide atoms
        # the linkers are actually bonded to, found from the bond graph rather than
        # assumed.  This has to be derived and not hard-coded, because the two builds
        # attach differently: the corrected chimera puts both linkers on ONE strand
        # (the 5' P of nucleotide 1 and the 3' O3' of nucleotide N), while the earlier
        # two-arm build put one on each strand's 5' P.  Reading the graph gives the
        # right pair either way, and both are the same physical quantity -- the distance
        # across the duplex between the points where the protein holds it -- which is
        # the x that enters Zocchi's f(x).
        linkers = {"DNL", "DNP"}
        nuc = {r.index for _, res in strands for r in res}
        attach = []
        for b in top.bonds:
            ra, rb = b[0].residue, b[1].residue
            for near, far in ((b[0], b[1]), (b[1], b[0])):
                if far.residue.name in linkers and near.residue.index in nuc:
                    attach.append((near.index, near.residue.index, near.name))
        attach = sorted(set(attach), key=lambda x: x[1])
        if len(attach) != 2:
            print(f"  WARNING: found {len(attach)} linker-to-DNA bonds, expected 2; "
                  "falling back to the two strands' 5' phosphorus")
            attach = []
            for _, res in strands:
                first = res[0]
                p = next((a for a in first.atoms if a.name == "P"), None)
                a = p if p is not None else next(x for x in first.atoms
                                                 if x.name == "C5'")
                attach.append((a.index, first.index, a.name))
        one_strand = len({top.atom(i).residue.chain.index for i, _, _ in attach}) == 1
        summary["attachment"] = {
            "atoms": [f"{top.atom(i).residue.name}{top.atom(i).residue.resSeq}.{nm}"
                      for i, _, nm in attach],
            "same_strand": one_strand,
            "topology": "closed covalent loop" if one_strand else "two separate arms"}
        print(f"  attachment: {summary['attachment']['atoms']} -> "
              f"{summary['attachment']['topology']}")
        # THIS is the one the convention actually mattered for.  The spring's ends are
        # 60-75 A apart in a box whose minimum image is only good to ~42 A, so the default
        # gave a value that flipped between periodic images from frame to frame: for
        # S3_spring27 it read 53.8 A against a true 70.0, and for the closed loop 45.4
        # against 71.2.  Every force this project quoted from Zocchi's f(x) was evaluated
        # at one of those.  Both ends are on the same covalently connected molecule, so
        # the raw coordinates are contiguous and are the honest answer.
        d, spring_pbc = geometry.long_distance(
            t, attach[0][0], attach[1][0], "spring end-to-end", mols)
        summary["pbc"]["spring_end_to_end"] = spring_pbc
        rows["spring_end_to_end_A"] = d
        n_bp = len(strands[0][1])
        tau = TAU_C_NICKED if nicked else TAU_C_CONTINUOUS
        sp = Spring(n_bp, tau_c=tau)
        f = np.array([sp.force(x / 10.0) if x / 10.0 < sp.L else np.nan for x in d])
        e = np.array([sp.energy(x / 10.0) / KT if x / 10.0 < sp.L else np.nan for x in d])
        rows["spring_5P_5P_A"] = d          # kept under the old name for the
        #                                      earlier runs' CSVs to stay readable
        rows["spring_force_pN"] = f
        rows["spring_energy_kT"] = e
        summary["spring"] = {
            "n_bp": n_bp, "tau_c_pN_nm": tau, "nicked": nicked,
            "contour_nm": round(sp.L, 3),
            "x_mean_A": round(float(d.mean()), 3), "x_sd_A": round(float(d.std()), 3),
            "force_pN_mean": round(float(np.nanmean(f)), 3),
            "force_pN_sd": round(float(np.nanstd(f)), 3),
            "energy_kT_mean": round(float(np.nanmean(e)), 3),
            "regime": "kinked" if d.mean() / 10.0 < sp.x_c else "smoothly bent"}
        print(f"  spring 5'P-5'P: {d.mean():.2f} +/- {d.std():.2f} A -> "
              f"f = {np.nanmean(f):.2f} +/- {np.nanstd(f):.2f} pN, "
              f"E = {np.nanmean(e):.2f} kT ({summary['spring']['regime']})")

        # --- 3. where the duplex gives ---
        a_res, b_res = strands[0][1], strands[1][1]
        n = min(len(a_res), len(b_res))
        c1 = []
        for k in range(n):
            ai = next(x.index for x in a_res[k].atoms if x.name == "C1'")
            bi = next(x.index for x in b_res[n - 1 - k].atoms if x.name == "C1'")
            c1.append((ai, bi))
        pos = t.xyz[:, np.array(c1).ravel(), :].reshape(t.n_frames, n, 2, 3)
        centres = pos.mean(axis=2) * 10.0

        # Base-pair centres spiral around the helical axis at ~1.9 A, so the direction
        # between *consecutive* centres swings by 20-30 deg per step from the helix
        # itself, with nothing to do with bending -- measured naively this duplex
        # reported 27 deg/bp and 687 deg of total bend.  Averaging the centres over one
        # full turn cancels the spiral and leaves the axis.
        w = int(min(11, max(3, n // 3)))
        if w % 2 == 0:
            w -= 1
        kern = np.ones(w) / w
        smooth = np.stack([np.apply_along_axis(
            lambda v: np.convolve(v, kern, mode="valid"), 1, centres[:, :, k])
            for k in range(3)], axis=2)
        tang = np.diff(smooth, axis=1)
        tang /= np.linalg.norm(tang, axis=2, keepdims=True)
        cosang = np.clip((tang[:, :-1] * tang[:, 1:]).sum(axis=2), -1, 1)
        bend = np.degrees(np.arccos(cosang))          # per axis step
        total = np.degrees(np.arccos(np.clip(
            (tang[:, 0] * tang[:, -1]).sum(axis=1), -1, 1)))
        offset = 1 + w // 2                           # smoothed index -> bp number
        rows["bend_max_deg"] = bend.max(axis=1)
        rows["bend_argmax_bp"] = bend.argmax(axis=1) + offset + 1.0
        rows["bend_total_deg"] = total
        summary["bend"] = {
            "axis_smoothing_window_bp": w,
            "per_step_mean_deg": round(float(bend.mean()), 3),
            "per_step_max_mean_deg": round(float(bend.max(axis=1).mean()), 3),
            "hottest_bp": int(np.bincount(bend.argmax(axis=1) + offset + 1,
                                          minlength=n + 2).argmax()),
            "total_bend_deg_mean": round(float(total.mean()), 2),
            "total_bend_deg_sd": round(float(total.std()), 2),
            "span_bp": [1 + w // 2, n - w // 2],
            "note": ("total_bend covers only the interior span the smoothing window "
                     "leaves, so it is smaller than the arc the duplex was built with "
                     "by roughly the fraction of base pairs trimmed")}
        print(f"  axis bend: {bend.mean():.2f} deg/step (window {w} bp), "
              f"total {total.mean():.1f} +/- {total.std():.1f} deg, "
              f"sharpest near bp {summary['bend']['hottest_bp']}")

        # per-step bend profile along the duplex, averaged over frames: this is where
        # a kink would announce itself as a localised spike
        prof = pd.DataFrame({
            "bp": np.arange(bend.shape[1]) + offset + 1,
            "bend_deg_mean": bend.mean(axis=0).round(4),
            "bend_deg_sd": bend.std(axis=0).round(4)})
        prof.to_csv(f"{out_prefix}_bend_profile.csv", index=False)

        # base-pair opening
        pairs, labels = [], []
        for k in range(n):
            ra, rb = a_res[k], b_res[n - 1 - k]
            key = (ra.name.rstrip("35"), rb.name.rstrip("35"))
            for na, nb in WC.get(key, []):
                try:
                    pairs.append([next(x.index for x in ra.atoms if x.name == na),
                                  next(x.index for x in rb.atoms if x.name == nb)])
                    labels.append(k + 1)
                except StopIteration:
                    pass
        if pairs:
            dd = md.compute_distances(t, np.array(pairs)) * 10.0
            lab = np.array(labels)
            open_frac = []
            for k in range(1, n + 1):
                sel = dd[:, lab == k]
                if sel.size:
                    open_frac.append(float((sel.mean(axis=1) > 4.0).mean()))
            rows["bp_open_count"] = np.array(
                [sum(dd[fr, lab == k].mean() > 4.0 for k in range(1, n + 1))
                 for fr in range(t.n_frames)], float)
            summary["base_pair_opening"] = {
                "criterion": "mean WC donor-acceptor distance > 4.0 A",
                "fraction_open_per_bp": [round(v, 4) for v in open_frac],
                "most_open_bp": int(np.argmax(open_frac) + 1),
                "mean_open_count": round(float(rows["bp_open_count"].mean()), 3)}
            print(f"  base pairs open (>4 A): "
                  f"{rows['bp_open_count'].mean():.2f} of {n} on average; "
                  f"most labile bp {summary['base_pair_opening']['most_open_bp']}")

    df = pd.DataFrame(rows)
    df.to_csv(f"{out_prefix}_mechanics.csv", index=False)
    Path(f"{out_prefix}_mechanics.json").write_text(json.dumps(summary, indent=2))
    print(f"  wrote {out_prefix}_mechanics.csv/.json")
    return summary


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--traj", required=True)
    ap.add_argument("--top", required=True)
    ap.add_argument("--resmap", required=True)
    ap.add_argument("--sites", nargs=2, type=int, default=[133, 149])
    ap.add_argument("--out-prefix", required=True)
    ap.add_argument("--stride", type=int, default=1)
    ap.add_argument("--skip", type=int, default=0,
                    help="frames to drop from the start, after striding")
    ap.add_argument("--clamp-pN", type=float, default=None)
    ap.add_argument("--nicked", action="store_true")
    a = ap.parse_args()
    analyse(a.traj, a.top, a.resmap, a.sites, a.out_prefix, a.stride,
            a.clamp_pN, a.nicked)
