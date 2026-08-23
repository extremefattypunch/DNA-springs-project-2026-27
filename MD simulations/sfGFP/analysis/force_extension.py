#!/usr/bin/env python3
"""Force-extension of the bare DNA springs, measured from the force field.

Why this exists
---------------
Doeke asked, of the forces quoted in the report, whether they come from a model of the
DNA's elasticity or from the simulation itself.  They came from the model: Zocchi's f(x)
evaluated at the end-to-end distance the chimera settled at.  That is a fair thing to do
but it is not a measurement, and the honest answer is to make the measurement.

So: take the duplex on its own, no protein, solvated; hold its two ends at a fixed
separation x with a harmonic restraint of known stiffness; and let it equilibrate.  At
equilibrium the restraint's mean pull balances whatever the molecule is pushing with, so

    f_DNA(x) = k (<r> - x0)

with no elasticity model anywhere in the derivation.  Scanning x0 traces out f(x)
directly, and it can then be laid over Zocchi's curve to see whether his numbers are
the right ones for our construct.

The precision of a mean force taken from equilibrium fluctuations is sqrt(k kT / N) for
N independent samples, however the restraint is applied, so the error bars here are set
by sampling and are reported rather than glossed.

The same runs answer the kinking question (point 8), which is why the bend profile and
base-pair integrity are measured here too.  The tight windows compress the duplex well
past where it sits in the chimera, so if a kink is going to nucleate in this force field
this is where it should appear.  What counts as a kink is stated explicitly: a local
bend far above the smooth background, or a Watson-Crick pair whose N1-N3 hydrogen bond
has opened -- the type II kink of Lankas et al. and of Curuksu et al. breaks exactly one
central pair, so a per-pair measure is the right one.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import mdtraj as md
import numpy as np
import pandas as pd

# 1 pN = 0.60221408 kJ/mol/nm, the same conversion the engine's force clamp uses.
KJ_PER_MOL_NM_TO_PN = 1.0 / 0.60221408
PURINES = {"DA", "DG"}
WC_OPEN_A = 4.0          # N1(purine)-N3(pyrimidine) beyond this: the pair has opened
TURN_BP = 10.5           # one helical turn, for smoothing the axis


def load_map(path):
    rows = json.loads(Path(path).read_text())
    return {r["index"]: (r["orig_resnum"], r["orig_resname"], r.get("orig_chain"))
            for r in rows}


def dna_strands(top, resmap):
    """[(chain, [residues in 5'->3' order])] for the nucleic acid, from the map.

    The chain has to come from the residue map and not from the topology: a prmtop
    carries no chain field, so mdtraj reads the whole solute as one chain and any
    per-strand analysis silently does nothing.
    """
    by = {}
    for idx, (num, name, ch) in resmap.items():
        base = name.rstrip("35")
        if base not in ("DA", "DT", "DG", "DC"):
            continue
        by.setdefault(ch, []).append((num, idx))
    return [(ch, [top.residue(i) for _, i in sorted(v)]) for ch, v in sorted(by.items())]


def axis_bend(t, strands):
    """Bend per base-pair step, in degrees, from base-pair centres.

    Base-pair centres spiral about the helical axis by ~1.9 A, so the angle between raw
    consecutive centre-to-centre vectors is dominated by that spiral -- measured that
    way an ordinary B-form duplex reads 27 deg per step and several hundred degrees in
    total, which says nothing about bending.  Smoothing the centres over one full turn
    removes the spiral and leaves the axis.
    """
    (ca, ra), (cb, rb) = strands[0], strands[1]
    n = min(len(ra), len(rb))
    pairs = []
    for i in range(n):
        a = next((x.index for x in ra[i].atoms if x.name == "C1'"), None)
        b = next((x.index for x in rb[n - 1 - i].atoms if x.name == "C1'"), None)
        if a is not None and b is not None:
            pairs.append((a, b))
    if len(pairs) < 12:
        return None, None
    xyz = t.xyz  # (frames, atoms, 3), nm
    ctr = 0.5 * (xyz[:, [p[0] for p in pairs], :] + xyz[:, [p[1] for p in pairs], :])
    w = int(round(TURN_BP))
    ker = np.ones(w) / w
    sm = np.stack([np.apply_along_axis(lambda v: np.convolve(v, ker, mode="valid"),
                                      1, ctr[:, :, k]) for k in range(3)], axis=-1)
    tang = np.diff(sm, axis=1)
    tang /= np.linalg.norm(tang, axis=-1, keepdims=True)
    cosang = np.clip((tang[:, :-1, :] * tang[:, 1:, :]).sum(axis=-1), -1.0, 1.0)
    return np.degrees(np.arccos(cosang)), w


def wc_opening(t, strands):
    """N1-N3 distance for every Watson-Crick pair, in angstrom."""
    (ca, ra), (cb, rb) = strands[0], strands[1]
    n = min(len(ra), len(rb))
    pairs, labels = [], []
    for i in range(n):
        r1, r2 = ra[i], rb[n - 1 - i]
        n1 = "N1" if r1.name.rstrip("35") in PURINES else "N3"
        n2 = "N1" if r2.name.rstrip("35") in PURINES else "N3"
        a = next((x.index for x in r1.atoms if x.name == n1), None)
        b = next((x.index for x in r2.atoms if x.name == n2), None)
        if a is not None and b is not None:
            pairs.append([a, b])
            labels.append(f"{r1.name.rstrip('35')}{i + 1}")
    if not pairs:
        return None, None
    return md.compute_distances(t, np.array(pairs)) * 10.0, labels


def analyse(traj, top_path, resmap_path, progress, out_prefix, stride=1, skip_ns=10.0,
            report_ps=20.0):
    prog = json.loads(Path(progress).read_text())
    x0 = prog.get("restraint_x0_A")
    k = prog.get("restraint_k")
    if x0 is None or k is None:
        raise ValueError("this run has no harmonic restraint recorded")
    i, j = (int(v) for v in prog["clamp_atoms"].split(","))

    t = md.load(str(traj), top=str(top_path), stride=stride)
    skip = int(round(skip_ns * 1000.0 / report_ps / stride))
    if len(t) <= skip + 10:
        raise ValueError(f"only {len(t)} frames after a {skip_ns} ns discard")
    t = t[skip:]
    top = t.topology

    r = md.compute_distances(t, np.array([[i, j]]))[:, 0] * 10.0     # angstrom
    # Correlation time from the lag-1 autocorrelation, so the error bar reflects the
    # number of INDEPENDENT samples rather than the number of frames.
    d = r - r.mean()
    rho = float((d[:-1] * d[1:]).sum() / max((d * d).sum(), 1e-12))
    rho = min(max(rho, 0.0), 0.99)
    n_ind = len(r) * (1.0 - rho) / (1.0 + rho)
    sem = float(r.std() / np.sqrt(max(n_ind, 1.0)))

    f_kjnm = k * (r.mean() - x0) / 10.0            # (A -> nm)
    f_pN = f_kjnm * KJ_PER_MOL_NM_TO_PN
    f_sem_pN = k * sem / 10.0 * KJ_PER_MOL_NM_TO_PN

    resmap = load_map(resmap_path)
    strands = dna_strands(top, resmap)
    out = {"x0_A": x0, "k_kJ_mol_nm2": k, "frames": len(t), "skip_ns": skip_ns,
           "mean_r_A": round(float(r.mean()), 3), "sd_r_A": round(float(r.std()), 3),
           "sem_r_A": round(sem, 4), "lag1_rho": round(rho, 4),
           "independent_samples": round(float(n_ind), 1),
           "force_pN": round(float(f_pN), 3), "force_sem_pN": round(float(f_sem_pN), 3),
           "n_strands": len(strands)}
    rows = {"r_A": r}

    if len(strands) == 2:
        bend, w = axis_bend(t, strands)
        if bend is not None:
            per_step = bend.mean(axis=0)
            out["bend"] = {
                "smoothing_bp": w,
                "mean_per_step_deg": round(float(per_step.mean()), 3),
                "max_step_deg": round(float(per_step.max()), 3),
                "max_step_index": int(per_step.argmax()) + 1,
                # A kink is a LOCAL spike, so what matters is how far the worst step
                # stands above the smooth background, in units of the spread across
                # steps -- not its absolute value, which grows with overall curvature.
                "max_step_z": round(float((per_step.max() - per_step.mean())
                                          / max(per_step.std(), 1e-9)), 2),
                "per_step_deg": [round(float(v), 3) for v in per_step]}
            rows["max_bend_deg"] = bend.max(axis=1)

        wc, labels = wc_opening(t, strands)
        if wc is not None:
            open_frac = (wc > WC_OPEN_A).mean(axis=0)
            out["wc"] = {
                "criterion_A": WC_OPEN_A,
                "pairs": len(labels),
                "mean_N1_N3_A": round(float(wc.mean()), 3),
                "worst_pair": labels[int(open_frac.argmax())],
                "worst_pair_open_fraction": round(float(open_frac.max()), 4),
                "pairs_open_over_10pc": int((open_frac > 0.10).sum()),
                "frames_with_any_open": round(
                    float(((wc > WC_OPEN_A).any(axis=1)).mean()), 4)}
            rows["n_open_pairs"] = (wc > WC_OPEN_A).sum(axis=1).astype(float)

    if out_prefix:
        pd.DataFrame(rows).to_csv(f"{out_prefix}_fx.csv", index=False)
        Path(f"{out_prefix}_fx.json").write_text(json.dumps(out, indent=2))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--runs", default="data/runs")
    ap.add_argument("--systems", default="build/systems")
    ap.add_argument("--out", default="analysis/force_extension.csv")
    ap.add_argument("--stride", type=int, default=2)
    ap.add_argument("--skip-ns", type=float, default=10.0)
    ap.add_argument("--min-ns", type=float, default=15.0)
    a = ap.parse_args()

    rows = []
    # The trailing _k<stiffness> is optional, so the glob must not require the name to
    # end at the window: a soft-restraint run is D_27bp_x50A_k33.
    for traj in sorted(Path(a.runs).glob("*_x*A*/rep*/03_production/traj.dcd")):
        tag = traj.parent.parent.parent.name
        base = tag.rsplit("_x", 1)[0]
        sysdir = Path(a.systems) / base
        prog = traj.parent / "progress.json"
        if not (sysdir / "system.prmtop").exists() or not prog.exists():
            print(f"  skip {tag}: no topology or progress"); continue
        ns = json.loads(prog.read_text()).get("ns_done", 0.0)
        if ns < a.min_ns:
            print(f"  skip {tag}/{traj.parent.parent.name}: only {ns} ns"); continue
        try:
            r = analyse(traj, sysdir / "system.prmtop", sysdir / "residue_map.json",
                        prog, str(traj.parent / "analysis"), a.stride, a.skip_ns)
        except Exception as e:                                    # noqa: BLE001
            print(f"  {tag}/{traj.parent.parent.name} failed: "
                  f"{type(e).__name__}: {e}")
            continue
        flat = {"system": base, "tag": tag, "replicate": traj.parent.parent.name,
                "ns": ns}
        for kk, vv in r.items():
            if isinstance(vv, dict):
                for k2, v2 in vv.items():
                    if not isinstance(v2, list):
                        flat[f"{kk}.{k2}"] = v2
            else:
                flat[kk] = vv
        rows.append(flat)
        print(f"  {tag}/{traj.parent.parent.name}: x0 {r['x0_A']:.0f} A, "
              f"<r> {r['mean_r_A']:.2f}, f = {r['force_pN']:+.2f} "
              f"+/- {r['force_sem_pN']:.2f} pN"
              + (f", max step {r['bend']['max_step_deg']:.1f} deg "
                 f"(z {r['bend']['max_step_z']:.1f})" if "bend" in r else "")
              + (f", pairs open {r['wc']['pairs_open_over_10pc']}"
                 if "wc" in r else ""))
    if not rows:
        raise SystemExit("no force-extension windows with enough production yet")
    df = pd.DataFrame(rows).sort_values(["system", "x0_A", "replicate"])
    Path(a.out).parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(a.out, index=False)
    print(f"\nwrote {a.out}: {len(df)} windows")


if __name__ == "__main__":
    main()
