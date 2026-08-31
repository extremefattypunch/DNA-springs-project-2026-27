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
import re
import sys
from pathlib import Path

import mdtraj as md
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import geometry  # noqa: E402

# 1 pN = 0.60221408 kJ/mol/nm, the same conversion the engine's force clamp uses.
KJ_PER_MOL_NM_TO_PN = 1.0 / 0.60221408
PURINES = {"DA", "DG"}
WC_OPEN_A = 4.0          # N1(purine)-N3(pyrimidine) beyond this: the pair has opened
# A kink breaks an INTERNAL pair.  The terminal two at each end fray in any duplex at any
# curvature, so counting them would report kinking in a perfectly straight molecule; and
# the occupancy threshold keeps a single fluctuating frame from counting as a kink.  Both
# numbers are fixed here rather than chosen after looking at the data.
TERMINAL_MARGIN = 2
OPEN_OCCUPANCY = 0.10
# A window whose target separation approaches the box's perpendicular width is measuring a
# duplex that can see its own periodic image, and its force is not the molecule's alone.
# The solvent boxes here were padded around the BENT duplex, which is 61.7 A long against
# the straight molecule's 93, so the extended windows quietly ran in a box too small for
# them -- the 27-mer's 66-88 A windows and the 40-mers' 82-128 A ones.  Flagged rather
# than silently dropped, so the curve shows where it stops being trustworthy.
BOX_FRACTION_MAX = 0.85
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
            report_ps=20.0, clamp_atoms=None, seeded_bp=None):
    """Bend profile, base-pair opening and (where there is a restraint) the force.

    Two kinds of run come through here.  A harmonic *window* records restraint_x0_A and
    restraint_k in progress.json, and the mean restraint force is the measurement.  A
    *free* run -- the seeded-kink duplexes -- has no restraint at all, and the question
    is instead whether the kink it was built with survives.  The geometry analysis is
    identical either way, which is the reason this shares one code path rather than
    getting a second module: the Watson-Crick criterion, the helical-axis smoothing and
    the minimum-image handling are all things it took several attempts to get right, and
    a parallel implementation would drift from them.
    """
    prog = json.loads(Path(progress).read_text())
    x0 = prog.get("restraint_x0_A")
    k = prog.get("restraint_k")
    restrained = x0 is not None and k is not None
    ca = prog.get("clamp_atoms")
    if ca is None and clamp_atoms is not None:
        ca = Path(clamp_atoms).read_text().strip() if Path(str(clamp_atoms)).exists() \
            else str(clamp_atoms)
    if ca is None:
        raise ValueError("no clamp atoms recorded and none supplied")
    i, j = (int(v) for v in str(ca).split(","))

    t = md.load(str(traj), top=str(top_path), stride=stride)
    skip = int(round(skip_ns * 1000.0 / report_ps / stride))
    if len(t) <= skip + 10:
        raise ValueError(f"only {len(t)} frames after a {skip_ns} ns discard")
    t = t[skip:]
    top = t.topology

    # The restraint OpenMM applied used the true, unwrapped separation
    # (setUsesPeriodicBoundaryConditions(False) on unwrapped internal coordinates), so the
    # simulation is right whatever this analysis does.  But the windows run from 50 to
    # 128 A in a box whose minimum image is only good to ~38 A, and mdtraj minimum-images
    # by default -- which is what produced the nonsensical -89 and -178 pN in the first
    # pass.  Both clamp atoms are on the same strand, hence the same molecule, so the raw
    # coordinates are contiguous and give the separation the restraint actually saw.
    r, pbc = geometry.long_distance(
        t, i, j, f"restraint pair (window x0={x0} A)" if restrained
        else "end-to-end pair (free run, no restraint)")
    # Correlation time from the lag-1 autocorrelation, so the error bar reflects the
    # number of INDEPENDENT samples rather than the number of frames.
    d = r - r.mean()
    rho = float((d[:-1] * d[1:]).sum() / max((d * d).sum(), 1e-12))
    rho = min(max(rho, 0.0), 0.99)
    n_ind = len(r) * (1.0 - rho) / (1.0 + rho)
    sem = float(r.std() / np.sqrt(max(n_ind, 1.0)))

    if restrained:
        f_kjnm = k * (r.mean() - x0) / 10.0        # (A -> nm)
        f_pN = f_kjnm * KJ_PER_MOL_NM_TO_PN
        f_sem_pN = k * sem / 10.0 * KJ_PER_MOL_NM_TO_PN
    else:
        # A free duplex exerts no force on anything, and reporting 0 would put a
        # meaningless point on the force-extension curve. None keeps the row out of it.
        f_pN = f_sem_pN = None

    resmap = load_map(resmap_path)
    strands = dna_strands(top, resmap)
    box_width = 2.0 * geometry.min_image_limit(t)
    # Judge the box on where the duplex ACTUALLY sits, not only on where we aimed.  A
    # soft restraint is not a constraint: at k = 33 kJ/mol/nm2 the duplex settles 2-7 A
    # beyond its target, so a window that passes on x0 can still spend its whole
    # trajectory close enough to its own image to matter.  Both tests must pass.
    r_mean = float(r.mean())
    # A free run has no target, so only the separation it actually adopted can be judged
    # against the box; ok_x0 is vacuously true there rather than being skipped, so the
    # two columns still mean what their names say.
    ok_x0 = bool(x0 < BOX_FRACTION_MAX * box_width) if restrained else True
    ok_r = bool(r_mean < BOX_FRACTION_MAX * box_width)
    box_ok = ok_x0 and ok_r
    if not box_ok:
        which = ("x0 = %.0f A" % x0) if not ok_x0 else ("the duplex settled at %.0f A"
                                                        % r_mean)
        print(f"  WARNING: {which} is "
              f"{max(x0 or 0.0, r_mean) / box_width:.2f} of the box's perpendicular "
              f"width ({box_width:.0f} A) -- the duplex can reach its own periodic "
              "image and this force is not its alone")
    out = {"x0_A": x0, "k_kJ_mol_nm2": k, "restrained": restrained,
           "frames": len(t), "skip_ns": skip_ns,
           "pbc": pbc, "box_width_A": round(box_width, 1),
           "x0_over_box": round(float(x0 / box_width), 3) if restrained else None,
           "r_over_box": round(float(r_mean / box_width), 3),
           "box_ok": box_ok, "box_ok_on_x0": ok_x0, "box_ok_on_measured_r": ok_r,
           "mean_r_A": round(float(r.mean()), 3), "sd_r_A": round(float(r.std()), 3),
           "sem_r_A": round(sem, 4), "lag1_rho": round(rho, 4),
           "independent_samples": round(float(n_ind), 1),
           "force_pN": round(float(f_pN), 3) if f_pN is not None else None,
           "force_sem_pN": round(float(f_sem_pN), 3) if f_sem_pN is not None else None,
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
            n_pairs = len(labels)
            internal = [(i + 1, float(f)) for i, f in enumerate(open_frac)
                        if TERMINAL_MARGIN < i + 1 <= n_pairs - TERMINAL_MARGIN
                        and f > OPEN_OCCUPANCY]
            out["wc"] = {
                "criterion_A": WC_OPEN_A,
                "pairs": n_pairs,
                "mean_N1_N3_A": round(float(wc.mean()), 3),
                "worst_pair": labels[int(open_frac.argmax())],
                "worst_pair_open_fraction": round(float(open_frac.max()), 4),
                "pairs_open_over_10pc": int((open_frac > OPEN_OCCUPANCY).sum()),
                "frames_with_any_open": round(
                    float(((wc > WC_OPEN_A).any(axis=1)).mean()), 4),
                # This is the one that decides whether a kink has formed.  Terminal
                # fraying is not a kink and is excluded by construction.
                "internal_open": [{"bp": b, "occupancy": round(o, 3)}
                                  for b, o in internal],
                "n_internal_open": len(internal),
                # the headline number for the threshold plot: how much of the time the
                # most persistent internal opening is actually open
                "max_internal_occupancy": round(
                    max((o for _, o in internal), default=0.0), 3),
                "most_open_internal_bp": (max(internal, key=lambda x: x[1])[0]
                                          if internal else None),
                "kinked": bool(internal)}
            rows["n_open_pairs"] = (wc > WC_OPEN_A).sum(axis=1).astype(float)

            # ---- did a kink we PUT there survive? ----
            # The seeded runs start with one interior pair deliberately broken and 114
            # degrees of local bend.  Left free, the duplex either keeps that -- which
            # would make the kink metastable, and Zocchi's softened branch reachable --
            # or heals it, which says the kinked state is not a minimum at zero load.
            # The distinction needs the pair's history, not its average, because a pair
            # that shuts after 2 ns and a pair that shuts after 25 ns mean different
            # things and have similar means.
            if seeded_bp is not None and 1 <= seeded_bp <= n_pairs:
                col = seeded_bp - 1
                d_seed = wc[:, col]
                openness = d_seed > WC_OPEN_A
                # first frame after which it is shut for the rest of the trajectory
                healed_at = None
                if not openness[-1]:
                    idx = np.where(openness)[0]
                    last_open = int(idx[-1]) if len(idx) else -1
                    healed_at = (last_open + 1) * report_ps * stride / 1000.0
                flank = [c for c in (col - 1, col + 1) if 0 <= c < n_pairs]
                entry = {
                    "bp": seeded_bp,
                    "label": labels[col],
                    "N1_N3_mean_A": round(float(d_seed.mean()), 3),
                    "N1_N3_first_A": round(float(d_seed[0]), 3),
                    "N1_N3_last_A": round(float(d_seed[-1]), 3),
                    "open_fraction": round(float(openness.mean()), 4),
                    "open_at_end": bool(openness[-1]),
                    "healed_after_ns": (round(healed_at, 2)
                                        if healed_at is not None else None),
                    "flanking_N1_N3_A": [round(float(wc[:, c].mean()), 3)
                                         for c in flank],
                    # a kink that moves rather than heals is a third outcome, and the
                    # most interesting one, so it is named rather than left implicit
                    "most_open_bp_now": int(open_frac.argmax()) + 1,
                    "migrated": bool(internal
                                     and max(internal, key=lambda x: x[1])[0]
                                     != seeded_bp)}
                if openness.mean() > 0.5:
                    entry["verdict"] = "persisted"
                elif openness.mean() > 0.05:
                    entry["verdict"] = "intermittent"
                elif entry["migrated"]:
                    entry["verdict"] = "migrated"
                else:
                    entry["verdict"] = "healed"
                out["seeded"] = entry
                rows["seeded_N1_N3_A"] = d_seed

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
    ap.add_argument("--free", default="K27_*",
                    help="glob for FREE bare-duplex systems with no restraint -- the "
                         "seeded kinks. Their geometry is analysed the same way; they "
                         "simply have no force to report.")
    a = ap.parse_args()

    # Harmonic windows carry the target in the directory name; free runs do not, so they
    # are matched separately rather than by loosening the window glob -- which would also
    # have swept up every protein chimera.
    targets = [(t_, True) for t_ in
               sorted(Path(a.runs).glob("*_x*A*/rep*/03_production/traj.dcd"))]
    if a.free:
        targets += [(t_, False) for t_ in
                    sorted(Path(a.runs).glob(f"{a.free}/rep*/03_production/traj.dcd"))]

    rows = []
    for traj, is_window in targets:
        tag = traj.parent.parent.parent.name
        base = tag.rsplit("_x", 1)[0] if is_window else tag
        sysdir = Path(a.systems) / base
        # K27_bp7 -> seeded pair 7. Taken from the name and cross-checked against what
        # the seeder actually wrote, so a rename cannot silently point this at the wrong
        # base pair.
        seeded = None
        if not is_window:
            m = re.search(r"_bp(\d+)$", base)
            if m:
                seeded = int(m.group(1))
                kj = Path("build/03_dna") / f"kink27_bp{seeded}.json"
                if kj.exists():
                    rec = json.loads(kj.read_text()).get("step")
                    if rec is not None and int(rec) != seeded:
                        print(f"  skip {tag}: name says bp{seeded} but the seeder "
                              f"recorded step {rec}")
                        continue
        prog = traj.parent / "progress.json"
        if not (sysdir / "system.prmtop").exists() or not prog.exists():
            print(f"  skip {tag}: no topology or progress"); continue
        ns = json.loads(prog.read_text()).get("ns_done", 0.0)
        if ns < a.min_ns:
            print(f"  skip {tag}/{traj.parent.parent.name}: only {ns} ns"); continue
        try:
            r = analyse(traj, sysdir / "system.prmtop", sysdir / "residue_map.json",
                        prog, str(traj.parent / "analysis"), a.stride, a.skip_ns,
                        clamp_atoms=sysdir / "clamp_atoms.txt", seeded_bp=seeded)
        except Exception as e:                                    # noqa: BLE001
            print(f"  {tag}/{traj.parent.parent.name} failed: "
                  f"{type(e).__name__}: {e}")
            continue
        # `duplex` is the PHYSICAL molecule, which is not the same thing as the
        # system name: D_27bp and B_27bp are the identical 27 bp duplex differing only
        # in how much water surrounds it (78.3 A across against 85.5).  The force curve
        # is a property of the molecule, so it is grouped by this, and the two boxes are
        # distinguished by box_width_A -- which makes the overlapping windows a
        # cross-check between boxes rather than two unrelated series.
        nick = "nick" in base.lower()
        npair = (r.get("wc") or {}).get("pairs")
        duplex = (f"{npair} bp{', nicked' if nick else ''}" if npair
                  else base.split("_", 1)[-1])
        flat = {"system": base, "duplex": duplex, "box_A": r.get("box_width_A"),
                "tag": tag, "replicate": traj.parent.parent.name, "ns": ns}
        for kk, vv in r.items():
            if isinstance(vv, dict):
                for k2, v2 in vv.items():
                    if not isinstance(v2, list):
                        flat[f"{kk}.{k2}"] = v2
            else:
                flat[kk] = vv
        rows.append(flat)
        if not r.get("restrained"):
            sd = r.get("seeded") or {}
            print(f"  {tag}/{traj.parent.parent.name}: FREE, <r> {r['mean_r_A']:.2f} A, "
                  f"seeded bp{sd.get('bp')} {sd.get('verdict', '?')} "
                  f"(open {100 * sd.get('open_fraction', 0):.0f}% of frames, "
                  f"N1-N3 {sd.get('N1_N3_first_A', 0):.1f} -> "
                  f"{sd.get('N1_N3_last_A', 0):.1f} A)")
            continue                      # already appended above; the print differs only
        print(f"  {tag}/{traj.parent.parent.name}: x0 {r['x0_A']:.0f} A, "
              f"<r> {r['mean_r_A']:.2f}, f = {r['force_pN']:+.2f} "
              f"+/- {r['force_sem_pN']:.2f} pN"
              + (f", max step {r['bend']['max_step_deg']:.1f} deg "
                 f"(z {r['bend']['max_step_z']:.1f})" if "bend" in r else "")
              + (f", internal pairs open {r['wc']['n_internal_open']}"
                 + (f" (bp {[x['bp'] for x in r['wc']['internal_open']]})"
                    if r['wc']['internal_open'] else "")
                 if "wc" in r else ""))
    if not rows:
        raise SystemExit("no harmonic windows and no free duplex runs with "
                         "enough production yet")
    df = pd.DataFrame(rows).sort_values(["system", "x0_A", "replicate"])
    Path(a.out).parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(a.out, index=False)
    print(f"\nwrote {a.out}: {len(df)} windows")


if __name__ == "__main__":
    main()
