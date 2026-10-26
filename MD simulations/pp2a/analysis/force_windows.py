#!/usr/bin/env python3
"""The measured force-extension curve of the bare 50 nt strand.

Runs the sfGFP window analysis (force_extension.analyse: mean restraint force
f = k(<r> - x0) with an autocorrelation-corrected error and the box-validity guard)
over every D_ss50_x<X0>A_k<K> window, and sets each point against the Marko-Siggia
tension the repository's SsSpring predicts at the same <r>.  This is the comparison the
sfGFP work never made for its single strand.
"""
from __future__ import annotations

import argparse
import json
import os
import sys
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from common import SCRATCH, dedup_mdtraj_load  # noqa: E402
import ssdna  # noqa: E402
import force_extension as fx  # noqa: E402  (sfGFP)
from spring_model import Spring, SsSpring, KT, XI_SS  # noqa: E402


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--runs", default=str(SCRATCH / "runs"))
    ap.add_argument("--systems", default=str(HERE.parent / "build" / "systems"))
    ap.add_argument("--out", default=str(HERE / "force_extension.csv"))
    ap.add_argument("--stride", type=int, default=2)
    ap.add_argument("--skip-ns", type=float, default=5.0)
    ap.add_argument("--min-ns", type=float, default=10.0)
    a = ap.parse_args()
    rows = []
    for traj in sorted(Path(a.runs).glob("D_*_x*A_k*/rep*/03_production/traj.dcd")):
        tag = traj.parents[2].name
        base = tag.rsplit("_x", 1)[0]
        sysdir = Path(a.systems) / base
        prog = traj.parent / "progress.json"
        if not prog.exists():
            continue
        ns = json.loads(prog.read_text()).get("ns_done", 0.0)
        if ns < a.min_ns:
            print(f"  skip {tag}: only {ns} ns"); continue
        try:
            with dedup_mdtraj_load(20.0):
                r = fx.analyse(traj, sysdir / "system.prmtop", sysdir / "residue_map.json", prog,
                               str(traj.parent / "analysis"), a.stride, a.skip_ns,
                               clamp_atoms=sysdir / "clamp_atoms.txt")
        except Exception as e:                                    # noqa: BLE001
            print(f"  {tag} failed: {type(e).__name__}: {e}"); continue
        # the strand's own shape in this window: does stretching it straighten the kinks?
        kink = {}
        try:
            k = ssdna.analyse(traj, sysdir / "system.prmtop", sysdir / "residue_map.json",
                              str(traj.parent / "analysis"), a.stride,
                              int(round(a.skip_ns * 1000.0 / 20.0 / a.stride)), 20.0).get("kink", {})
            kink = {"kink.axis_curvature_deg_per_nt": k.get("axis_curvature_deg_per_nt", {}).get("mean"),
                    "kink.sharp_turns_per_frame": k.get("sharp_turns_per_frame", {}).get("mean"),
                    "kink.frames_with_any_sharp_turn": k.get("frames_with_any_sharp_turn"),
                    "kink.hottest_nt": k.get("hottest_nt"), "kink.hottest_step": k.get("hottest_step"),
                    "kink.bii_fraction": k.get("bii_fraction", {}).get("mean")}
        except Exception as e:                                    # noqa: BLE001
            print(f"  {tag}: strand analysis failed: {type(e).__name__}: {e}")
        n_nt = sum(1 for row in json.loads((sysdir / "residue_map.json").read_text())
                   if row["orig_resname"].rstrip("35") in ("DA", "DT", "DG", "DC"))
        sp = SsSpring(n_nt)
        x = r["mean_r_A"] / 10.0
        s = min(x / sp.L0, 0.999)
        flat = {"system": base, "tag": tag, "replicate": traj.parents[1].name, "ns": ns, "n_nt": n_nt,
                "ms_tension_pN": sp.tension(x),
                "ms_tension_unclipped_pN": (KT / XI_SS) * (1.0 / (4.0 * (1.0 - s) ** 2) - 0.25 + s),
                "ideal_chain_pN": 3.0 * KT * x / (2.0 * XI_SS * sp.L0),
                "ms_state": sp.state(x), **kink}
        if r.get("n_strands") == 2:                      # a duplex: Zocchi's rod at the same <r>
            zs = Spring(n_nt)
            flat["zocchi_force_pN"] = round(zs.force(x), 3) if zs.is_feasible(x) else None
            flat["zocchi_regime"] = ("kinked" if x < zs.x_c else "bent") if zs.is_feasible(x) else "infeasible"
            flat["zocchi_unkinked_pN"] = round(zs.unkinked_force(x), 3) if zs.is_feasible(x) else None
        for kk, vv in r.items():
            if isinstance(vv, dict):
                for k2, v2 in vv.items():
                    if not isinstance(v2, list):
                        flat[f"{kk}.{k2}"] = v2
            else:
                flat[kk] = vv
        rows.append(flat)
        # sign: force_extension reports k(<r>-x0) > 0 when the molecule PUSHES its ends
        # apart (duplex convention); a strand that pulls gives a negative number.
        print(f"  {tag}: x0 {r['x0_A']:.0f} A, <r> {r['mean_r_A']:.1f} A, measured "
              f"{r['force_pN']:+.2f} +/- {r['force_sem_pN']:.2f} pN (push +), "
              f"Marko-Siggia pull {sp.tension(x):.2f} pN ({sp.state(x)}), box_ok {r['box_ok']}")
    if not rows:
        sys.exit("no finished windows yet")
    df = pd.DataFrame(rows).sort_values("x0_A")
    df.to_csv(a.out, index=False)
    print(f"wrote {a.out}: {len(df)} windows")


if __name__ == "__main__":
    main()
