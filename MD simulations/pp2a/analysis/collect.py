#!/usr/bin/env python3
"""Run every analysis over every run and aggregate to analysis/summary.csv.

One row per run.  The spring systems are discovered under $PP2A_MD_SCRATCH/runs the way
the sfGFP collector does it (<TAG>/rep*/03_production/traj.dcd, topology under
build/systems/<base>); the apo baseline is the cryo-EM reproduction's Aa_rep1 (flat
layout, 100 ps frames, ff19SB/OPC) and is added as system ``P0_apo`` with its own frame
spacing.  Every module writes its per-frame CSV next to the trajectory under
``analysis_*``; figures.py reads those and this summary, never the trajectories.

Modules: rmsd (rigidity), heat (solenoid bend/twist/aperture), mechanics (sfGFP: the
Cb-Cb anchor distance and the spring end-to-end with the periodic-image fix), ssdna
(the strand), contacts (sfGFP: protein-DNA contacts), tet_sites (phi/psi + ring
neighbours).  Benchmarks (``*_bench``) and restraint windows (``*_x*A_k*``, handled by
force_extension) are skipped.
"""
from __future__ import annotations

import argparse
import json
import os
import re
import sys
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from common import SITES, SFGFP_ANALYSIS, dedup_mdtraj_load  # noqa: E402,F401
import heat, rmsd, ssdna, tet_sites       # noqa: E402
import contacts, mechanics                 # noqa: E402  (sfGFP)

SCRATCH = Path(os.environ.get("PP2A_MD_SCRATCH", "/n/netscratch/hekstra_lab/Lab/ian_poon/pp2a-md"))
APO = {"system": "P0_apo",
       "traj": Path(os.environ.get("PP2A_DATA", "/n/netscratch/hekstra_lab/Lab/ian_poon/pp2a_day2025"))
       / "md/runs/Aa_rep1/traj.dcd",
       "top": Path(os.environ.get("PP2A_DATA", "/n/netscratch/hekstra_lab/Lab/ian_poon/pp2a_day2025"))
       / "md/systems/Aa/system.prmtop",
       "report_ps": 100.0, "note": "ff19SB/OPC, Day et al. protocol, cryo-EM reproduction"}
WINDOW_RE = re.compile(r"_x[\d.]+A_k[\d.]+$")


def flat(prefix, d, out):
    for k, v in (d or {}).items():
        key = f"{prefix}{k}"
        if isinstance(v, dict):
            flat(f"{key}.", v, out)
        elif isinstance(v, (int, float, str, bool)) or v is None:
            out[key] = v
    return out


def run_all(row, traj, top, resmap, pre, stride, skip, sites, report_ps=20.0):
    def step(name, fn):
        try:
            flat(f"{name}.", fn(), row)
        except Exception as e:                                   # noqa: BLE001
            print(f"  {name} failed: {type(e).__name__}: {e}")
    step("rmsd", lambda: rmsd.analyse(traj, top, resmap, pre, stride, skip, sites, report_ps))
    step("heat", lambda: heat.analyse(traj, top, resmap, pre, stride, skip, report_ps))
    step("site", lambda: tet_sites.analyse(traj, top, resmap, pre, stride, skip, sites, report_ps))
    if resmap is not None:
        step("mech", lambda: mechanics.analyse(traj, top, resmap, list(sites), pre, stride, None, False, skip))
        step("ssdna", lambda: ssdna.analyse(traj, top, resmap, pre, stride, skip, report_ps))

        def ct():
            c = contacts.analyse(traj, top, resmap, pre, stride, skip)
            if not c.get("has_dna"):
                return {}
            out = {k: v for k, v in c.items() if not isinstance(v, list)}
            out["residues_over_half"] = ";".join(c["residues_over_half_the_time"])
            return out
        step("contact", ct)
    else:
        # the apo run has no residue map: give mechanics an identity map on disk
        ident = Path(pre + "_identity_resmap.json")
        if not ident.exists():
            import mdtraj as md
            from common import PROTEIN_SKIP
            t = md.load_prmtop(str(top))
            ident.write_text(json.dumps([{"index": r.index, "orig_chain": "A", "orig_resnum": r.index + 1,
                                          "orig_resname": r.name, "prmtop_resname": r.name}
                                         for r in t.residues if r.name not in PROTEIN_SKIP]))
        step("mech", lambda: mechanics.analyse(traj, top, ident, list(sites), pre, stride, None, False, skip))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--runs", default=str(SCRATCH / "runs"))
    ap.add_argument("--systems", default=str(HERE.parent / "build" / "systems"))
    ap.add_argument("--out", default=str(HERE))
    ap.add_argument("--stride", type=int, default=2)
    ap.add_argument("--min-ns", type=float, default=1.0)
    ap.add_argument("--skip-ns", type=float, default=5.0)
    ap.add_argument("--report-ps", type=float, default=20.0)
    ap.add_argument("--only", nargs="+", default=None)
    ap.add_argument("--no-apo", action="store_true")
    args = ap.parse_args()

    jobs = []
    for traj in sorted(Path(args.runs).glob("*/rep*/03_production/traj.dcd")):
        tag = traj.parents[2].name
        if tag.endswith("_bench") or WINDOW_RE.search(tag):
            continue
        sysdir = Path(args.systems) / tag
        if not (sysdir / "system.prmtop").exists():
            print(f"  skip {tag}: no topology")
            continue
        prog = traj.parent / "progress.json"
        ns = json.loads(prog.read_text())["ns_done"] if prog.exists() else 0.0
        jobs.append({"system": tag, "replicate": traj.parents[1].name, "ns": ns, "traj": traj,
                     "top": sysdir / "system.prmtop", "resmap": sysdir / "residue_map.json",
                     "report_ps": args.report_ps, "pre": str(traj.parent / "analysis")})
    if not args.no_apo and APO["traj"].exists():
        prog = APO["traj"].parent / "progress.json"
        ns = json.loads(prog.read_text())["ns_done"] if prog.exists() else 0.0
        apo_out = SCRATCH / "apo_analysis"
        apo_out.mkdir(parents=True, exist_ok=True)
        jobs.append({"system": APO["system"], "replicate": "rep1", "ns": ns, "traj": APO["traj"],
                     "top": APO["top"], "resmap": None, "report_ps": APO["report_ps"],
                     "pre": str(apo_out / "analysis"), "note": APO["note"]})

    rows = []
    for j in jobs:
        if args.only and j["system"] not in args.only:
            continue
        if j["ns"] < args.min_ns:
            print(f"  skip {j['system']}: only {j['ns']} ns")
            continue
        skip = int(round(args.skip_ns * 1000.0 / j["report_ps"] / args.stride))
        print(f"\n=== {j['system']} {j['replicate']}  ({j['ns']} ns, first {args.skip_ns} ns discarded, "
              f"stride {args.stride} x {j['report_ps']} ps) ===")
        row = {"system": j["system"], "replicate": j["replicate"], "ns": j["ns"],
               "skip_ns": args.skip_ns, "stride": args.stride, "report_ps": j["report_ps"],
               "sites": "/".join(map(str, SITES)), "note": j.get("note", "")}
        with dedup_mdtraj_load(j["report_ps"]):      # the sfGFP modules load through md.load
            run_all(row, j["traj"], j["top"], j["resmap"], j["pre"], args.stride, skip, SITES, j["report_ps"])
        rows.append(row)
    if not rows:
        sys.exit("nothing to collect yet")
    df = pd.DataFrame(rows)
    out = Path(args.out) / "summary.csv"
    if args.only and out.exists():
        old = pd.read_csv(out)
        df = pd.concat([old[~old.system.isin(args.only)], df], ignore_index=True, sort=False)
    df = df.sort_values(["system", "replicate"])
    df.to_csv(out, index=False)
    print(f"\nwrote {out}: {len(df)} runs, {len(df.columns)} columns")


if __name__ == "__main__":
    main()
