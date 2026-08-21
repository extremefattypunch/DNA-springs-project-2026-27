#!/usr/bin/env python3
"""Run every analysis over every replicate and aggregate to tidy tables.

One row per (system, clamp force, replicate) in ``analysis/summary.csv``, plus the
per-frame CSVs each module writes next to the trajectory.  Everything downstream --
figures, the report -- reads these, so no plot recomputes anything.

A run directory named ``S6_clamp_7pN`` means system S6_clamp at 7 pN; the topology
lives under the base system name, since the whole point of the clamp ladder is that
six forces share one topology.
"""
from __future__ import annotations

import argparse
import json
import re
import sys
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import barrel  # noqa: E402
import chromophore  # noqa: E402
import mechanics  # noqa: E402

CLAMP_RE = re.compile(r"^(?P<base>.+?)_(?P<f>[\d.]+)pN$")


def split_tag(tag: str):
    m = CLAMP_RE.match(tag)
    if m:
        return m.group("base"), float(m.group("f"))
    return tag, None


def flat(prefix, d, out):
    for k, v in (d or {}).items():
        key = f"{prefix}{k}"
        if isinstance(v, dict):
            flat(f"{key}.", v, out)
        elif isinstance(v, (int, float, str, bool)) or v is None:
            out[key] = v
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--runs", default="data/runs")
    ap.add_argument("--systems", default="build/systems")
    ap.add_argument("--out", default="analysis")
    ap.add_argument("--stride", type=int, default=2)
    ap.add_argument("--sites", nargs=2, type=int, default=[133, 149])
    ap.add_argument("--min-ns", type=float, default=1.0,
                    help="skip replicates with less than this much production")
    ap.add_argument("--skip-ns", type=float, default=5.0,
                    help="production time to discard as settling before averaging")
    ap.add_argument("--report-ps", type=float, default=20.0,
                    help="trajectory frame spacing, used to convert --skip-ns")
    args = ap.parse_args()

    rows = []
    for traj in sorted(Path(args.runs).glob("*/rep*/03_production/traj.dcd")):
        rep_dir = traj.parent.parent
        tag = rep_dir.parent.name
        if tag.startswith("_"):
            continue
        base, clamp = split_tag(tag)
        sysdir = Path(args.systems) / base
        top, rmap = sysdir / "system.prmtop", sysdir / "residue_map.json"
        if not top.exists() or not rmap.exists():
            print(f"  skip {tag}/{rep_dir.name}: no topology at {sysdir}")
            continue
        prog = traj.parent / "progress.json"
        ns = json.loads(prog.read_text())["ns_done"] if prog.exists() else 0.0
        if ns < args.min_ns:
            print(f"  skip {tag}/{rep_dir.name}: only {ns} ns so far")
            continue
        skip = int(round(args.skip_ns * 1000.0 / args.report_ps / args.stride))
        pre = str(traj.parent / "analysis")
        print(f"\n=== {tag} {rep_dir.name}  ({ns} ns, first {args.skip_ns} ns discarded) ===")
        row = {"system": base, "tag": tag, "clamp_pN": clamp,
               "replicate": rep_dir.name, "ns": ns, "skip_ns": args.skip_ns}
        try:
            c = chromophore.analyse(traj, top, rmap, pre, args.stride, skip)
            flat("chromo.", c, row)
        except Exception as e:                              # noqa: BLE001
            print(f"  chromophore failed: {type(e).__name__}: {e}")
        try:
            m = mechanics.analyse(traj, top, rmap, args.sites, pre, args.stride,
                                  clamp, nicked="nick" in base, skip=skip)
            flat("mech.", m, row)
        except Exception as e:                              # noqa: BLE001
            print(f"  mechanics failed: {type(e).__name__}: {e}")
        try:
            b = barrel.analyse(traj, top, rmap, args.sites, pre, args.stride, skip)
            flat("barrel.", b, row)
        except Exception as e:                              # noqa: BLE001
            print(f"  barrel failed: {type(e).__name__}: {e}")
        rows.append(row)

    if not rows:
        sys.exit("no replicates with enough production yet")
    df = pd.DataFrame(rows).sort_values(["system", "clamp_pN", "replicate"])
    out = Path(args.out)
    df.to_csv(out / "summary.csv", index=False)
    print(f"\nwrote {out / 'summary.csv'}: {len(df)} replicates, {len(df.columns)} columns")
    for s, g in df.groupby("system"):
        print(f"  {s:<16} {len(g)} replicate(s), {g['ns'].sum():.0f} ns total")


if __name__ == "__main__":
    main()
