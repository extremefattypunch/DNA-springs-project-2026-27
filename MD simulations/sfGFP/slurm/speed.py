#!/usr/bin/env python3
"""Current throughput per run, refusing to quote a stale log.

Why the freshness check is the whole point
------------------------------------------
Run directories are REUSED.  A job that is cancelled and resubmitted writes into the
same ``02_equilibrate/npt.log`` and ``03_production/production.log``, and the engine
appends rather than truncating -- which is exactly right for requeue, and a trap for
anyone reading throughput out of those files.

On 2026-08-31 I cancelled a 4-way pack running on an A40, resubmitted it constrained to
an A100, then read npt.log to see what the faster card was giving.  It read 18-19 ns/day,
the same as the A40, and I nearly reported that the A100 constraint had bought nothing.
The file had last been written at 09:49:58; the A40 job was cancelled at 09:56 and the
A100 job had not yet reached NPT.  I was reading the dead job's numbers.

So this reports the age of every number it prints, and marks anything older than
``--stale-after`` seconds as STALE rather than quoting it as current.  It also names the
GPU each job actually holds, because a throughput figure means nothing without it:
gpu_requeue bills 418.25 whether it hands you an A100-SXM4-80GB at 2039 GB/s or an A40
at 696, and comparing across cards is what made the packing measurement wrong too.

    python slurm/speed.py                 # every run with a log
    python slurm/speed.py --stale-after 120
"""
from __future__ import annotations

import argparse
import glob
import os
import re
import subprocess
import time
from collections import defaultdict

STAGE = {"03_production/production.log": "prod", "02_equilibrate/npt.log": "npt",
         "02_equilibrate/nvt.log": "nvt"}


def last_speed(path):
    """(ns done, ns/day) from the last numeric line, or None."""
    try:
        last = None
        with open(path) as fh:
            for line in fh:
                if line[:1].isdigit():
                    last = line
        if not last:
            return None
        f = last.split()
        # OpenMM's StateDataReporter: step, time(ps), ..., speed(ns/day) last-but-two
        return float(f[1]) / 1000.0, float(f[-2] if len(f) > 8 else f[-1])
    except (OSError, IndexError, ValueError):
        return None


def gpu_of_running_jobs():
    """job id -> GPU model, for anything of mine in the queue."""
    out = subprocess.run(["squeue", "-h", "-u", os.environ.get("USER", ""),
                          "-o", "%i|%N|%j"], stdout=subprocess.PIPE,
                         stderr=subprocess.DEVNULL, text=True).stdout
    jobs = {}
    for line in out.strip().splitlines():
        p = line.split("|")
        if len(p) < 3:
            continue
        jid = p[0]
        tres = subprocess.run(["sacct", "-n", "-X", "-j", jid, "-o", "AllocTRES%90"],
                              stdout=subprocess.PIPE, stderr=subprocess.DEVNULL,
                              text=True).stdout
        m = re.search(r"gres/gpu:([a-zA-Z0-9_.-]+)", tres)
        jobs[jid] = (p[2], m.group(1) if m else "?", p[1])
    return jobs


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--runs", default="data/runs")
    ap.add_argument("--stale-after", type=float, default=300.0,
                    help="seconds; a log untouched for longer is reported STALE")
    a = ap.parse_args()

    now = time.time()
    rows = []
    for sub, tag in STAGE.items():
        for path in glob.glob(os.path.join(a.runs, "*", "rep*", *sub.split("/"))):
            r = last_speed(path)
            if not r:
                continue
            ns, speed = r
            age = now - os.path.getmtime(path)
            run = os.sep.join(path.split(os.sep)[-4:-2])
            rows.append((run, tag, ns, speed, age))

    # one row per run: the freshest stage wins, since that is the one still advancing
    best = {}
    for run, tag, ns, speed, age in rows:
        if run not in best or age < best[run][3]:
            best[run] = (tag, ns, speed, age)

    jobs = gpu_of_running_jobs()
    if jobs:
        print("\n  running jobs and the card each actually holds:")
        for jid, (name, card, node) in sorted(jobs.items()):
            print(f"    {jid}  {name:18s} {card:26s} {node}")

    print(f"\n  {'run':34s} {'stage':>6} {'ns':>8} {'ns/day':>8} {'log age':>9}")
    live = defaultdict(list)
    for run in sorted(best):
        tag, ns, speed, age = best[run]
        stale = age > a.stale_after
        mark = "  STALE" if stale else ""
        print(f"  {run:34s} {tag:>6} {ns:8.2f} {speed:8.0f} "
              f"{age:8.0f}s{mark}")
        if not stale:
            live[tag].append(speed)
    for tag, sp in sorted(live.items()):
        print(f"\n  {tag}: {len(sp)} live runs, {sum(sp):.0f} ns/day aggregate, "
              f"{sum(sp)/len(sp):.0f} each")
    if any(best[r][3] > a.stale_after for r in best):
        print("\n  STALE rows are logs from a cancelled or finished job in a reused run")
        print("  directory. They are not what is running now -- do not quote them.")


if __name__ == "__main__":
    main()
