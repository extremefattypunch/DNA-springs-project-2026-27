#!/usr/bin/env python3
"""What a job will cost the lab's fairshare, before it is submitted.

Why this exists
---------------
On 2026-08-23 I submitted 120 GPU jobs and spent **194,084 billing-hours** of
hekstra_lab's fairshare in a day, which exhausted the allocation and stopped everyone's
work.  Nothing in ``sbatch`` output, ``squeue`` or ``scontrol`` tells you a job is
expensive; the cost is invisible until it is gone.  So this prints it first.

What the numbers actually are, measured on this cluster rather than assumed
--------------------------------------------------------------------------
``sacct`` reports ``billing=842`` for the jobs I ran.  That decomposes exactly as

    Gres/gpu 836.5 x 1  = 836.50    99.34 %
    CPU        0.9  x 4 =   3.60     0.43 %
    Mem        0.06 x 32=   1.92     0.23 %
                          ------
                          842.02   -> sacct's 842, so billing is the SUM of the
                                      weighted TRES (this cluster sets no MAX_TRES)

Three consequences, and they are not the intuitive ones:

1. **The GPU is the entire bill.**  Right-sizing 32 GB -> 4 GB and 4 cores -> 1 saves
   0.5 %.  Worth doing because it is free, but it is not the lever.  Requesting less
   *VRAM* saves nothing at all: ``Gres/gpu=1`` is charged for a whole device however
   little of it you touch, and this cluster publishes no per-GPU-type billing weight.

2. **Billing accrues on ELAPSED time, not on the time you request.**  My jobs asked for
   2,880 GPU-hours and used 230.5; the bill came to 194,084 billing-hours, which is
   842 x 230.5 exactly.  The unused 92 % cost nothing.  An accurate ``--time`` still
   matters a great deal -- it is what lets the backfill scheduler slot the job into a
   gap, so it starts sooner -- but it is a scheduling lever, not a cost lever.

3. **The partition is the lever.**  Identical work, identical hardware:

       gpu_test          0        free; 2 jobs at a time, 12 h, a 3g.20gb MIG slice
       gpu_requeue     418.25     half of `gpu`, and it can land an H100 or H200
       gpu             836.5      what I used
       seas_gpu       1727.8      2x `gpu`
       gpu_h200       2648.8      6.3x `gpu_requeue` for the same H200

   ``gpu_h200`` charges 2,648.8 for hardware ``gpu_requeue`` will give you for 418.25.
   The only thing the premium buys is not being preempted -- and ``openmm_engine.py``
   checkpoints every 10 reporter intervals and resumes by appending, so a preemption
   costs a few hundred picoseconds of wall-clock and nothing else.

Cost per nanosecond is the figure that actually matters, and this script computes it
from the throughput *observed in our own run logs* rather than from a guess.  A GPU type
we have never run on is reported as unmeasured instead of being given an invented
speed-up.

    python slurm/cost.py                          # price the partitions, cheapest first
    python slurm/cost.py --ns 50 --reps 3 --sys L3_loop27
    python slurm/cost.py --plan slurm/planned.txt # price a whole planned batch
"""
from __future__ import annotations

import argparse
import glob
import os
import re
import subprocess
import sys
from collections import defaultdict

# Partitions that bill nothing.  Kept as a set rather than inferred from a zero weight so
# that a config change shows up as a surprise rather than being silently absorbed.
FREE = {"gpu_test"}


def partitions():
    """Live TRESBillingWeights per partition -> {name: {"gpu":w, "cpu":w, "mem":w, ...}}."""
    out = subprocess.run(["scontrol", "show", "partition"],
                         stdout=subprocess.PIPE, stderr=subprocess.DEVNULL,
                         text=True).stdout
    parts = {}
    for block in out.split("\n\n"):
        nm = re.search(r"PartitionName=(\S+)", block)
        w = re.search(r"TRESBillingWeights=(\S+)", block)
        if not nm:
            continue
        d = {"gpu": 0.0, "cpu": 0.0, "mem": 0.0, "maxtime": None, "gputypes": []}
        mt = re.search(r"MaxTime=(\S+)", block)
        if mt:
            d["maxtime"] = mt.group(1)
        if w:
            for tok in w.group(1).split(","):
                k, _, v = tok.partition("=")
                v = v.rstrip("G")
                try:
                    fv = float(v)
                except ValueError:
                    continue
                kl = k.lower()
                if kl.startswith("gres/gpu"):
                    d["gpu"] = fv
                elif kl == "cpu":
                    d["cpu"] = fv
                elif kl == "mem":
                    d["mem"] = fv
        parts[nm.group(1)] = d
    return parts


def gpu_types(part):
    out = subprocess.run(["sinfo", "-h", "-p", part, "-o", "%G"],
                         stdout=subprocess.PIPE, stderr=subprocess.DEVNULL,
                         text=True).stdout
    return sorted({m for m in re.findall(r"gpu:([a-zA-Z0-9_.-]+)", out)})


def observed_speed(runs="data/runs"):
    """ns/day actually achieved, per system, from our own production logs.

    Column 8 of OpenMM's StateDataReporter output is the running speed estimate.  Taken
    from the last line, which is the estimate over the longest window.
    """
    by_sys = defaultdict(list)
    for log in glob.glob(os.path.join(runs, "*", "rep*", "03_production", "production.log")):
        sysname = log.split(os.sep)[-4]
        try:
            with open(log) as fh:
                last = None
                for line in fh:
                    if line[:1].isdigit():
                        last = line
            if last:
                by_sys[sysname].append(float(last.split()[7]))
        except (OSError, IndexError, ValueError):
            continue
    return {k: sum(v) / len(v) for k, v in by_sys.items() if v}


def bill(p, gpus, cpus, mem_gb):
    """Billing units per hour of ELAPSED time."""
    return p["gpu"] * gpus + p["cpu"] * cpus + p["mem"] * mem_gb


def fairshare():
    out = subprocess.run(["sshare", "-U", "-u", os.environ.get("USER", ""),
                          "-o", "Account,FairShare,RawUsage", "-n", "-P"],
                         stdout=subprocess.PIPE, stderr=subprocess.DEVNULL,
                         text=True).stdout
    for line in out.strip().splitlines():
        f = line.split("|")
        if len(f) >= 3:
            return f[0].strip(), f[1].strip(), f[2].strip()
    return None


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--ns", type=float, default=50.0, help="nanoseconds per replicate")
    ap.add_argument("--reps", type=int, default=3)
    ap.add_argument("--sys", default=None,
                    help="system name, to use its own observed throughput")
    ap.add_argument("--gpus", type=int, default=1)
    ap.add_argument("--cpus", type=int, default=2)
    ap.add_argument("--mem", type=float, default=8.0, help="GB")
    ap.add_argument("--plan", default=None,
                    help="file of 'SYSTEM ns reps' lines to price as one batch")
    ap.add_argument("--runs", default="data/runs")
    a = ap.parse_args()

    parts = partitions()
    gpu_parts = {k: v for k, v in parts.items() if v["gpu"] or k in FREE}
    if not gpu_parts:
        sys.exit("no GPU partitions visible -- is this a login node?")
    speeds = observed_speed(a.runs)

    fs = fairshare()
    if fs:
        print(f"\n  account {fs[0]}   FairShare {fs[1]}   RawUsage {int(float(fs[2])):,}")
        if float(fs[1]) < 0.1:
            print("  FairShare is near zero: the allocation is spent. Submit nothing "
                  "billed.")

    # ---- the plan, and its total ----
    if a.plan:
        jobs = []
        for line in open(a.plan):
            line = line.split("#")[0].strip()
            if not line:
                continue
            f = line.split()
            jobs.append((f[0], float(f[1]) if len(f) > 1 else a.ns,
                         int(f[2]) if len(f) > 2 else a.reps))
    else:
        jobs = [(a.sys or "(unnamed)", a.ns, a.reps)]

    def hours(sysname, ns):
        sp = speeds.get(sysname)
        return (ns / sp * 24.0, True) if sp else (None, False)

    print(f"\n  {'system':16s} {'ns':>6} {'reps':>5} {'ns/day':>8} {'GPU-hours':>10}")
    total_h, unknown = 0.0, []
    for sysname, ns, reps in jobs:
        h, known = hours(sysname, ns)
        sp = speeds.get(sysname)
        if known:
            total_h += h * reps
            print(f"  {sysname:16s} {ns:6.0f} {reps:5d} {sp:8.0f} {h * reps:10.1f}")
        else:
            unknown.append(sysname)
            print(f"  {sysname:16s} {ns:6.0f} {reps:5d} {'—':>8} {'unmeasured':>10}")
    if unknown:
        print(f"\n  no throughput logged yet for: {', '.join(sorted(set(unknown)))}")
        print("  run one replicate on gpu_test (free) first; the cost of everything "
              "else follows from its ns/day.")
    if not total_h:
        return

    print(f"\n  {total_h:.1f} GPU-hours of elapsed time. What each partition charges "
          f"for it:\n")
    print(f"  {'partition':16s} {'per hour':>9} {'batch cost':>14} {'vs cheapest':>12}  "
          f"{'limit':>12}  hardware")
    rows = []
    for name, p in gpu_parts.items():
        per = bill(p, a.gpus, a.cpus, a.mem)
        rows.append((per * total_h, per, name, p))
    rows.sort()
    best = next((c for c, *_ in rows if c > 0), None)
    for cost, per, name, p in rows:
        rel = "free" if cost == 0 else (f"{cost / best:.1f}x" if best else "—")
        types = gpu_types(name)
        hw = ", ".join(t.replace("nvidia_", "") for t in types[:3]) or "—"
        if len(types) > 3:
            hw += f", +{len(types) - 3} more"
        star = "  <-- " if cost == 0 else ""
        print(f"  {name:16s} {per:9.1f} {cost:14,.0f} {rel:>12}  "
              f"{str(p['maxtime']):>12}  {hw}{star}")

    print(f"\n  For reference, the day that exhausted the allocation cost 194,084 "
          f"billing-hours.")
    cheap_billed = next((r for r in rows if r[0] > 0), None)
    dear = next((r for r in rows if r[2] == "gpu"), None)
    if cheap_billed and dear and dear[0] > cheap_billed[0]:
        print(f"  Choosing {cheap_billed[2]} over gpu saves "
              f"{dear[0] - cheap_billed[0]:,.0f} billing-hours on this batch alone "
              f"({dear[0] / cheap_billed[0]:.1f}x).")
    print("  Requesting less memory or fewer cores changes the total by well under 1 %: "
          "the GPU is 99.3 % of it.")
    print("  An accurate --time does not reduce the bill (that is charged on elapsed "
          "time) but it does get the job scheduled sooner.\n")


if __name__ == "__main__":
    main()
