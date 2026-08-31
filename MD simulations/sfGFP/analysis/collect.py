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
import contacts  # noqa: E402
import sites as sites_mod  # noqa: E402

CLAMP_RE = re.compile(r"^(?P<base>.+?)_(?P<f>[\d.]+)pN$")


def split_tag(tag: str):
    m = CLAMP_RE.match(tag)
    if m:
        return m.group("base"), float(m.group("f"))
    return tag, None


_HAS_PROTEIN = {}


def has_protein(sysdir: Path) -> bool:
    """Does this system contain a protein at all?

    The bare-duplex systems -- the force-extension windows and the seeded kinks -- have
    no protein, so every analysis in this file is meaningless for them: there is no
    chromophore, no attachment sites, no barrel and no protein-DNA contact.  They would
    still be globbed (K27_bp7/rep1 matches `*/rep*/03_production/traj.dcd` and has a
    build directory), and would land in summary.csv as rows of failed analyses plus
    three spurious entries in the system count.  force_extension.py is where they belong.

    Decided from the topology rather than from a name pattern, so a new bare-duplex
    system is handled without anyone remembering to add a prefix here.
    """
    key = str(sysdir)
    if key not in _HAS_PROTEIN:
        try:
            import mdtraj as md
            top = md.load_prmtop(str(sysdir / "system.prmtop"))
            _HAS_PROTEIN[key] = bool(top.select("protein").size)
        except Exception:                                        # noqa: BLE001
            _HAS_PROTEIN[key] = True        # if unsure, analyse it and let it fail loudly
    return _HAS_PROTEIN[key]


def flat(prefix, d, out):
    for k, v in (d or {}).items():
        key = f"{prefix}{k}"
        if isinstance(v, dict):
            flat(f"{key}.", v, out)
        elif isinstance(v, (int, float, str, bool)) or v is None:
            out[key] = v
    return out


def sites_for(sysdir: Path, default):
    """The two modified residues of THIS system, in its own numbering.

    Not a constant.  sfGFP is modified at 133/149 and ubiquitin at 9/63, and passing
    sfGFP's pair to a ubiquitin trajectory does not fail -- it silently produces no
    anchor distance at all, which is how every ubiquitin system came to have a blank
    Cb-Cb column.  The system builder already writes the two clamp atom indices, so the
    residue numbers are recoverable from them and the residue map with no rebuild.
    """
    ca = sysdir / "clamp_atoms.txt"
    rm = sysdir / "residue_map.json"
    if not (ca.exists() and rm.exists()):
        return default
    try:
        idx = [int(v) for v in ca.read_text().strip().split(",")]
        rows = json.loads(rm.read_text())
        import mdtraj as md
        top = md.load_prmtop(str(sysdir / "system.prmtop"))
        by_index = {r["index"]: r["orig_resnum"] for r in rows}
        out = []
        for a in idx:
            ri = top.atom(a).residue.index
            if ri in by_index:
                out.append(by_index[ri])
        return out if len(out) == 2 else default
    except Exception:                                        # noqa: BLE001
        return default


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--runs", default="data/runs")
    ap.add_argument("--systems", default="build/systems")
    ap.add_argument("--out", default="analysis")
    ap.add_argument("--stride", type=int, default=2)
    ap.add_argument("--sites", nargs=2, type=int, default=[133, 149],
                    help="override the modified residues; by default each system's own "
                         "pair is recovered from its clamp_atoms.txt")
    ap.add_argument("--only", nargs="+", default=None, metavar="SYSTEM",
                    help="re-analyse just these systems and MERGE the result into the "
                         "existing summary, leaving every other row untouched.  A full "
                         "pass costs the better part of an hour, and most code fixes "
                         "touch one code path on a handful of systems.")
    ap.add_argument("--min-ns", type=float, default=1.0,
                    help="skip replicates with less than this much production")
    ap.add_argument("--skip-ns", type=float, default=5.0,
                    help="production time to discard as settling before averaging")
    ap.add_argument("--report-ps", type=float, default=20.0,
                    help="trajectory frame spacing, used to convert --skip-ns")
    args = ap.parse_args()
    args.sites_given = any(a.startswith("--sites") for a in sys.argv)

    rows = []
    for traj in sorted(Path(args.runs).glob("*/rep*/03_production/traj.dcd")):
        rep_dir = traj.parent.parent
        tag = rep_dir.parent.name
        if tag.startswith("_"):
            continue
        base, clamp = split_tag(tag)
        if args.only and base not in args.only:
            continue
        sysdir = Path(args.systems) / base
        top, rmap = sysdir / "system.prmtop", sysdir / "residue_map.json"
        if not top.exists() or not rmap.exists():
            print(f"  skip {tag}/{rep_dir.name}: no topology at {sysdir}")
            continue
        if not has_protein(sysdir):
            print(f"  skip {tag}/{rep_dir.name}: protein-free duplex -- "
                  f"force_extension.py handles these")
            continue
        prog = traj.parent / "progress.json"
        ns = json.loads(prog.read_text())["ns_done"] if prog.exists() else 0.0
        if ns < args.min_ns:
            print(f"  skip {tag}/{rep_dir.name}: only {ns} ns so far")
            continue
        skip = int(round(args.skip_ns * 1000.0 / args.report_ps / args.stride))
        pre = str(traj.parent / "analysis")
        print(f"\n=== {tag} {rep_dir.name}  ({ns} ns, first {args.skip_ns} ns discarded) ===")
        sites = args.sites if args.sites_given else sites_for(sysdir, args.sites)
        row = {"system": base, "tag": tag, "clamp_pN": clamp,
               "replicate": rep_dir.name, "ns": ns, "skip_ns": args.skip_ns,
               "sites": "/".join(str(s) for s in sites)}
        try:
            c = chromophore.analyse(traj, top, rmap, pre, args.stride, skip)
            flat("chromo.", c, row)
        except Exception as e:                              # noqa: BLE001
            print(f"  chromophore failed: {type(e).__name__}: {e}")
        try:
            m = mechanics.analyse(traj, top, rmap, sites, pre, args.stride,
                                  clamp, nicked="nick" in base, skip=skip)
            flat("mech.", m, row)
        except Exception as e:                              # noqa: BLE001
            print(f"  mechanics failed: {type(e).__name__}: {e}")
        try:
            b = barrel.analyse(traj, top, rmap, sites, pre, args.stride, skip)
            flat("barrel.", b, row)
        except Exception as e:                              # noqa: BLE001
            print(f"  barrel failed: {type(e).__name__}: {e}")
        try:
            ct = contacts.analyse(traj, top, rmap, pre, args.stride, skip)
            if ct.get("has_dna"):
                flat("contact.", {k: v for k, v in ct.items()
                                  if not isinstance(v, list)}, row)
                row["contact.residues_over_half"] = ";".join(
                    ct["residues_over_half_the_time"])
        except Exception as e:                              # noqa: BLE001
            print(f"  contacts failed: {type(e).__name__}: {e}")
        try:
            sc = sites_mod.analyse(traj, top, rmap, sites, pre, args.stride, skip)
            # DSSP is a per-residue dict of fractions; keep it in the per-run JSON and
            # out of the flat CSV, which would otherwise gain 84 columns of mostly zeros.
            flat("site.", {k: v for k, v in sc.items() if k != "dssp"}, row)
        except Exception as e:                              # noqa: BLE001
            print(f"  sites failed: {type(e).__name__}: {e}")
        rows.append(row)

    if not rows:
        sys.exit("no replicates with enough production yet")
    df = pd.DataFrame(rows)
    out = Path(args.out)
    prev = out / "summary.csv"
    if args.only and prev.exists():
        # Keep every row this pass did not recompute.  Concat rather than update in place,
        # because a re-analysed system can legitimately gain or lose columns -- that is
        # usually WHY it is being re-analysed -- and pandas aligns them for us.
        old = pd.read_csv(prev)
        keep = old[~old.system.isin(args.only)]
        df = pd.concat([keep, df], ignore_index=True, sort=False)
        print(f"  merged {len(rows)} recomputed rows into {len(keep)} kept rows")
    df = df.sort_values(["system", "clamp_pN", "replicate"])
    df.to_csv(prev, index=False)
    print(f"\nwrote {out / 'summary.csv'}: {len(df)} replicates, {len(df.columns)} columns")
    for s, g in df.groupby("system"):
        print(f"  {s:<16} {len(g)} replicate(s), {g['ns'].sum():.0f} ns total")


if __name__ == "__main__":
    main()
