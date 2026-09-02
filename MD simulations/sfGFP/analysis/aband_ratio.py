#!/usr/bin/env python3
"""The A-band shoulder in the 260828 ds/ss/nicked fluorometry, ordered by spring force.

Why this is worth its own script
--------------------------------
The whole question of whether a DNA spring can change sfGFP's colour points at one
observable: the ratio of the ~400 nm neutral-chromophore excitation band to the ~490 nm
anionic band.  That ratio has never been measured here, because the excitation scans
START at 400 nm -- exactly the A-band maximum -- so the band is cut in half and its
amplitude cannot be fitted.

What CAN be extracted from the existing scans is the 400 nm shoulder relative to the
peak.  It moves monotonically with modelled spring force.  That is the predicted
signature and it is already in the data, which makes it worth stating carefully -- and
worth being explicit that three a/b/c scans are technical replicates of ONE preparation,
so the small scatter is instrument repeatability and not evidence about the constructs.

    python analysis/aband_ratio.py
"""
from __future__ import annotations

import csv
import glob
import json
import os
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[3]
RUN = ROOT / "flurometer" / "260828_sfGFPDNA(dsVSssVSnicked)_att1"
CTRL = ROOT / "flurometer" / "260828_sfGFPDNA(dsVSssVSnickedCONTROLS)_att1"
# Group -> construct and its modelled spring force. The ds/ss assignment is the corrected
# one: an earlier version of the viewer had 1 and 2 the wrong way round.
GROUPS = {"1": ("ss", 0.0), "2": ("ds", 3.1), "3": ("nicked", 2.5)}


def load(p):
    w, s, r = [], [], []
    for row in csv.reader(open(p)):
        try:
            w.append(float(row[0])); s.append(float(row[1])); r.append(float(row[3]))
        except (ValueError, IndexError):
            continue
    return np.array(w), np.array(s), np.array(r)


def main():
    out = {"run": RUN.name, "note": "a/b/c are technical replicates of one preparation",
           "groups": {}}
    print(f"  {RUN.name}\n")
    first = sorted(glob.glob(str(RUN / "Ex1?_Data.csv")))[0]
    w0, _, _ = load(first)
    print(f"  excitation scans span {w0.min():.0f}-{w0.max():.0f} nm at "
          f"{np.diff(w0)[0]:.1f} nm")
    print(f"  the A band peaks at ~400 nm, so it is TRUNCATED at the first point and the")
    print(f"  400/490 ratio -- the observable the whole question turns on -- is unmeasured.\n")
    out["scan_range_nm"] = [float(w0.min()), float(w0.max())]
    out["scan_step_nm"] = float(np.diff(w0)[0])
    out["a_band_truncated"] = bool(abs(w0.min() - 400.0) < 1.0)

    print(f"  {'construct':10s} {'F (pN)':>7} {'I(400)/I(peak)':>26} {'mean':>8} {'sd':>8}")
    for g, (name, force) in GROUPS.items():
        vals, blanks = [], []
        for rep in "abc":
            p = RUN / f"Ex{g}{rep}_Data.csv"
            if not p.exists():
                continue
            w, s, r = load(p)
            sc = s / r                                    # lamp-ratio corrected
            vals.append(float(sc[np.argmin(abs(w - 400))] / sc.max()))
            b = CTRL / f"Ex{g}{rep}_Data.csv"
            if b.exists():
                wb, sb, rb = load(b)
                blanks.append(float((sb / rb)[np.argmin(abs(wb - 400))]
                                    / sc[np.argmin(abs(w - 400))]))
        if not vals:
            continue
        v = np.array(vals)
        out["groups"][name] = {
            "force_pN": force, "n_scans": len(vals),
            "ratio": [round(x, 5) for x in vals],
            "mean": round(float(v.mean()), 5),
            "sd": round(float(v.std(ddof=1)), 5),
            "blank_fraction_at_400nm": (round(float(np.mean(blanks)), 4)
                                        if blanks else None)}
        print(f"  {name:10s} {force:7.1f} {'  '.join(f'{x:.4f}' for x in vals):>26} "
              f"{v.mean():8.4f} {v.std(ddof=1):8.4f}")

    gs = out["groups"]
    if len(gs) == 3:
        by_f = sorted(gs.values(), key=lambda e: e["force_pN"])
        seq = [e["mean"] for e in by_f]
        mono = all(a > b for a, b in zip(seq, seq[1:]))
        rel = 100.0 * (seq[-1] - seq[0]) / seq[0]
        out["monotone_in_force"] = bool(mono)
        out["relative_change_pct"] = round(rel, 1)
        print(f"\n  ordered by modelled force: "
              + " -> ".join(f"{x:.4f}" for x in seq)
              + ("   MONOTONE" if mono else "   not monotone"))
        print(f"  relative change, zero-force to highest: {rel:+.1f} %")
        bf = [e["blank_fraction_at_400nm"] for e in gs.values()
              if e["blank_fraction_at_400nm"]]
        if bf:
            print(f"  DNA-only blanks are {100*min(bf):.2f}-{100*max(bf):.2f} % of the "
                  f"signal at 400 nm, so this is not scatter")
        print("\n  BUT: a/b/c are technical replicates of ONE preparation per condition,")
        print("  so the sd above is instrument repeatability. It does not license a claim")
        print("  about the constructs. Independent preparations, matched pH and matched")
        print("  absorbance are what would turn this into evidence.")

    p = Path(__file__).resolve().parent / "aband_ratio.json"
    p.write_text(json.dumps(out, indent=2))
    print(f"\n  wrote {p.name}")


if __name__ == "__main__":
    main()
