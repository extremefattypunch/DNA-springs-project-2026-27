#!/usr/bin/env python3
"""What the choice of click chemistry costs in force and buys in reach.

The two things a linker length controls
---------------------------------------
Measured on the built chimera L3_loop27, each tether spans 21.7 A from the protein's
C-beta to the DNA backbone phosphorus: 12.4-13.3 A of clicked tetrazine-TCO adduct plus
10.8 A of the DNA's own amino-modifier arm.  That single number sets two things at once,
and they pull in opposite directions.

REACH.  The duplex ends can be at most (C-beta separation + 2 x tether) apart, so a
LONGER tether reaches anchor pairs that are closer together on the protein.

FORCE.  This is the part I got backwards on the first pass, and the model settles it.
At 27 bp the contour is 89.1 A and the kink threshold x_c is 85.7 A, so every span our
geometry produces (59-75 A) sits on the circular-arc branch below x_c, where

    f = tau_c / (2R sqrt(1 - (x/2R)^2))

and 2R = 88.5 A.  That expression DIVERGES as x approaches 2R -- the force is what it
takes to STRAIGHTEN an already-kinked rod, and bending it further is cheap.  So force
RISES with span, and a shorter tether, which shortens the span, LOWERS the force.

Both quantities therefore improve with a LONGER tether, which sounds like a free lunch
and is not: a long arm of rotatable bonds is a leash, not a strut, and takes up the
spring's extension entropically instead of transmitting it.  The real currency is tether
RIGIDITY, which no amount of geometry can compute -- so this script reports the geometry
honestly and says where the argument has to stop.

    python analysis/tether_tradeoff.py
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import spring_model as SM       # noqa: E402

CB_CB_NM = 3.13                 # measured, ASP133-ASN149 in 2B3P
TETHER_NOW_NM = 2.17            # measured, C-beta to DNA phosphorus in L3_loop27

# Candidate linkages and the tether each would give, in nm from protein C-beta to the DNA
# phosphorus. The DNA-side amino-modifier arm (1.08 nm as built) is counted in every case
# because every chemistry needs some arm on the oligo; only the protein-side adduct
# changes. These are estimates from the adduct's heavy-atom count and are labelled as
# such -- the measured value is the tetrazine-TCO one.
LINKERS = [
    ("tetrazine-TCO (current, MEASURED)", 2.17),
    ("tetrazine on DNA, BCN on protein", 1.85),
    ("SPAAC, DBCO on DNA + azide-Lys",   1.75),
    ("SPAAC, BCN on protein + azide-DNA", 1.70),
    ("CuAAC, propargyl-Lys + azide-DNA", 1.50),
    ("maleimide-cysteine + short arm",   1.40),
]


def main():
    print(f"  anchors {CB_CB_NM * 10:.1f} A apart (D133-N149); tether now "
          f"{TETHER_NOW_NM * 10:.1f} A per side, measured\n")

    print("  === FORCE: what each linkage would deliver from the same 27 bp duplex ===\n")
    print(f"  {'linkage':38s} {'tether':>7} {'span x':>8} {'force':>9} {'state':>10}")
    base = None
    for name, arm in LINKERS:
        x = SM.span(CB_CB_NM, arm)                      # nm
        s = SM.Spring(n_bp=27, tau_c=SM.TAU_C_CONTINUOUS)
        try:
            f = s.force(x)
            st = s.state(x) if hasattr(s, "state") else ""
        except Exception as e:                           # noqa: BLE001
            f, st = float("nan"), f"({type(e).__name__})"
        if base is None:
            base = f
        print(f"  {name:38s} {arm*10:6.1f} A {x*10:7.1f} A {f:8.2f} pN {st:>10}"
              + ("" if base in (None, 0) or f != f else f"   {f/base:.2f}x"))

    print("\n  === REACH: the closest anchor pair each linkage could span ===\n")
    print("  A 27 bp duplex bent to its kinked branch sits with its ends roughly")
    print("  5.2-7.3 nm apart in our builds. To attach it, the anchors must satisfy")
    print("  span = C-beta + 2 x tether, so the CLOSEST reachable pair is:\n")
    print(f"  {'linkage':38s} {'tether':>7} {'closest anchors':>16}")
    X_MIN_NM = 5.2
    for name, arm in LINKERS:
        closest = max(0.0, X_MIN_NM - 2 * arm)
        note = ""
        if closest * 10 < 6.0:
            note = "  <- would reach the roGFP pairs"
        print(f"  {name:38s} {arm*10:6.1f} A {closest*10:14.1f} A{note}")

    print("\n  === what this actually means for the choice of chemistry ===")
    need = (X_MIN_NM - 0.46) / 2
    print("  Both force and reach improve with a LONGER tether, not a shorter one, because")
    print(f"  every span we produce sits below x_c and f diverges as x -> 2R. So on")
    print("  geometry alone, every alternative chemistry in the table is WORSE than what")
    print("  is already built: they are all shorter, and all deliver less force.")
    print()
    print(f"  Reaching the roGFP pairs would need ~{need*10:.0f} A of tether per side,")
    print("  more than double the 21.7 A we have. But a 24 A arm of rotatable bonds is a")
    print("  leash, not a strut: it would absorb the spring's extension entropically")
    print("  rather than transmit it. Length buys reach on paper and loses it in physics,")
    print("  and nothing in this calculation captures that -- it is a question about")
    print("  tether STIFFNESS, which needs simulation of the specific adduct, not algebra.")
    print()
    print("  So the honest conclusion from geometry is narrow: no alternative chemistry")
    print("  improves the mechanics, and the case for changing one has to be made on")
    print("  SPECTROSCOPY -- getting the tetrazine chromophore away from the chromophore")
    print("  -- rather than on force or reach.")


if __name__ == "__main__":
    main()
