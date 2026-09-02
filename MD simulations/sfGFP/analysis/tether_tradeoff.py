#!/usr/bin/env python3
"""What the choice of click chemistry costs in force, and what it buys in reach.

The question this answers
------------------------
"Instead of tetrazine-TCO, what if we used a different click linker?"  The chemistry is
argued elsewhere (rate, selectivity, spectroscopy); this script settles only the
mechanics, because the mechanics turn out to settle the question.

The one measured number
-----------------------
On the built chimera L3_loop27 each tether spans 21.7 A from the protein's C-beta to the
DNA backbone phosphorus: 12.4-13.3 A of clicked tetrazine-TCO adduct plus 10.8 A of the
DNA's own amino-modifier arm.  That single number sets both reach and force.

FORCE rises with span, which is the part that is easy to get backwards.  At 27 bp the
contour is 89.1 A and the kink threshold x_c is 85.7 A, so every span our geometry
produces (36-75 A) sits on the circular-arc branch below x_c, where

    f = tau_c / (2R sqrt(1 - (x/2R)^2)),   2R = 88.5 A

which DIVERGES as x -> 2R.  The force is what it takes to STRAIGHTEN an already-kinked
rod; bending it further is cheap.  So a shorter tether shortens the span and LOWERS the
force, and every alternative chemistry is shorter than what is already built.

REACH is not a constraint at all.  A compression spring does not have to be stretched to
fit: if the anchors are closer together, the duplex simply bends further.  The old
"maximum span = C-beta + 2 x tether" line in this script treated the tether as a rigid
strut and asked whether it was long enough to reach.  That was the wrong question.

Where the tether is compliant, and what that costs
--------------------------------------------------
The clicked adduct is a fused polycyclic -- dihydropyridazine fused to the sTCO bicycle,
no sp3 chain bonds -- so it is a strut.  The DNA-side amino-modifier arm is an alkyl
chain, and it is not.  The load path runs through both in series, so the equilibrium span
is a FIXED POINT: the duplex pushes the anchors apart, the compliant arm extends until its
tension matches, and the span that results sets the force.  Marko-Siggia closes the loop
without any new simulation:

    f lp / kT = z + 1/(4(1-z)^2) - 1/4,     z = extension / contour

At lp = 0.4 nm an alkyl arm is 12 % extended at 2 pN and 35 % at 7 pN.  Inverted: holding
such an arm at 90 % of contour would take 263 pN and at 99 %, 25,600 pN.  A few-pN spring
cannot pull a flexible arm taut, full stop -- so added flexible length is added contour
the spring never sees.

The result that answers the question
------------------------------------
Reading the measured 21.7 A as a rigid strut gives 7.2 pN at the built pair; solving the
compliant fixed point gives 5.2 pN.  That 2.0 pN modelling gap is the same size as the
ENTIRE spread across the six candidate chemistries (2.0 pN).  Chemistry choice sits
inside the error bar of the tether-stiffness model, and on our own clamp stiffness
(140 pN/A) 2 pN moves the anchors 0.014 A.  So mechanics cannot pick a chemistry -- but
it does say what to change instead: stiffen the DNA-side arm at constant length.  That is
worth more than any chemistry on the list.

    python analysis/tether_tradeoff.py [--json analysis/tether_tradeoff.json]
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import spring_model as SM       # noqa: E402

CB_CB_NM = 3.13                 # measured, ASP133-ASN149 in 2B3P
TETHER_NOW_NM = 2.17            # measured, C-beta to DNA phosphorus in L3_loop27
ADDUCT_NM = 1.30                # rigid, fused-ring clicked tetrazine-TCO (12.4-13.3 A)
ARM_NM = 1.08                   # compliant, DNA-side amino-modifier
LP_ALKYL_NM = 0.40              # persistence length of an alkyl/PEG chain
KT = 4.1                        # pN nm
CLAMP_STIFFNESS_PN_PER_A = 140.0   # measured, figures/fig1_force_response.csv S6_clamp

# Candidate linkages and the tether each would give, in nm from protein C-beta to the DNA
# phosphorus. The DNA-side amino-modifier arm (1.08 nm as built) is counted in every case
# because every chemistry needs some arm on the oligo; only the protein-side adduct
# changes. These are estimates from the adduct's heavy-atom count and are labelled as
# such -- the measured value is the tetrazine-TCO one.
LINKERS = [
    ("tetrazine-sTCO (current, MEASURED)", 2.17),
    ("inverted: tetrazine-DNA + TCO*-Lys", 1.85),
    ("SPAAC, DBCO on DNA + azide-Lys",     1.75),
    ("SPAAC, BCN on protein + azide-DNA",  1.70),
    ("CuAAC, propargyl-Lys + azide-DNA",   1.50),
    ("maleimide-cysteine + short arm",     1.40),
    ("oxime, pAcF + aminooxy-DNA",         1.33),
]

# Candidate anchor pairs, C-beta to C-beta (nm), measured on 2B3P.
PAIRS = [
    ("SER147-GLN204", 0.46),
    ("ASN149-SER202", 0.59),
    ("TYR151-SER208", 2.80),
    ("ASP133-ASN149 (built)", 3.13),
]


def ext_frac(f: float, lp: float = LP_ALKYL_NM) -> float:
    """Marko-Siggia inverted: fractional extension of a WLC at tension ``f`` (pN)."""
    if f <= 0.0:
        return 0.0
    lo, hi, target = 0.0, 1.0 - 1e-12, f * lp / KT
    for _ in range(200):                     # bisection; the forward form is monotone
        z = 0.5 * (lo + hi)
        if z + 1.0 / (4.0 * (1.0 - z) ** 2) - 0.25 < target:
            lo = z
        else:
            hi = z
    return 0.5 * (lo + hi)


def wlc_force(z: float, lp: float = LP_ALKYL_NM) -> float:
    """Forward Marko-Siggia: tension (pN) needed to hold a WLC at extension ``z``."""
    return (KT / lp) * (z + 1.0 / (4.0 * (1.0 - z) ** 2) - 0.25)


def solve(cb: float, rigid: float, arm: float, *, n_bp: int = 27,
          tau_c: float = SM.TAU_C_CONTINUOUS) -> tuple[float, float, float] | tuple[None, None, None]:
    """Self-consistent span and force for a rigid strut plus a compliant arm in series.

    Returns ``(force_pN, span_nm, tether_nm)``, or a triple of ``None`` if the geometry
    demands a span the duplex cannot reach even fully straight.
    """
    s = SM.Spring(n_bp=n_bp, tau_c=tau_c)
    f = 5.0
    for _ in range(500):                     # damped fixed point; converges monotonically
        tether = rigid + arm * ext_frac(f)
        x = cb + 2.0 * tether
        if not s.is_feasible(x):
            return None, None, None
        f = 0.5 * f + 0.5 * s.force(x)
    return f, x, tether


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--json", default="analysis/tether_tradeoff.json")
    a = ap.parse_args()
    s27 = SM.Spring(n_bp=27, tau_c=SM.TAU_C_CONTINUOUS)
    out: dict = {"cb_cb_A": CB_CB_NM * 10, "tether_measured_A": TETHER_NOW_NM * 10,
                 "adduct_rigid_A": ADDUCT_NM * 10, "arm_contour_A": ARM_NM * 10,
                 "lp_alkyl_A": LP_ALKYL_NM * 10, "n_bp": 27,
                 "x_c_A": s27.x_c * 10, "contour_A": s27.L * 10}

    print(f"  anchors {CB_CB_NM * 10:.1f} A apart (D133-N149); tether "
          f"{TETHER_NOW_NM * 10:.1f} A per side, measured on L3_loop27")
    print(f"  27 bp duplex: contour {s27.L*10:.1f} A, kink threshold x_c "
          f"{s27.x_c*10:.1f} A -- every span below is on the arc branch\n")

    print("  === how far a few pN can pull a flexible arm (Marko-Siggia, lp = 4 A) ===\n")
    out["wlc"] = []
    for f in (1, 2, 4, 7, 12, 20):
        z = ext_frac(f)
        out["wlc"].append({"force_pN": f, "ext_frac": z})
        print(f"  {f:5.1f} pN  ->  {z*100:5.1f} % of contour")
    print()
    for z in (0.90, 0.99):
        print(f"  holding it at {z*100:.0f} % of contour would take "
              f"{wlc_force(z):,.0f} pN")
    out["force_at_90pct_pN"] = wlc_force(0.90)
    out["force_at_99pct_pN"] = wlc_force(0.99)
    print("\n  -> a 2-7 pN spring cannot pull a flexible arm taut. Added flexible length")
    print("     is added contour the spring never sees.\n")

    print("  === FORCE: what each linkage delivers from the same 27 bp duplex ===\n")
    print(f"  {'linkage':36s} {'tether':>7} {'rigid read':>11} {'compliant':>11}")
    out["linkers"] = []
    rigid_f, comp_f = [], []
    for name, tether in LINKERS:
        x = SM.span(CB_CB_NM, tether)
        f_rigid = s27.force(x) if s27.is_feasible(x) else float("nan")
        f_comp, x_comp, t_comp = solve(CB_CB_NM, max(tether - ARM_NM, 0.0), ARM_NM)
        rigid_f.append(f_rigid)
        comp_f.append(f_comp)
        out["linkers"].append({"name": name, "tether_A": tether * 10,
                               "span_rigid_A": x * 10, "force_rigid_pN": f_rigid,
                               "span_compliant_A": None if x_comp is None else x_comp * 10,
                               "force_compliant_pN": f_comp})
        print(f"  {name:36s} {tether*10:5.1f} A {f_rigid:8.2f} pN {f_comp:8.2f} pN")

    spread_r = max(rigid_f) - min(rigid_f)
    spread_c = max(comp_f) - min(comp_f)
    gap = rigid_f[0] - comp_f[0]
    out.update({"spread_rigid_pN": spread_r, "spread_compliant_pN": spread_c,
                "model_gap_pN": gap})
    print(f"\n  spread across all {len(LINKERS)} chemistries : {spread_r:.2f} pN (rigid read), "
          f"{spread_c:.2f} pN (compliant)")
    print(f"  rigid-vs-compliant gap, current chemistry : {gap:.2f} pN")
    print(f"  -> the MODELLING gap is as large as the whole chemistry spread.")
    print(f"     On our measured clamp stiffness ({CLAMP_STIFFNESS_PN_PER_A:.0f} pN/A), "
          f"{spread_r:.1f} pN moves the anchors "
          f"{spread_r/CLAMP_STIFFNESS_PN_PER_A:.3f} A.\n")
    out["anchor_motion_from_spread_A"] = spread_r / CLAMP_STIFFNESS_PN_PER_A

    print("  === REACH is not the constraint: a compression spring just bends further ===\n")
    print(f"  {'anchor pair':24s} {'Cb-Cb':>7} {'tether':>8} {'span':>8} {'arm ext':>8} {'force':>9}")
    out["pairs"] = []
    for name, cb in PAIRS:
        f, x, t = solve(cb, ADDUCT_NM, ARM_NM)
        z = ext_frac(f)
        out["pairs"].append({"pair": name, "cb_cb_A": cb * 10, "tether_A": t * 10,
                             "span_A": x * 10, "arm_ext_frac": z, "force_pN": f})
        print(f"  {name:24s} {cb*10:5.1f} A {t*10:6.2f} A {x*10:6.2f} A "
              f"{z*100:6.1f} % {f:7.2f} pN")
    print("\n  Closer anchors give a shorter span and therefore LESS force, monotonically.")
    print("  Nothing is unreachable; the close pairs are simply weaker.\n")

    # The measured tether and the WLC tether disagree, and the gap is the honest error bar.
    z_measured = (TETHER_NOW_NM - ADDUCT_NM) / ARM_NM
    out["arm_ext_frac_measured"] = z_measured
    out["force_to_hold_measured_arm_pN"] = wlc_force(min(z_measured, 0.999999))
    print("  === where the measured geometry and this model disagree ===\n")
    print(f"  The MD tether is {TETHER_NOW_NM*10:.1f} A, i.e. the arm sits at "
          f"{z_measured*100:.0f} % of its {ARM_NM*10:.1f} A contour.")
    print(f"  Entropic elasticity alone would need "
          f"{wlc_force(min(z_measured, 0.999999)):,.0f} pN to hold it there, which the")
    print("  spring does not have. So either the 50 ns snapshot is not at mechanical")
    print("  equilibrium, or the arm is extended by excluded volume and solvation rather")
    print("  than by tension -- an alkyl chain pinned between a charged backbone and a")
    print("  protein surface is not a free WLC. Both readings are reported and neither is")
    print(f"  discarded: the force at the built pair is {comp_f[0]:.1f}-{rigid_f[0]:.1f} pN,")
    print("  and the WLC estimate is the lower bound.\n")


    print("  === what actually moves the force ===\n")
    f0, _, _ = solve(CB_CB_NM, ADDUCT_NM, ARM_NM)
    f_stiff, _, _ = solve(CB_CB_NM, ADDUCT_NM + ARM_NM, 0.0)
    out["levers"] = [{"lever": "as built (compliant arm)", "force_pN": f0}]
    print(f"  {'lever':38s} {'force':>9} {'change':>9}")
    print(f"  {'as built, compliant arm':38s} {f0:7.2f} pN {'--':>9}")
    print(f"  {'STIFFEN the arm at constant contour':38s} {f_stiff:7.2f} pN "
          f"{100*(f_stiff/f0-1):+8.0f} %")
    out["levers"].append({"lever": "stiffen the arm at constant contour",
                          "force_pN": f_stiff, "pct": 100 * (f_stiff / f0 - 1)})
    for n in (25, 30, 40):
        fn, _, _ = solve(CB_CB_NM, ADDUCT_NM, ARM_NM, n_bp=n)
        if fn is None:
            continue
        print(f"  {f'27 -> {n} bp':38s} {fn:7.2f} pN {100*(fn/f0-1):+8.0f} %")
        out["levers"].append({"lever": f"27 -> {n} bp", "force_pN": fn,
                              "pct": 100 * (fn / f0 - 1)})
    print(f"  {'best -> worst chemistry swap':38s} {'--':>9} "
          f"{-100*spread_c/comp_f[0]:+8.0f} %")

    print("\n  === the conclusion mechanics can actually support ===")
    print("  No alternative chemistry improves the mechanics: all are shorter, all")
    print("  deliver less force, and the whole spread sits inside the tether-stiffness")
    print("  model's own uncertainty. RIGIDITY, not identity, is the variable -- and the")
    print("  one compliant segment in the linkage is the DNA-side amino-modifier arm,")
    print("  which is also the cheapest component in the construct to change.")

    Path(a.json).write_text(json.dumps(out, indent=2))
    print(f"\n  wrote {a.json}")


if __name__ == "__main__":
    main()
