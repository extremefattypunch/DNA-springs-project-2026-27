#!/usr/bin/env python3
"""How much wavelength shift can 2-7 pN of DNA spring actually buy?

The question the wiki page has to answer
----------------------------------------
A memo proposed re-siting the spring to move sfGFP's excitation or emission wavelength.
Before choosing between candidate geometries it is worth asking whether ANY geometry can
deliver a 10-20 nm shift at the force a bent duplex provides, because if the answer is no
then the ranking of sites is beside the point.

Two independent routes, from measured quantities in the primary literature.  They use
completely different inputs and should not agree unless the physics is right.

Route A -- barrel compliance
    Dietz & Rief (2006) PNAS, doi:10.1073/pnas.0602995103, measured GFP's directional
    spring constants by mechanical triangulation: 1 to 17 N/m depending on the pulling
    axis.  A force F therefore displaces the barrel by F/k.
    Barstow et al. (2008) PNAS, doi:10.1073/pnas.0802252105, measured the spectral
    consequence of a known chromophore deformation: hydrostatic pressure blue-shifts
    Citrine's emission 527 -> 510 nm, 17 nm, while deforming the chromophore "up to
    0.8 A at an applied pressure of 500 MPa".  That is the nm-per-angstrom conversion.

Route B -- pressure equivalence
    Take the pressure that produced the 17 nm shift and ask what force it exerts across
    the barrel's cross-section.  Brejc et al. (1997) PNAS give the barrel as "a nearly
    perfect cylinder 42 A long and 24 A in diameter".

Everything printed as "measured" is verbatim from those papers.  Everything printed as
"derived" is arithmetic here, and is labelled as such -- the numbers below are not in any
paper and should not be cited as if they were.
"""
from __future__ import annotations

import math

# ---- measured inputs, all verbatim from the primary literature ----
K_SOFT_N_PER_M = 1.0        # Dietz & Rief 2006: "fracture forces widely varying", k = 1-17 N/m
K_STIFF_N_PER_M = 17.0
CITRINE_SHIFT_NM = 17.0     # Barstow 2008: 527 nm ambient -> 510 nm at 360 MPa
CITRINE_DEFORM_A = 0.8      # Barstow 2008: "progressive deformation of up to 0.8 A" at 500 MPa
PRESSURE_MPA = 360.0        # the pressure at which the 17 nm shift is complete
BARREL_DIAM_A = 24.0        # Brejc 1997: "24 A in diameter"
SPRING_PN = (2.0, 7.0)      # this project's own measured/modelled spring range
WANTED_NM = (10.0, 20.0)    # the shift the experiment would need to be worth doing


def main():
    print("  MEASURED INPUTS")
    print(f"    GFP barrel spring constants        {K_SOFT_N_PER_M:.0f}-{K_STIFF_N_PER_M:.0f} N/m"
          "   Dietz & Rief 2006 PNAS")
    print(f"    Citrine emission shift             {CITRINE_SHIFT_NM:.0f} nm "
          f"(527 -> 510)      Barstow 2008 PNAS")
    print(f"    for a chromophore deformation of   {CITRINE_DEFORM_A:.1f} A"
          "                Barstow 2008 PNAS")
    print(f"    the spring this project delivers   {SPRING_PN[0]:.0f}-{SPRING_PN[1]:.0f} pN"
          "                 this project")
    print()

    sens = CITRINE_SHIFT_NM / CITRINE_DEFORM_A
    print("  ROUTE A -- barrel compliance x spectral sensitivity   [derived]")
    print(f"    spectral sensitivity = {CITRINE_SHIFT_NM:.0f} nm / {CITRINE_DEFORM_A:.1f} A"
          f" = {sens:.0f} nm per A of chromophore deformation")
    for name, k in (("softest axis", K_SOFT_N_PER_M), ("stiffest axis", K_STIFF_N_PER_M)):
        # 1 N/m = 1000 pN/nm
        k_pn_per_nm = k * 1000.0
        for F in SPRING_PN:
            d_nm = F / k_pn_per_nm
            d_A = d_nm * 10.0
            print(f"      {name:13s} k={k:4.0f} N/m   F={F:.0f} pN  ->  "
                  f"{d_A:7.4f} A  ->  {d_A * sens:6.2f} nm")
    best = SPRING_PN[1] / (K_SOFT_N_PER_M * 1000.0) * 10.0 * sens
    print(f"    best case, softest axis at {SPRING_PN[1]:.0f} pN: {best:.1f} nm")
    print()

    r_m = BARREL_DIAM_A * 1e-10 / 2.0
    area = math.pi * r_m ** 2
    force_N = PRESSURE_MPA * 1e6 * area
    force_pN = force_N * 1e12
    print("  ROUTE B -- what force the pressure experiment actually applied   [derived]")
    print(f"    barrel cross-section  pi x ({BARREL_DIAM_A/2:.0f} A)^2 = {area:.3e} m^2")
    print(f"    {PRESSURE_MPA:.0f} MPa across it = {force_pN:,.0f} pN")
    print(f"    ratio to our spring: {force_pN / SPRING_PN[1]:,.0f}x at "
          f"{SPRING_PN[1]:.0f} pN")
    print(f"    so scaling the 17 nm linearly: "
          f"{CITRINE_SHIFT_NM / (force_pN / SPRING_PN[1]):.3f} nm")
    print()

    print("  THE TWO ROUTES, AND WHAT THEY MEAN   [derived]")
    print(f"    Route A (compliance):   {best:.1f} nm at the most favourable axis and force")
    print(f"    Route B (pressure):     {CITRINE_SHIFT_NM / (force_pN / SPRING_PN[1]):.2f} nm")
    print(f"    wanted:                 {WANTED_NM[0]:.0f}-{WANTED_NM[1]:.0f} nm")
    lo = min(best, CITRINE_SHIFT_NM / (force_pN / SPRING_PN[1]))
    print(f"    short by a factor of {WANTED_NM[0] / best:.0f} on the optimistic route "
          f"and {WANTED_NM[0] / lo:.0f} on the conservative one.")
    print()
    print("    The two routes use different measurements and disagree by an order of")
    print("    magnitude with each other, which is itself informative: Route A takes the")
    print("    softest of five measured axes and applies the full spectral sensitivity to")
    print("    it, so it is a ceiling. Route B assumes the deformation scales linearly")
    print("    with load, which for a barrel loaded at two points rather than isotropically")
    print("    is a floor. The honest statement is that the answer lies between them, and")
    print("    both ends are far below what the experiment would need.")


if __name__ == "__main__":
    main()
