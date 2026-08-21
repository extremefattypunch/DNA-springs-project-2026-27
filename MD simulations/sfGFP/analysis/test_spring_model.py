#!/usr/bin/env python3
"""Regression tests: reproduce Zocchi's own published forces and energies.

Run with `python analysis/test_spring_model.py` (no pytest needed).
"""
import sys
from spring_model import KT, Spring, choose_n_bp, span

FAILS = []


def check(label, got, want, tol, unit=""):
    ok = abs(got - want) <= tol
    print(f"{'PASS' if ok else 'FAIL'}  {label:<58} got {got:8.3f} {unit:<6} "
          f"want {want:.3f} +/- {tol}")
    if not ok:
        FAILS.append(label)


# --- Renilla luciferase chimera, Tseng & Zocchi JACS 2013 / Book pp.99-101 ---
# "s = (1.9 + 2 x 2.1) nm = 6.1 nm"; nicked spring, tau_c = 27 pN nm.
s_rluc = span(1.9, 2.1)
check("RLuc span s = d(Cys161,Cys188) + 2*crosslinker", s_rluc, 6.1, 1e-9, "nm")
# JACS p.11884: "we find from eq 6 that the force on the enzyme is f = 1.5 pN"
check("RLuc 60-mer nicked, f (JACS Eq.6 = 1.5 pN)",
      Spring(60, tau_c=27.0).force(s_rluc), 1.5, 0.05, "pN")
# Book p.101: "a force f = 2.4 pN slows down the reaction by a factor 0.7"
check("RLuc 40-mer nicked, f (Book p.101 = 2.4 pN)",
      Spring(40, tau_c=27.0).force(s_rluc), 2.4, 0.1, "pN")
# JACS p.11884: the unkinked (WLC) treatment "gives f = 5 pN, significantly
# overestimating the force" for the 60-mer.  Our linear-elasticity branch is the
# same physical statement (kinking suppressed) and should land in that ballpark.
f_unk = Spring(60, tau_c=27.0).unkinked_force(s_rluc)
check("RLuc 60-mer, kinking suppressed (paper: ~5 pN, i.e. ~3x too big)",
      f_unk, 5.0, 1.5, "pN")
assert f_unk > Spring(60).force(s_rluc), "suppressing kinking must raise the force"

# --- Guanylate kinase chimera, Book p.105 / Methods Enzymol 2021 p.277 ---
# "If the enzyme was not deformed, the end-to-end distance of the DNA would be
#  s = 10 nm.  The elastic energy of just the DNA ... is 6.6 kT."
check("GK 60-bp spring at s = 10 nm, E_DNA (Book p.105 = 6.6 kT)",
      Spring(60, tau_c=27.0).energy(10.0) / KT, 6.6, 0.7, "kT")

# --- Book p.73/p.80 parameter self-consistency ---
# The book quotes B two ways -- "200 pN nm^2" (p.73) and "50 kT nm" (p.80) -- which
# differ by its own rounding of kT: 50 x 4.2 = 210, not 200.  Tolerance covers that.
check("B in kT units (Book p.73 = 200 pN nm^2 vs p.80 = 50 kT nm)",
      200.0 / KT, 50.0, 3.0, "kT nm")
# gamma = 1 exactly at L = 2B/tau_c = 13.33 nm = 40.4 bp for tau_c = 30.  Book p.75
# rounds this up to "about 15 nm or 45 bp"; test the exact boundary, not the rounding.
L_gamma1 = 2.0 * 200.0 / 30.0
check("L where gamma = 1 for tau_c = 30 (Book p.75 rounds to '~15 nm')",
      L_gamma1, 13.33, 0.02, "nm")
check("  ... in base pairs (Book p.75 rounds to '~45 bp')",
      L_gamma1 / 0.33, 40.4, 0.2, "bp")

# --- Structural sanity, and the counterintuitive sign of the softening transition ---
# In the kinked branch f = tau_c / (2R sqrt(1 - (x/2R)^2)), so the force *falls* as the
# spring is compressed further, even though the stored energy rises.  A kinked rod is a
# constant-torque hinge: f = tau_c / (lever arm), and compressing it lengthens the lever
# arm.  This is exactly Zocchi's statement that the softening transition "limits the
# force it can deliver to a few pN" (Methods Enzymol 2021 p.278) -- you cannot get a
# harder spring by squeezing a long one harder, only by using a shorter one.
for n in (25, 40, 60):
    sp = Spring(n)
    assert sp.x0 < sp.L, "thermal bending must shorten the zero-force EED"
    assert sp.x_c < sp.x0, f"{n} bp: kinking must set in below the relaxed EED"
    f_squeezed, f_relaxed = sp.force(sp.x_c * 0.80), sp.force(sp.x_c * 0.999)
    e_squeezed, e_relaxed = sp.energy(sp.x_c * 0.80), sp.energy(sp.x_c * 0.999)
    assert f_squeezed < f_relaxed, f"{n} bp: kinked force must fall with compression"
    assert e_squeezed > e_relaxed, f"{n} bp: stored energy must rise with compression"
print("PASS  x_c < x0 < L, and kinked force falls / energy rises with compression")

# At a fixed span, a longer spring is a *softer* spring -- the practical consequence.
forces = [Spring(n).force(x_) for n, x_ in ((26, 6.73), (40, 6.73), (60, 6.73))]
assert forces[0] > forces[1] > forces[2], forces
print(f"PASS  at fixed span 6.73 nm, force falls with length: "
      f"26 bp {forces[0]:.2f} > 40 bp {forces[1]:.2f} > 60 bp {forces[2]:.2f} pN")

# --- Feasibility guard: a 20 bp spring cannot span the sfGFP geometry ---
x_sfgfp = span(3.13, 2.1)          # measured Asp133-Asn149 Cb-Cb + 2.1 nm arms
try:
    Spring(20).force(x_sfgfp)
    FAILS.append("20 bp should have been rejected as infeasible")
    print("FAIL  20 bp spring wrongly accepted at the sfGFP span")
except ValueError:
    print(f"PASS  20 bp spring rejected at sfGFP span x = {x_sfgfp:.2f} nm "
          f"(contour only {Spring(20).L:.2f} nm)")

n_strong = choose_n_bp(x_sfgfp, 7.0)
print(f"PASS  strong-spring selection at sfGFP span: {n_strong} bp -> "
      f"{Spring(n_strong).force(x_sfgfp):.2f} pN nicked, "
      f"{Spring(n_strong, tau_c=34.0).force(x_sfgfp):.2f} pN continuous")

print()
if FAILS:
    print(f"{len(FAILS)} FAILURE(S): " + "; ".join(FAILS))
    sys.exit(1)
print("all checks passed")
