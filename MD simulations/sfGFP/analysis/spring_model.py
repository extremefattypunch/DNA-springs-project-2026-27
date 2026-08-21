#!/usr/bin/env python3
"""Zocchi's kinked-DNA spring model, implemented verbatim from the primary sources.

A dsDNA rod attached at both ends to two points on a protein surface must bend,
because its persistence length (~50 nm) dwarfs the protein.  Below a critical
end-to-end distance the rod stops bending smoothly and *kinks*, which caps the
force it can deliver at a few pN.  Both branches are implemented here.

Sources
-------
Zocchi, *Molecular Machines: A Materials-Science Approach* (Princeton, 2018):
    Eq. 2.188  gamma = L*tau_c / (2B)
    Eq. 2.193  E0    = (1/6)(L/B) tau_c**2
    Eq. 2.194  x_c   = L (1 - gamma**2/15)
    Eq. 2.197  R     = (L/2)(1 - gamma**2/90)
    Eq. 2.205  x0    = L (1 - T*L/(10 B))
    Eq. 2.206  E(x)  (= Eq. 3.51)
    Eq. 3.52   f(x)
Tseng & Zocchi, *JACS* 135:11879 (2013), Eq. 6 and p. 11884 for the parameters.
    That paper writes 2L = 0.33*N_bp and R = L(1 - 2 gamma**2/45); its L is a
    half-length, so the two forms are algebraically identical.  We use the book's.

Parameter values (Book p.73, p.80; JACS p.11884)
-----------------------------------------------
    B      = 200 pN nm**2  (= 50 kT nm)   bending modulus
    tau_c  = 27 pN nm                     critical bending torque, *nicked* rod
    tau_c  = 31-36 pN nm                  *continuous* rod (measured indirectly)
    rise   = 0.33 nm/bp
    T      = 4.2 pN nm                    kT at room temperature, Zocchi's unit
    valid  for gamma < 1, i.e. L < ~15 nm (~45 bp)

Deliberately NOT used: the F ~ 10 pN / W ~ 25 kT figure of Choi & Zocchi,
*Biophys J* 92:1651 (2007) p.1657.  That paper labels it "an upper bound" from a
worm-like-chain treatment that ignores kinking; the kinked model above supersedes
it and gives 1.5-2.4 pN for the same constructs.

All lengths in nm, forces in pN, energies in pN*nm (divide by T for kT).
"""
from __future__ import annotations

import math
from dataclasses import dataclass

B_DEFAULT = 200.0        # pN nm^2
TAU_C_NICKED = 27.0      # pN nm
TAU_C_CONTINUOUS = 34.0  # pN nm, midpoint of the measured 31-36 range
RISE = 0.33              # nm per base pair
KT = 4.2                 # pN nm, Zocchi's room-temperature thermal energy


@dataclass(frozen=True)
class Spring:
    """A dsDNA leaf spring of ``n_bp`` base pairs."""

    n_bp: int
    tau_c: float = TAU_C_NICKED
    B: float = B_DEFAULT
    T: float = KT

    @property
    def L(self) -> float:
        """Contour length (nm).  Book p.25: rise 0.33 nm/bp."""
        return RISE * self.n_bp

    @property
    def gamma(self) -> float:
        """Book Eq. 2.188.  The model's own validity condition is gamma < 1."""
        return self.L * self.tau_c / (2.0 * self.B)

    @property
    def R(self) -> float:
        """Radius of curvature in the kinked state, Book Eq. 2.197."""
        return 0.5 * self.L * (1.0 - self.gamma**2 / 90.0)

    @property
    def x_c(self) -> float:
        """End-to-end distance at the kinking (softening) transition, Eq. 2.194."""
        return self.L * (1.0 - self.gamma**2 / 15.0)

    @property
    def x0(self) -> float:
        """Zero-force end-to-end distance, Eq. 2.205.  Below L from thermal bending."""
        return self.L * (1.0 - self.T * self.L / (10.0 * self.B))

    @property
    def E0(self) -> float:
        """Elastic energy stored at the transition, Eq. 2.193 (pN nm)."""
        return (self.L / self.B) * self.tau_c**2 / 6.0

    def is_feasible(self, x: float) -> bool:
        """A spring can only be *compressed*: its contour must exceed the span."""
        return x < self.L

    def force(self, x: float) -> float:
        """Force pushing the attachment points apart, Book Eq. 3.52 (pN).

        ``x`` is the end-to-end distance the protein plus linkers impose.
        """
        if not self.is_feasible(x):
            raise ValueError(
                f"{self.n_bp} bp spring has contour length {self.L:.2f} nm but must "
                f"span x = {x:.2f} nm: it would be stretched, not bent, and the "
                f"linkers rather than the DNA would carry the strain."
            )
        if x >= self.x_c:
            return 10.0 * self.B / self.L**2 - self.T / (self.L - x)
        two_R = 2.0 * self.R
        return self.tau_c / (two_R * math.sqrt(1.0 - (x / two_R) ** 2))

    def energy(self, x: float) -> float:
        """Elastic energy stored in the bent/kinked rod, Book Eq. 2.206 (pN nm)."""
        if not self.is_feasible(x):
            raise ValueError(f"{self.n_bp} bp spring cannot span x = {x:.2f} nm")
        if x >= self.x_c:
            return (-10.0 * self.B / self.L**2) * (x - self.x0) - self.T * math.log(
                (self.L - x) / (self.L - self.x0)
            )
        return self.tau_c * math.acos(x / (2.0 * self.R))

    def unkinked_force(self, x: float) -> float:
        """Linear-elasticity force with kinking suppressed, for comparison only.

        This is the x >= x_c branch extrapolated below x_c.  It plateaus at the
        Euler buckling force 10B/L^2 (Book Eq. 2.169, within 2% of pi^2 B/L^2) and
        is what overestimates the force by ~3-5x -- the error in the 2007 paper.
        """
        return 10.0 * self.B / self.L**2 - self.T / (self.L - x)


def span(d_cb_cb: float, arm: float) -> float:
    """End-to-end distance the spring must span (nm).

    Book p.99: ``s = (distance between attachment residues) + 2 x (crosslinker length)``.
    For the RLuc chimera Zocchi used s = 1.9 + 2 x 2.1 = 6.1 nm.
    """
    return d_cb_cb + 2.0 * arm


def choose_n_bp(x: float, target_force: float, *, tau_c: float = TAU_C_NICKED,
                slack: float = 1.5, n_max: int = 80) -> int:
    """Smallest feasible spring whose force does not exceed ``target_force``.

    ``slack`` (nm) keeps the contour length safely above the span so the rod is
    genuinely bent rather than teetering at the divergence of f(x) as x -> 2R.
    Returns the n_bp whose force is closest to the target from above.
    """
    best, best_err = None, math.inf
    for n in range(4, n_max + 1):
        s = Spring(n, tau_c=tau_c)
        if s.L < x + slack:
            continue
        err = abs(s.force(x) - target_force)
        if err < best_err:
            best, best_err = n, err
    if best is None:
        raise ValueError(f"no spring up to {n_max} bp can span x = {x:.2f} nm")
    return best


def ladder(x: float, n_bps=(20, 22, 24, 26, 28, 30, 35, 40, 50, 60),
           tau_c: float = TAU_C_NICKED) -> str:
    """Human-readable force table at a given span."""
    rows = [f"span x = {x:.2f} nm,  tau_c = {tau_c:.1f} pN nm,  B = {B_DEFAULT:.0f} pN nm^2",
            f"{'n_bp':>5} {'L/nm':>7} {'gamma':>6} {'x_c/nm':>7} {'regime':>9} "
            f"{'f/pN':>8} {'E/kT':>7} {'f_unkinked':>11}"]
    for n in n_bps:
        s = Spring(n, tau_c=tau_c)
        if not s.is_feasible(x):
            rows.append(f"{n:>5} {s.L:7.2f} {s.gamma:6.3f} {s.x_c:7.2f} "
                        f"{'INFEASIBLE':>9} {'--':>8} {'--':>7} {'--':>11}")
            continue
        regime = "kinked" if x < s.x_c else "bent"
        flag = " *" if s.gamma >= 1.0 else ""
        rows.append(f"{n:>5} {s.L:7.2f} {s.gamma:6.3f} {s.x_c:7.2f} {regime:>9} "
                    f"{s.force(x):8.2f} {s.energy(x)/KT:7.2f} {s.unkinked_force(x):11.2f}{flag}")
    rows.append("  * gamma >= 1: outside the model's stated validity range (L > ~45 bp)")
    return "\n".join(rows)


if __name__ == "__main__":
    import argparse

    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--d-cb-cb", type=float, default=3.13,
                   help="Cbeta-Cbeta distance between attachment sites, nm "
                        "(default: 3.13, measured Asp133-Asn149 in 2B3P)")
    p.add_argument("--arm", type=float, nargs="+", default=[1.8, 2.1, 2.4],
                   help="one-sided linker length(s), nm")
    p.add_argument("--target-force", type=float, default=7.0)
    args = p.parse_args()

    for arm in args.arm:
        x = span(args.d_cb_cb, arm)
        print(f"\n=== arm = {arm:.2f} nm  ->  x = {x:.2f} nm ===")
        for tau in (TAU_C_NICKED, TAU_C_CONTINUOUS):
            print(ladder(x, tau_c=tau))
            print()
        n = choose_n_bp(x, args.target_force)
        print(f"strong spring targeting {args.target_force:.1f} pN: "
              f"{n} bp -> {Spring(n).force(x):.2f} pN (nicked)")
