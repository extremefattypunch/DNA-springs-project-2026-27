#!/usr/bin/env python3
"""Zocchi's model, evaluated for this construct: the spring, the kink, and the protein.

The sfGFP page sets every measurement against the theory the experiment was designed
on -- Zocchi's kinked-rod spring (Molecular Machines 2018, eqs 2.188-2.206 and 3.51-3.52;
Tseng & Zocchi JACS 2013) -- and this page does the same, using the repository's own
implementation (``sfGFP/analysis/spring_model.py``, which reproduces the book's RLuc and
GK numbers).  Three things are computed at the span the simulations actually settle at:

1. **The spring.**  f(x), the regime (bent or kinked), the softening point x_c, the
   stored energy and the Euler load, for a duplex of 27, 40 and 50 bp -- the duplex that
   *would* push at this span.  The chimera carries a single strand, which is Zocchi's
   "ss state: essentially no stress on the enzyme"; the Marko-Siggia tension is quoted
   beside it so the page shows, at the same x, what each model says.
2. **The kink.**  The critical bending torque tau_c sets the curvature at which a duplex
   kinks, tau_c / B = 27/200 nm^-1 = 2.6 deg/bp (3.2 deg/bp at tau_c = 34), against which
   the strand's own curvature and the sfGFP windows' empirical opening thresholds
   (5.5 deg/bp onset, 8.7 deg/bp certain) are set.
3. **The protein spring.**  Zocchi's two-spring picture (book pp. 117-119) gives the
   protein a stiffness kappa ~ 100 kT/nm^2 along the attachment coordinate and a 4 pN
   yield; sfGFP measured 230-350 kT/nm^2 with a force clamp.  Without a clamp ladder,
   equipartition on the Cb-Cb fluctuation, kappa = kT / var(x), gives the horseshoe's
   stiffness along the same coordinate from each run, and with it the displacement the
   modelled forces would produce.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np
import pandas as pd

from common import SYSTEMS, analysis_prefix
from spring_model import (B_DEFAULT, KT, RISE, TAU_C_CONTINUOUS, TAU_C_NICKED, Spring, SsSpring,
                          XI_SS)

HERE = Path(__file__).resolve().parent
KT_PN_NM = KT                      # 4.2 pN nm, Zocchi's room-temperature unit
SFGFP_THRESHOLDS = {"onset_deg_per_bp": 5.5, "certain_deg_per_bp": 8.7,
                    "chimera_deg_per_bp": "4.8-5.2", "curuksu_typeII_deg_per_bp": "8.6-9.3",
                    "source": "sfGFP page, bare-duplex restraint windows (docs/wiki/Home.md)"}


def critical_curvature_deg_per_bp(tau_c, B=B_DEFAULT):
    """The curvature at which a rod of bending modulus B yields at torque tau_c: kappa = tau_c/B."""
    return math.degrees(tau_c / B * RISE)


def spring_table(x_nm):
    rows = []
    for n in (27, 40, 50):
        for tau, kind in ((TAU_C_NICKED, "nicked"), (TAU_C_CONTINUOUS, "continuous")):
            sp = Spring(n, tau_c=tau)
            row = {"n_bp": n, "tau_c_pN_nm": tau, "kind": kind, "contour_nm": round(sp.L, 2),
                   "gamma": round(sp.gamma, 3), "in_validity_range": bool(sp.gamma < 1.0),
                   "x_c_nm": round(sp.x_c, 2), "euler_pN": round(10.0 * sp.B / sp.L ** 2, 2),
                   "feasible": sp.is_feasible(x_nm)}
            if sp.is_feasible(x_nm):
                row.update({"regime": "kinked" if x_nm < sp.x_c else "bent",
                            "force_pN": round(sp.force(x_nm), 2),
                            "unkinked_force_pN": round(sp.unkinked_force(x_nm), 2),
                            "energy_kT": round(sp.energy(x_nm) / KT, 2)})
            rows.append(row)
    return rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=str(HERE / "zocchi.json"))
    a = ap.parse_args()
    out = {"sources": {
        "spring": "Zocchi, Molecular Machines (2018) eqs 2.188-2.206, 3.51-3.52; Tseng & Zocchi JACS 135:11879 (2013) eq 6",
        "parameters": f"B = {B_DEFAULT:.0f} pN nm^2, tau_c = {TAU_C_NICKED:.0f} (nicked) / {TAU_C_CONTINUOUS:.0f} (continuous) pN nm, rise {RISE} nm/bp, kT = {KT} pN nm",
        "ssDNA": f"Marko-Siggia with xi = {XI_SS} nm (Zocchi 2021 chapter), 0.59 nm/nt (Smith, Cui & Bustamante 1996)",
        "protein": "book pp. 117-119: kappa ~ 100 kT/nm^2, yield ~4 pN at 3 A (eq. 3.100); sfGFP clamp ladder 230-350 kT/nm^2",
        "implementation": "sfGFP/analysis/spring_model.py (reproduces RLuc 1.5 pN, 40-mer 2.4 pN, GK 6.6 kT)"}}

    # --- the span: the DNA end-to-end the chimera settles at, else the build target ---
    span = {}
    m = Path(str(analysis_prefix("P3_ss50")) + "_mechanics.csv")
    if m.exists():
        d = pd.read_csv(m)
        if "spring_end_to_end_A" in d:
            x = d["spring_end_to_end_A"].to_numpy()
            span = {"source": "P3_ss50 DNA end-to-end (5'P to 3'O3'), simulated", "A": float(x.mean()),
                    "sd_A": float(x.std()), "p5_A": float(np.percentile(x, 5)), "p95_A": float(np.percentile(x, 95))}
    if not span:
        try:
            att = json.loads((HERE.parent / "build/04_assemble/pp2a_ss50.json").read_text())
            span = {"source": "build target (assemble_chimera.py anchor span); the chimera has not run",
                    "A": float(att["duplex_anchor_span_A"]), "sd_A": 0.0}
        except Exception:                                         # noqa: BLE001
            span = {"source": "default", "A": 90.0, "sd_A": 0.0}
    x_nm = span["A"] / 10.0
    out["span"] = span

    # --- 1. the spring at this span ---
    out["duplex_springs_at_span"] = spring_table(x_nm)
    ss = SsSpring(50)
    s = min(x_nm / ss.L0, 0.999)
    out["ssDNA_at_span"] = {"n_nt": 50, "contour_nm": round(ss.L0, 2), "natural_rms_nm": round(ss.rms_end_to_end, 3),
                            "state": ss.state(x_nm), "tension_pN_clipped": round(ss.tension(x_nm), 3),
                            "tension_pN_unclipped": round((KT / XI_SS) * (1 / (4 * (1 - s) ** 2) - 0.25 + s), 3),
                            "zocchi_2021": "Fig. 1C: with the DNA spring in the ss state there is essentially no stress on the enzyme"}
    # the force curve for the figure
    xs = np.linspace(2.0, 17.0, 301)
    curves = {"x_A": (xs * 10).round(2).tolist()}
    for n in (27, 40, 50):
        for tau, kind in ((TAU_C_NICKED, "nicked"), (TAU_C_CONTINUOUS, "continuous")):
            sp = Spring(n, tau_c=tau)
            curves[f"f_{n}bp_{kind}_pN"] = [round(sp.force(v), 3) if sp.is_feasible(v) else None for v in xs]
            curves[f"f_{n}bp_{kind}_unkinked_pN"] = [round(sp.unkinked_force(v), 3) if sp.is_feasible(v) else None for v in xs]
    curves["ss50_pull_pN"] = [round(ss.tension(v), 3) if v < ss.L0 else None for v in xs]
    pd.DataFrame(curves).to_csv(HERE / "zocchi_curves.csv", index=False)

    # --- 2. the kink ---
    kink = {"critical_curvature_deg_per_bp": {"nicked_27": round(critical_curvature_deg_per_bp(TAU_C_NICKED), 2),
                                              "continuous_34": round(critical_curvature_deg_per_bp(TAU_C_CONTINUOUS), 2)},
            "softening_points": {f"{n} bp": {"x_c_A": round(Spring(n).x_c * 10, 1),
                                             "half_rod_chord_A": round(Spring(n).R * 10, 1)   # book eqs 2.189/2.197: the kink-to-anchor chord, ~L/2}
                                 for n in (27, 40, 50)},
            "sfGFP_empirical": SFGFP_THRESHOLDS}
    for n in (27, 40, 50):
        sp = Spring(n)
        # a uniformly bent rod of contour L with chord x has arc angle theta: x = (L/theta)*2 sin(theta/2)
        def chord(theta, L=sp.L):
            return L if theta < 1e-9 else (L / theta) * 2 * math.sin(theta / 2)
        lo, hi = 1e-6, 2 * math.pi - 1e-6
        for _ in range(80):
            mid = 0.5 * (lo + hi)
            if chord(mid) > sp.x_c:
                lo = mid
            else:
                hi = mid
        kink["softening_points"][f"{n} bp"]["curvature_at_x_c_deg_per_bp"] = round(math.degrees(lo / n), 2)
        theta = None
        if x_nm < sp.L:
            lo2, hi2 = 1e-6, 2 * math.pi - 1e-6
            for _ in range(80):
                mid = 0.5 * (lo2 + hi2)
                if chord(mid) > x_nm:
                    lo2 = mid
                else:
                    hi2 = mid
            theta = lo2
        kink["softening_points"][f"{n} bp"]["curvature_at_span_deg_per_bp"] = (round(math.degrees(theta / n), 2) if theta else None)
    ssj = Path(str(analysis_prefix("P3_ss50")) + "_ssdna.json")
    if ssj.exists():
        k = json.loads(ssj.read_text()).get("kink", {})
        kink["strand_measured"] = {"axis_curvature_deg_per_nt": k.get("axis_curvature_deg_per_nt", {}).get("mean"),
                                   "sharp_turns_per_frame": k.get("sharp_turns_per_frame", {}).get("mean"),
                                   "hottest_nt": k.get("hottest_nt"), "hottest_step": k.get("hottest_step"),
                                   "bii_fraction": k.get("bii_fraction", {}).get("mean")}
    out["kink"] = kink

    # --- 3. the protein spring, by equipartition on Cb-Cb ---
    prot = {"zocchi_generic_kT_nm2": 100, "sfGFP_measured_kT_nm2": "230-350",
            "method": "kappa = kT / var(Cb-Cb), from each run's own fluctuation (equipartition; slow drift inflates the variance, so this is "
                      "a lower bound on the local stiffness; read it against the block SEM)"}
    for s in SYSTEMS:
        m = Path(str(analysis_prefix(s)) + "_mechanics.csv")
        if not m.exists():
            continue
        x = pd.read_csv(m)["anchor_CB_CB_A"].to_numpy() / 10.0
        var = float(x.var())
        kappa = 1.0 / var                                   # kT / nm^2
        prot[s] = {"mean_A": round(float(x.mean()) * 10, 2), "sd_A": round(float(x.std()) * 10, 2),
                   "kappa_kT_nm2": round(kappa, 2), "kappa_pN_per_nm": round(kappa * KT_PN_NM, 2),
                   "displacement_per_pN_A": round(10.0 / (kappa * KT_PN_NM), 2),
                   "displacement_at_2pN_A": round(2.0 * 10.0 / (kappa * KT_PN_NM), 2)}
    out["protein_spring"] = prot

    Path(a.out).write_text(json.dumps(out, indent=2))
    print(f"wrote {a.out}: span {span['A']:.1f} A ({span['source']}); "
          f"50 bp nicked would push {next((r['force_pN'] for r in out['duplex_springs_at_span'] if r['n_bp'] == 50 and r['kind'] == 'nicked' and r.get('force_pN')), 'n/a')} pN; "
          f"ssDNA {out['ssDNA_at_span']['state']} {out['ssDNA_at_span']['tension_pN_clipped']} pN; "
          f"kink at {kink['critical_curvature_deg_per_bp']['nicked_27']}-{kink['critical_curvature_deg_per_bp']['continuous_34']} deg/bp; "
          + ", ".join(f"{s}: kappa {v['kappa_kT_nm2']} kT/nm2" for s, v in prot.items() if isinstance(v, dict)))


if __name__ == "__main__":
    main()
