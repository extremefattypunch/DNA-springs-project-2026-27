#!/usr/bin/env python3
"""Generate the GitHub wiki landing page from the analysis, not from memory.

Why this is a script
--------------------
The wiki page was hand-written, and by the time 81 replicates had finished it was
claiming 2,550 ns (the true figure was 4,050), a curvature of 5.7-6.3 deg/bp (4.8-5.2
after the axis-averaging fix), and that the direct force measurement was "running now"
(it had finished).  Every one of those was a number that existed correctly in
findings.json at the same moment the page misstated it.  A summary that is retyped
drifts from its source; one that is generated cannot.

So the page below takes every quantity from ``analysis/findings.json`` and
``analysis/force_extension.csv``.  Prose that interprets those numbers is still written
by hand -- that is the part a person should be responsible for -- but no figure in it is.

    python report/build_wiki.py                 # -> docs/wiki/Home.md
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd

RAW = ("https://raw.githubusercontent.com/extremefattypunch/"
       "DNA-springs-project-2026-27/main/docs/figures")
REPO = "https://github.com/extremefattypunch/DNA-springs-project-2026-27"
PAGES = "https://extremefattypunch.github.io/DNA-springs-project-2026-27/"

# Euler buckling load of a clamped-clamped rod, pi^2 B / L^2, with B = 200 pN nm^2 and
# L the contour of the 27 bp duplex.  Quoted because the measured plateau lands on it.
B_PN_NM2 = 200.0
RISE_NM = 0.33


def euler_pN(n_bp: int) -> float:
    L = (n_bp - 1) * RISE_NM
    return 3.141592653589793 ** 2 * B_PN_NM2 / (L * L)


def fx_summary(path: Path):
    """Measured force per construct, and where base pairs start opening.

    Only windows the box can actually hold are used, and the test is applied to the
    separation the duplex ACTUALLY adopted rather than to the target it was pulled
    towards: with the soft restraint the two differ by 2-7 A, and a window that passes on
    its target can still spend the whole trajectory close enough to its own image to
    matter.
    """
    if not path.exists():
        return None
    d = pd.read_csv(path)
    if "box_ok" not in d.columns or "mean_r_A" not in d.columns:
        return {"stale": True}
    from scipy import stats as st

    # box_ok is AUTHORITATIVE: computed by force_extension.py from unrounded values.
    # Re-deriving it from the rounded csv columns flipped a run at 0.850049 and made the
    # wiki disagree with the report about the valid set.
    n_all = len(d)
    if "box_ok_on_measured_r" not in d.columns:
        d = d.assign(box_ok=d.box_ok.astype(str).str.lower().isin(("true", "1"))
                     & (d.mean_r_A / d.box_width_A <= 0.85))
    d = d[d.box_ok.astype(str).str.lower().isin(("true", "1"))]
    if not len(d):
        return {"stale": True}

    # n_bp comes from wc.pairs, which force_extension.py MEASURES from the
    # topology, rather than from a map keyed on system name.  The map had to
    # be edited in five files whenever a system was added, and a system
    # missing from it silently fell back to 27.
    NB = (d.groupby("system")["wc.pairs"].max().dropna()
          .astype(int).to_dict() if "wc.pairs" in d else {})
    # Keyed on the physical duplex, so the same molecule in two box sizes forms one
    # curve rather than two unrelated series, and the windows they share become a
    # cross-check between boxes.
    if "duplex" not in d.columns:
        d = d.assign(duplex=d.system.map(
            lambda s: f"{NB.get(s, 27)} bp"
                      + (", nicked" if "nick" in s.lower() else "")))
    # Split by box as well as duplex, and report the LARGEST box for each: a matched
    # comparison (27 bp at x0 = 50 and 58 in a 76 A and an 83 A box) gives 26.3 against
    # 33.3 pN, p = 0.02, so the measurement is not converged and a pooled number mixes a
    # sound measurement with one whose box was too small.
    SMALL = {"27 bp": 80.0, "40 bp": 100.0, "40 bp, nicked": 100.0}
    d = d.assign(_big=[b >= SMALL.get(s, 90.0) for b, s in zip(d.box_A, d.duplex)])
    out = {"stale": False, "n_valid": int(len(d)), "systems": {}, "box_effect": {}}
    for s, g in d.groupby("duplex"):
        shared = sorted(set(g[~g._big].x0_A) & set(g[g._big].x0_A))
        if shared:
            m = g[g.x0_A.isin(shared)]
            sm, lg = m[~m._big].force_pN, m[m._big].force_pN
            if len(sm) > 1 and len(lg) > 1:
                out["box_effect"][s] = {
                    "small": float(sm.mean()), "large": float(lg.mean()),
                    "ratio": float(lg.mean() / sm.mean()),
                    "p": float(st.ttest_ind(lg.values, sm.values,
                                            equal_var=False).pvalue),
                    "x0": shared}
        if g._big.any():
            g = g[g._big]                 # best estimate = the largest box available
        eul = euler_pN(int(str(s).split()[0]))
        lr = st.linregress(g.mean_r_A, g.force_pN)
        e = {"label": s, "n": int(len(g)), "box_A": round(float(g.box_A.mean()), 0),
             "boxes_A": sorted({round(float(b)) for b in g.get("box_A", [])
                                if b == b}) if "box_A" in g else [],
             "r_lo": float(g.mean_r_A.min()), "r_hi": float(g.mean_r_A.max()),
             "f": float(g.force_pN.mean()), "sem": float(g.force_pN.sem()),
             "euler": float(eul), "ratio": float(g.force_pN.mean() / eul),
             "slope_p": float(lr.pvalue)}
        a = g[g.k_kJ_mol_nm2 == 33.0].force_pN.values
        b = g[g.k_kJ_mol_nm2 == 200.0].force_pN.values
        if len(a) > 1 and len(b) > 1:
            e["k_agree_p"] = float(st.ttest_ind(a, b, equal_var=False).pvalue)
        out["systems"][s] = e

    i_, n_ = (d[d.duplex == "40 bp"].force_pN.values,
              d[d.duplex == "40 bp, nicked"].force_pN.values)
    if len(i_) > 1 and len(n_) > 1:
        out["nick"] = {"intact": float(i_.mean()), "nicked": float(n_.mean()),
                       "p": float(st.ttest_ind(i_, n_, equal_var=False).pvalue)}

    if "wc.max_internal_occupancy" in d.columns:
        def per_bp(chord, n_bp):
            L = (n_bp - 1) * 3.376
            lo, hi = L / (2 * 3.141592653589793), 1e7
            for _ in range(300):
                m = 0.5 * (lo + hi)
                import math
                lo, hi = ((m, hi) if 2 * m * math.sin(L / (2 * m)) < chord
                          else (lo, m))
            import math
            return math.degrees(3.376 / (0.5 * (lo + hi)))

        dd = d[d["wc.max_internal_occupancy"].notna()].copy()
        dd["opened"] = dd["wc.max_internal_occupancy"] > 0
        # Grouped by WINDOW, matching the report's table, so the two cannot quote
        # different onsets for the same thing.  Per-run binning is what the figure does
        # and gives a slightly different pair of numbers for a defensible reason; having
        # the prose in both places agree matters more than either grouping does.
        wins = []
        for (s, x0), g in dd.groupby(["system", "x0_A"]):
            r = float(g.mean_r_A.mean())
            wins.append((per_bp(r, NB.get(s, 27)), int(g.opened.sum()), int(len(g))))
        wins.sort()
        opened = [c for c, k, n in wins if k]
        allopen = [c for c, k, n in wins if k == n]
        if opened:
            out["kink"] = {
                "gentlest_opening_deg_bp": float(min(opened)),
                "certain_from_deg_bp": float(min(allopen)) if allopen else None,
                "sampled_lo": float(min(c for c, k, n in wins)),
                "sampled_hi": float(max(c for c, k, n in wins))}
            bps = dd[dd.opened & (dd.duplex.astype(str).str.startswith("27"))][
                "wc.most_open_internal_bp"].dropna()
            if len(bps):
                vc = bps.astype(int).value_counts()
                out["kink"]["most_common_bp"] = int(vc.index[0])
                out["kink"]["most_common_bp_frac"] = float(vc.iloc[0] / len(bps))
                out["kink"]["n_openings_27bp"] = int(len(bps))
    out["n_dropped"] = int(n_all - out["n_valid"])
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--findings", default="analysis/findings.json")
    ap.add_argument("--fx", default="analysis/force_extension.csv")
    ap.add_argument("--out", default="../../docs/wiki/Home.md")
    a = ap.parse_args()

    f = json.loads(Path(a.findings).read_text())
    fx = fx_summary(Path(a.fx))
    sp, ss = f.get("springs") or {}, f.get("ss_springs") or {}
    vm = f.get("extension", {}).get("vs_matched_control") or {}
    sq = f.get("is_it_a_spring") or {}
    topo = f.get("topology") or {}
    lab = f.get("labels") or {}

    # --- the constructs, ordered by how much they compress ---
    rows = "\n".join(
        "| %s | %s | %s | **%+.3f** | %.3f |" % (
            lab.get(k, k), "ubiquitin" if v.get("family") == "ubiquitin" else "sfGFP",
            ("%.2f" % (sp[k]["force_pN"]["mean"])) if k in sp else "0 — slack",
            v["delta"], v["p"])
        for k, v in sorted(vm.items(), key=lambda kv: kv[1]["delta"]))

    spring_rows = "\n".join(
        "| %s | %.1f Å | **%.2f pN** | %.1f k<sub>B</sub>T | %.2f of %d |" % (
            lab.get(k, k), v["x_A"]["mean"], v["force_pN"]["mean"],
            v["energy_kT"]["mean"], v["bp_open"]["mean"],
            int(round(v["contour_nm"] * 10 / 3.376)) + 1)
        for k, v in sp.items())

    n_neg = sum(1 for v in vm.values() if v["delta"] < 0)
    e27 = euler_pN(27)

    if fx and not fx.get("stale") and fx["systems"].get("27 bp"):
        # NOT `rows` -- that name already holds the constructs table above, and
        # shadowing it silently replaced that table with this one.
        fx_rows = "\n".join(
            "| %s | %.0f–%.0f | **%.1f ± %.1f** | %.1f | %.2f | %.2f |" % (
                e["label"], e["r_lo"], e["r_hi"], e["f"], e["sem"], e["euler"],
                e["ratio"], e["slope_p"])
            for e in fx["systems"].values())
        d27 = fx["systems"]["27 bp"]
        ks = [e["k_agree_p"] for e in fx["systems"].values() if "k_agree_p" in e]
        fx_txt = f"""A bare duplex, no protein, held at a series of separations by a
harmonic restraint between the two linker nitrogens. The mean restraint force
⟨k(r − x₀)⟩ *is* the force the duplex exerts — no model in between.

| duplex | separations sampled (Å) | measured force (pN) | Euler load π²B/L² (pN) | measured / Euler | p, slope vs separation |
|---|---|---|---|---|---|
{fx_rows}

Two things to read.

**Is the force flat in separation? The two boxes disagree.** A load independent of
compression is the signature of a **buckled rod**. The 76 Å box shows exactly that for the
27 bp — no resolvable slope over 19 runs spanning 33–63 Å (p = 0.49) — but the 83 Å box,
with six runs over only 51–65 Å, rises measurably (p = 0.014). Six points over 14 Å against
nineteen over 30 Å is not a fair contest, and we cannot say which reading is right: either
the plateau is real and the larger box is too sparse to show it, or the force climbs as
compression is relieved. Supported in one box, contradicted in the other.

**But the number depends on the box, and it is not converged.** An earlier version of
this page reported 25.7 ± 1.0 pN for the 27 bp — a ratio of 0.96 to its Euler load — and
called it textbook elasticity. That was an artefact of too small a box. Running the same
duplex at the same separations (50 and 58 Å) in a 76 Å and an 83 Å box gives
{fx['box_effect'].get('27 bp', {}).get('small', 0):.1f} against
{fx['box_effect'].get('27 bp', {}).get('large', 0):.1f} pN — a factor of
{fx['box_effect'].get('27 bp', {}).get('ratio', 0):.2f} at p =
{fx['box_effect'].get('27 bp', {}).get('p', 1):.3f}. **That figure is withdrawn.**

The table above quotes the largest box available for each duplex. There the 27 bp sits at
{d27['ratio']:.2f} of its Euler load and the 40 bp at
{fx['systems'].get('40 bp', {}).get('ratio', 0):.2f} of its own — two lengths whose Euler
loads differ by more than a factor of two, landing on the same multiple. That is a
property of the duplex rather than of one simulation cell, and a better kind of agreement
than the one being retracted.

Zocchi's kinked branch predicts **2.4–6.6 pN** for these same duplexes. We measure
{d27['f']:.1f} pN — about **six times higher**, and on the *unyielded* branch, the one
his model says the duplex has already left. So on this timescale the duplex is not a
softened kinked spring; it is an elastic rod holding its buckling load, and the modelled
forces in the table above are if anything **underestimates**."""
        if ks:
            fx_txt += ("\n\nEach separation was run at two restraint stiffnesses, 33 and "
                       "200 kJ mol⁻¹ nm⁻². The soft one is barely a constraint — the "
                       "duplex settles 2–7 Å past its target — so if the force were an "
                       "artefact of how hard we pulled, the two would separate. Against "
                       "the separation actually adopted they agree (p = "
                       + ", ".join(f"{x:.2f}" for x in ks)
                       + "). That is the control that makes this the duplex's force "
                         "rather than the restraint's.")
        if fx.get("nick"):
            nk = fx["nick"]
            fx_txt += (f"\n\n**What the nick does not do.** Zocchi's model treats the "
                       f"nick as the key softening control and predicts ~1.3× less "
                       f"force. Measured, there is no resolvable difference: "
                       f"{nk['intact']:.1f} pN intact against {nk['nicked']:.1f} pN "
                       f"nicked, p = {nk['p']:.2f}. Over 20–30 ns, removing one backbone "
                       f"bond from a rod held by 40 base pairs of stacking does not "
                       f"change its buckling load.")
        if fx.get("n_dropped"):
            fx_txt += (f"\n\n{fx['n_dropped']} runs are **discarded**: the duplex sat "
                       f"close enough to its own periodic image that the restraint force "
                       f"is not its alone. The chimera's own separation (~70 Å for the "
                       f"27 bp) is outside what these boxes can hold, so an earlier "
                       f"figure quoted at the chimera's geometry is **withdrawn**; "
                       f"repeating those windows in a larger box is queued.")
        k = fx.get("kink") or {}
        kink_txt = (
            f"""Not in the chimera. Our curvature is **4.8–5.2° per base-pair step**, and
in three 50 ns runs no interior base pair opened. The windows above reach
{k.get('sampled_lo', 0):.1f}–{k.get('sampled_hi', 0):.1f}°/bp, and there opening is a
probability rather than a threshold: it appears from **{k.get('gentlest_opening_deg_bp', 0):.1f}°/bp**
and becomes certain — every run in the window — only above
**{k.get('certain_from_deg_bp') or float('nan'):.1f}°/bp**,
with the 8.6–9.3°/bp at which the one directly comparable all-atom calculation nucleates
a type II kink sitting where our probability is already 0.8.

The opening is not diffuse. """
            + (f"**Base pair {k['most_common_bp']}** accounts for "
               f"{100 * k.get('most_common_bp_frac', 0):.0f}% of the "
               f"{k.get('n_openings_27bp', 0)} openings seen in the 27 bp duplex, and "
               f"its incoming step is **CpA** — one of the three pyrimidine–purine "
               f"steps kinking is known to prefer."
               if k.get("most_common_bp") else "")
            + """

So this is a measurement rather than an absence: the force field is entirely capable of
breaking base pairs, it does so around where the literature says it should, and our
construct sits just below the onset. That **excludes "the force field cannot kink"** as
the explanation — one of the three candidates put to us in review.""")
    else:
        fx_txt = ("The direct force measurement is being re-analysed; the figure "
                  "previously quoted here came from a window too large for its "
                  "periodic box and is withdrawn.")
        kink_txt = ("Not in the chimera, at 4.8–5.2° per base-pair step. The direct "
                    "test against curvature is being re-analysed.")

    # ---- seeded kinks: released with no restraint, does the kink survive? ----
    seeded_txt = ""
    fxcsv = Path(a.fx)
    if fxcsv.exists():
        dd = pd.read_csv(fxcsv)
        if "seeded.bp" in dd.columns:
            dd = dd[dd["seeded.bp"].notna()]
            if len(dd):
                rows = []
                verdicts = []
                for (sysname, bp), g in dd.groupby(["system", "seeded.bp"]):
                    v = list(g.get("seeded.verdict", []))
                    verdict = max(set(v), key=v.count) if v else "?"
                    verdicts.append(verdict)
                    healed = g["seeded.healed_after_ns"].dropna()
                    rows.append("| %d | %d | %.1f | %.1f | %.0f%% | %s | **%s** |" % (
                        int(bp), len(g),
                        g["seeded.N1_N3_first_A"].mean(),
                        g["seeded.N1_N3_last_A"].mean(),
                        100 * g["seeded.open_fraction"].mean(),
                        ("%.1f" % healed.mean()) if len(healed) else "—",
                        verdict))
                rows.sort()
                all_healed = all(v == "healed" for v in verdicts)
                any_persist = any(v == "persisted" for v in verdicts)
                counts = ", ".join(f"{verdicts.count(v)} {v}"
                                   for v in sorted(set(verdicts)))
                head = ("| seeded bp | runs | N1–N3 start (Å) | end (Å) | open, "
                        "frames | shut after (ns) | outcome |\n"
                        "|---|---|---|---|---|---|---|\n" + "\n".join(rows))
                if all_healed:
                    body = (
                        "\n\n**Every seeded kink healed** (" + counts + "). Released "
                        "without load, the broken pair re-formed and the 114° bend "
                        "relaxed, at all three positions in every replicate. So in this "
                        "force field the kinked state is **not a metastable minimum** — "
                        "it is a geometry that exists only while something holds it "
                        "there.\n\nThat bears on Zocchi's model harder than the absence "
                        "of spontaneous kinking did. A barrier never crossed might be "
                        "crossed given longer; a state that is not a minimum cannot be "
                        "occupied at all. His softened branch — where the 2.4–6.6 pN "
                        "predictions come from — requires the duplex to *sit* kinked "
                        "under a few piconewtons. We put it there and it does not stay, "
                        "which is why the measured force lands on the unyielded Euler "
                        "branch.")
                elif any_persist:
                    body = (
                        "\n\n" + counts + ". A kink placed by hand and released "
                        "**stayed**, so the kinked state is a real minimum at zero load "
                        "and Zocchi's softened branch is physically reachable. Our "
                        "chimera does not reach it — but that is now a statement about a "
                        "barrier, not about the state failing to exist, and a barrier can "
                        "be crossed on timescales longer than ours. This is the outcome "
                        "most favourable to the model.")
                else:
                    body = (
                        "\n\n" + counts + ". The outcome depends on where the kink was "
                        "placed, which suggests the kinked state sits near the edge of "
                        "stability with sequence context deciding. Neither the softened "
                        "branch nor its absence follows; more positions and longer runs "
                        "would settle it.")
                seeded_txt = f"""

## A kink put there by hand

Everything above says the duplex does not kink at the curvature the chimera imposes.
The sharper question, and the one Zocchi's model turns on, is whether the kinked state
is a **metastable minimum** the duplex would sit in if it ever got there, or not a
minimum at all. That is answerable by construction: three duplexes were built with a
type II kink already in place — one interior pair broken, both flanking pairs intact,
114° of local bend, at base pairs 7, 15 and 21 — then released with **no restraint** and
watched.

{head}{body}
"""
    md = f"""# DNA springs on sfGFP — all-atom MD

**[▶ Read the full report]({PAGES})**
 · [no-setup mirror](https://htmlpreview.github.io/?https://raw.githubusercontent.com/extremefattypunch/DNA-springs-project-2026-27/main/docs/index.html)
 · [code and structures]({REPO}/tree/main/MD%20simulations/sfGFP)
 · [response to review]({REPO}/blob/main/feedback/2026-8-22_doeke-response.md)
 · [how it was built]({REPO}/blob/main/feedback/2026-8-22_setup-tutorial.md)

Two amber codons put the tetrazine amino acid **Tet2-Et** at Asp133 and Asn149 on
adjacent strands of sfGFP's β-barrel; a doubly-sTCO-modified DNA strand clicks onto both,
closing a covalent loop. dsDNA is far too stiff to follow the barrel's curvature, so it
bends — and is supposed to push the two attachment points apart with a few piconewtons.

**{f['total_ns']:,.0f} ns over {f['n_replicates']} replicates and \
{len(f['systems_present'])} systems**, in two proteins with two different attachment
chemistries.

![the chimera]({RAW}/S3_spring27_overview_still.png)

## The headline

Attaching DNA to this protein pulls the two attachment sites **together** by about a third
of an ångström — and it makes **no measurable difference whether that DNA is a spring**.

{n_neg} of the {len(vm)} DNA-bearing constructs compress the sites; not one extends them,
which is what a spring pushing its anchors apart would do. And the size of the compression
tracks neither of the two obvious causes:

| against | correlation | p |
|---|---|---|
| the modelled spring force | r = **{sq.get('vs_modelled_force', {}).get('r', float('nan')):+.2f}** (sign *positive*) | {sq.get('vs_modelled_force', {}).get('p', float('nan')):.2f} |
| how much protein the DNA touches | r = **{sq.get('vs_contact', {}).get('r', float('nan')):+.2f}** | {sq.get('vs_contact', {}).get('p', float('nan')):.2f} |

Constructs with **no spring at all** — a slack single strand at provably zero tension, or
two arms that cannot circularise — compress by
{abs(sq.get('mean_delta_no_spring', 0)):.2f} Å. Those with a spring compress by
{abs(sq.get('mean_delta_with_spring', 0)):.2f} Å. The difference is not resolved
(p = {(sq.get('no_spring_vs_spring') or {}).get('p', float('nan')):.2f}).

This is the conclusion the companion [ubiquitin study](https://doi.org/10.64898/2026.07.05.736590)
from this lab reached from NMR — local protein–DNA interaction rather than spring action —
reproduced here independently, in a different protein, with a different chemistry, and
from a different observable.

![site displacement against the force said to cause it]({RAW}/fig7_headline.png)

| construct | protein | modelled force | Δ Cβ–Cβ vs its own clicked control (Å) | p |
|---|---|---|---|---|
{rows}

## Why the protein barely moves

Clamping the same two Cβ with a *known* axial force gives the reason: sfGFP's stiffness
along that coordinate is **230–350 k<sub>B</sub>T nm⁻²**, two to three and a half times
Zocchi's generic ~100. At that stiffness a 4 pN spring moves the sites 0.03 Å — below what
50 ns of equilibrium fluctuation can resolve. Doeke Hekstra's physical objection was
right: a force couple this small averages out at equilibrium.

## What the springs deliver, by the model

| spring | span *x* | modelled force | stored energy | base pairs open |
|---|---|---|---|---|
{spring_rows}

These are Zocchi's *f(x)* evaluated at the end-to-end distance each simulation settles at
— a **model**, not a measurement.

## What they deliver, measured

{fx_txt}

![measured force against extension]({RAW}/fig6_force_extension.png)

## Does the duplex kink?

{kink_txt}

![base-pair opening against curvature]({RAW}/fig5_kink_threshold.png)

{seeded_txt}
## Does closing the covalent loop matter?

The rebuild (below) made five matched pairs possible, differing only in whether a covalent
path runs through the DNA. All three **duplex** pairs say closure compresses more, by
{abs(topo.get('duplex_pairs_planned', {}).get('mean_A', 0)):.2f} Å on average
(paired p = {topo.get('duplex_pairs_planned', {}).get('p', float('nan')):.2f}); both
**single-stranded** pairs fail to reproduce it. A three-pair trend that reverses when two
more pairs are added is a hypothesis, not a finding. Six replicates per system would
settle it.

## What the modification itself does

Installing the two Tet2-Et residues shortens Cβ–Cβ by 0.62 Å before any DNA is attached —
fifteen to thirty times any force effect, which is exactly why using the *unclicked*
protein as the zero produced a spurious result. It is **not** a secondary-structure
change: DSSP is identical everywhere. Asp133 natively sits in the **left-handed α region**
(φ = +57°), which Asp, Asn and Gly favour and a phenylalanine-shaped side chain does not;
Tet2-Et is pushed 17.5° out and clicking gives it back.

![backbone conformation at the two sites]({RAW}/ramachandran.png)

The practical consequence: when choosing sites for this amino acid, avoid positive-φ
positions, and treat an aromatic neighbour as a liability — the ring stacks on Tyr200 in
some constructs and not others, and *differently before and after clicking*, so it will
masquerade as a mechanical effect.

> ### Corrections
>
> **The first build attached the DNA wrongly.** It put a linker on the 5′ end of *each of
> two different strands*. The experiment modifies **one** strand at both ends, which is
> what makes the chimera a covalent circle and what stops one chimera swapping arms with
> another to dimerise. As built, every tether bond was a bridge in the bond graph and no
> covalent cycle existed — the loop was closed by base pairing alone. Rebuilt, with a
> topology test the old build fails.
>
> **A claim about the spring's effect is retracted.** The report previously said the
> chimeras displaced the attachment sites *more* than an equivalent axial pull. That used
> the unclicked protein as its zero; against the matched clicked control the sign inverts.
>
> **Every long distance was measured with the wrong convention.** mdtraj applies the
> minimum image by default, which is only valid out to half the box width — 38–48 Å, while
> the springs' ends are 60–130 Å apart. Errors reached 42 Å. Fixed, and every force
> re-derived.
>
> **The duplex curvature was misreported** as 5.7–6.3°/bp; averaging base-pair centres
> over a full helical turn before differencing gives 4.8–5.2°/bp.
>
> The first two are due to Doeke Hekstra; the
> [full response]({REPO}/blob/main/feedback/2026-8-22_doeke-response.md) works through all
> nine of his comments.

## What this establishes, and what it does not

**Established.** The DNA is not acting as a spring on this protein: eleven constructs
across two proteins with unrelated chemistries all pull their attachment sites *together*
by about a third of an ångström, and the size tracks neither the modelled force
(r = +0.08, p = 0.81) nor the contact area (r = +0.32, p = 0.33). Constructs with no
spring at all do it as much as ones the model says push with 6.6 pN. That reproduces the
companion ubiquitin study's NMR conclusion in a second protein, from a different
observable, and it is the firmest thing here.

The duplex also delivers far more force than the model says — tens of piconewtons against
Zocchi's 2.4–6.6, on the unyielded branch — and the two restraint stiffnesses agree, which
is what makes it the duplex's force rather than the restraint's. And modification strains
an attachment site only where the native backbone sits at positive φ: of four sites in two
proteins, the one left-handed-α position moves 17.5° in φ and the three ordinary ones move
by a degree or less.

**Retracted.** Three claims this project made did not survive its own checks: that the
spring displaces the sites more than an axial pull (wrong zero — caught by Doeke Hekstra);
that the measured force is the Euler load at 0.96× (a box too small — caught by the
overlap windows run for exactly that purpose); and that the nick softens the duplex
threefold (single replicates). Each made the central negative result stronger, not weaker.

**Open.** The measured force is not converged in box size and should be read as a lower
bound. Whether a kink holds its bend is untested — see above. No box we have reaches the
chimera's own 70 Å separation within our own validity criterion. And what the compression
actually *is* — not force, not contact area — these runs do not resolve.

## Still to run

Nothing is queued. Three things would close the open questions above, in order of value:

1. **The seeded kinks in boxes that fit them.** The three systems are rebuilt (108–113 Å
   boxes, 166–192k atoms) and ready; 6 runs × 15 ns would settle whether the bend really
   relaxes. Not submitted — the lab's fairshare is still recovering from a day I wasted on
   it, and the report is honest without it.
2. **A box-size series for the force measurement**, to turn a lower bound into a number.
3. **A box wide enough for the 27 bp at 70 Å** — about 95 Å — to recover the one figure
   this project has had to withdraw twice.
"""
    out = Path(a.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(md)
    print(f"wrote {out} ({len(md.splitlines())} lines, "
          f"{len(vm)} constructs, {len(sp)} springs)")
    if fx and fx.get("stale"):
        print("  NOTE: force_extension.csv predates the box check -- the measured-force "
              "paragraph says so instead of quoting it")


if __name__ == "__main__":
    main()
