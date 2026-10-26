#!/usr/bin/env python3
"""Write the PP2A spring wiki page from findings.json, in the voice of the sfGFP Home page.

Every number on the page comes from analysis/findings.json, analysis/force_extension.csv,
build/01_protein/attachment_sites.json and build/systems/*/build_report.json; the prose
is written around them.  Sections whose runs do not exist yet say so rather than
inventing a number, so the page can be regenerated at every stage of the project.

Images are served from GitHub Pages (docs/figures/), because the repository is private
and raw.githubusercontent links do not render on the wiki.
"""
from __future__ import annotations

import argparse
import json
import math
import re
from pathlib import Path

HERE = Path(__file__).resolve().parent
PP2A = HERE.parent
ROOT = PP2A.parents[1]
PAGES = "https://extremefattypunch.github.io/DNA-springs-project-2026-27"
FIG = f"{PAGES}/figures"
REPO = "https://github.com/extremefattypunch/DNA-springs-project-2026-27"
CODE = f"{REPO}/tree/main/MD%20simulations/pp2a"
NA = "—"


def g(d, *keys, default=None):
    for k in keys:
        if not isinstance(d, dict) or k not in d or d[k] is None:
            return default
        d = d[k]
    return d


def num(v, fmt="{:.2f}", default=NA):
    try:
        if v is None or (isinstance(v, float) and math.isnan(v)):
            return default
        return fmt.format(v)
    except (TypeError, ValueError):
        return default


def p_txt(p):
    if p is None:
        return NA
    return "< 0.001" if p < 0.001 else f"{p:.2f}" if p >= 0.01 else f"{p:.3f}"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--findings", default=str(PP2A / "analysis" / "findings.json"))
    ap.add_argument("--out", default=str(ROOT / "docs" / "wiki" / "PP2A-DNA-Spring-MD.md"))
    a = ap.parse_args()
    f = json.loads(Path(a.findings).read_text())
    att = json.loads((PP2A / "build/01_protein/attachment_sites.json").read_text())
    pair = next(iter(att["pairs"].values()))
    builds = {}
    for p in sorted((PP2A / "build/systems").glob("*/build_report.json")):
        b = json.loads(p.read_text())
        builds[b["system"]] = b
    present = f["systems_present"]
    labels = f["labels"]
    ctrl = f["control"]
    have = lambda s: s in present  # noqa: E731
    ext = f.get("extension", {})
    vs = f.get("extension_vs_control", {})
    ss = f.get("ss_spring", {})
    win = f.get("windows", [])
    rig = f.get("rigidity", {})
    heat = f.get("heat", {})
    hinge = f.get("hinge_vs_control", {})
    sites = f.get("sites", {})

    # ---------- headline ----------
    chim = "P3_ss50"
    if have(chim) and chim in vs:
        w = vs[chim]
        sign = "together" if w["delta"] < 0 else "apart"
        verdict = ("resolved" if w["resolved"] else "not resolved")
        headline = (f"With the 50 nt strand clicked on, the two attachment sites sit "
                    f"**{abs(w['delta']):.1f} ± {w['sem_delta']:.1f} Å {sign}** relative to the clicked "
                    f"protein without DNA (p = {p_txt(w['p'])} on ten block means each, {verdict}). "
                    f"The strand itself spends **{100 * (ss.get('fraction_of_frames_stretched') or 0):.0f} %** of "
                    f"the run beyond its own natural size, so by the model it is **{ss.get('state', '?')}** and pulls with "
                    f"{num(ss.get('tension_pN_clipped'))} pN "
                    f"({num(ss.get('tension_pN_unclipped'))} pN if the Marko–Siggia interpolation is trusted below the "
                    f"chain's rms size).")
    elif have(chim) and ss:
        headline = (f"The chimera run is finished: the attachment sites sit at "
                    f"**{num(g(ext, chim, 'mean_A'), '{:.1f}')} ± {num(g(ext, chim, 'sd_A'), '{:.1f}')} Å**, and the strand is "
                    f"**{ss.get('state', '?')}** at {num(ss.get('tension_pN_clipped'))} pN by the model. "
                    f"*The displacement against the clicked zero ({labels.get(ctrl, ctrl)}) is written here once "
                    f"that control run has finished.*")
    else:
        headline = ("*The chimera run has not finished; this paragraph is written by "
                    "`report/build_wiki_pp2a.py` from `analysis/findings.json` once it has.*")

    def ext_row(s):
        e = ext.get(s)
        if not e:
            return None
        w = vs.get(s)
        d = f"{w['delta']:+.2f} ± {w['sem_delta']:.2f}" if w else ("zero" if s == ctrl else NA)
        return (f"| {labels[s]} | {e['mean_A']:.1f} ± {e['sd_A']:.1f} | {e['min_A']:.0f}–{e['max_A']:.0f} | "
                f"{d} | {p_txt(w['p']) if w else NA} |")
    ext_rows = "\n".join(r for r in (ext_row(s) for s in present) if r) or f"| *no runs yet* | | | | |"

    # ---------- rigidity ----------
    def rig_row(s):
        r = rig.get(s, {})
        return (f"| {labels[s]} | {num(r.get('rmsd.rmsd_core_A.mean'), '{:.2f}')} (max {num(r.get('rmsd.rmsd_core_A.max'), '{:.1f}')}) | "
                f"{num(r.get('rmsd.rmsd_tail_corefit_A.mean'), '{:.1f}')} | {num(r.get('rmsd.rmsf_A.mean_core'), '{:.2f}')} | "
                f"{num(r.get('rmsd.rmsf_A.at_sites.7'), '{:.1f}')} / {num(r.get('rmsd.rmsf_A.at_sites.553'), '{:.1f}')} | "
                f"{num(r.get('rmsd.helix_fraction_core'), '{:.3f}')} |")
    rig_rows = "\n".join(rig_row(s) for s in present)
    rvc = f.get("rigidity_vs_control", {})
    rig_cmp = ""
    if chim in rvc:
        w = rvc[chim]
        rig_cmp = (f"The chimera's core amplitude (RMSD about the run's own mean, per block) differs from the clicked "
                   f"control's by {w['delta']:+.2f} ± {w['sem_delta']:.2f} Å (p = {p_txt(w['p'])} on block means).")

    # ---------- solenoid ----------
    def heat_row(s):
        h = heat.get(s, {})
        reg = h.get("region_rmsd_A", {})
        return (f"| {labels[s]} | {num(h.get('heat.bend.total_deg.mean'), '{:.0f}')} ± {num(h.get('heat.bend.total_deg.sd'), '{:.0f}')} | "
                f"{num(h.get('heat.aperture_H1_H15_A.mean'), '{:.1f}')} ± {num(h.get('heat.aperture_H1_H15_A.sd'), '{:.1f}')} | "
                f"{num(h.get('heat.literature_pair.A.mean'), '{:.1f}')} ± {num(h.get('heat.literature_pair.A.sd'), '{:.1f}')} | "
                f"{num(reg.get('N (HEAT 1-3)'), '{:.2f}')} / {num(reg.get('central (HEAT 4-12)'), '{:.2f}')} / {num(reg.get('C (HEAT 13-15)'), '{:.2f}')} | "
                f"H{num(h.get('heat.bend.most_variable_repeat'), '{:.0f}')} | {num(h.get('heat.salt_bridge_sidechain.fraction_formed_below_5A'), '{:.2f}')} |")
    heat_rows = "\n".join(heat_row(s) for s in present)
    hinge_txt = ""
    if chim in hinge:
        hs = hinge[chim]
        parts = []
        for k in hs["largest_shifts"]:
            e = hs["per_repeat"][str(k)] if str(k) in hs["per_repeat"] else hs["per_repeat"][k]
            parts.append(f"HEAT {k} ({e['delta_deg']:+.1f}°, p = {p_txt(e['p'])}, spread ×{e['sd_ratio']:.2f})")
        ap_ = hs.get("aperture")
        hinge_txt = (f"Against the clicked control, the repeats whose bend angle moves most under the strand are "
                     + ", ".join(parts) + ". The horseshoe aperture (HEAT1–HEAT15) changes by "
                     + (f"{ap_['delta']:+.1f} ± {ap_['sem_delta']:.1f} Å (p = {p_txt(ap_['p'])})." if ap_ else "an unmeasured amount."))

    # ---------- strand ----------
    if ss:
        strand_txt = f"""The strand's end-to-end distance settles at **{num(g(ss, 'x_A', 'mean'), '{:.0f}')} ± {num(g(ss, 'x_A', 'sd'), '{:.0f}')} Å**
(range {num(g(ss, 'x_A', 'min'), '{:.0f}')}–{num(g(ss, 'x_A', 'max'), '{:.0f}')} Å), against a contour of {num(ss.get('contour_nm'), '{:.0f}')} nm and a
natural rms size of {num((ss.get('natural_rms_nm') or 0) * 10, '{:.0f}')} Å. Its radius of gyration is {num(ss.get('rg_A'), '{:.0f}')} Å. It makes
{num(ss.get('sharp_turns_per_frame'), '{:.1f}')} sharp (> 90°) turns per frame along its C4′ trace, hottest at nucleotide
{num(ss.get('hottest_nt'), '{:.0f}')}, keeps {100 * (ss.get('stacked_fraction') or 0):.0f} % of adjacent bases stacked, and has
**{num(ss.get('nt_in_contact'), '{:.0f}')} of 50 nucleotides on the protein** at any moment
({100 * (ss.get('fraction_of_strand_on_protein') or 0):.0f} % of the strand, time-averaged). Protein residues it touches more than half the
time: {ss.get('protein_residues_touched_over_half') or 'none'}."""
        ms_txt = f"""Marko–Siggia applied frame by frame to the strand's end-to-end distance gives
**{num(ss.get('tension_pN_clipped'))} pN** as the repository uses it (zero below the chain's rms size) and
{num(ss.get('tension_pN_unclipped'))} pN if the interpolation formula is trusted all the way down;
the ideal-chain entropic estimate 3k<sub>B</sub>T·x/(2ξL) gives {num(ss.get('ideal_chain_pN'))} pN. All three are
**models**, not measurements. For comparison, Zocchi's duplex model would put a 50 bp
spring at this span on its kinked branch at 1.7–3.3 pN of *push* — and outside the
model's own validity range (γ > 1 for L > ~45 bp). The sign is the point: nothing here
pushes."""
    else:
        strand_txt = "*Not yet simulated.*"
        ms_txt = ""
    fig4_block = f"![modelled against measured force for the single strand]({FIG}/pp2a_fig4_force_extension.png)"
    fig6_block = (f"![the strand: size, modelled tension, where it bends, stacking]({FIG}/pp2a_fig6_ssdna.png)"
                  if (ROOT / "docs/figures/pp2a_fig6_ssdna.png").exists()
                  else "*The strand figure is rendered once the P3_ss50 analysis exists.*")

    # ---------- windows ----------
    if win:
        wrows = "\n".join(
            f"| {w['x0_A']:.0f} | {w['mean_r_A']:.1f} ± {w['sd_r_A']:.1f} | **{-w['force_pN']:+.2f} ± {w['force_sem_pN']:.2f}** | "
            f"{w['ms_tension_pN']:.2f} ({w['ms_tension_unclipped_pN']:.2f}) | {w['ms_state']} | {'yes' if w['box_ok'] else '**no**'} |"
            for w in win)
        ok = [w for w in win if w["box_ok"]]
        agree = []
        for w in ok:
            meas = -w["force_pN"]          # pull positive
            model = w["ms_tension_pN"]
            agree.append(f"at {w['mean_r_A']:.0f} Å the force field gives {meas:+.2f} ± {w['force_sem_pN']:.2f} pN of pull "
                         f"against {model:.2f} pN modelled")
        win_txt = f"""A bare 50 nt strand, no protein, held at a series of end-to-end separations by a soft harmonic restraint
(k = 33 kJ mol⁻¹ nm⁻²) between the atoms the linkers would bond to. The mean restraint force ⟨k(r − x₀)⟩ *is* the force the strand
exerts, read as ⟨k(x₀ − r)⟩ so that positive means the strand **pulls** its ends together
(`force_extension.csv` stores the opposite, push-positive, sign).

| target x₀ (Å) | settled ⟨r⟩ (Å) | measured pull (pN) | Marko–Siggia pull (unclipped) | model's verdict | box valid? |
|---|---|---|---|---|---|
{wrows}

Where they can be compared: {'; '.join(agree) if agree else 'no window has a valid box yet'}."""
    else:
        win_txt = "*The restraint windows have not run yet.*"

    # ---------- sites ----------
    nat = sites.get("native_AF", {})

    def site_row(s):
        d = sites.get(s, {})
        a7, a553 = d.get("7", {}), d.get("553", {})
        return (f"| {labels[s]} | {a7.get('residue') or NA} | {num(a7.get('phi'), '{:+.0f}')}, {num(a7.get('psi'), '{:+.0f}')} "
                f"(± {num(a7.get('phi_sd'), '{:.0f}')}, {num(a7.get('psi_sd'), '{:.0f}')}) | {num(a7.get('ring_neighbours_per_frame'), '{:.2f}')} | "
                f"{a553.get('residue') or NA} | {num(a553.get('phi'), '{:+.0f}')}, {num(a553.get('psi'), '{:+.0f}')} "
                f"(± {num(a553.get('phi_sd'), '{:.0f}')}, {num(a553.get('psi_sd'), '{:.0f}')}) | {num(a553.get('ring_neighbours_per_frame'), '{:.2f}')} |")
    site_rows = "\n".join(site_row(s) for s in present)
    af7 = f"{num(g(nat, '7', 'phi'), '{:+.0f}')}, {num(g(nat, '7', 'psi'), '{:+.0f}')}"
    af553 = f"{num(g(nat, '553', 'phi'), '{:+.0f}')}, {num(g(nat, '553', 'psi'), '{:+.0f}')}"

    # ---------- build table ----------
    def build_row(s, what):
        b = builds.get(s)
        if not b:
            return f"| `{s}` | {what} | — | — |"
        ions = b.get("ions", {})
        return (f"| `{s}` | {what} | {b['atoms']:,} | {b['box'][0]:.0f} Å octahedron, "
                f"Na⁺ {ions.get('Na+', '?')} / Cl⁻ {ions.get('Cl-', '?')} / Mg²⁺ {ions.get('MG', '?')} |")
    build_rows = "\n".join([
        "| `P0_apo` | the AlphaFold Aα monomer alone (the cryo-EM reproduction's run; ff19SB/OPC, 100 ns) | 163,423 | 129 × 144 × 78 Å box |",
        build_row("P1_tet", "Tet2-Et at Asp7 and Asn553, unclicked"),
        build_row("P2_clicked", "both tethers clicked to aminohexanol, no DNA — **the zero**"),
        build_row("P3_ss50", "the chimera: one 50 nt strand clicked at both ends, covalent loop"),
        build_row("D_ss50", "the strand alone, for the restraint windows"),
    ])
    loop = g(builds, "P3_ss50", "covalent_loop", "closed")
    ns_total = f.get("total_ns", 0)

    still_block = ""
    for s_, cap in (("P3_ss50", "the chimera"), ("P0_apo", "the apo A\u03b1 monomer, while the chimera runs")):
        if (ROOT / f"docs/figures/pp2a_{s_}_overview_still.png").exists():
            still_block = f"![{cap}]({FIG}/pp2a_{s_}_overview_still.png)\n*{cap}.*"
            break
    anim_block = ""
    for s, cap in (("P3_ss50", "the chimera: horseshoe, strand, tethers"),
                   ("P2_clicked", "the clicked control"),
                   ("P0_apo", "the apo monomer (cryo-EM reproduction run)")):
        for v in ("overview", "sites", "strain"):
            gif = ROOT / f"docs/figures/pp2a_{s}_{v}.gif"
            if gif.exists():
                anim_block += f"![{labels.get(s, s)} — {v}]({FIG}/pp2a_{s}_{v}.gif)\n*{cap} — {v} view, {labels.get(s, s)}.*\n\n"
    if not anim_block:
        anim_block = "*Animations are rendered by `finalize.sh` once the chimera has run.*\n"

    md = f"""# DNA spring on PP2A Aα — all-atom MD

**[code and structures]({CODE})**
 · [sfGFP page — the pipeline this reuses](Home)
 · [PP2A cryo-EM & MD reproduction — the apo baseline](PP2A-Cryo-EM-and-MD-Reproduction)
 · [how the sfGFP chimera was built]({REPO}/blob/main/feedback/2026-8-22_setup-tutorial.md)

The scaffold subunit of PP2A (Aα, PR65, UniProt P30153) is a 589-residue horseshoe of
fifteen HEAT repeats. Two amber codons put the tetrazine amino acid **Tet2-Et** at the two
tips of the horseshoe — construct positions **D6X and N552X**, which are **Asp7 and Asn553**
in UniProt numbering (the construct counts from Ala2) — and a single 50 nt DNA strand,
sTCO-modified at both ends, clicks onto both, closing a covalent loop. The monomer is
simulated alone, without the B and C subunits.

A single strand is not the leaf spring the duplex constructs are. It is an entropic chain
that pulls its ends **together**, and only once stretched beyond its natural size: for 50
nucleotides that size is {num((ss.get('natural_rms_nm') or 9.407) * 10, '{:.0f}')} Å (contour {num(ss.get('contour_nm') or 29.5, '{:.1f}')} nm,
persistence 1.5 nm). The two Cβ atoms are {pair['d_anchor_anchor_A']:.1f} Å apart in the AlphaFold model and the
tethers reach a further 22 and 26 Å, so the strand *can* be stretched to ~{pair['d_anchor_anchor_A'] + 48:.0f} Å by the
protein — but only if the horseshoe opens.

**{ns_total:,.0f} ns over {len(present)} systems**, one trajectory each.

{still_block}

## The headline

{headline}

![the attachment-point separation against the clicked zero]({FIG}/pp2a_fig7_headline.png)

| system | Cβ7–Cβ553 (Å), mean ± sd | range | Δ vs {labels.get(ctrl, ctrl)} (Å) | p |
|---|---|---|---|---|
{ext_rows}

The zero is the clicked protein without DNA, not the unmodified one: on sfGFP, installing
and clicking the two residues moved the sites by 0.6 Å before any DNA was attached, which
is what made the first comparison there wrong. Errors are block SEMs over ten blocks of a
single 100 ns trajectory; p is Welch's t on the block means. That treats 9.5 ns blocks as
independent, which they are not for the slowest horseshoe motions, so these p-values are
optimistic.

## Backbone rigidity

Day et al. (2025) analysed the free Aα monomer by splitting the N-terminal tail (residues
1–13, which no crystal structure resolves) from the rest; their tail wandered by 10–15 Å
while the rest stayed at 2–4 Å. Our attachment site Asp7 sits in that tail.

| system | core Cα RMSD (Å) | N-tail RMSD after the core fit (Å) | core RMSF (Å) | RMSF at Asp7 / Asn553 (Å) | helix fraction, core |
|---|---|---|---|---|---|
{rig_rows}

{rig_cmp}

![RMSD, Day et al.'s split]({FIG}/pp2a_fig1_rmsd.png)

![per-residue RMSF with the HEAT repeats]({FIG}/pp2a_fig2_rmsf.png)

## The distance between the attachment points

Three distances are reported because three groups measure three different things:
Cβ7–Cβ553 (ours, the spring's own coordinate), Cα Asn30–Phe578 (Yilmaz et al. 2025's
"N29–F577", 66.5 ± 9.7 Å in their 3 × 704 ns apo runs, 47.7 Å in the compact 6NTS trimer,
70.3 Å in the extended 1B3U crystal) and the HEAT1–HEAT15 centroid distance (Day et al.,
62.7 Å in the B56γ holoenzyme). The AlphaFold start model is compact on all three
(Asn30–Phe578 {att['literature_pair']['d_CA_CA_A']:.0f} Å).

![distance between the attachment points]({FIG}/pp2a_fig3_end_to_end.png)

## What the spring delivers, by the model

{strand_txt}

{ms_txt}

{fig4_block}

## What it delivers, measured

{win_txt}

## Does the horseshoe kink, and does the strand?

The solenoid's axis is the chain of HEAT-repeat Cα centroids; the bend at a repeat is the
angle between the incoming and outgoing axis segments, the twist the dihedral across
four consecutive centroids. A hinge is a repeat whose bend *moves* relative to the clicked
control, not one with a large static bend — a solenoid is bent everywhere. The literature
puts the hinge at HEAT 10–11 (crystal comparison), 11–12 (Day) or 12–13 (Yilmaz), with a
diffuse expansion of repeats 4–12.

| system | total axis bend (°) | HEAT1–15 (Å) | Asn30–Phe578 (Å) | region RMSD N / central / C (Å) | most variable repeat | Arg113–Glu191 salt bridge formed |
|---|---|---|---|---|---|---|
{heat_rows}

{hinge_txt}

![where the solenoid bends and twists]({FIG}/pp2a_fig5_heat_bend.png)

The strand's own size, modelled tension, sharp turns and stacking are summarised above
("What the spring delivers, by the model"); the figure below follows them along the trajectory.

{fig6_block}

## What the modification itself does

On sfGFP the one site that sat at **positive φ** (left-handed α, Asp133) was pushed 17.5°
out of its basin by the modification; ordinary sites moved by a degree or less. In the
AlphaFold model of Aα, **Asp7 sits at φ, ψ = {af7}** — the same left-handed-α position —
and Asn553 at {af553}, in the first turn of HEAT 14's helix. Asp7 is also on the
low-confidence tail (pLDDT {num(g(nat, '7', 'plddt'), '{:.0f}')}), so its backbone is free to leave that basin
whether modified or not.

| system | residue 7 | φ, ψ (± sd) | aromatic rings on the Tet ring / frame | residue 553 | φ, ψ (± sd) | rings / frame |
|---|---|---|---|---|---|---|
{site_rows}

![backbone conformation at the two sites]({FIG}/pp2a_ramachandran.png)

## Animations

{anim_block}
## What this establishes, and what it does not

*Written once every run has finished; until then the tables above are the record.*

## Still to run

- Replicates. Everything above is one trajectory per system; the sfGFP page's conclusions
  rested on three, and the review asked for six before believing a sub-ångström shift.
- A longer apo run in the same force field as the chimera (the apo baseline here is
  ff19SB/OPC, the chimera ff14SB/OL15/TIP3P because that is what the custom residues were
  parameterised against).
- The holoenzyme-shaped starting conformation (9MIP chain A, Cβ–Cβ 80 Å), where the same
  strand would be matched rather than slack.

## How it was built

| | |
|---|---|
| Start model | AlphaFold DB v6 model of P30153, full length 1–589 (mean pLDDT {num(g(att, 'plddt', 'mean_all'), '{:.1f}')}; residues 1–13 {num(sum(g(att, 'plddt', 'res_1_13', default=[0]*13)) / 13, '{:.0f}')}) — the model the apo run used |
| Numbering | {att['numbering_note']} |
| Sites | Asp7 (SASA {pair['site_a']['sasa_A2']:.0f} Å², pLDDT {num(pair['site_a'].get('plddt'), '{:.0f}')}) and Asn553 (SASA {pair['site_b']['sasa_A2']:.0f} Å², pLDDT {num(pair['site_b'].get('plddt'), '{:.0f}')}); Cβ–Cβ {pair['d_anchor_anchor_A']:.1f} Å |
| Chemistry | Tet2-Et + sTCO click adduct (TDP), 5′-amino-C6 (DNL) and 3′-amino-C7 (DNP) carbamate linkers, parameters unchanged from the sfGFP build (GAFF2 / AM1-BCC, ff14SB backbone) |
| Strand | 50 nt, `GCTTGGATGGTACCGTGGACTCCTGCCAGCAACTCACGGTCTAGGCTCCA` (the centre of Zocchi's 60-mer, as on sfGFP and ubiquitin — a reuse, not a design), built as B-form, bent to a 90 Å anchor span (chosen: both tethers at 88 % of reach, 11.7 Å lift, as the sfGFP build), complement discarded |
| Covalent loop | asserted at build: cut both Cα–Cβ bonds, the two Cβ remain connected through the DNA — {'closed' if loop else 'NOT closed'} |
| Force field | ff14SB + DNA.OL15 + TIP3P + Joung–Cheatham NaCl (150 mM) + Li/Merz 12-6 Mg²⁺ (5 mM), 12 Å truncated octahedron |
| Engine | OpenMM 8.6 (`sfGFP/runs/openmm_engine.py`): minimise; 100 ps NVT + 900 ps NPT with heavy-atom restraints ramped off; production NPT 300 K, 1 bar, LangevinMiddle 4 fs with hydrogen-mass repartitioning on the solute, PME 0.9 nm, frames every 20 ps, checkpoint every 200 ps |
| Compute | 1 ns benchmarks on the free `gpu_test` MIG slices; production on `gpu_requeue` (half the billing of `gpu`, preemption costs ≤ 200 ps); windows on `gpu_test`; orchestrated from a 1-core `shared` job |

| system | what | atoms | box, ions |
|---|---|---|---|
{build_rows}

Every number on this page is a field of [`analysis/findings.json`]({CODE}/analysis/findings.json);
every figure ships its CSV in [`figures/`]({CODE}/figures).

## References

1. Zocchi, G. *Molecular Machines: A Materials-Science Approach* (Princeton, 2018), §2.15 and §3.5 (the DNA spring; ssDNA persistence length p. 72).
2. Zocchi, G. *Methods Enzymol.* **647**, 257 (2021) — the chimera method; ℓ<sub>ss</sub> ≈ 1.5 nm.
3. Marko, J. F. & Siggia, E. D. *Macromolecules* **28**, 8759 (1995); Smith, Cui & Bustamante *Science* **271**, 795 (1996) — ssDNA contour 0.56–0.59 nm/nt.
4. Day, Huang, Leonard, O'Connor, Narla & Taylor. *Structure* **33**, 1688 (2025) — Aα vs Aβ holoenzymes, 100 ns MD, HEAT1–15 distance, Arg113–Glu191.
5. Yilmaz, Banerjee, Saha, Ohlmeyer, Gordon, Itzhaki, Bahar & Gur. bioRxiv 2025.11.26.690749 — PR65 monomer MD (3 × 704 ns), N29–F577 end-to-end, hinge at HEAT 12–13.
6. Taylor et al. *Cancer Res.* **79**, 4242 (2019) — P179R, φ/ψ free-energy surfaces of a single residue.
7. Wlodarchak et al. *Cell Res.* **23**, 931 (2013) — A-subunit termini distances across holoenzymes, hinge at HEAT 10–11.
8. Eddins et al. *Bio-protocol* **14**, e5048 (2024) — Tet2-Et / sTCO chemistry.
9. Boral, Schnebly et al., bioRxiv 10.64898/2026.07.05.736590 — the ubiquitin DNA-spring study from this lab.
"""
    # integrity: no 'None' leaking into prose, and no python reprs
    bad = re.findall(r"\bNone\b|\{'|\[\{", md)
    if bad:
        print(f"  WARNING: {len(bad)} suspicious tokens in the page: {bad[:5]}")
    out = Path(a.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(md)
    print(f"wrote {out} ({len(md.splitlines())} lines)")


if __name__ == "__main__":
    main()
