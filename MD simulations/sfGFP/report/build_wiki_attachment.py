#!/usr/bin/env python3
"""Generate the wiki page on choosing DNA-spring attachment sites on sfGFP.

Generated rather than written by hand for the same reason Home.md is: every geometric
number in it comes from analysis/attachment_candidates.json, measured on 2B3P, so the
page cannot drift from the structure it describes.  The prose is the argument; the numbers
are the data.

    python report/build_wiki_attachment.py
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

RAW = ("https://raw.githubusercontent.com/extremefattypunch/"
       "DNA-springs-project-2026-27/main/docs/figures")
REPO = "https://github.com/extremefattypunch/DNA-springs-project-2026-27"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--geom", default="analysis/attachment_candidates.json")
    ap.add_argument("--out", default="../../docs/wiki/Choosing-Attachment-Sites.md")
    a = ap.parse_args()
    g = json.loads(Path(a.geom).read_text())
    ab = json.loads(Path("analysis/aband_ratio.json").read_text())
    AB = ab["groups"]
    P, R = g["pairs"], g["residues"]
    ds = g.get("design_search", {})

    def pair_row(k, note):
        e = P[k]
        return (f"| **{e['residues'].replace('-', '–')}** | {e['cb_cb_A']:.1f} Å | "
                f"{e['sasa_A2'][0]:.0f} / {e['sasa_A2'][1]:.0f} Å² | {note} |")

    md = f"""# Choosing attachment sites on sfGFP

Where to click the DNA spring if the aim is to move the fluorescence **wavelength**
rather than its brightness. This page exists to answer a specific proposal — put the
spring where roGFP puts its redox disulfide — and the short answer is that the target is
well chosen and the delivery is impossible.

Everything geometric below is measured on **2B3P** (sfGFP, 1.40 Å, Pédelacq 2006) by
[`analysis/attachment_candidates.py`]({REPO}/blob/main/MD%20simulations/sfGFP/analysis/attachment_candidates.py).

---

## The answer first

**No attachment site delivers a 10–20 nm wavelength shift at the force a bent duplex
provides.** Three independent routes converge on **0.03–0.3 nm, ceiling ~1 nm** — short
of the target by a factor of 30 or more, where choosing the best site buys a factor of
2–3.

| route | measured input | predicted at 2–7 pN |
|---|---|---|
| eGFP's own hydrostatic response | eGFP emission blue-shifts **0.7 nm per 50 MPa**, room temperature, in solution; a 5 pN spring over the barrel cross-section is ≈1 MPa equivalent | **≈0.014 nm** |
| barrel compliance × spectral sensitivity | spring constants 1–17 N/m; **0.5 Å of chromophore displacement plus ~10° tilt buys ~1 nm** (glycine-insertion structures, room temperature) | **0.01 – 0.44 nm** |
| electrochromic (Stark) | in-situ Δμ gives **1.56 nm per MV/cm** for EGFP; a perfectly coupled 5 pN spring moves the mean field ≈0.3 MV/cm | **≲0.5 nm** as a ceiling |

> **A correction to an earlier version of this page.** It quoted Citrine's 17 nm
> pressure shift and derived a sensitivity of 21 nm per Å. That 17 nm is a
> **cryo-trapped** number — the deformation frozen in. The same protein at room
> temperature moves **≈1 nm per 100 MPa**, so the sensitivity was ~3.4× too generous and
> the ceiling too high. The corrected figures are above.

Two direct experiments exist on a mechanically loaded GFP, and **neither reports a colour
change**. Pulled at 30–60 pN — five to twenty times our spring — the chromophore goes
**dark 3.5 ms before any length change**, with no fluorescence from any intermediate:
*"complete structural integrity is necessary to observe single-molecule fluorescence of
GFP"* (Ganim & Rief 2017). In cells, reversible **on/off** switching is seen at
single-digit pN (Shoyer 2024), whose authors also state plainly that *"forced-induced
changes in the excitation or emission wavelengths of FPs have not been described."*

**GFP's chromophore behaves as a threshold device in the mechanical coordinate: intact
and green, or broken and dark. No state has been observed that is bright and a different
colour.**

### And the energy budget is not close

Our own clamp measures the barrel at 230–350 k<sub>B</sub>T nm⁻² along D133–N149, i.e.
946–1440 pN nm⁻¹. A 5 pN load stores ½Fx ≈ **0.003 k<sub>B</sub>T**. Granting an
implausible 3 Å of coupled compliance still gives only 0.5 k<sub>B</sub>T. The barriers
gating the chromophore's soft coordinates are **1.2 and 7.7 k<sub>B</sub>T** — 2.4× and
15× the spring's entire budget. A leaf spring cannot buy a hydrogen-bond rearrangement.

The noise floor is the protein, not the instrument: thermal field fluctuations at the
chromophore are **δF = 5 ± 0.3 MV/cm**, which is 7–10 nm of instantaneous per-molecule
spectral wander and accounts for GFP's inhomogeneous linewidth. A perfectly coupled 5 pN
spring shifts the *mean* by ~6 % of that.

**The honest restatement of the goal**: shift a two-state population and read it as a
**band-amplitude ratio at fixed emission wavelength** — not as a wavelength shift. There
is already a signal of exactly that shape in your own data, [below](#there-is-already-a-force-ordered-signal-in-your-own-data).

---

## Why the roGFP positions cannot take a spring

The proposal was to anchor at **S147↔Q204** or **N149↔S202**, the pairs roGFP uses,
because constraining those two β-strands is known to swing the 400/490 nm excitation
equilibrium. The reasoning is right: residues 147–149 sit on β-strand 7 and 202–205 on
β-strand 10, the two strands are **adjacent and antiparallel**, and they carry the
chromophore's own hydrogen bonds — in 2B3P's own `SHEET` record, His148 pairs directly
with Thr203.

But roGFP puts *cysteines* there, and a disulfide is a 2 Å bond. Measured Cβ–Cβ:

| pair | separation | solvent exposure | |
|---|---|---|---|
{pair_row('S147-Q204', 'roGFP1/2 disulfide — **too close for a spring**')}
{pair_row('N149-S202', 'roGFP3/4 disulfide — **too close for a spring**')}
{pair_row('N149-Q204', 'the mixed diagonal')}
{pair_row('D133-N149', '**the pair we build**')}

A DNA spring needs its anchors ~30 Å apart. Each tether — the Tet2-Et side chain, the
sTCO adduct, the DNA's own terminal arm — is 10–15 Å on its own, and the built geometry
holds the duplex ends 62–76 Å apart across a 31 Å Cβ–Cβ separation. Asked to span 4.6 Å,
the duplex would simply lie alongside the protein; it could not load those two residues
against each other at all.

![where the sites are and why the geometry fails]({RAW}/fig8_attachment_sites.png)

### Nor can any buildable pair on those strands pry them apart

Searching both chromophore-supporting strands for exposed pairs in the buildable range
(24–42 Å) returns **{len(ds.get('buildable', []))} candidates**, and every one of them
pushes *across* the His148→Thr203 axis rather than along it. The best alignment is
{max((b['alignment_with_H148_T203_axis'] for b in ds.get('buildable', [])), default=0):.2f}
(1.00 would be a spring pushing exactly along the axis that separates the two H-bond
donors, which sit {ds.get('H148_T203_cb_cb_A', 0):.1f} Å apart).

The reason is unavoidable: two points far enough apart to span with a spring must be
separated *along* the strands, and the strands run the height of the barrel — nearly
perpendicular to the 5 Å direction that separates them. Such a spring shears or twists
the strands; it does not pry them open.

---

## Corrections to the proposal's supporting literature

Five sources were cited. **One is the wrong paper, two are reviews that do not contain
the numbers attributed to them**, and two check out. Six further factual points needed
fixing — set out here because the underlying claims are mostly *true*, just mis-sourced,
and someone ordering constructs needs the right numbers.

| the claim | status |
|---|---|
| roGFP1 = S147C/Q204C, roGFP2 = N149C/S202C | **wrong.** roGFP1 *and* roGFP2 both carry S147C/Q204C (roGFP2 = roGFP1 + S65T); roGFP3 *and* roGFP4 both carry N149C/S202C. roGFP5/6 carry **both** disulfides. All carry C48S and Q80R |
| 6.1-fold for 147–204, 4.3-fold for 149–202 | **numbers right, source wrong and comparison unsafe.** Both are in Hanson 2004 Table I (roGFP1 δ=6.1, roGFP3 δ=4.3) but measured at 475/400 nm, while the S65T pair (roGFP2 δ=5.8, roGFP4 δ=2.6) was measured at 490/400 nm. Not apples-to-apples. roGFP1 also has a *second* row, δ=4.5, from a lipoic-acid rather than DTT titration |
| the 149/202 geometry is roGFP's | **priority is split.** The 149–202 disulfide was first built in YFP as rxYFP by Østergaard 2001, whose 2.2-fold figure is an *intensity* change, not an excitation ratio |
| residues 147–149 lie on β-strand 7 | **147 does not.** β7 begins at His148 in every structure checked; Ser147 is the preceding residue and sits across from Gln204 as a non-hydrogen-bonded wide pair |
| Gln204 couples to the chromophore through a conserved water | **only via the backbone, and the geometry is marginal.** Across ten deposited structures the conserved water W1 sits 2.5–2.8 Å from the chromophore phenol and from Ser205 OG, but **3.5–3.6 Å from Gln204 N** — too long for a good hydrogen bond. At 0.90 Å resolution the Gln204 amide hydrogen is spent on the cross-strand β-sheet bond to Phe223. The **side chain** is not water-coupled to the chromophore, and no source says it is |
| a big excitation shift does not imply a big emission shift | **true but conditional.** It holds only while an excited-state proton acceptor survives. The deGFP series is the counterexample *and* it is an H148 series — see below |
| T203Y stacks with the chromophore and red-shifts emission | **true; the cited paper is about a copepod protein.** The real sources are Ormö 1996 (GFP S65T 489/511 nm → T203Y variant 513/527 nm) and Wachter 1998 ("nearly coplanar π–π stacking", chromophore "displaced by nearly 1 Å") |

---

## What our own simulations say

The proposal argued that our existing spring already transmits stress into the
chromophore neighbourhood, citing *"Ser205 → Glu222: 29% WT → 0% with the 40-bp spring"*.
That number does not support the claim.

Ser205–Glu222 is **bistable** — `findings.json` already flags it — and per replicate it is
bimodal: each 50 ns run either forms the bond and holds it (~0.95 occupancy) or never
forms it (0.00). Its "mean occupancy" is therefore the *fraction of replicates* that
landed in the formed state, a Bernoulli proportion:

| construct | replicates forming it | Fisher exact vs control |
|---|---|---|
| clicked control, no DNA | 3 of 9 | — |
| 40 bp, two-arm *(superseded build)* | **0 of 6** | p = 0.23 |
| 40 bp, closed loop *(the correct molecule)* | **2 of 3** | p = 0.52 |

The quoted "29% → 0%" is six runs landing tails, taken from the two-arm topology that has
since been [retracted](Home) — and in the corrected closed-loop build it moves the other
way. Across all nine constructs the only resolved changes are Gln94 (−7%) and an Arg96
"difference" of 1.000 versus 0.999 that is an artefact of near-zero variance.

**Our own data provides no evidence that a 2–7 pN spring reaches the chromophore's
hydrogen-bond network** — which is what the estimate at the top of this page predicts.

---

## There is already a force-ordered signal in your own data

The 260828 ds/ss/nicked run contains the predicted signature. Taking the 400 nm shoulder
relative to the excitation peak, lamp-ratio corrected:

| construct | modelled force | I(400)/I(peak) | |
|---|---|---|---|
| **ss** | 0 pN | **{AB['ss']['mean']:.4f}** ± {AB['ss']['sd']:.4f} | slack strand, no load |
| **nicked** | ~2.5 pN | **{AB['nicked']['mean']:.4f}** ± {AB['nicked']['sd']:.4f} | |
| **ds** | ~3.1 pN | **{AB['ds']['mean']:.4f}** ± {AB['ds']['sd']:.4f} | |

**Monotone in modelled spring force, {abs(ab['relative_change_pct']):.1f} % relative
decrease**, with the emission peak unchanged. DNA-only blanks are
{100*min(e['blank_fraction_at_400nm'] for e in AB.values()):.2f}–{100*max(e['blank_fraction_at_400nm'] for e in AB.values()):.2f} %
of the signal at 400 nm, so it is not scatter. "Band ratio moves, emission λ fixed" is
exactly the signature the physics predicts.

**Do not believe it yet.** The three a/b/c scans per condition are *technical* replicates
of one preparation, so the quoted scatter is instrument repeatability and says nothing
about the constructs. Absolute brightness differs across the samples, and pH,
concentration and inner-filter effects are not controlled at the 10⁻³ level this needs.
Tet2-Et also **quenches sfGFP at N150, with clicking relieving the quench** — a
first-order confound for any amplitude measurement on these constructs.

### The one cheap measurement that would settle it

**Extend the excitation scan.** It currently runs {ab['scan_range_nm'][0]:.0f}–{ab['scan_range_nm'][1]:.0f} nm
at {ab['scan_step_nm']:.1f} nm, which **starts exactly at the A-band maximum** — the band
is cut in half, its amplitude cannot be fitted, and the 400/490 ratio that this whole
question turns on has therefore never actually been measured.

Re-run **340–510 nm at 0.2 nm, lamp-corrected**, on the ss / nicked / ds triplet with
matched pH and matched absorbance, from independent preparations. No new construct, no
new chemistry. Fit the bands rather than taking `argmax` — on the existing 1 nm grid a
parabola fit already bounds any emission shift at **|Δλ| ≲ 0.2 nm**, which is consistent
with the prediction above and excludes the ambitious version of the goal.

---

## What to do instead

Ranked by how well supported each is, not by ambition. The first is better
supported than any change of site.

**1. Change the force, not the site.** Every number above is linear in force. Our own
restrained-duplex measurement puts the bare 27 bp duplex at **tens of piconewtons on the
unyielded branch**, about six times Zocchi's model value
([see the report]({REPO}/blob/main/docs/index.html)). Going from 5 to 35 pN moves the
prediction from ~0.05 to ~0.35 nm — across, rather than below, the detection threshold.
A shorter, stiffer duplex costs one oligonucleotide. Measuring the force the *chimera*
actually delivers is cheaper than building four constructs, and it gates everything else.

**2. Read a band ratio, not a peak position.** Follow the signal above: fit the 400 and
490 nm bands and report their ratio at fixed emission wavelength. This is the observable
the physics supports, the one the existing data already hints at, and the one a
two-state population shift would move.

**3. If you do re-site, Tyr151 ↔ Ser208.** At
{next((b['cb_cb_A'] for b in ds.get('buildable', []) if b['pair'].startswith('TYR151')), 0):.1f} Å
it is the only buildable pair on the two chromophore-bearing strands that nearly fits the
existing 27 bp geometry, and it loads β7 against β10 in shear. Its justification is
better than a geometric one: a purely *topological* perturbation to β7 — cutting the
chain at 146 rather than 140, with no chemistry near the chromophore — moves the site
coupling by ~4 pKa units and rewrites the proton-transfer dynamics. **β7 mechanics is a
demonstrated lever into the phenolate equilibrium; bulk barrel compression is not.**

**4. Keep D133 ↔ N149 and build the isogenic zero-force twin.** The one experiment nobody
has done is a *paired* comparison at identical chemistry — same clicks, same adducts,
spring versus slack — which is exactly what separates a mechanical effect from the
Tet2-Et quenching. That value does not depend on picking the best site.

### Two things the earlier version of this page got wrong

**T203Y is not the amplifier it looks like.** Its red shift comes mainly from the **loss
of a hydrogen bond**, not from π-stacking geometry: the aliphatic substitutions T203V and
T203I produce a similar shift (Drobizhev 2015). Mechanically modulating a stacking
distance therefore would not recover it, and an earlier version of this page recommended
a **T203Y-sensitised background on that mistaken basis**.

**"Measure intensity instead" is confounded here, not free.** Tet2-Et quenches sfGFP at
N150 and clicking relieves the quench, so amplitude on these constructs already carries a
chemistry signal of the same sign and larger size. That is the same trap as the retracted
point-3 claim, which used the unclicked protein as its zero. Amplitude is only usable
against the clicked, spring-free twin in recommendation 4.

**One point mutation may be worth more than any site change: T65S.** In sfGFP's Thr65
background the proton wire is severed and the protein is >98 % anionic with a dark A
state — which is *why* the A band here is a shoulder rather than a peak. T65S restores
internal proton transfer, a bright A state, and an A:B ratio near unity, where the
ratiometric gain is roughly an order of magnitude higher. Nobody has built T65S sfGFP
with Tet2-Et at 133/149, and it will cost brightness and folding robustness.

---

## References

Verified against the primary source. Where a claim above rests on a secondary source or
on arithmetic of ours, it says so at that point.

**The spring, and what force it delivers**
1. Dietz H, Rief M (2006) Protein structure by mechanical triangulation. *PNAS* — GFP's directional spring constants, 1–17 N/m; fracture 356–548 pN. [10.1073/pnas.0602995103](https://doi.org/10.1073/pnas.0602995103)
2. Dietz H, Rief M (2004) Exploring the energy landscape of GFP by single-molecule mechanical experiments. *PNAS* 101:16192 — first yield at ~35 pN, N-terminal helix detachment. [10.1073/pnas.0404549101](https://doi.org/10.1073/pnas.0404549101)

**Mechanically loaded fluorescent proteins — the direct precedents**
3. Ganim Z, Rief M (2017) Mechanically switching single-molecule fluorescence of GFP by unfolding and refolding. *PNAS* 114:11052 — optical tweezers, 30–60 pN; fluorescence lost 3.5 ms *before* any length change; no intermediate fluoresces. [10.1073/pnas.1704937114](https://doi.org/10.1073/pnas.1704937114)
4. Shoyer TC *et al.* (2024) *Cell Rep Methods* 4:100815 — reversible on/off FP switching under load in cells, modelled 0–10 pN; states outright that force-induced wavelength changes have not been described. [10.1016/j.crmeth.2024.100815](https://doi.org/10.1016/j.crmeth.2024.100815)

**Pressure as the calibrated analogue — and why the room-temperature number is the one to use**
5. Morikawa TJ, Nishiyama M, Yoshizawa K, Fujita H, Watanabe TM (2021) Glycine insertion modulates the fluorescence properties of *Aequorea victoria* GFP and its variants. *Biophys Physicobiol* 18:145 — **eGFP emission blue-shifts 0.7 nm from 0.1 to 50 MPa; eYFP red-shifts 0.5 nm**, same apparatus, room temperature. [10.2142/biophysico.bppb-v18.016](https://doi.org/10.2142/biophysico.bppb-v18.016)
6. Watanabe TM *et al.* (2013) Glycine insertion makes YFP sensitive to hydrostatic pressure. *PLoS ONE* 8:e73212 — one glycine in β-strand 7 displaces the chromophore **~0.5 Å with ~10° tilt for a 1 nm blue shift**; PDB 3VGQ/3VGR. [10.1371/journal.pone.0073212](https://doi.org/10.1371/journal.pone.0073212)

**The electrochromic route**
7. Drobizhev M, Callis PR, Nifosì R, Wicks G, Stoltzfus C, Barnett L, Hughes TE, Sullivan P, Rebane A (2015) Long- and short-range electrostatic fields in the GFP family. *Sci Rep* 5:13223 — 26 mutants spanning ~65 nm on one chromophore; Δμ = 4.47 ± 0.16 D; **δF = 5 ± 0.3 MV/cm** thermal field noise; T203V and T203I give a similar red shift to T203Y, so **T203Y acts by hydrogen-bond loss, not π-stacking distance**. [10.1038/srep13223](https://doi.org/10.1038/srep13223)
8. Slocum JD, First JT, Webb LJ (2017) Orthogonal electric field measurements near the GFP chromophore. *J Phys Chem B* 121:6799 — Stark shifts and chromophore pKa shifts measured for the same mutations, in agreement. [10.1021/acs.jpcb.7b03935](https://doi.org/10.1021/acs.jpcb.7b03935)

**Mechanical and pressure effects on fluorescent-protein spectra**
3. Barstow B, Ando N, Kim CU, Gruner SM (2008) Alteration of citrine structure by hydrostatic pressure explains the accompanying spectral shift. *PNAS* 105:13362 — emission 527 → 510 nm; chromophore deformation up to 0.8 Å. [10.1073/pnas.0802252105](https://doi.org/10.1073/pnas.0802252105)
4. Barstow B, Ando N, Kim CU, Gruner SM (2009) Coupling of pressure-induced structural shifts to spectral changes in a yellow fluorescent protein. *Biophys J* 97:1719 — differential motion of two β-barrel atom clusters. [10.1016/j.bpj.2009.06.039](https://doi.org/10.1016/j.bpj.2009.06.039)
5. *Cell Reports Methods* (2024) — "forced-induced changes in the excitation or emission wavelengths of FPs have not been described". [10.1016/j.crmeth.2024.100815](https://doi.org/10.1016/j.crmeth.2024.100815)

**roGFP and rxYFP: the disulfide precedent**
6. Hanson GT, Aggeler R, Oglesbee D, Cannon M, Capaldi RA, Tsien RY, Remington SJ (2004) Investigating mitochondrial redox potential with redox-sensitive GFP indicators. *J Biol Chem* 279:13044 — roGFP1–6, Table I dynamic ranges. [10.1074/jbc.M312846200](https://doi.org/10.1074/jbc.M312846200)
7. Dooley CT, Dore TM, Hanson GT, Jackson WC, Remington SJ, Tsien RY (2004) Imaging dynamic redox changes in mammalian cells with GFP indicators. *J Biol Chem* 279:22284 — "up to 6-fold in vitro and up to ~3-fold in cytosol". [10.1074/jbc.M312847200](https://doi.org/10.1074/jbc.M312847200)
8. Østergaard H, Henriksen A, Hansen FG, Winther JR (2001) Shedding light on disulfide bond formation: engineering a redox switch in green fluorescent protein. *EMBO J* 20:5853 — rxYFP[149-202], 2.2-fold intensity change. [10.1093/emboj/20.21.5853](https://doi.org/10.1093/emboj/20.21.5853)

**The chromophore, its hydrogen bonds and its colour**
9. Ormö M, Cubitt AB, Kallio K, Gross LA, Tsien RY, Remington SJ (1996) Crystal structure of the *Aequorea victoria* green fluorescent protein. *Science* 273:1392 — 11-stranded barrel; T203Y red-shift. [10.1126/science.273.5280.1392](https://doi.org/10.1126/science.273.5280.1392)
10. Wachter RM, Elsliger M-A, Kallio K, Hanson GT, Remington SJ (1998) Structural basis of spectral shifts in the yellow-emission variants of GFP. *Structure* 6:1267 — π–π stacking; chromophore displaced ~1 Å. [10.1016/S0969-2126(98)00127-0](https://doi.org/10.1016/S0969-2126%2898%2900127-0)
11. Brejc K, Sixma TK, Kitts PA, Kain SR, Tsien RY, Ormö M, Remington SJ (1997) Structural basis for dual excitation and photoisomerization of GFP. *PNAS* 94:2306 — barrel 42 Å × 24 Å; Glu222/Ser205 water. [10.1073/pnas.94.6.2306](https://doi.org/10.1073/pnas.94.6.2306)
12. Hanson GT *et al.* (2002) GFP variants as ratiometric dual emission pH sensors, 1. *Biochemistry* 41:15477 — the deGFP series, green↔blue emission switching. [10.1021/bi026609p](https://doi.org/10.1021/bi026609p)
13. McAnaney TB, Park ES, Hanson GT, Remington SJ, Boxer SG (2002) …2. Excited-state dynamics. *Biochemistry* 41:15489. [10.1021/bi026610o](https://doi.org/10.1021/bi026610o)
14. Ahmed RD *et al.* (2025) MD-guided identification of a brighter sfGFP variant. *Commun Chem* 8:174 — sfGFP H148A/H148C: excitation to ~400 nm, emission stays ~510 nm. [10.1038/s42004-025-01573-4](https://doi.org/10.1038/s42004-025-01573-4)
15. Hartley AM, Worthy HL, Reddington SC, Rizkallah PJ, Jones DD (2016) Molecular basis for functional switching of GFP. *Chem Sci* 7:6484 — sfGFP148azF shifts excitation ~90 nm, emission retained at ~511 nm. [10.1039/c6sc00944a](https://doi.org/10.1039/c6sc00944a)
16. Arpino JAJ, Rizkallah PJ, Jones DD (2012) Crystal structure of EGFP to 1.35 Å. *PLoS ONE* 7:e47132 — the conserved waters W1–W6. [10.1371/journal.pone.0047132](https://doi.org/10.1371/journal.pone.0047132)

**sfGFP, its topology, and where it tolerates modification**
17. Pédelacq J-D, Cabantous S, Tran T, Terwilliger TC, Waldo GS (2006) Engineering and characterization of a superfolder GFP. *Nat Biotechnol* 24:79 — sfGFP; **PDB 2B3P** at 1.40 Å. [10.1038/nbt1172](https://doi.org/10.1038/nbt1172)
18. Huang YM *et al.* (2015) Toward computationally designed self-reporting biosensors. *Biochemistry* 54:6263 — strand order 3-2-1-6-5-4-9-8-7-10-11, all pairings antiparallel except 1–6. [10.1021/acs.biochem.5b00786](https://doi.org/10.1021/acs.biochem.5b00786)
19. Baird GS, Zacharias DA, Tsien RY (1999) Circular permutation and receptor insertion within GFPs. *PNAS* 96:11241 — 10 tolerant permutation sites; camgaroo. [10.1073/pnas.96.20.11241](https://doi.org/10.1073/pnas.96.20.11241)
20. Pope JR *et al.* (2021) Association of fluorescent protein pairs… *Adv Sci* 8:2003167 — Q204 surface-exposed in sfGFP; PDB 5NI3. [10.1002/advs.202003167](https://doi.org/10.1002/advs.202003167)
21. Reddington SC, Tippmann EM, Jones DD (2012) Residue choice defines efficiency of bioorthogonal modification. *Chem Commun* 48:8419 — a dye clicked at position 204. [10.1039/c2cc31887c](https://doi.org/10.1039/c2cc31887c)
22. Xiong Q *et al.* (2022) Expanding protein functionality with genetically encoded DBTD. *Chem Sci* 13:3571 — N149 and Q204 as solvent-exposed sfGFP encoding sites. [10.1039/d1sc05710c](https://doi.org/10.1039/d1sc05710c)

**Cited in the original proposal, and not usable for the claim made**
- Suto K *et al.* (2009) *Genes Cells* 14:727 — **wrong protein.** CpYGFP from the copepod *Chiridius poppei*, cited for T203Y in *Aequorea* GFP. [10.1111/j.1365-2443.2009.01305.x](https://doi.org/10.1111/j.1365-2443.2009.01305.x)
- Wang X *et al.* (2013) *J Mol Med* 91:917 — a **review** of ROS imaging; contains neither the 6.1-fold nor the 4.3-fold figure, and does not mention N149C/S202C, His148 or Thr203. [10.1007/s00109-013-1067-4](https://doi.org/10.1007/s00109-013-1067-4)
- Marchetti M *et al.* (2023) *Sensors* 23:8517 — a **review**; states the roGFP1 genotype correctly but is not where it was measured. [10.3390/s23208517](https://doi.org/10.3390/s23208517)
"""
    out = Path(a.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(md)
    print(f"wrote {out} ({len(md.splitlines())} lines, "
          f"{md.count('https://doi.org/')} resolvable DOIs)")


if __name__ == "__main__":
    main()
