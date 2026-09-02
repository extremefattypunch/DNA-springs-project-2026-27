# Choosing attachment sites on sfGFP

Where to click the DNA spring if the aim is to move the fluorescence **wavelength**
rather than its brightness. This page exists to answer a specific proposal — put the
spring where roGFP puts its redox disulfide — and the short answer is that the target is
well chosen and the delivery is impossible.

It also answers the follow-up question, *"what if we used a different click linker instead
of tetrazine?"* — [that section is here](#would-a-different-click-chemistry-help), and the
answer is no, for reasons that turn out to be neither geometric nor spectroscopic.

Everything geometric below is measured on **2B3P** (sfGFP, 1.40 Å, Pédelacq 2006) by
[`analysis/attachment_candidates.py`](https://github.com/extremefattypunch/DNA-springs-project-2026-27/blob/main/MD%20simulations/sfGFP/analysis/attachment_candidates.py).

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

## Why the roGFP positions are attachable, and still not worth using

> **A retraction.** An earlier version of this page said a DNA spring **could not be
> attached** at these positions at all, on the grounds that a spring "needs its anchors
> ~30 Å apart" while roGFP's pairs are 4.6–5.9 Å apart. **That argument was wrong.** A
> compression spring does not need to be *stretched* to fit — if the anchors are closer
> together the duplex simply bends further, and reach is never the binding constraint.
> Anchors 4.6 Å apart hold a 27 bp duplex at a 4.2 pN working point against
> 5.5 pN at the built pair. The real objection is narrower and is set out
> below. The figure's original caption was wrong for the same reason.

The proposal was to anchor at **S147↔Q204** or **N149↔S202**, the pairs roGFP uses,
because constraining those two β-strands is known to swing the 400/490 nm excitation
equilibrium. The reasoning is right: residues 147–149 sit on β-strand 7 and 202–205 on
β-strand 10 (Ser147 itself falls just outside β7, which begins at His148), the two
strands are **adjacent and antiparallel**, and they carry the chromophore's own hydrogen
bonds — in 2B3P's own `SHEET` record, His148 pairs directly with Thr203.

roGFP puts *cysteines* there, and a disulfide is a 2 Å bond. Measured Cβ–Cβ:

| pair | separation | solvent exposure | |
|---|---|---|---|
| **SER147–GLN204** | 4.6 Å | 66 / 75 Å² | roGFP1/2 disulfide — **too close for a spring** |
| **ASN149–SER202** | 5.8 Å | 61 / 32 Å² | roGFP3/4 disulfide — **too close for a spring** |
| **ASN149–GLN204** | 10.6 Å | 61 / 75 Å² | the mixed diagonal |
| **ASP133–ASN149** | 31.3 Å | 108 / 61 Å² | **the pair we build** |

### What the tether actually does

Measured on the built chimera, each tether spans **21.7 Å** from the protein Cβ to the
DNA backbone phosphorus: 12.4–13.3 Å of clicked tetrazine–TCO adduct plus 10.8 Å of the
DNA's own amino-modifier arm.

The obvious next move is to ask whether the arms are *long enough to reach* a given pair,
and that is the wrong question. **A compression spring does not have to be stretched to
fit.** If the anchors sit closer together, the duplex simply bends further. Reach is
never the binding constraint; the span that results is, because the span sets the force.

The clicked adduct is a strut — dihydropyridazine fused to the sTCO bicycle, no sp³ chain
bonds anywhere in its 13 Å. The DNA-side amino-modifier arm is an alkyl chain, and is
not. They carry the load in series, so the equilibrium span is a **fixed point**: the
duplex pushes, the compliant arm extends until its tension matches, and the span that
results sets the force. Marko–Siggia closes the loop
([`analysis/tether_tradeoff.py`](https://github.com/extremefattypunch/DNA-springs-project-2026-27/blob/main/MD%20simulations/sfGFP/analysis/tether_tradeoff.py)):

| anchor pair | Cβ–Cβ | equilibrium span | arm extension | 27 bp force |
|---|---|---|---|---|
| S147–Q204 | 4.6 Å | 35.6 Å | 23 % | 4.2 pN |
| N149–S202 | 5.9 Å | 37.0 Å | 24 % | 4.2 pN |
| Tyr151–Ser208 | 28.0 Å | 60.0 Å | 28 % | 5.2 pN |
| **D133–N149 (built)** | 31.3 Å | 63.6 Å | 29 % | 5.5 pN |

Every span is on the circular-arc branch below *x*<sub>c</sub> = 85.7 Å, where
*f* = τ<sub>c</sub>/(2R√(1−(x/2R)²)) **diverges** as *x* → 2R. So closer anchors give a
shorter span and, monotonically, **less force**. That is the real answer to "should we
re-site to the chromophore strands": **you can, and it buys you nothing** — 4.2 pN
against 5.5 pN.

> **A second correction to this section.** An earlier version said the close pairs could
> be reached "only with both arms at 100 % of contour — a knife-edge". That was still too
> generous to the idea. A few-pN spring cannot pull a flexible arm anywhere near taut:
> holding a 4 Å-persistence-length chain at 90 % of contour takes **263 pN**, and
> at 99 %, **25,633 pN**. At 2 pN such an arm is 12 % extended and at 7 pN
> 35 %. There is no knife-edge because there is no reaching: added flexible
> length is added contour the spring never sees.

**Where this model and the measured geometry disagree, and by how much.** The MD tether
is 21.7 Å, which puts the arm at 81 % of its contour — entropic elasticity alone
would need 73 pN to hold it there, which the spring does not have. Either the
50 ns snapshot is not at mechanical equilibrium, or the arm is extended by excluded volume
and solvation rather than by tension: an alkyl chain pinned between a charged backbone and
a protein surface is not a free chain. **Both readings are kept.** The force at the built
pair is **5.5–7.2 pN**, and everything downstream on this page is
linear in it.

![the chromophore network, and the candidate anchors to scale](https://raw.githubusercontent.com/extremefattypunch/DNA-springs-project-2026-27/main/docs/figures/fig8_attachment_sites.png)

### Nor can any buildable pair on those strands pry them apart

Searching both chromophore-supporting strands for exposed pairs in the buildable range
(24–42 Å) returns **6 candidates**, and every one of them
pushes *across* the His148→Thr203 axis rather than along it. The best alignment is
0.39
(1.00 would be a spring pushing exactly along the axis that separates the two H-bond
donors, which sit 5.1 Å apart).

The reason is unavoidable: two points far enough apart to span with a spring must be
separated *along* the strands, and the strands run the height of the barrel — nearly
perpendicular to the 5 Å direction that separates them. Such a spring shears or twists
the strands; it does not pry them open.

---

## Would a different click chemistry help?

**No — and both premises behind the question fail.** The mechanics say every alternative
is weaker, and the spectroscopic case for switching, which is the one that looked strong,
does not survive contact with the photophysics.

Seven linkages, priced on the same 27 bp duplex at the built anchor pair. Force is given
in both readings of the tether (rigid strut / compliant fixed point) because the gap
between them is the whole story:

| linkage | tether | *f* rigid | *f* compliant | *k*₂ (M⁻¹s⁻¹) |
|---|---|---|---|---|
| **tetrazine–sTCO (current, measured)** | 21.7 Å | 7.16 pN | 5.16 pN | **7.2 × 10⁴** |
| inverted: tetrazine–DNA + TCO\*-Lys | 18.5 Å | 6.04 pN | 4.76 pN | 5.2 × 10³ – 1.7 × 10⁴ |
| SPAAC, DBCO–DNA + AzK | 17.5 Å | 5.80 pN | 4.66 pN | 0.8 – 0.9 |
| SPAAC, BCN–protein + azide–DNA | 17.0 Å | 5.69 pN | 4.62 pN | 0.75 (aryl azide) |
| CuAAC, propargyl-Lys + azide–DNA | 15.0 Å | 5.32 pN | 4.46 pN | 10 – 100 (per 10–100 µM Cu) |
| maleimide–cysteine, short arm | 14.0 Å | 5.17 pN | 4.39 pN | 10¹ – 10³ |
| oxime, pAcF + aminooxy–DNA | 13.3 Å | 5.08 pN | 4.35 pN | 10⁻³ – 10⁻¹ (aniline-catalysed) |

**The spread across all seven chemistries is 2.1 pN.** The gap between the two
ways of reading the *same* tether is **2.0 pN**. Chemistry choice therefore sits
inside the error bar of the tether-stiffness model — and on our own measured clamp
stiffness (140 pN Å⁻¹, `figures/fig1_force_response.csv`), 2.1 pN moves the
anchors **0.015 Å**, which on the sensitivity established above is ~0.09 nm of
wavelength. Mechanics cannot pick a chemistry, in either direction.

### What *can* decide it: completeness

The rate column spans five orders of magnitude, and that is the axis that matters,
because **unreacted protein is spectroscopically indistinguishable from unloaded
protein**. Tetrazine–sTCO is measured at **72,500 ± 1660 M⁻¹s⁻¹** on GFP-bound tetrazine
in cellulo (Blizzard 2015) and delivers *">95 % encoding and labeling fidelity"* with
*"quantitative"* reaction *"within 5 min at room temperature and physiological pH"*
(Eddins 2024). At 10 µM protein against 100 µM duplex that is *t*<sub>95</sub> ≈ 0.4 s.

The same arithmetic on the alternatives, from primary kinetics rather than supplier
literature — Dommerholt's Table 2 gives DBCO/aliphatic azide 0.8–0.9 M⁻¹s⁻¹ and
BCN/aryl azide 0.75 M⁻¹s⁻¹ in aqueous acetonitrile:

| linkage | *t*<sub>95</sub> at 10 µM + 100 µM |
|---|---|
| tetrazine–sTCO | **0.4 s** |
| CuAAC (50 M⁻¹s⁻¹ at 100 µM Cu) | 10 min |
| maleimide–cysteine (10² M⁻¹s⁻¹) | 5 min |
| SPAAC, DBCO + AzK | **9.8 h** |
| oxime, aniline-catalysed (10⁻² M⁻¹s⁻¹) | 35 days |

SPAAC is genuinely the cleanest chemistry on the list — azide, alkyne and triazole are
all deep-UV absorbers, and a triazole-linked GFP heterodimer keeps normal photophysics
(Pope 2021 route; Hartley 2016 for azF at 148). It is not disqualified by cleanliness but
by the clock: an overnight incubation at ten-fold excess of a doubly-modified duplex,
during which the reagent is also free to cap one end and stop. **One pairing detail worth
having: DBCO pairs with alkyl azides and BCN with aryl azides. Getting it the wrong way
round costs an order of magnitude** (Dommerholt 2016, Table 2).

### The spectroscopic case for switching collapses

The reason to move away from tetrazine is that it absorbs 450–570 nm, peaking ~515 nm,
right on sfGFP's emission band, with *R*<sub>0</sub> = **2.6 nm** and **89 % FRET** at
site 150 — measured, with FQY falling 0.65 → 0.070 (Tang 2022). That is a 9.1-fold
quench, and it is real.

**But it is a property of the *unreacted* state only.** Losing N₂ destroys the tetrazine's
6π azadiene: *"reaction with sTCO disrupts the tetrazine conjugation, changes the S₀–S₁
electronic transition to ππ\*, increases the oscillator strength, and blue-shifts the
reddest absorption peak to UV region … outside the GFP chromophore emission range"*
(Tang 2022). The dihydropyridazine is colourless, and the quantum yields come back:

| construct | FQY unreacted | FQY clicked | turn-on |
|---|---|---|---|
| GFP134-Tet (Asp133) | 0.53 ± 0.02 | 0.62 ± 0.02 | 1.2 ± 0.1 |
| GFP150-Tet (Asn149) | 0.070 ± 0.002 | 0.64 ± 0.02 | 9.1 ± 0.4 |
| sfGFP reference | — | 0.65 | — |

**Our measurement is made on the clicked construct, where the quencher no longer exists.**
0.62 and 0.64 against 0.65 is full recovery. So switching chemistry to escape the quench
trades a solved problem for unsolved ones — and this retracts what I said earlier, that
the case for a different linker had to be made on spectroscopy. It cannot be made there
either.

The residual photoinduced-electron-transfer quench of the dihydropyridazine product is
real (Pinto-Pacheco 2020 — and it is the one argument for pairing tetrazine with a
strained *alkyne* instead, which gives the aromatic pyridazine and *"better … maximum
contrast ratio"*), but PET needs orbital contact and sfGFP's chromophore is buried 15.7 Å
from a surface adduct. Geometrically excluded here, which is exactly why the FQYs recover.

### The one inversion that would move the quencher off the protein — and why not to

Put the tetrazine on the DNA and a genetically encoded strained alkene on the protein, and
the absorber ends up on the partner that is consumed, added in excess, and washed away.
This is not speculative: **TCO\*-Lys is encoded at sfGFP position 150** by an engineered
PylRS (*"Y271A, L274M, C313A"*, "TCORS") in *E. coli* and in HEK-293 (Lang 2012), it is
*"stable for several days under physiological conditions"*, and rates stay usable —
Beliu 2019 measured TCO against two of their tetrazine conjugates at
**5,235 ± 258 and 17,248 ± 3132 M⁻¹s⁻¹**, 4–14× below 72,500 but still *t*½ ≈ 12 s at
5 µM, so completeness survives.

**It is disqualified by the bond it makes.** TCO\*-Lys attaches its cyclooctene to the
lysine ε-amine through a **carbamate**, and tetrazine + TCO-carbamate is the field's
standard *click-to-release* linkage. Beliu 2019 watched it happen to TCO\*-Lys itself —
*"a slight but highly reproducible decrease in fluorescence intensity after the initial
strong increase … attributed most likely to tautomerization of the TCO\*-tetrazine
conjugate resulting in the release of lysine by decarboxylative elimination"* — and
linkers deliberately engineered around that motif reach *"up to 96 % elimination after
2 hours"* (Keppel 2022) — an upper bound on how fast such a bond *can* break, not a
prediction for TCO\*-Lys, but the relevant bound when the same bond is asked to hold a
load for hours. **Building a load-bearing anchor out of the bond the field uses when it
wants a bond to break** is not a trade worth making for a load-bearing anchor.

Two further costs, for completeness. It adds ~4 rotatable sp³ bonds where the current
adduct has none — the wrong direction, given the section below. And it moves the
quencher onto the **moving** part of the construct: a tetrazine on the DNA arm sits at a
chromophore distance that *changes* when the duplex hybridises or the spring loads,
manufacturing a force-dependent FRET signal indistinguishable from the effect being
measured. The current quench's one virtue is that it is fixed.

### What would actually help: stiffen the arm, don't change the click

| lever | force | change |
|---|---|---|
| as built (compliant arm) | 5.52 pN | — |
| **stiffen the DNA-side arm at constant contour** | **8.47 pN** | **+53 %** |
| 27 → 25 bp | 6.72 pN | +22 % |
| best → worst chemistry swap | — | −0.8 pN |

The one genuinely compliant segment in the whole linkage is the DNA-side amino-modifier
arm; the protein half is fused-ring rigid. Replacing a flexible C6/C12 amino-modifier with
a rigid aromatic or alkynyl spacer of the same contour is worth more than any chemistry on
the list, and it is the cheapest component in the construct to change — one oligo order.

**The constraint to respect while doing it.** Spampinato 2023 compared TCO- and
BCN-modified nucleotides against tetrazine dyes and found *"the longer linker is crucial
for efficient labeling"* — so **shortening** the arm to raise the force costs reaction
efficiency, coupling force directly against completeness. But *stiffening at constant
contour is not shortening*, and that option appears unexplored.

### Hazards specific to this chemistry that the experiment must control

1. **The tetrazine's redox state mimics a successful click.** Reduction to the
   dihydrotetrazine both destroys reactivity *and* relieves the quench, so a chemically
   dead protein looks spectroscopically identical to a clicked one:
   *"when it is reduced, the sfGFP<sup>Tet2-Et</sup> will be fluorescent green"*
   (Eddins 2024). Oxidation back is fast after desalting — *"a few minutes (up to 1 h at
   most)"* — but **fluorescence recovery alone cannot be the completeness readout.** Use
   the orthogonal gel shift with sTCO-PEG5000, as Eddins 2023 does — that paper's own
   finding is that the dominant cause of unreactive protein is **near-cognate suppression
   of the amber codon**, not the chemistry at all.
2. **Do not run long 400–405 nm excitation on unclicked samples.** Tetrazine photolyses
   there, producing a spurious turn-on that mimics clicking — at exactly the wavelength
   the A-band measurement needs.
3. **A doubly-sTCO duplex is a double-headed sTCO reagent, and that class misbehaves on
   this exact protein.** Jang 2020: dh-sTCO *"only create protein dimers if the dh-sTCO is
   the limiting reagent, because an excessive concentration … results in each Tet-protein
   with a singly reacted dh-sTCO molecule."* Intermolecular crosslinks when limiting,
   single-end caps when in excess — neither is the intended 1:1 loop, and no titration
   window is published. **The completeness number is the yield of correctly looped monomer
   by SEC or native gel, not tetrazine consumption.**
4. **Tet2-Et at Asn149 perturbs the very observable we are reading.** Incorporation at
   site 150 *"reduces the chromophore pKa (suppression of protonated state at same buffer
   pH) and blue-shifts the deprotonated state absorption peak by 4 nm. Reaction with sTCO
   partially mitigates these effects"* (Tang 2022). Our readout is the A-band ratio, which
   is a direct measure of the protonated-state population. A 4 nm shift and a pKa change,
   only partly reversed on clicking, are far larger than any mechanical effect — which is
   the same trap as using unclicked protein as the zero, and the reason recommendation 4
   below is the isogenic *clicked*, spring-free twin. At Asp133 the perturbation is
   negligible: *"the absorption profile of GFP134-Tet almost matches sfGFP."*
5. **Do not chase the S202Y fluorogenicity mutant** (31-fold turn-on) or the Y151A/Y200A
   variants. The Y→A openings admit solvent to the chromophore — on-state FQY falls to
   0.56–0.57 with a new ps solvent-quenching channel — and S202Y is a mutation in the
   chromophore's own pocket, which corrupts the electrostatics the excitation ratio reads.

### Gaps I could not close

No molar extinction coefficient is published anywhere for the Tet-v2.0/Tet2-Et **visible**
band; two independent derivations from Tang 2022's own reported quantities bracket it at
**ε(515 nm) ≈ 4–8 × 10² M⁻¹cm⁻¹**, which makes the tetrazine a negligible inner filter but
an efficient intramolecular FRET acceptor. No rate constant exists for Tet2-Et (the ethyl
variant) separately from Tet-v2.0-Me. No crystal structure of a tetrazine-encoded sfGFP
exists, so our 12.4–13.3 Å adduct span has nothing to be validated against — Tang 2022 say
so and fall back on 50 ns MD themselves. And **no sTCO-oligonucleotide has been coupled to
a genetically encoded tetrazine protein anywhere in the literature**: TCO-on-DNA is
established, but always with an ordinary TCO or BCN reacting with a small tetrazine dye,
i.e. the reverse polarity from this construct.

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
| **ss** | 0 pN | **0.1278** ± 0.0014 | slack strand, no load |
| **nicked** | ~2.5 pN | **0.1195** ± 0.0003 | |
| **ds** | ~3.1 pN | **0.1159** ± 0.0005 | |

**Monotone in modelled spring force, 9.3 % relative
decrease**, with the emission peak unchanged. DNA-only blanks are
0.35–0.46 %
of the signal at 400 nm, so it is not scatter. "Band ratio moves, emission λ fixed" is
exactly the signature the physics predicts.

**Do not believe it yet.** The three a/b/c scans per condition are *technical* replicates
of one preparation, so the quoted scatter is instrument repeatability and says nothing
about the constructs. Absolute brightness differs across the samples, and pH,
concentration and inner-filter effects are not controlled at the 10⁻³ level this needs.
Tet2-Et also **quenches sfGFP at N150, with clicking relieving the quench** — a
first-order confound for any amplitude measurement on these constructs.

### The one cheap measurement that would settle it

**Extend the excitation scan.** It currently runs 400–500 nm
at 1.0 nm, which **starts exactly at the A-band maximum** — the band
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
([see the report](https://github.com/extremefattypunch/DNA-springs-project-2026-27/blob/main/docs/index.html)). Going from 5 to 35 pN moves the
prediction from ~0.05 to ~0.35 nm — across, rather than below, the detection threshold.
A shorter, stiffer duplex costs one oligonucleotide. Measuring the force the *chimera*
actually delivers is cheaper than building four constructs, and it gates everything else.

**The cheapest single change available is the DNA-side arm.** Replacing the flexible
C6/C12 amino-modifier with a rigid spacer of the same contour is worth
**+53 %** in force (5.5 → 8.5 pN), which is more than
any change of click chemistry can offer and more than going 27 → 25 bp. It is also one
oligo order. [The argument is below.](#what-would-actually-help-stiffen-the-arm-dont-change-the-click)

**2. Read a band ratio, not a peak position.** Follow the signal above: fit the 400 and
490 nm bands and report their ratio at fixed emission wavelength. This is the observable
the physics supports, the one the existing data already hints at, and the one a
two-state population shift would move.

**3. If you do re-site, Tyr151 ↔ Ser208.** At
27.9 Å
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

### What earlier versions of this page got wrong

**T203Y is not the amplifier it looks like.** Its red shift comes mainly from the **loss
of a hydrogen bond**, not from π-stacking geometry: the aliphatic substitutions T203V and
T203I produce a similar shift (Drobizhev 2015). Mechanically modulating a stacking
distance therefore would not recover it, and an earlier version of this page recommended
a **T203Y-sensitised background on that mistaken basis**.

**"If you change the chemistry, do it for spectroscopy."** I said that in answer to
the question about alternative click linkers, on the grounds that tetrazine absorbs where
sfGFP emits. **The quench is a property of the unreacted tetrazine only** — the clicked
dihydropyridazine is colourless and the quantum yield returns to 0.62–0.64 against
sfGFP's 0.65. Since every measurement is made on the clicked construct, there is no
spectroscopic case for switching either. Both premises of that question fail, and the
[section above](#would-a-different-click-chemistry-help) is the corrected answer.

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
on arithmetic of ours, it says so at that point. References 23–34 were located and their
bibliographic details confirmed through **PubMed**; every DOI below resolves.

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

**The click chemistry: rates, encoding, and photophysics**
23. Blizzard RJ, Backus DR, Brown W, Bazewicz CG, Li Y, Mehl RA (2015) Ideal bioorthogonal reactions using a site-specifically encoded tetrazine amino acid. *J Am Chem Soc* 137(32):10044 — Tet-v2.0 + sTCO at **72,500 ± 1660 M⁻¹s⁻¹**, measured on GFP-bound tetrazine in cellulo. [10.1021/jacs.5b03275](https://doi.org/10.1021/jacs.5b03275)
24. Tang L, Bednar RM, Rozanov ND, Hemshorn ML, Mehl RA, Fang C (2022) Rational design for high bioorthogonal fluorogenicity of tetrazine-encoded green fluorescent proteins. *Nat Sci (Weinh)* 2:e20220028 — the full photophysics: tetrazine absorbs 450–570 nm, *R*₀ = 2.6 nm, 89 % FRET at site 150, FQY 0.070 → 0.64 on clicking, and the product's *"reddest absorption peak"* blue-shifted *"to UV region"*. Also the 4 nm shift and pKa change at site 150. [10.1002/ntls.20220028](https://doi.org/10.1002/ntls.20220028)
25. Eddins AJ, Pung AH, Cooley RB, Mehl RA (2024) Tetrazine amino acid encoding for rapid and complete protein bioconjugation. *Bio-protocol* 14(16):e5048 — Tet2-Et, the pAJE-E7 machinery, >95 % encoding and labeling fidelity, ~200 mg L⁻¹; the redox equilibrium that makes a reduced protein look clicked. [10.21769/BioProtoc.5048](https://doi.org/10.21769/BioProtoc.5048)
26. Eddins AJ, Bednar RM, Jana S, Pung AH, Mbengi L, Meyer K, Perona JJ, Cooley RB, Karplus PA, Mehl RA (2023) Truncation-free genetic code expansion with tetrazine amino acids for quantitative protein ligations. *Bioconjug Chem* 34(12):2243 — near-cognate suppression as the real cause of unreactive protein, and the sTCO-PEG5000 gel shift as the orthogonal completeness assay. [10.1021/acs.bioconjchem.3c00380](https://doi.org/10.1021/acs.bioconjchem.3c00380)
27. Jang HS, Jana S, Blizzard RJ, Meeuwsen JC, Mehl RA (2020) Access to faster eukaryotic cell labeling with encoded tetrazine amino acids. *J Am Chem Soc* 142(16):7245 — Tet-v3.0; and the **double-headed sTCO** result: dimers when limiting, single-end caps in excess. [10.1021/jacs.9b11520](https://doi.org/10.1021/jacs.9b11520)
28. Lang K, Davis L, Wallace S, Mahesh M, Cox DJ, Blackman ML, Fox JM, Chin JW (2012) Genetic encoding of bicyclononynes and *trans*-cyclooctenes for site-specific protein labeling. *J Am Chem Soc* 134(25):10317 — TCO\*-Lys encoded at **sfGFP position 150** by TCORS (Y271A/L274M/C313A). [10.1021/ja302832g](https://doi.org/10.1021/ja302832g)
29. Beliu G, Kurz AJ, Kuhlemann AC, Behringer-Pliess L, Meub M, Wolf N, Seibel J, Shi Z-D, Schnermann M, Grimm JB, Lavis LD, Doose S, Sauer M (2019) Bioorthogonal labeling with tetrazine-dyes for super-resolution microscopy. *Commun Biol* 2:261 — TCO\*-Lys *"stable for several days"*; tetrazine absorption *"peaking at ~515 nm"*; and the **decarboxylative elimination** that releases the lysine. [10.1038/s42003-019-0518-z](https://doi.org/10.1038/s42003-019-0518-z)
30. Keppel P, Sohr B, Kuba W, Goldeck M, Skrinjar P, Carlson JCT, Mikula H (2022) Tetrazine-triggered bioorthogonal cleavage of *trans*-cyclooctene-caged phenols using a minimal self-immolative linker strategy. *ChemBioChem* 23(20):e202200363 — *"up to 96 % elimination after 2 hours"*. **Read as an upper bound, not a prediction:** this is a linker deliberately engineered to eliminate. It bounds how fast a TCO-carbamate *can* break, which is the relevant number when the same bond is asked to hold a load. [10.1002/cbic.202200363](https://doi.org/10.1002/cbic.202200363)
31. Pinto-Pacheco B, Carbery WP, Khan S, Turner DB, Buccella D (2020) Fluorescence quenching effects of tetrazines and their Diels-Alder products: mechanistic insight toward fluorogenic efficiency. *Angew Chem Int Ed* 59(49):22140 — the dihydropyridazine's residual **PET** quench, and why a strained *alkyne* gives *"maximum contrast ratio"*. [10.1002/anie.202008757](https://doi.org/10.1002/anie.202008757)
32. Dommerholt J, Rutjes FPJT, van Delft FL (2016) Strain-promoted 1,3-dipolar cycloaddition of cycloalkynes and organic azides. *Top Curr Chem (Cham)* 374:16 — a **review**, so its rate tables are a compilation rather than one measurement — the primary SPAAC rate tables: BCN 0.14, DIBAC/DBCO 0.31 M⁻¹s⁻¹ with benzyl azide, and the DBCO↔alkyl / BCN↔aryl pairing rule. [10.1007/s41061-016-0016-4](https://doi.org/10.1007/s41061-016-0016-4)
33. Spampinato A, Kuzmová E, Pohl R, Sýkorová V, Vrábel M, Kraus T, Hocek M (2023) *trans*-Cyclooctene- and bicyclononyne-linked nucleotides for click modification of DNA with fluorogenic tetrazines. *Bioconjug Chem* 34(4):772 — *"the longer linker is crucial for efficient labeling"*: shortening the DNA arm costs completeness. [10.1021/acs.bioconjchem.3c00064](https://doi.org/10.1021/acs.bioconjchem.3c00064)
34. Uttamapinant C, Tangpeerachaikul A, Grecian S, Clarke S, Singh U, Slade P, Gee KR, Ting AY (2012) Fast, cell-compatible click chemistry with copper-chelating azides for biomolecular labeling. *Angew Chem Int Ed* 51(24):5852 — CuAAC at *"10–100 M⁻¹s⁻¹ per 10–100 µM Cu"*, and why a defensible protein-scale *k*₂ does not exist for it. [10.1002/anie.201108181](https://doi.org/10.1002/anie.201108181)

**Cited in the original proposal, and not usable for the claim made**
- Suto K *et al.* (2009) *Genes Cells* 14:727 — **wrong protein.** CpYGFP from the copepod *Chiridius poppei*, cited for T203Y in *Aequorea* GFP. [10.1111/j.1365-2443.2009.01305.x](https://doi.org/10.1111/j.1365-2443.2009.01305.x)
- Wang X *et al.* (2013) *J Mol Med* 91:917 — a **review** of ROS imaging; contains neither the 6.1-fold nor the 4.3-fold figure, and does not mention N149C/S202C, His148 or Thr203. [10.1007/s00109-013-1067-4](https://doi.org/10.1007/s00109-013-1067-4)
- Marchetti M *et al.* (2023) *Sensors* 23:8517 — a **review**; states the roGFP1 genotype correctly but is not where it was measured. [10.3390/s23208517](https://doi.org/10.3390/s23208517)
