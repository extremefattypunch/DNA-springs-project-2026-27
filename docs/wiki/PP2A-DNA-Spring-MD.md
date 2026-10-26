# DNA spring on PP2A Aα — all-atom MD

**[code and structures](https://github.com/extremefattypunch/DNA-springs-project-2026-27/tree/main/MD%20simulations/pp2a)**
 · [sfGFP page — the pipeline this reuses](Home)
 · [PP2A cryo-EM & MD reproduction — the apo baseline](PP2A-Cryo-EM-and-MD-Reproduction)
 · [how the sfGFP chimera was built](https://github.com/extremefattypunch/DNA-springs-project-2026-27/blob/main/feedback/2026-8-22_setup-tutorial.md)

The scaffold subunit of PP2A (Aα, PR65, UniProt P30153) is a 589-residue horseshoe of
fifteen HEAT repeats. Two amber codons put the tetrazine amino acid **Tet2-Et** at the two
tips of the horseshoe — construct positions **D6X and N552X**, which are **Asp7 and Asn553**
in UniProt numbering (the construct counts from Ala2) — and a single 50 nt DNA strand,
sTCO-modified at both ends, clicks onto both, closing a covalent loop. The monomer is
simulated alone, without the B and C subunits.

A single strand is not the leaf spring the duplex constructs are. It is an entropic chain
that pulls its ends **together** at any separation — with a tension 3k<sub>B</sub>T·x/(2ξL) of about
1.3 pN at the 90 Å build span (Marko–Siggia 1.6 pN) that climbs steeply only as the span approaches
the contour. For 50 nucleotides the free chain's rms size is 94 Å
(contour 29.5 nm, persistence 1.5 nm); it says whether the strand is slack or taut,
not whether it pulls. The sfGFP page quotes the tension as zero below that size — a clipping
convention of that analysis, kept here only as a second column. The two Cβ atoms are 54.4 Å apart in the AlphaFold model and the
tethers reach a further 22 and 26 Å, so the strand *can* be stretched to ~102 Å by the
protein — but only if the horseshoe opens.

**400 ns over 4 systems**, one trajectory each.

![the 50 nt single-strand chimera](https://extremefattypunch.github.io/DNA-springs-project-2026-27/figures/pp2a_P3_ss50_overview_still.png)
*the 50 nt single-strand chimera.*

## The headline

With the 50 nt strand clicked on, the two attachment sites sit **3.1 ± 3.5 Å together** relative to the clicked protein without DNA (p = 0.39 on ten block means each, not resolved). The strand itself spends **52 %** of the run beyond its own natural size, so by the model it is **matched** and pulls with **1.74 pN** (Marko–Siggia, Zocchi eq. 2.207; ideal chain 1.36 pN; 0.98 pN under the sfGFP page's zero-below-rms convention).

Zocchi's construct is the *pair*: the same strand clicked at both ends, unhybridised (an entropic
coil that can only pull) and hybridised (a leaf spring that pushes its anchors apart — the book's
point that the duplex, with its 50 nm persistence length, delivers far more force than the coil).
Both are run at 50 and 80 nt, so the two axes of the design can be read separately:

| chimera | DNA end-to-end (Å) | Cβ7–Cβ553 (Å) | Δ vs clicked zero (Å) | modelled force (pN; + push, − pull) | Zocchi f(x) a duplex of this length would give at this span | nucleotides on the protein | kinking |
|---|---|---|---|---|---|---|---|
| 50 nt ssDNA, closed loop | 95 ± 9 | 77.8 | -3.1 ± 3.5 (p = 0.39) | −1.74 (Marko–Siggia, matched) | 2.05 (kinked) | 3 of 50 | 0.18 turns/frame, 11.5°/nt |



![the attachment-point separation against the clicked zero](https://extremefattypunch.github.io/DNA-springs-project-2026-27/figures/pp2a_fig7_headline.png)

| system | Cβ7–Cβ553 (Å), mean ± sd | range | Δ vs clicked, no DNA (Å) | p |
|---|---|---|---|---|
| apo Aα (ff19SB/OPC) | 56.4 ± 7.1 | 40–75 | -24.40 ± 2.99 | < 0.001 |
| 2× Tet2-Et, unclicked | 78.6 ± 7.6 | 56–94 | -2.32 ± 3.11 | 0.47 |
| clicked, no DNA | 80.9 ± 8.9 | 52–102 | zero | — |
| 50 nt ssDNA, closed loop | 77.8 ± 9.1 | 54–100 | -3.08 ± 3.51 | 0.39 |

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
| apo Aα (ff19SB/OPC) | 2.72 (max 6.1) | 11.6 | 2.11 | 6.8 / 3.9 | 0.780 |
| 2× Tet2-Et, unclicked | 4.88 (max 8.8) | 12.8 | 2.57 | 6.1 / 4.4 | 0.775 |
| clicked, no DNA | 5.99 (max 11.0) | 17.2 | 2.41 | 7.5 / 4.2 | 0.772 |
| 50 nt ssDNA, closed loop | 4.81 (max 8.5) | 10.4 | 2.26 | 6.4 / 4.6 | 0.773 |

The chimera's core amplitude (RMSD about the run's own mean, per block) differs from the clicked control's by -0.12 ± 0.26 Å (p = 0.65 on block means).

![RMSD, Day et al.'s split](https://extremefattypunch.github.io/DNA-springs-project-2026-27/figures/pp2a_fig1_rmsd.png)

![per-residue RMSF with the HEAT repeats](https://extremefattypunch.github.io/DNA-springs-project-2026-27/figures/pp2a_fig2_rmsf.png)

**The horseshoe as a whole.** The solenoid's axis is the chain of HEAT-repeat Cα centroids; the
bend at a repeat is the angle between the incoming and outgoing axis segments, the twist the
dihedral across four consecutive centroids. The literature puts PR65's hinge at HEAT 10–11
(crystal comparison), 11–12 (Day) or 12–13 (Yilmaz), with a diffuse expansion of repeats 4–12.

| system | total axis bend (°) | HEAT1–15 (Å) | Asn30–Phe578 (Å) | region RMSD N / central / C (Å) | most variable repeat | Arg113–Glu191 salt bridge formed |
|---|---|---|---|---|---|---|
| apo Aα (ff19SB/OPC) | 150 ± 9 | 53.5 ± 6.7 | 42.0 ± 5.6 | 1.19 / 1.48 / 0.73 | H12 | 0.82 |
| 2× Tet2-Et, unclicked | 159 ± 7 | 73.5 ± 7.6 | 60.5 ± 7.7 | 1.44 / 2.08 / 0.81 | H12 | 0.94 |
| clicked, no DNA | 165 ± 8 | 73.0 ± 7.5 | 59.2 ± 7.4 | 1.20 / 2.96 / 0.76 | H12 | 0.29 |
| 50 nt ssDNA, closed loop | 159 ± 7 | 73.9 ± 7.5 | 61.0 ± 7.1 | 1.31 / 2.23 / 0.72 | H12 | 0.97 |

Against the clicked control, the repeats whose bend angle moves most under the strand are HEAT 5 (-1.6°, p = < 0.001, spread ×0.94), HEAT 13 (-2.6°, p = 0.002, spread ×0.89), HEAT 2 (+3.7°, p = 0.009, spread ×1.23). The horseshoe aperture (HEAT1–HEAT15) changes by +0.8 ± 2.8 Å (p = 0.77).

![where the solenoid bends and twists](https://extremefattypunch.github.io/DNA-springs-project-2026-27/figures/pp2a_fig5_heat_bend.png)

## The distance between the attachment points

Three distances are reported because three groups measure three different things:
Cβ7–Cβ553 (ours, the spring's own coordinate), Cα Asn30–Phe578 (Yilmaz et al. 2025's
"N29–F577", 66.5 ± 9.7 Å in their 3 × 704 ns apo runs, 47.7 Å in the compact 6NTS trimer,
70.3 Å in the extended 1B3U crystal) and the HEAT1–HEAT15 centroid distance (Day et al.,
62.7 Å in the B56γ holoenzyme). The AlphaFold start model is compact on all three
(Asn30–Phe578 39 Å).

![distance between the attachment points](https://extremefattypunch.github.io/DNA-springs-project-2026-27/figures/pp2a_fig3_end_to_end.png)

## How stiff the protein is along the spring's coordinate

Zocchi's two-spring picture (book pp. 117–119, eq. 3.100) gives the protein a stiffness κ ≈ 100 k<sub>B</sub>T nm⁻²
along the line between the attachment points on its linear branch, softening at a strain x<sub>e</sub> ≈ 3 Å
(κx<sub>e</sub> ≈ 30 k<sub>B</sub>T nm⁻¹ ≈ 125 pN) to a constant restoring force f<sub>e</sub> ≈ 1 k<sub>B</sub>T nm⁻¹ ≈ 4 pN;
his criterion for a spring to matter is that it injects ≳ 1 k<sub>B</sub>T into the enzyme (p. 117). sfGFP's
force-clamp ladder measured 230–350 k<sub>B</sub>T nm⁻². No clamp ladder was run here, but equipartition on
the Cβ7–Cβ553 fluctuation of each run, κ = k<sub>B</sub>T / var(x), gives the stiffness of that coordinate, N-terminal
tail included; the HEAT1–HEAT15 aperture alone gives nearly the same figure, so the softness is the
horseshoe's. It is a single-run estimate whose bias has no fixed sign: the slow opening and closing
of the horseshoe *is* the compliance, and 100 ns samples it only partly.

| system | Cβ7–Cβ553 (Å) | κ (k<sub>B</sub>T nm⁻²) | κ from the HEAT1–15 aperture | κ (pN nm⁻¹) | displacement per pN (Å) | energy a 2 pN spring would store (k<sub>B</sub>T) |
|---|---|---|---|---|---|---|
| apo Aα (ff19SB/OPC) | 56.4 ± 7.0 | 2.0 | 2.2 | 8.4 | 1.2 | 0.06 |
| 2× Tet2-Et, unclicked | 78.5 ± 7.6 | 1.8 | 1.7 | 7.3 | 1.4 | 0.07 |
| clicked, no DNA | 80.9 ± 8.9 | 1.2 | 1.8 | 5.3 | 1.9 | 0.09 |
| 50 nt ssDNA, closed loop | 77.8 ± 9.1 | 1.2 | 1.8 | 5.0 | 2.0 | 0.10 |

This is the number that decides whether a spring of a few piconewtons is visible at all. On sfGFP
(κ ≈ 300) a 4 pN spring moved the sites 0.03 Å and stored 0.002 k<sub>B</sub>T, below what 50 ns can
resolve; a HEAT horseshoe is two orders of magnitude softer along this coordinate, so even the ~2 pN
a single strand could deliver would move the tips by ångströms and store a fair fraction of
k<sub>B</sub>T — if it were under tension.

## What a spring would deliver, by Zocchi's model

Zocchi's kinked-rod model (eqs 2.188–2.206 and 3.51–3.52 of the book; Tseng & Zocchi 2013, eq. 6;
B = 200 pN nm², τ<sub>c</sub> = 27 pN nm nicked / 34 continuous, 0.33 nm per bp, k<sub>B</sub>T = 4.2 pN nm),
evaluated with the repository's implementation at the DNA end-to-end distance the chimera settles at, 95 ± 9 Å. A duplex pushes its ends apart on the
bent branch, f = 10B/L² − k<sub>B</sub>T/(L − x), until the span drops below x<sub>c</sub>, where it
kinks and the force collapses to τ<sub>c</sub> / (2R √(1 − (x/2R)²)):

| duplex | contour (nm) | regime at this span | f(x) (pN) | unkinked rod (pN) | Euler load 10B/L² (pN) | stored energy (k<sub>B</sub>T) | x<sub>c</sub> (Å) | γ = Lτ<sub>c</sub>/2B |
|---|---|---|---|---|---|---|---|---|
| 27 bp, nicked | 8.9 | *infeasible: contour shorter than the span* | — | — | 25.2 | — | 87 | 0.60 |
| 27 bp, continuous | 8.9 | *infeasible: contour shorter than the span* | — | — | 25.2 | — | 86 | 0.76 |
| 40 bp, nicked τc = 27 pN nm | 13.2 | kinked | **3.01** | 10.3 | 11.5 | 4.8 | 125 | 0.89 |
| 40 bp, continuous τc = 34 pN nm | 13.2 | kinked | **3.84** | 10.3 | 11.5 | 6.1 | 121 | 1.12 * |
| 50 bp, nicked τc = 27 pN nm | 16.5 | kinked | **2.05** | 6.7 | 7.3 | 6.1 | 151 | 1.11 * |
| 50 bp, continuous τc = 34 pN nm | 16.5 | kinked | **2.61** | 6.7 | 7.3 | 7.6 | 143 | 1.40 * |
| 80 bp, nicked τc = 27 pN nm | 26.4 | kinked | **1.14** | 2.6 | 2.9 | 7.6 | 208 | 1.78 * |
| 80 bp, continuous τc = 34 pN nm | 26.4 | kinked | **1.48** | 2.6 | 2.9 | 9.5 | 175 | 2.24 * |

\* γ ≥ 1: beyond the small-γ condition the book states after eq. 2.188 (≈ 45 bp at τ<sub>c</sub> ≈ 30 pN nm).
Zocchi applies eq. 3.52 in this regime himself — his nicked 60-mer (eq. 3.53, 1.5 pN) sits at γ = 1.34 —
so the star marks an extrapolation, not a disqualification. x<sub>c</sub> is eq. 2.194, the small-γ
estimate, as on the sfGFP page.

The **unkinked rod** column is the bent-branch force continued below x<sub>c</sub>, and it is there
because the force field did not follow the model on sfGFP: its bare-duplex restraint windows measured
**35.3 ± 1.7 pN** for the 27 bp and **16.4 ± 1.8 pN** for the 40 bp — about 1.3× their Euler loads, on
the *unyielded* branch the model says the duplex has left. On MD timescales the bold kinked-branch
f(x) is therefore an underestimate by roughly 3–6×.

The chimera carries none of these. Its 50 nt **single strand** is Zocchi's "ss state", of which his
2021 chapter says *"with the DNA spring in the ss state there is essentially no stress on the
enzyme"* — true in his geometry, where a 60-mer bridges 6.1 nm, a fraction of the chain's rms size.
Here the strand is a worm-like chain with ξ = 1.5 nm and a contour of 29.5 nm whose rms size is
94 Å, and the span is a comparable 95 Å: **matched**, pulling with
**1.73 pN** by Marko–Siggia (eq. 2.207; 1.73 pN under the
sfGFP page's zero-below-rms convention). With the book's own ℓ<sub>s</sub> ≈ 0.8 nm instead of the 2021 chapter's
1.5 nm the rms size shrinks to 69 Å and the pull rises to
3.25 pN — the verdict "slack" rests on which of Zocchi's two persistence lengths one takes.
Every number in this section is a **model**, not a measurement.

![Zocchi's spring model at this span](https://extremefattypunch.github.io/DNA-springs-project-2026-27/figures/pp2a_fig9_zocchi_model.png)

## What it delivers, measured

**The DNA in the chimeras.** The strand's end-to-end distance settles at **95 ± 9 Å**
(range 75–122 Å), against a contour of 30 nm and a
natural rms size of 94 Å. Its radius of gyration is 38 Å. It makes
0.2 sharp (> 90°) turns per frame along its C4′ trace, hottest at nucleotide
3, keeps 78 % of adjacent bases stacked, and has
**3 of 50 nucleotides on the protein** at any moment
(6 % of the strand, time-averaged). Protein residues it touches more than half the
time: TDP553;MET1;TDP7;ALA2;ASP5;ALA3;ALA4;SER554;GLN557.

Marko–Siggia (Zocchi eq. 2.207) applied frame by frame to the strand's end-to-end distance gives
**1.74 pN** of pull; the Gaussian-chain limit 3k<sub>B</sub>T·x/(2ξL) (eq. 2.173)
gives 1.36 pN; the sfGFP page's convention of calling the tension zero below the
chain's rms size gives 0.98 pN. All three are **models**, not measurements. For comparison, Zocchi's duplex model would put a 50 bp
spring at this span on its kinked branch at 2.0–2.6 pN of *push* (γ = 1.11–1.40, beyond the book's stated
γ < 1 condition but the regime Zocchi's own 60-mer sits in). The sign is the point: nothing here
pushes.

![the strand's size and modelled tension over time](https://extremefattypunch.github.io/DNA-springs-project-2026-27/figures/pp2a_fig8_strand.png)

**The bare strand under restraint.** Each bare spring — the 50 and 80 nt strands and the 50 and 80 bp duplexes, no protein — held at a series of
end-to-end separations by a soft harmonic restraint (k = 33 kJ mol⁻¹ nm⁻²) between the atoms the linkers
would bond to, exactly as the sfGFP duplexes were. The mean restraint force ⟨k(r − x₀)⟩ *is* the force the strand
exerts, read as ⟨k(x₀ − r)⟩ so that positive means the strand **pulls** its ends together
(`force_extension.csv` stores the opposite, push-positive, sign).

| spring | target x₀ (Å) | settled ⟨r⟩ (Å) | measured force (pN; + push, − pull) | model at ⟨r⟩ (pN) | kinking in the window | box valid? |
|---|---|---|---|---|---|---|
| ss50 | 60 | 60.3 ± 2.8 | **+1.52 ± 1.51** | −0.98 (Marko–Siggia, slack) | 0.02 turns/frame | yes |
| ss50 | 80 | 80.1 ± 2.7 | **+0.50 ± 1.42** | −1.38 (Marko–Siggia, matched) | 0.01 turns/frame | yes |

Where they can be compared: ss50 at 60 Å: force field +1.52 ± 1.51 pN against -0.98 pN by Marko–Siggia (within the error); ss50 at 80 Å: force field +0.50 ± 1.42 pN against -1.38 pN by Marko–Siggia (within the error).

![modelled against measured force for the single strand](https://extremefattypunch.github.io/DNA-springs-project-2026-27/figures/pp2a_fig4_force_extension.png)

## Does the DNA kink?

**By the model.** A rod of bending modulus B yields at the critical torque τ<sub>c</sub>, i.e. at a
curvature τ<sub>c</sub>/B = **2.55° per bp** (nicked) to
**3.21° per bp** (continuous) — Zocchi's kink. That is the curvature at the
rod's centre, where the torque peaks (the book's cubic shape, eq. 2.183, has zero curvature at the ends);
a uniform arc with the same chord bends only ≈ 0.63 τ<sub>c</sub>/B ≈ 1.6°/bp, the same for every length.
Below the softening span x<sub>c</sub> the duplex is two nearly straight arms meeting at a constant-torque
kink (book fig. 2.31); R is the chord of each half-rod from the kink to its anchor, R ≈ L/2 (eqs 2.189,
2.197), not a radius of curvature. The two curvature columns are chord-equivalent uniform arcs, for scale:

| duplex | x<sub>c</sub> (Å) | half-rod chord R (Å) | chord-equivalent uniform curvature at x<sub>c</sub> (°/bp) | …at this span (°/bp) |
|---|---|---|---|---|
| 27 bp | 87 | 44 | 1.62 | stretched |
| 40 bp | 125 | 65 | 1.63 | 3.87 |
| 50 bp | 151 | 81 | 1.64 | 3.93 |
| 80 bp | 208 | 127 | 1.67 | 3.17 |

The force field disagrees with the model on *where* a duplex kinks: on sfGFP's bare duplexes base pairs
began to open at 5.5°/bp and opened in every run only above 8.7°/bp
(the chimera there sat at 4.8-5.2°/bp and never kinked; the one directly comparable all-atom
calculation nucleates a type II kink at 8.6-9.3°/bp). A single strand has no base pairs to
open, so for it "kink" means a sharp turn of the backbone: a C4′-trace angle above 90°, a BI → BII
switch of the ε/ζ torsions, or a step whose bases have come apart.

**Measured on the strand — Barely.** The strand makes a > 90° backbone turn in
**18 % of frames** (0.18 per frame); its hottest nucleotide,
3 (CpT-TpT steps), turns in 12 % of them — a terminus bonded to the tether, not an interior kink. A sharp turn coincides with an
unstacked step 100 % of the time. The smoothed-axis curvature is
**11.5° per nucleotide** (2.6× what the same construction gives on an ideal straight B-form strand, 4.4°/nt; the rest is the thermal bending of a chain with a 1.5 nm persistence length). 15 % of the backbone sits in the
BII (ε − ζ > 0) substate and 78 % of adjacent bases stay stacked. Whether stretching the
strand straightens it is read off the restraint windows (the squares in the last panel).

![does the DNA kink: bend, sharp turns, torsions, stacking along the strand](https://extremefattypunch.github.io/DNA-springs-project-2026-27/figures/pp2a_fig6_dna_kink.png)

## Ramachandran plots at the two Tet2-Et sites

On sfGFP the one site that sat at **positive φ** (left-handed α, Asp133) was pushed 17.5°
out of its basin by the modification; ordinary sites moved by a degree or less. In the
AlphaFold model of Aα, **Asp7 sits at φ, ψ = +65, +57** — the same left-handed-α position —
and Asn553 at -62, -27, in the first turn of HEAT 14's helix. Asp7 is also on the
low-confidence tail (pLDDT 52), so its backbone is free to leave that basin
whether modified or not.

| system | residue 7 | φ, ψ (± sd) | aromatic rings on the Tet ring / frame | residue 553 | φ, ψ (± sd) | rings / frame |
|---|---|---|---|---|---|---|
| apo Aα (ff19SB/OPC) | ASP | -89, +5 (± 43, 75) | — | ASN | -57, -37 (± 9, 9) | — |
| 2× Tet2-Et, unclicked | TET | -53, -46 (± 12, 14) | 0.43 | TET | -52, -46 (± 9, 10) | 0.00 |
| clicked, no DNA | TDP | +35, +53 (± 77, 58) | 0.00 | TDP | -57, -36 (± 10, 10) | 0.00 |
| 50 nt ssDNA, closed loop | TDP | +49, +33 (± 12, 15) | 0.39 | TDP | -59, -35 (± 9, 8) | 0.09 |

![backbone conformation at the two sites](https://extremefattypunch.github.io/DNA-springs-project-2026-27/figures/pp2a_ramachandran.png)

## Animations

![50 nt ssDNA, closed loop — overview](https://extremefattypunch.github.io/DNA-springs-project-2026-27/figures/pp2a_P3_ss50_overview.gif)
*the 50 nt coil: horseshoe, strand, tethers — overview view, 50 nt ssDNA, closed loop.*

![50 nt ssDNA, closed loop — sites](https://extremefattypunch.github.io/DNA-springs-project-2026-27/figures/pp2a_P3_ss50_sites.gif)
*the 50 nt coil: horseshoe, strand, tethers — sites view, 50 nt ssDNA, closed loop.*

![50 nt ssDNA, closed loop — strain](https://extremefattypunch.github.io/DNA-springs-project-2026-27/figures/pp2a_P3_ss50_strain.gif)
*the 50 nt coil: horseshoe, strand, tethers — strain view, 50 nt ssDNA, closed loop.*

![clicked, no DNA — overview](https://extremefattypunch.github.io/DNA-springs-project-2026-27/figures/pp2a_P2_clicked_overview.gif)
*the clicked control — overview view, clicked, no DNA.*

![clicked, no DNA — sites](https://extremefattypunch.github.io/DNA-springs-project-2026-27/figures/pp2a_P2_clicked_sites.gif)
*the clicked control — sites view, clicked, no DNA.*

![apo Aα (ff19SB/OPC) — overview](https://extremefattypunch.github.io/DNA-springs-project-2026-27/figures/pp2a_P0_apo_overview.gif)
*the apo monomer (cryo-EM reproduction run) — overview view, apo Aα (ff19SB/OPC).*

![apo Aα (ff19SB/OPC) — strain](https://extremefattypunch.github.io/DNA-springs-project-2026-27/figures/pp2a_P0_apo_strain.gif)
*the apo monomer (cryo-EM reproduction run) — strain view, apo Aα (ff19SB/OPC).*


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
| Start model | AlphaFold DB v6 model of P30153, full length 1–589 (mean pLDDT 94.9; residues 1–13 55) — the model the apo run used |
| Numbering | File numbering is UniProt P30153 (Met1 = 1).  The construct labels D6X/N552X count from Ala2 (construct n = UniProt n+1): D6X = Asp7, N552X = Asn553.  UniProt 6 is Gly (no CB) and 552 is Asp, so the literal reading is impossible.  Confirmed with Ian 2026-10-26.  Yilmaz et al.'s 'N29-F577' is the same convention = Asn30-Phe578 here. |
| Sites | Asp7 (SASA 85 Å², pLDDT 52) and Asn553 (SASA 106 Å², pLDDT 94); Cβ–Cβ 54.4 Å |
| Chemistry | Tet2-Et + sTCO click adduct (TDP), 5′-amino-C6 (DNL) and 3′-amino-C7 (DNP) carbamate linkers, parameters unchanged from the sfGFP build (GAFF2 / AM1-BCC, ff14SB backbone) |
| Strand | 50 nt, `GCTTGGATGGTACCGTGGACTCCTGCCAGCAACTCACGGTCTAGGCTCCA` (the centre of Zocchi's 60-mer, as on sfGFP and ubiquitin — a reuse, not a design), built as B-form, bent to a 90 Å anchor span (chosen: both tethers at 88 % of reach, 11.7 Å lift, as the sfGFP build), complement discarded |
| Covalent loop | asserted at build: cut both Cα–Cβ bonds, the two Cβ remain connected through the DNA — closed |
| Force field | ff14SB + DNA.OL15 + TIP3P + Joung–Cheatham NaCl (150 mM) + Li/Merz 12-6 Mg²⁺ (5 mM), 12 Å truncated octahedron |
| Engine | OpenMM 8.6 (`sfGFP/runs/openmm_engine.py`): minimise; 100 ps NVT + 900 ps NPT with heavy-atom restraints ramped off; production NPT 300 K, 1 bar, LangevinMiddle 4 fs with hydrogen-mass repartitioning on the solute, PME 0.9 nm, frames every 20 ps, checkpoint every 200 ps |
| Compute | 1 ns benchmarks on the free `gpu_test` MIG slices; production on `gpu_requeue` (half the billing of `gpu`, preemption costs ≤ 200 ps); windows on `gpu_test`; orchestrated from a 1-core `shared` job |

| system | what | atoms | box, ions |
|---|---|---|---|
| `P0_apo` | the AlphaFold Aα monomer alone (the cryo-EM reproduction's run; ff19SB/OPC, 100 ns) | 163,423 | 129 × 144 × 78 Å box |
| `P1_tet` | Tet2-Et at Asp7 and Asn553, unclicked | 121,808 | 120 Å octahedron, Na⁺ 125 / Cl⁻ 108 / Mg²⁺ 3 |
| `P2_clicked` | both tethers clicked to aminohexanol, no DNA — **the zero** | 125,624 | 122 Å octahedron, Na⁺ 128 / Cl⁻ 113 / Mg²⁺ 4 |
| `P3_ss50` | 50 nt single strand clicked at both ends, covalent loop (the coil) | 259,465 | 154 Å octahedron, Na⁺ 299 / Cl⁻ 239 / Mg²⁺ 7 |
| `P4_ds50` | the same strand with its complement annealed: 50 bp duplex, covalent loop (the spring) | 259,015 | 154 Å octahedron, Na⁺ 346 / Cl⁻ 237 / Mg²⁺ 7 |
| `P5_ss80` | 80 nt single strand, covalent loop | 445,318 | 184 Å octahedron, Na⁺ 496 / Cl⁻ 418 / Mg²⁺ 13 |
| `P6_ds80` | 80 bp duplex, covalent loop | 456,923 | 185 Å octahedron, Na⁺ 583 / Cl⁻ 426 / Mg²⁺ 13 |
| `D_ss50` | bare 50 nt strand, for the restraint windows | 277,063 | 157 Å octahedron, Na⁺ 298 / Cl⁻ 265 / Mg²⁺ 8 |
| `D_ds50` | bare 50 bp duplex, for the restraint windows | 303,807 | 162 Å octahedron, Na⁺ 370 / Cl⁻ 290 / Mg²⁺ 9 |
| `D_ss80` | bare 80 nt strand | 371,816 | 173 Å octahedron, Na⁺ 413 / Cl⁻ 356 / Mg²⁺ 11 |
| `D_ds80` | bare 80 bp duplex | 381,304 | 174 Å octahedron, Na⁺ 499 / Cl⁻ 363 / Mg²⁺ 11 |

Every number on this page is a field of [`analysis/findings.json`](https://github.com/extremefattypunch/DNA-springs-project-2026-27/tree/main/MD%20simulations/pp2a/analysis/findings.json);
every figure ships its CSV in [`figures/`](https://github.com/extremefattypunch/DNA-springs-project-2026-27/tree/main/MD%20simulations/pp2a/figures).

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
