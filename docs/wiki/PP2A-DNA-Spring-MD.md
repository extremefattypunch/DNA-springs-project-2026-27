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
that pulls its ends **together**, and only once stretched beyond its natural size: for 50
nucleotides that size is 94 Å (contour 29.5 nm,
persistence 1.5 nm). The two Cβ atoms are 54.4 Å apart in the AlphaFold model and the
tethers reach a further 22 and 26 Å, so the strand *can* be stretched to ~102 Å by the
protein — but only if the horseshoe opens.

**100 ns over 1 systems**, one trajectory each.

![the apo Aα monomer, while the chimera runs](https://extremefattypunch.github.io/DNA-springs-project-2026-27/figures/pp2a_P0_apo_overview_still.png)
*the apo Aα monomer, while the chimera runs.*

## The headline

*The chimera run has not finished; this paragraph is written by `report/build_wiki_pp2a.py` from `analysis/findings.json` once it has.*

![the attachment-point separation against the clicked zero](https://extremefattypunch.github.io/DNA-springs-project-2026-27/figures/pp2a_fig7_headline.png)

| system | Cβ7–Cβ553 (Å), mean ± sd | range | Δ vs P2_clicked (Å) | p |
|---|---|---|---|---|
| apo Aα (ff19SB/OPC) | 56.4 ± 7.1 | 40–75 | — | — |

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



![RMSD, Day et al.'s split](https://extremefattypunch.github.io/DNA-springs-project-2026-27/figures/pp2a_fig1_rmsd.png)

![per-residue RMSF with the HEAT repeats](https://extremefattypunch.github.io/DNA-springs-project-2026-27/figures/pp2a_fig2_rmsf.png)

**The horseshoe as a whole.** The solenoid's axis is the chain of HEAT-repeat Cα centroids; the
bend at a repeat is the angle between the incoming and outgoing axis segments, the twist the
dihedral across four consecutive centroids. The literature puts PR65's hinge at HEAT 10–11
(crystal comparison), 11–12 (Day) or 12–13 (Yilmaz), with a diffuse expansion of repeats 4–12.

| system | total axis bend (°) | HEAT1–15 (Å) | Asn30–Phe578 (Å) | region RMSD N / central / C (Å) | most variable repeat | Arg113–Glu191 salt bridge formed |
|---|---|---|---|---|---|---|
| apo Aα (ff19SB/OPC) | 150 ± 9 | 53.5 ± 6.7 | 42.0 ± 5.6 | 1.19 / 1.48 / 0.73 | H12 | 0.82 |



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

Zocchi's two-spring picture (book pp. 117–119) gives the protein a stiffness κ ≈ 100 k<sub>B</sub>T nm⁻²
along the line between the attachment points, yielding at ~4 pN; sfGFP's force-clamp ladder measured
230–350. No clamp ladder was run here, but equipartition on the Cβ7–Cβ553 fluctuation of each run,
κ = k<sub>B</sub>T / var(x), gives the horseshoe's stiffness along the same coordinate (a lower bound —
slow opening and closing of the horseshoe inflates the variance):

| system | Cβ7–Cβ553 (Å) | κ (k<sub>B</sub>T nm⁻²) | κ (pN nm⁻¹) | displacement per pN of spring force (Å) |
|---|---|---|---|---|
| apo Aα (ff19SB/OPC) | 56.4 ± 7.0 | 2.0 | 8.4 | 1.2 |

This is the number that decides whether a spring of a few piconewtons is visible at all. On sfGFP
(κ ≈ 300) a 4 pN spring moved the sites 0.03 Å, below what 50 ns can resolve; a HEAT horseshoe is
two orders of magnitude softer along this coordinate, so even the ~2 pN a single strand could
deliver would move the tips by ångströms — if it were under tension.

## What a spring would deliver, by Zocchi's model

Zocchi's kinked-rod model (eqs 2.188–2.206 and 3.51–3.52 of the book; Tseng & Zocchi 2013, eq. 6;
B = 200 pN nm², τ<sub>c</sub> = 27 pN nm nicked / 34 continuous, 0.33 nm per bp, k<sub>B</sub>T = 4.2 pN nm),
evaluated with the repository's implementation at the build's anchor span of 90 Å (the chimera has not run yet). A duplex pushes its ends apart on the
bent branch, f = 10B/L² − k<sub>B</sub>T/(L − x), until the span drops below x<sub>c</sub>, where it
kinks and the force collapses to τ<sub>c</sub> / (2R √(1 − (x/2R)²)):

| duplex | contour (nm) | regime at this span | f(x) (pN) | unkinked rod (pN) | Euler load 10B/L² (pN) | stored energy (k<sub>B</sub>T) | x<sub>c</sub> (Å) | γ = Lτ<sub>c</sub>/2B |
|---|---|---|---|---|---|---|---|---|
| 27 bp, nicked | 8.9 | *infeasible: contour shorter than the span* | — | — | 25.2 | — | 87 | 0.60 |
| 27 bp, continuous | 8.9 | *infeasible: contour shorter than the span* | — | — | 25.2 | — | 86 | 0.76 |
| 40 bp, nicked τc = 27 pN nm | 13.2 | kinked | **2.84** | 10.5 | 11.5 | 5.2 | 125 | 0.89 |
| 40 bp, continuous τc = 34 pN nm | 13.2 | kinked | **3.62** | 10.5 | 11.5 | 6.5 | 121 | 1.12 * |
| 50 bp, nicked τc = 27 pN nm | 16.5 | kinked | **1.99** | 6.8 | 7.3 | 6.3 | 151 | 1.11 * |
| 50 bp, continuous τc = 34 pN nm | 16.5 | kinked | **2.54** | 6.8 | 7.3 | 7.9 | 143 | 1.40 * |

\* γ ≥ 1: outside the model's stated validity range (L ≳ 45 bp).

The chimera carries none of these. Its 50 nt **single strand** is Zocchi's "ss state", of which his
2021 chapter says *"with the DNA spring in the ss state there is essentially no stress on the
enzyme"*: a worm-like chain with ξ = 1.5 nm and a contour of 29.5 nm whose natural rms size is
94 Å. At this span it is **matched**, pulling with
0.00 pN as the model is used (1.60 pN if the Marko–Siggia
interpolation is trusted below the rms size). Every number in this section is a **model**, not a measurement.

![Zocchi's spring model at this span](https://extremefattypunch.github.io/DNA-springs-project-2026-27/figures/pp2a_fig9_zocchi_model.png)

*Not yet simulated.*





## What it delivers, measured

*The restraint windows have not run yet.*

![modelled against measured force for the single strand](https://extremefattypunch.github.io/DNA-springs-project-2026-27/figures/pp2a_fig4_force_extension.png)

## Does the DNA kink?

**By the model.** A rod of bending modulus B yields at the critical torque τ<sub>c</sub>, i.e. at a
curvature τ<sub>c</sub>/B = **2.55° per bp** (nicked) to
**3.21° per bp** (continuous) — Zocchi's kink. Below the softening
span x<sub>c</sub> the duplex sits in the kinked state with a hinge of radius R:

| duplex | x<sub>c</sub> (Å) | R of the kinked state (Å) | uniform curvature at x<sub>c</sub> (°/bp) | uniform curvature at this span (°/bp) |
|---|---|---|---|---|
| 27 bp | 87 | 44 | 1.62 | stretched |
| 40 bp | 125 | 65 | 1.63 | 4.17 |
| 50 bp | 151 | 81 | 1.64 | 4.10 |

The force field disagrees with the model on *where* a duplex kinks: on sfGFP's bare duplexes base pairs
began to open at 5.5°/bp and opened in every run only above 8.7°/bp
(the chimera there sat at 4.8-5.2°/bp and never kinked; the one directly comparable all-atom
calculation nucleates a type II kink at 8.6-9.3°/bp). A single strand has no base pairs to
open, so for it "kink" means a sharp turn of the backbone: a C4′-trace angle above 90°, a BI → BII
switch of the ε/ζ torsions, or a step whose bases have come apart.



*The strand's kink figure is rendered once the P3_ss50 analysis exists.*

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

![backbone conformation at the two sites](https://extremefattypunch.github.io/DNA-springs-project-2026-27/figures/pp2a_ramachandran.png)

## Animations

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
| `P3_ss50` | the chimera: one 50 nt strand clicked at both ends, covalent loop | 259,465 | 154 Å octahedron, Na⁺ 299 / Cl⁻ 239 / Mg²⁺ 7 |
| `D_ss50` | the strand alone, for the restraint windows | 277,063 | 157 Å octahedron, Na⁺ 298 / Cl⁻ 265 / Mg²⁺ 8 |

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
