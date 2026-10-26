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
| apo Aα (ff19SB/OPC) | 56.5 ± 7.2 | 40–79 | — | — |

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
| apo Aα (ff19SB/OPC) | 2.73 (max 6.1) | 11.6 | 2.11 | 6.9 / 3.9 | 0.780 |



![RMSD, Day et al.'s split](https://extremefattypunch.github.io/DNA-springs-project-2026-27/figures/pp2a_fig1_rmsd.png)

![per-residue RMSF with the HEAT repeats](https://extremefattypunch.github.io/DNA-springs-project-2026-27/figures/pp2a_fig2_rmsf.png)

## The distance between the attachment points

Three distances are reported because three groups measure three different things:
Cβ7–Cβ553 (ours, the spring's own coordinate), Cα Asn30–Phe578 (Yilmaz et al. 2025's
"N29–F577", 66.5 ± 9.7 Å in their 3 × 704 ns apo runs, 47.7 Å in the compact 6NTS trimer,
70.3 Å in the extended 1B3U crystal) and the HEAT1–HEAT15 centroid distance (Day et al.,
62.7 Å in the B56γ holoenzyme). The AlphaFold start model is compact on all three
(Asn30–Phe578 39 Å).

![distance between the attachment points](https://extremefattypunch.github.io/DNA-springs-project-2026-27/figures/pp2a_fig3_end_to_end.png)

## What the spring delivers, by the model

*Not yet simulated.*

Marko–Siggia applied frame by frame to the strand's end-to-end distance gives
**— pN** as the repository uses it (zero below the chain's rms size) and
— pN if the interpolation formula is trusted all the way down;
the ideal-chain entropic estimate 3k<sub>B</sub>T·x/(2ξL) gives — pN. All three are
**models**, not measurements. For comparison, Zocchi's duplex model would put a 50 bp
spring at this span on its kinked branch at 1.7–3.3 pN of *push* — and outside the
model's own validity range (γ > 1 for L > ~45 bp). The sign is the point: nothing here
pushes.

![modelled against measured force for the single strand](https://extremefattypunch.github.io/DNA-springs-project-2026-27/figures/pp2a_fig4_force_extension.png)

## What it delivers, measured

*The restraint windows have not run yet.*

## Does the horseshoe kink, and does the strand?

The solenoid's axis is the chain of HEAT-repeat Cα centroids; the bend at a repeat is the
angle between the incoming and outgoing axis segments, the twist the dihedral across
four consecutive centroids. A hinge is a repeat whose bend *moves* relative to the clicked
control, not one with a large static bend — a solenoid is bent everywhere. The literature
puts the hinge at HEAT 10–11 (crystal comparison), 11–12 (Day) or 12–13 (Yilmaz), with a
diffuse expansion of repeats 4–12.

| system | total axis bend (°) | HEAT1–15 (Å) | Asn30–Phe578 (Å) | region RMSD N / central / C (Å) | most variable repeat | Arg113–Glu191 salt bridge formed |
|---|---|---|---|---|---|---|
| apo Aα (ff19SB/OPC) | 150 ± 9 | 53.5 ± 6.8 | 42.0 ± 5.6 | 1.19 / 1.48 / 0.73 | H12 | 0.82 |



![where the solenoid bends and twists](https://extremefattypunch.github.io/DNA-springs-project-2026-27/figures/pp2a_fig5_heat_bend.png)

The strand: *not yet simulated.*

![the strand: size, modelled tension, where it bends, stacking](https://extremefattypunch.github.io/DNA-springs-project-2026-27/figures/pp2a_fig6_ssdna.png)

## What the modification itself does

On sfGFP the one site that sat at **positive φ** (left-handed α, Asp133) was pushed 17.5°
out of its basin by the modification; ordinary sites moved by a degree or less. In the
AlphaFold model of Aα, **Asp7 sits at φ, ψ = +65, +57** — the same left-handed-α position —
and Asn553 at -62, -27, in the first turn of HEAT 14's helix. Asp7 is also on the
low-confidence tail (pLDDT 52), so its backbone is free to leave that basin
whether modified or not.

| system | residue 7 | φ, ψ (± sd) | aromatic rings on the Tet ring / frame | residue 553 | φ, ψ (± sd) | rings / frame |
|---|---|---|---|---|---|---|
| apo Aα (ff19SB/OPC) | ASP | -89, +5 (± 43, 76) | — | ASN | -57, -37 (± 9, 9) | — |

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
