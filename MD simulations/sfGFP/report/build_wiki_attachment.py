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

**No attachment geometry delivers a 10–20 nm wavelength shift at the force a bent duplex
provides.** The calculation is short and can be done two independent ways, and both fall
far short:

| route | input | predicted shift at 2–7 pN |
|---|---|---|
| barrel compliance × spectral sensitivity | GFP's spring constants are 1–17 N/m; Citrine's emission moves 17 nm per 0.8 Å of chromophore deformation | **0.03 – 1.5 nm** |
| pressure equivalence | 360 MPa across a 24 Å barrel is ≈1,600 pN, i.e. 230× our spring | **≈0.07 nm** |

Route A takes the *softest* of five measured axes and grants it the full spectral
sensitivity, so it is a ceiling. Route B assumes the deformation scales linearly with
load, which for a barrel pushed at two points rather than squeezed isotropically is a
floor. The honest reading is that the answer lies between them: **somewhere under
1.5 nm, and short of the target by between 7× and 137×.**

This is consistent with the state of the field. As of 2024 the leading laboratory on
fluorescent-protein mechanics states plainly that *"forced-induced changes in the
excitation or emission wavelengths of FPs have not been described"*
([Cell Rep Methods 2024](https://doi.org/10.1016/j.crmeth.2024.100815)). Nobody has done
this, and the reason appears to be that a few piconewtons is three orders of magnitude
below the ~1,600 pN equivalent that moves Citrine by 17 nm.

It is also consistent with **our own simulations**: across all nine sfGFP constructs, no
chromophore hydrogen bond differs from the clicked control (details
[below](#what-our-own-simulations-say)).

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

## What to do instead

Ranked by expected value, not by ambition.

**1. Measure what a spring can actually change: intensity and lifetime, not colour.**
Quantum yield and excited-state lifetime respond to the hydrogen-bond network far more
readily than the peak position does, and our His148 occupancy already ranges 0.54–0.85
across constructs. This is the observable the physics supports.

**2. If colour is the goal, shift an equilibrium — do not deform a chromophore.**
T203Y is the wrong amplifier: its 16 nm needs 0.8 Å of deformation, which is what the
estimate says we cannot supply. The **deGFP** variants are the right idea instead —
`S65T/H148G/T203C` and `C48S/S65T/H148C/T203C` switch emission between **~515 nm green and
~460 nm blue**, a 55 nm swing, because they sit at a two-state equilibrium with a pKa near
7 (Hanson 2002; McAnaney 2002). Tipping a poised equilibrium costs far less energy than
compressing a barrel. This is speculative — nobody has driven a deGFP mechanically — but
it is the only route where the amplification mechanism is large enough to matter.

**3. Settle the force before choosing a site.** The estimate above uses Zocchi's
2–7 pN. Our own direct measurement on bare duplexes gives tens of piconewtons on the
*unyielded* branch, and the value rises with box size and has not converged
([see the report]({REPO}/blob/main/docs/index.html)). If the real construct delivers 30 pN
rather than 5, every number on this page improves by that factor. Measuring the delivered
force in the chimera is cheaper than building four constructs.

**4. If attaching anyway, these are the buildable positions.** Both **N149** and **Q204**
are already proven unnatural-amino-acid and click sites in sfGFP — a dye has been clicked
onto 204 (Reddington 2012), a strained amino acid encoded there (Xiong 2022), and Pope
2021 built sfGFP204×2. Of pairs that a spring could actually span on the chromophore
strands, **Tyr151–Ser208 at {next((b['cb_cb_A'] for b in ds.get('buildable', []) if b['pair'].startswith('TYR151')), 0):.1f} Å**
is closest to our built 31.3 Å span, so the existing 27 bp spring geometry nearly fits.
One hard constraint: every circular permutant with new termini *inside* a β-strand was
non-chromophoric, so attach at exposed side chains, never by breaking the chain.

---

## References

Verified against the primary source. Where a claim above rests on a secondary source or
on arithmetic of ours, it says so at that point.

**The spring, and what force it delivers**
1. Dietz H, Rief M (2006) Protein structure by mechanical triangulation. *PNAS* — GFP's directional spring constants, 1–17 N/m; fracture 356–548 pN. [10.1073/pnas.0602995103](https://doi.org/10.1073/pnas.0602995103)
2. Dietz H, Rief M (2004) Exploring the energy landscape of GFP by single-molecule mechanical experiments. *PNAS* 101:16192 — first yield at ~35 pN, N-terminal helix detachment. [10.1073/pnas.0404549101](https://doi.org/10.1073/pnas.0404549101)

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
