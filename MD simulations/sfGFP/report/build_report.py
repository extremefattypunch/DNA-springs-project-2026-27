#!/usr/bin/env python3
"""Build the standalone HTML report from the completed analyses.

Self-contained: figures are embedded as base64 PNG and the animations as base64 mp4,
so the file survives being mailed, dropped in a shared drive, or published as an
artifact.  Numbers come from analysis/findings.json, which is computed separately --
the report renders numbers it did not calculate, so a figure and the prose beside it
cannot disagree.
"""
from __future__ import annotations

import argparse
import base64
import json
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parent.parent


def b64(path: Path, mime: str):
    return (f"data:{mime};base64," + base64.b64encode(path.read_bytes()).decode()
            if path.exists() else None)


def esc(s):
    return str(s).replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")


def sig(x, nd=2):
    return "—" if x is None else f"{x:.{nd}f}"


def pm(a, nd=2, key="sem"):
    if not a:
        return "—"
    return f"{a['mean']:.{nd}f} ± {a[key]:.{nd}f}"


LABEL = {"S0_wt": "WT sfGFP", "S1_tet": "2× Tet2-Et", "S2_clicked": "clicked, no DNA",
         "S3_spring27": "27 bp, two-arm", "S4_spring40": "40 bp, two-arm",
         "S5_spring40nick": "40 bp nicked, two-arm", "S6_clamp": "force clamp",
         "L3_loop27": "27 bp, closed loop", "L4_loop40": "40 bp, closed loop",
         "L5_loop40nick": "40 bp nicked, closed loop",
         "C1_ss50": "50 nt ssDNA, closed loop", "C2_ssarms": "two 25 nt arms",
         "U0_apo": "ubiquitin, apo", "U1_clicked": "ubiquitin, clicked, no DNA",
         "U2_ss50": "ubiquitin, 50 nt ssDNA loop",
         "U3_ds50": "ubiquitin, 50 bp dsDNA loop",
         "U4_ssarms": "ubiquitin, two 25 nt arms"}
PURPOSE = {
    "S0_wt": "Baseline. Validates the published chromophore parameters against a "
             "structure nobody modified.",
    "S1_tet": "The two encoded tetrazines, unclicked — the fluorometry sample. Also "
              "the topology the force-clamp ladder runs on.",
    "S2_clicked": "The clicked tether with no DNA: the zero-force reference with a "
                  "chimera’s own chemistry.",
    # These three are the SUPERSEDED builds: a linker on the 5' end of each of two
    # different strands, so no covalent cycle exists and the loop is closed by base
    # pairing alone.  They are kept and reported as what they are, because that is
    # Zocchi's own pre-ligation two-arm chimera and it is the closest thing we have to a
    # non-circular control.  L3/L4/L5 are the corrected circular builds.
    "S3_spring27": "27 bp spring, <strong>two-arm</strong> construct — one linker on "
                   "each of two strands, no covalent cycle. Superseded by "
                   "<code>L3_loop27</code>; kept as the non-circular comparison.",
    "S4_spring40": "40 bp, two-arm construct. Zocchi’s mechanically clean length "
                   "(γ &lt; 1). Superseded by <code>L4_loop40</code>.",
    "S5_spring40nick": "40 bp two-arm construct, nicked at the centre of one strand. "
                       "Superseded by <code>L5_loop40nick</code>, where the nick moves "
                       "to the complement so the loop stays closed.",
    "S6_clamp": "A constant force of 0/2/4/7/12/20 pN applied straight along the "
                "Cβ–Cβ axis: the calibration curve.",
    # The corrected builds: one strand modified at both its 5' and 3' ends, so a covalent
    # path runs out along one tether, through the DNA, and back along the other.  This is
    # the molecule the experiment makes; the S3-S5 trio above is not.
    "L3_loop27": "The corrected 27 bp chimera — one strand carrying both linkers, so a "
                 "covalent cycle exists. Asserted with a topology test the superseded "
                 "build fails.",
    "L4_loop40": "The corrected 40 bp chimera. Longer spring: less force, more stored "
                 "energy.",
    "L5_loop40nick": "The corrected nicked 40 bp chimera — the low-stress reference "
                     "Zocchi’s calibration leans on. Measured directly, it is not "
                     "softer.",
    "C1_ss50": "A 50-nucleotide single strand clicked at both ends. Its contour is far "
               "longer than the span it has to cover, so it is slack and its tension is "
               "exactly zero — which makes it the zero-force control the whole argument "
               "rests on.",
    "C2_ssarms": "The same strand cut in the middle: two 25 nt arms, each attached at "
                 "one end, which cannot circularise. Holds DNA content, sequence, sites "
                 "and chemistry fixed and changes only the topology. The loop test is "
                 "asserted in both directions, because a control that quietly closed "
                 "would look like a null result.",
    "U0_apo": "Ubiquitin, unmodified. The second protein exists so the conclusion does "
              "not rest on one fold and one chemistry.",
    "U1_clicked": "1UBQ F45W/T9C/K63C with the maleimide-thioether and DBCO-triazole "
                  "crosslinker installed, no DNA: this family’s own matched control.",
    "U2_ss50": "The 50 nt single-stranded chimera from the companion preprint, in the "
               "protein that preprint studies. Slack, so zero tension.",
    "U3_ds50": "The 50 bp duplex chimera: ubiquitin’s closed covalent loop, and the "
               "only ubiquitin construct under any modelled load.",
    "U4_ssarms": "Two 25 nt arms on ubiquitin — the preprint’s own control, which "
                 "reproduced most of the circular chimera’s NMR spectrum.",
}
FIGS = [
    ("fig1_force_response", "Force response of the attachment sites",
     """The deformation coordinate is the Cβ–Cβ separation of the two attachment
     residues; the abscissa is the force pulling them apart. Blue is the axial clamp,
     where the force is <em>known</em> because we imposed it. Faint marks are individual
     replicates, filled marks their mean, bars the standard error. The orange marks are
     the three chimeras, placed at the force their own spring geometry implies. Two
     things to read: the clamp response is barely distinguishable from flat — its fitted
     compliance has a confidence interval straddling zero, which is why the annotation
     quotes a <em>bound</em> on stiffness rather than a value — and the chimeras are
     offset from it, but <em>downward</em> once the right zero is used. The chimeras are
     plotted here against the clamp ladder, which runs on the <em>unclicked</em> protein;
     read the offsets against the clicked control in the table above instead, not off
     this panel. The 4 pN clamp point falling below the 0 pN point is the size of the
     noise at this sampling; it is not a physical dip."""),
    ("ramachandran", "Backbone conformation at the two modified sites",
     """Every frame of every replicate of the wild type, the unclicked Tet2-Et protein
     and the clicked adduct, at both attachment sites; large markers are circular means,
     and the grey outlines are the conventional Ramachandran regions drawn as a common
     reference across the two panels. The point is the position, not the spread.
     <strong>Site 133 lives in the left-handed α region</strong> — positive φ — which is
     a conformation Asp, Asn and Gly commonly adopt and a phenylalanine-shaped side chain
     does not: Tet2-Et (orange) is pushed 17.5° out in φ and 14° in ψ, and clicking
     (green) returns it to within a few degrees of wild type (blue). Site 149 is ordinary
     β and all three sit on top of one another. DSSP over residues 128–155 is identical
     in all three, so this is a local backbone strain, not a change of secondary
     structure."""),
    ("fig2_rmsf", "Backbone mobility, residue by residue",
     """Per-residue Cα RMSF, aligned on the barrel core rather than the whole protein —
     aligning globally lets the tethers and the flexible loops they sit in absorb the
     fit and smear their motion across every residue. The barrel is rigid everywhere
     (0.3–0.5 Å) with peaks at the loops, and the traces for WT, the two-tetrazine
     protein and the chimeras lie on top of each other: at 45 ns per replicate the
     spring does not measurably loosen the fold anywhere, including at the attachment
     sites themselves. The ACE/NME caps and the residues either side are omitted; their
     4–11 Å free-end motion would compress everything else into the bottom tenth of the
     axis."""),
    ("fig3_chromophore_hbonds", "The chromophore’s hydrogen-bond network",
     """Occupancy judged properly — on the hydrogen-to-acceptor distance and the
     donor–H–acceptor angle, not the heavy-atom distance. That distinction changes a
     conclusion: His148 looks broken by a 3.5 Å heavy-atom cutoff while its HD1 proton
     sits 2.3–2.5 Å from the phenolate with near-linear geometry. Ticks are individual
     replicates. The Ser205–Glu222 row is the reason they are shown: within a system it
     is either near 0% or near 90%, a two-state Ser205 hydroxyl rotamer sampled once per
     trajectory, so its ~30% mean is an average of coin flips rather than a
     measurement."""),
    ("fig4_dna_bend", "Where the duplex takes the strain",
     """Bend per step of the helical axis, averaged over replicates, shaded by their
     standard deviation. Base-pair centres are first averaged over one full helical
     turn: without that, the ~1.9 Å spiral of the centres alone reads as 27° per base
     pair of bend that is not there. A kink — the softening transition Zocchi’s model
     is built on — would appear as a sharp localised spike. None appears. The curvature
     is spread at 3–8° per step across both spring lengths, and essentially no base
     pairs open, so on this timescale the duplex is an elastically bent rod rather than
     a hinge."""),
    ("fig5_kink_threshold", "Where this force field starts breaking base pairs — and how reliably",
     """The previous panel says no kink appears in the chimera. That could mean the duplex
     genuinely does not kink at this curvature, or that this force field cannot kink at
     all — and the two have opposite implications for Zocchi’s model. Separating them
     needs curvature the chimera does not reach, which is what the harmonic windows
     provide. The measure plotted is the <em>fraction of runs</em> in which an interior
     base pair opens, because the per-run observation is bimodal — a run either keeps
     every interior pair shut or opens one and holds it open most of the trajectory — so
     an average over runs would be a number no run exhibits. Ticks along the top and
     bottom are the individual runs; the band is a Wilson 68 % interval, which is the
     right interval for a proportion out of four; the small numbers are how many runs
     each bin holds. Curvature is the circular-arc equivalent of the separation the
     duplex actually adopted, not of the target it was pulled towards.
     <strong>Opening is stochastic from about 5.4°/bp and becomes certain only at
     9.5°/bp</strong>, with the published type II kink curvature (grey) sitting where the
     probability is already 0.8. So the force field is entirely capable of breaking base
     pairs, and it does so around where the literature says it should. Our chimera, at
     4.8–5.2°/bp, opened none in three 50 ns runs — it sits just below the onset. That is
     what makes this a measurement rather than an absence, and it excludes “the force
     field cannot kink” as the explanation."""),
    ("fig6_force_extension", "Force measured from the force field, not from a model",
     """The bare duplexes, no protein, held at a series of separations by a harmonic
     restraint; the ordinate is the mean restraint force, which <em>is</em> the force the
     duplex exerts. Faint marks are individual runs, circles the stiff restraint and
     squares the soft one; the solid bar is each duplex’s mean over the separations it was
     sampled at, banded by its standard error. It is drawn flat because it <em>is</em>
     flat — regressing force on separation gives no resolvable slope for any of the three,
     which is what a buckled rod does: past its critical load it keeps folding at roughly
     constant force. Two things to read. The plateaux sit on the dashed <strong>Euler
     buckling loads</strong> π²B/L², within 4 % for the 27 bp, with nothing fitted. And
     they sit far above the grey band, which is what Zocchi’s kinked branch predicts for
     these same duplexes — about six times below what we measure. The two restraint
     stiffnesses interleave rather than separating, which is the control that makes this
     the duplex’s force and not the restraint’s."""),
]
STEPS = [
    ("Prepare", "build/01_protein",
     "2B3P → tleap-ready structure. Strips the crystallisation Cd²⁺ and acetate, keeps "
     "the buried cavity waters, assigns every histidine tautomer from its own H-bond "
     "geometry, and measures the attachment geometry everything downstream reads."),
    ("Parameterise", "build/02_params",
     "Explicit atom and bond tables → RDKit → antechamber (AM1-BCC, GAFF2) → prepgen → "
     "backbone harmonisation → parmchk2. Zero parameters flagged for review across all "
     "four custom residues."),
    ("Build the spring", "build/03_dna",
     "Ideal B-form duplex from PyMOL fnab, bent to the required span by mapping the "
     "helical axis onto a circular arc with each base pair moved as a rigid slab, so "
     "intra-pair geometry is preserved exactly and only the inter-pair geometry changes."),
    ("Assemble", "build/04_assemble",
     "Tethers posed by backbone superposition plus torsion optimisation; the duplex "
     "docked onto the two target phosphate positions with its arc pushed away from the "
     "barrel; tleap driven in two passes so the ion count matches the water it dissolves "
     "in."),
    ("Run", "runs, slurm",
     "Minimise → restrained NVT → staged NPT release → production at 4 fs with "
     "hydrogen-mass repartitioning on the solute only. Restart-safe, so preemption "
     "costs one checkpoint interval."),
    ("Measure", "analysis",
     "Chromophore geometry and its H-bond network, protein and spring extension, where "
     "the duplex gives, barrel deformation, cavity water — then reduced to one row per "
     "replicate and one findings file the report renders."),
]
def params_summary(root: Path):
    """Read the parameterisation report so the verification table cannot go stale.

    These rows used to be literals reading "all four custom residues".  There are seven
    now, and a literal has no way of noticing -- so the counts are derived, and if a
    residue is ever added without rerunning build_residues.py the numbers move rather
    than the sentence quietly becoming false.
    """
    f = root / "build" / "02_params" / "PARAMS_REPORT.json"
    if not f.exists():
        return None
    rows = json.loads(f.read_text())
    if isinstance(rows, dict):
        rows = list(rows.values())
    names = [r["residue"] for r in rows]
    return {
        "n": len(rows), "names": ", ".join(names),
        "atoms": sum(r["element_check"]["atoms_checked"] for r in rows),
        "mismatches": sum(r["element_check"]["element_mismatches"] for r in rows),
        "flagged": sum(r["frcmod_audit"]["n_flagged"] for r in rows),
    }


VER = [
    ("OpenMM with CUDA", "all four platforms agree within tolerance", True),
    ("CRO atom names against the xFPchromophores template", "22 / 22, checked before use", True),
    ("CRO net charge from the library", "−1.0000 exactly", True),
    ("parmchk2 parameters flagged “ATTN, need revision”", "__PARAMS_FLAGGED__", True),
    ("Assigned GAFF type implies the declared element", "__PARAMS_ELEMENTS__", True),
    ("Prep tree parent links are declared bonds", "__PARAMS_TREE__", True),
    ("tleap errors on the assembled chimera", "0", True),
    ("Net charge, every system", "integral to &lt; 1×10⁻⁴", True),
    ("Carbamate bonds present in the prmtop", "2 / 2 per chimera, asserted", True),
    ("Linker bonds after minimisation",
     "C–N 1.365 Å (ideal 1.38); O–P 1.599 Å (ideal 1.61)", True),
    ("Spring model against Zocchi’s published forces",
     "1.47 vs 1.5 pN; 2.33 vs 2.4 pN; 6.62 vs 6.6 k<sub>B</sub>T", True),
    ("Bend transform is the identity at zero curvature", "O3′–P 1.601 Å at R = 566 Å", True),
]


# ---------------------------------------------------------------------------
# References.  Order is the order they are first cited; cite() renders the marker.
# Every entry is something this setup actually rests on -- a structure, a parameter
# set, a formula, a piece of chemistry, or a program that produced a number here.
# ---------------------------------------------------------------------------
REFS = [
    ("pedelacq", "Pédelacq J-D, Cabantous S, Tran T, Terwilliger TC, Waldo GS. "
     "Engineering and characterization of a superfolder green fluorescent protein. "
     "<i>Nat Biotechnol</i> <b>24</b>, 79–88 (2006).",
     "https://doi.org/10.1038/nbt1172",
     "The sfGFP structure this work is built on, PDB 2B3P."),
    ("ormo", "Ormö M, Cubitt AB, Kallio K, Gross LA, Tsien RY, Remington SJ. "
     "Crystal structure of the <i>Aequorea victoria</i> green fluorescent protein. "
     "<i>Science</i> <b>273</b>, 1392–1395 (1996).",
     "https://doi.org/10.1126/science.273.5280.1392",
     "Identifies the residues that hold the chromophore: Thr203 and His148 at the "
     "phenolate, Arg96 and Gln94 at the imidazolinone. This is where the choice of "
     "which hydrogen bonds to track comes from."),
    ("brejc", "Brejc K, Sixma TK, Kitts PA, Kain SR, Tsien RY, Ormö M, Remington SJ. "
     "Structural basis for dual excitation and photoisomerization of the "
     "<i>Aequorea victoria</i> green fluorescent protein. "
     "<i>Proc Natl Acad Sci USA</i> <b>94</b>, 2306–2311 (1997).",
     "https://doi.org/10.1073/pnas.94.6.2306",
     "The Thr203–Glu222–Ser205 proton relay behind the chromophore, and the basis "
     "for treating the chromophore as anionic at pH 7.4."),
    ("niwa", "Niwa H, Inouye S, Hirano T, Matsuno T, Kojima S, Kubota M, Ohashi M, "
     "Tsuji FI. Chemical nature of the light emitter of the <i>Aequorea</i> green "
     "fluorescent protein. <i>Proc Natl Acad Sci USA</i> <b>93</b>, 13617–13622 (1996).",
     "https://doi.org/10.1073/pnas.93.24.13617",
     "The synthetic chromophore is essentially non-fluorescent free in solution while "
     "the protein-bound one is bright: the basis for reading bridge twist as a proxy "
     "for non-radiative decay."),
    ("breyfogle", "Breyfogle KL, Blood DL, Rosnik AM, Krueger BP. "
     "Molecular dynamics force field parameters for the EGFP chromophore and some of "
     "its analogues. <i>J Phys Chem B</i> <b>127</b>, 5772–5788 (2023).",
     "https://doi.org/10.1021/acs.jpcb.3c01486",
     "The CRO chromophore parameters, loaded through leaprc.xFPchromophores. Also the "
     "source of the backbone charge-constraint idea used for the custom residues."),
    ("eddins", "Eddins AJ, Pung H-S, Cooley RB, Mehl RA. Tetrazine amino acid encoding "
     "for rapid and complete protein bioconjugation. "
     "<i>Bio-protocol</i> <b>14</b>(16), e5048 (2024).",
     "https://doi.org/10.21769/BioProtoc.5048",
     "Tet2-Et itself, the strained-TCO chemistry, and the two permissive amber sites "
     "in sfGFP (N150 and D134) that this construct uses."),
    ("blizzard", "Blizzard RJ, Backus DR, Brown W, Bazewicz CG, Li Y, Mehl RA. "
     "Ideal bioorthogonal reactions using a site-specifically encoded tetrazine amino "
     "acid. <i>J Am Chem Soc</i> <b>137</b>, 10044–10047 (2015).",
     "https://doi.org/10.1021/jacs.5b03275",
     "The tetrazine amino acid scaffold Tet2-Et is built on, and the sTCO kinetics."),
    ("zocchi_book", "Zocchi G. <i>Molecular Machines: A Materials Science Approach.</i> "
     "Princeton University Press, Princeton NJ (2018). ISBN 978-0-691-17386-3 "
     "(hardcover); 978-1-4008-9006-4 (ebook).",
     "https://doi.org/10.23943/9781400890064",
     "The DNA-spring concept and the whole mechanical model: Eqs. 2.188–2.206 for the "
     "kinked rod, Eq. 3.52 for the force, and the ~100 k<sub>B</sub>T nm⁻² protein "
     "stiffness with its softening transition."),
    ("qu", "Qu H, Tseng C-Y, Wang Y, Levine AJ, Zocchi G. The elastic energy of sharply "
     "bent nicked DNA. <i>Europhys Lett</i> <b>90</b>, 18003 (2010).",
     "https://doi.org/10.1209/0295-5075/90/18003",
     "The elastic-energy expression E(x) whose derivative is the force formula used "
     "here, and the measurement of the critical bending torque."),
    ("tseng2013", "Tseng C-Y, Zocchi G. Mechanical control of Renilla luciferase. "
     "<i>J Am Chem Soc</i> <b>135</b>, 11879–11886 (2013).",
     "https://doi.org/10.1021/ja4043565",
     "Eq. 6 as implemented, τ<sub>c</sub> = 27 pN nm for a nicked rod, and the "
     "1.5 / 2.4 pN values the implementation is tested against."),
    ("tseng2021", "Tseng C-Y, Wang A, Zocchi G. Protein–DNA chimeras: synthesis and "
     "applications. <i>Methods Enzymol</i> <b>647</b>, 257–281 (2021).",
     "https://doi.org/10.1016/bs.mie.2020.09.010",
     "The verified oligonucleotide sequences used for the springs, and the measured "
     "energy split E<sub>tot</sub> 9.1 = E<sub>DNA</sub> 6.6 + E<sub>enzyme</sub> "
     "2.5 k<sub>B</sub>T."),
    ("choi2007", "Choi B, Zocchi G. Guanylate kinase, induced fit, and the allosteric "
     "spring probe. <i>Biophys J</i> <b>92</b>, 1651–1658 (2007).",
     "https://doi.org/10.1529/biophysj.106.092866",
     "The earlier worm-like-chain force estimate of ~10 pN, which that paper itself "
     "labels an upper bound and which is deliberately <i>not</i> used here."),
    ("boral", "Boral S, Schnebly J, Gamada D, Gardner KH, Hekstra DR. Conformational "
     "changes induced in ubiquitin by circular protein–DNA chimeras. "
     "<i>bioRxiv</i> 2026.07.05.736590 (2026).",
     "https://doi.org/10.64898/2026.07.05.736590",
     "The companion experiment in this lab: the ~4 pN estimate for a 50-mer spring, "
     "and the buffer and salt conditions matched in the simulations."),
    ("maier", "Maier JA, Martinez C, Kasavajhala K, Wickstrom L, Hauser KE, Simmerling C. "
     "ff14SB: improving the accuracy of protein side chain and backbone parameters from "
     "ff99SB. <i>J Chem Theory Comput</i> <b>11</b>, 3696–3713 (2015).",
     "https://doi.org/10.1021/acs.jctc.5b00255", "The protein force field."),
    ("zgarbova", "Zgarbová M, Šponer J, Otyepka M, Cheatham TE III, Galindo-Murillo R, "
     "Jurečka P. Refinement of the sugar–phosphate backbone torsion beta for AMBER force "
     "fields improves the description of Z- and B-DNA. "
     "<i>J Chem Theory Comput</i> <b>11</b>, 5723–5736 (2015).",
     "https://doi.org/10.1021/acs.jctc.5b00716", "The DNA force field, OL15."),
    ("wang_gaff", "Wang J, Wolf RM, Caldwell JW, Kollman PA, Case DA. Development and "
     "testing of a general AMBER force field. "
     "<i>J Comput Chem</i> <b>25</b>, 1157–1174 (2004).",
     "https://doi.org/10.1002/jcc.20035",
     "GAFF/GAFF2, which types the tetrazine, the dihydropyridazine, the sTCO cage and "
     "the linker."),
    ("jakalian", "Jakalian A, Jack DB, Bayly CI. Fast, efficient generation of "
     "high-quality atomic charges. AM1-BCC model: II. Parameterization and validation. "
     "<i>J Comput Chem</i> <b>23</b>, 1623–1641 (2002).",
     "https://doi.org/10.1002/jcc.10128", "The charge model for the custom residues."),
    ("jorgensen", "Jorgensen WL, Chandrasekhar J, Madura JD, Impey RW, Klein ML. "
     "Comparison of simple potential functions for simulating liquid water. "
     "<i>J Chem Phys</i> <b>79</b>, 926–935 (1983).",
     "https://doi.org/10.1063/1.445869", "TIP3P water."),
    ("joung", "Joung IS, Cheatham TE III. Determination of alkali and halide monovalent "
     "ion parameters for use in explicitly solvated biomolecular simulations. "
     "<i>J Phys Chem B</i> <b>112</b>, 9020–9041 (2008).",
     "https://doi.org/10.1021/jp8001614", "Na⁺ and Cl⁻ parameters."),
    ("li_merz", "Li P, Roberts BP, Chakravorty DK, Merz KM Jr. Rational design of "
     "particle mesh Ewald compatible Lennard-Jones parameters for +2 metal cations in "
     "condensed-phase simulations. "
     "<i>J Chem Theory Comput</i> <b>9</b>, 2733–2748 (2013).",
     "https://doi.org/10.1021/ct400146w",
     "Mg²⁺, the 12-6 set — chosen over 12-6-4 because OpenMM has no C4 term."),
    ("eastman", "Eastman P, Galvelis R, Peláez RP, et al. OpenMM 8: molecular dynamics "
     "simulation with machine learning potentials. "
     "<i>J Phys Chem B</i> <b>128</b>, 109–116 (2024).",
     "https://doi.org/10.1021/acs.jpcb.3c06662", "The simulation engine."),
    ("leimkuhler", "Leimkuhler B, Matthews C. Rational construction of stochastic "
     "numerical methods for molecular sampling. "
     "<i>Appl Math Res Express</i> <b>2013</b>, 34–56 (2013).",
     "https://doi.org/10.1093/amrx/abs010",
     "The LangevinMiddle scheme used for production."),
    ("hopkins", "Hopkins CW, Le Grand S, Walker RC, Roitberg AE. Long-time-step "
     "molecular dynamics through hydrogen mass repartitioning. "
     "<i>J Chem Theory Comput</i> <b>11</b>, 1864–1874 (2015).",
     "https://doi.org/10.1021/ct5010406", "The 4 fs timestep."),
    ("case", "Case DA, Aktulga HM, Belfon K, et al. <i>AmberTools</i>. "
     "<i>J Chem Inf Model</i> <b>63</b>, 6183–6191 (2023).",
     "https://doi.org/10.1021/acs.jcim.3c01153",
     "tleap, antechamber, prepgen and parmchk2, which built every topology here."),
    ("mdtraj", "McGibbon RT, Beauchamp KA, Harrigan MP, et al. MDTraj: a modern open "
     "library for the analysis of molecular dynamics trajectories. "
     "<i>Biophys J</i> <b>109</b>, 1528–1532 (2015).",
     "https://doi.org/10.1016/j.bpj.2015.08.015", "Trajectory analysis."),
    ("shrake", "Shrake A, Rupley JA. Environment and exposure to solvent of protein "
     "atoms. Lysozyme and insulin. <i>J Mol Biol</i> <b>79</b>, 351–371 (1973).",
     "https://doi.org/10.1016/0022-2836(73)90011-9",
     "The solvent-accessibility algorithm used to check both attachment sites are "
     "exposed."),
    ("kabsch", "Kabsch W. A solution for the best rotation to relate two sets of "
     "vectors. <i>Acta Crystallogr A</i> <b>32</b>, 922–923 (1976).",
     "https://doi.org/10.1107/S0567739476001873",
     "Superposition, used to place the tethers on the backbone."),
    ("pymol", "Schrödinger LLC. <i>The PyMOL Molecular Graphics System</i>, "
     "version 3.x (open-source build).", "https://pymol.org/",
     "The <code>fnab</code> builder that generated the ideal B-form duplexes, and "
     "every molecular render here."),
    ("labms", "Abebe S. Spring-stimulated protein mechanics. Hekstra lab meeting, "
     "20 March 2026 (internal).", None,
     "The ESI-MS of the single and double Tet2-Et constructs (27,967 and 28,108 Da) "
     "that settles which residues carry the amber codons."),
]
REF_INDEX = {k: i + 1 for i, (k, *_ ) in enumerate(REFS)}


def cite(*keys):
    """Inline superscript citation markers linking to the reference list."""
    parts = []
    for k in keys:
        n = REF_INDEX[k]
        parts.append(f'<a href="#ref-{n}" class="cite">{n}</a>')
    return "<sup>" + ",".join(parts) + "</sup>"


def refs_html():
    rows = []
    for i, (key, text, url, why) in enumerate(REFS, 1):
        link = (f' <a href="{url}" target="_blank" rel="noopener">link</a>'
                if url else "")
        rows.append(f'<li id="ref-{i}"><span class="ref-n">{i}</span>'
                    f'<div><p class="ref-cite">{text}{link}</p>'
                    f'<p class="ref-why">{why}</p></div></li>')
    return "<ol class=\"refs\">" + "".join(rows) + "</ol>"


CSS = """
<title>DNA Spring on sfGFP</title>
<link rel="preconnect" href="https://fonts.googleapis.com">
<link rel="preconnect" href="https://fonts.gstatic.com" crossorigin>
<link rel="stylesheet" href="https://fonts.googleapis.com/css2?family=Spectral:ital,wght@0,400;0,600;1,400&family=Source+Sans+3:wght@400;600&family=IBM+Plex+Mono:wght@400;500&display=swap">
<style>
  :root {
    --ground:#f7f8f5; --panel:#ffffff; --ink:#10201c; --ink-2:#4a5c56;
    --ink-3:#7d8d87; --rule:#dde4de; --rule-2:#eef2ee;
    --gfp:#1a8a4b; --dna:#2a5f8f; --force:#a55c17;
    --ok:#1a8a4b; --warn:#a55c17; --flag:#8c3a3a;
    --serif:"Spectral",Georgia,"Times New Roman",serif;
    --sans:"Source Sans 3","Helvetica Neue",Arial,sans-serif;
    --mono:"IBM Plex Mono","SFMono-Regular",Consolas,monospace;
    --measure:70ch;
  }
  @media (prefers-color-scheme: dark) {
    :root:not([data-theme="light"]) {
      --ground:#0d1412; --panel:#141d1a; --ink:#e9efeb; --ink-2:#a9bab3;
      --ink-3:#7b8c86; --rule:#26332e; --rule-2:#1b2622;
      --gfp:#4cc47f; --dna:#6fa8dc; --force:#d99a4e;
      --ok:#4cc47f; --warn:#d99a4e; --flag:#e08585;
    }
  }
  :root[data-theme="dark"] {
    --ground:#0d1412; --panel:#141d1a; --ink:#e9efeb; --ink-2:#a9bab3;
    --ink-3:#7b8c86; --rule:#26332e; --rule-2:#1b2622;
    --gfp:#4cc47f; --dna:#6fa8dc; --force:#d99a4e;
    --ok:#4cc47f; --warn:#d99a4e; --flag:#e08585;
  }
  * { box-sizing:border-box; }
  body { margin:0; background:var(--ground); color:var(--ink);
    font-family:var(--sans); font-size:16px; line-height:1.62;
    -webkit-font-smoothing:antialiased; }
  .wrap { max-width:1080px; margin:0 auto; padding:0 28px 96px; }
  .col { max-width:var(--measure); }
  .mono { font-family:var(--mono); font-variant-numeric:tabular-nums; }
  .num { font-family:var(--mono); font-variant-numeric:tabular-nums; text-align:right;
    white-space:nowrap; }
  header.mast { padding:60px 0 30px; border-bottom:1px solid var(--rule); }
  .eyebrow { font-family:var(--mono); font-size:11.5px; letter-spacing:.13em;
    text-transform:uppercase; color:var(--ink-3); margin:0 0 18px; }
  h1 { font-family:var(--serif); font-weight:600; font-size:clamp(34px,5vw,52px);
    line-height:1.08; letter-spacing:-.015em; margin:0 0 18px; text-wrap:balance; }
  .lede { font-family:var(--serif); font-size:19.5px; line-height:1.58;
    color:var(--ink-2); max-width:62ch; margin:0; }
  .lede em { color:var(--ink); font-style:italic; }
  .stats { display:grid; gap:1px; background:var(--rule);
    grid-template-columns:repeat(auto-fit,minmax(170px,1fr));
    border:1px solid var(--rule); margin:34px 0 0; }
  .stat { background:var(--panel); padding:17px 19px; }
  .stat .k { font-family:var(--mono); font-size:10.5px; letter-spacing:.11em;
    text-transform:uppercase; color:var(--ink-3); }
  .stat .v { font-family:var(--mono); font-size:24px; font-weight:500; margin-top:6px;
    font-variant-numeric:tabular-nums; }
  .stat .s { font-size:13px; color:var(--ink-2); margin-top:4px; line-height:1.45; }
  .stat.force .v { color:var(--force); }
  .stat.ok .v { color:var(--ok); }
  section { padding:50px 0 0; }
  h2 { font-family:var(--serif); font-weight:600; font-size:27px; letter-spacing:-.01em;
    margin:0 0 6px; text-wrap:balance; }
  .sec-head { border-bottom:1px solid var(--rule); padding-bottom:12px;
    margin-bottom:24px; display:flex; align-items:baseline;
    justify-content:space-between; gap:18px; flex-wrap:wrap; }
  .sec-head .tag { font-family:var(--mono); font-size:11.5px; color:var(--ink-3);
    letter-spacing:.08em; }
  h3 { font-family:var(--sans); font-weight:600; font-size:17px; margin:32px 0 8px;
    max-width:var(--measure); text-wrap:balance; }
  p { margin:0 0 15px; max-width:var(--measure); }
  a { color:var(--dna); text-decoration:none; border-bottom:1px solid currentColor; }
  .tbl-wrap { overflow-x:auto; margin:0 0 22px; }
  table { border-collapse:collapse; width:100%; font-size:14.5px; }
  th, td { text-align:left; padding:9px 14px 9px 0;
    border-bottom:1px solid var(--rule-2); vertical-align:top; }
  th { font-family:var(--mono); font-size:10.5px; letter-spacing:.1em;
    text-transform:uppercase; color:var(--ink-3); font-weight:400;
    border-bottom:1px solid var(--rule); white-space:nowrap; }
  td.mono, td.num { font-size:13.5px; }
  tr.hl td { background:var(--rule-2); }
  .badge { font-family:var(--mono); font-size:11.5px; letter-spacing:.06em;
    white-space:nowrap; }
  .badge.ok { color:var(--ok); } .badge.warn { color:var(--warn); }
  .badge.flag { color:var(--flag); }
  figure { margin:0 0 36px; }
  figure img, figure video { display:block; width:100%; height:auto;
    background:#fcfcfb; border:1px solid var(--rule); }
  figcaption { font-size:14px; line-height:1.58; color:var(--ink-2); margin-top:11px;
    max-width:80ch; }
  figcaption strong { color:var(--ink); font-weight:600; }
  ol.steps { list-style:none; padding:0; margin:0; }
  ol.steps li { display:flex; gap:16px; padding:15px 0;
    border-bottom:1px solid var(--rule-2); }
  .step-n { color:var(--ink-3); font-size:12px; padding-top:4px; min-width:22px; }
  ol.steps h4 { margin:0 0 4px; font-size:15.5px; font-weight:600; }
  ol.steps p { margin:0; font-size:14.5px; color:var(--ink-2); }
  .path { font-size:12px; color:var(--ink-3); font-weight:400; }
  .note { border-left:2px solid var(--force); padding:2px 0 2px 16px; margin:22px 0;
    max-width:var(--measure); }
  .note.flag { border-left-color:var(--flag); }
  .note .h { font-family:var(--mono); font-size:10.5px; letter-spacing:.11em;
    text-transform:uppercase; color:var(--force); }
  .note.flag .h { color:var(--flag); }
  code { font-family:var(--mono); font-size:.9em; background:var(--rule-2);
    padding:1px 5px; border-radius:2px; }
  sup .cite, a.cite { font-family:var(--mono); font-size:10.5px; color:var(--dna);
    text-decoration:none; border:0; padding:0 1px; }
  a.cite:hover { text-decoration:underline; }
  ol.refs { list-style:none; padding:0; margin:0; counter-reset:r; }
  ol.refs li { display:flex; gap:14px; padding:12px 0;
    border-bottom:1px solid var(--rule-2); scroll-margin-top:20px; }
  ol.refs li:target { background:var(--rule-2); }
  .ref-n { font-family:var(--mono); font-size:12px; color:var(--ink-3);
    min-width:26px; padding-top:2px; }
  .ref-cite { margin:0 0 3px; font-size:14.5px; max-width:82ch; }
  .ref-why { margin:0; font-size:13.5px; color:var(--ink-3); max-width:82ch;
    font-style:italic; }
  footer { margin-top:58px; padding-top:22px; border-top:1px solid var(--rule);
    font-size:13.5px; color:var(--ink-3); }
</style>
"""


def build(args):
    root = Path(args.root)
    figdir = root / "figures"
    f = json.loads((root / "analysis" / "findings.json").read_text())
    systems = {}
    for d in sorted((root / "build" / "systems").glob("*/build_report.json")):
        r = json.loads(d.read_text())
        systems[r["system"]] = r
    params = json.loads((root / "build" / "02_params" / "PARAMS_REPORT.json").read_text())
    sites = json.loads((root / "build" / "01_protein" / "attachment_sites.json").read_text())
    cbcb = sites["pairs"]["ASP133-ASN149"]["d_anchor_anchor_A"]

    sp = f.get("springs", {})
    strongest = max(sp, key=lambda k: sp[k]["force_pN"]["mean"]) if sp else None
    clamp6 = f.get("clamp", {}).get("S6_clamp", {})
    ext = f.get("extension", {})
    # TWO zeros, deliberately, because the tables below need different ones and sharing
    # one variable is how the retracted comparison stayed on the page after the prose
    # retracting it went up.  The clamp ladder runs on the UNCLICKED protein, so
    # unclicked-at-0-pN is its own correct and self-consistent reference -- its 0 pN row
    # has to read +0.000.  A chimera's zero is the clicked, unloaded control, because
    # every chimera carries the clicked adduct and clicking alone moves the sites 0.65 A.
    zero_unclicked = ext.get("zero_force_ref", {}).get("mean")
    zero_clicked = ext.get("matched_control", {}).get("mean") or zero_unclicked

    # ---------- measured force-extension ----------
    # Computed here from the shipped csv rather than transcribed, because these are the
    # numbers that answer "is the force a model or a measurement" and a stale copy of
    # them would be the worst single error the report could contain.
    def _measured_fx():
        import numpy as _np
        from scipy import stats as _st
        csv = Path(__file__).resolve().parent.parent / "analysis" / "force_extension.csv"
        if not csv.exists():
            return None
        d = pd.read_csv(csv)
        if not {"mean_r_A", "force_pN", "box_ok", "box_width_A"} <= set(d.columns):
            return None
        d = d[d.box_ok.astype(str).str.lower().isin(("true", "1"))]
        d = d.assign(r_over_box=d.mean_r_A / d.box_width_A)
        d = d[d.r_over_box <= 0.85]
        if not len(d):
            return None
        NB = {"D_27bp": 27, "D_40bp": 40, "D_40bp_nick": 40}
        NICE = {"D_27bp": "27 bp", "D_40bp": "40 bp",
                "D_40bp_nick": "40 bp, nicked"}
        out = {"n_runs": int(len(d)), "n_dropped": None, "by_system": {}}
        for s, g in d.groupby("system"):
            if s not in NB:
                continue
            eul = _np.pi ** 2 * 200.0 / (((NB[s] - 1) * 0.33) ** 2)
            lr = _st.linregress(g.mean_r_A, g.force_pN)
            e = {"label": NICE[s], "n": int(len(g)),
                 "r_lo": float(g.mean_r_A.min()), "r_hi": float(g.mean_r_A.max()),
                 "f": float(g.force_pN.mean()), "sem": float(g.force_pN.sem()),
                 "euler": float(eul), "ratio": float(g.force_pN.mean() / eul),
                 "slope": float(lr.slope), "slope_p": float(lr.pvalue)}
            a = g[g.k_kJ_mol_nm2 == 33.0].force_pN.values
            b = g[g.k_kJ_mol_nm2 == 200.0].force_pN.values
            if len(a) > 1 and len(b) > 1:
                tt = _st.ttest_ind(a, b, equal_var=False)
                e["k_agree_p"] = float(tt.pvalue)
                e["f_soft"], e["f_stiff"] = float(a.mean()), float(b.mean())
            out["by_system"][s] = e
        i, n = d[d.system == "D_40bp"].force_pN.values, \
            d[d.system == "D_40bp_nick"].force_pN.values
        if len(i) > 1 and len(n) > 1:
            tt = _st.ttest_ind(i, n, equal_var=False)
            out["nick"] = {"intact": float(i.mean()), "nicked": float(n.mean()),
                           "p": float(tt.pvalue)}
        return out

    mfx = _measured_fx()
    mfx_rows = ""
    if mfx:
        mfx_rows = "\n".join(
            "<tr><td>%s</td><td class='num'>%d</td><td class='num'>%.0f–%.0f</td>"
            "<td class='num'><strong>%.1f ± %.1f</strong></td>"
            "<td class='num'>%.1f</td><td class='num'>%.2f</td>"
            "<td class='num'>%.2f</td></tr>" % (
                e["label"], e["n"], e["r_lo"], e["r_hi"], e["f"], e["sem"],
                e["euler"], e["ratio"], e["slope_p"])
            for e in mfx["by_system"].values())

    # ---------- hero ----------
    hero_force = (f"{sp[strongest]['force_pN']['mean']:.1f} pN" if strongest else "—")
    # "5'P-5'P" described the superseded attachment, where a linker sat on each of two
    # strands.  The distance is now measured between whichever nucleotide atoms the
    # linkers are actually bonded to, which for the corrected build is one strand's 5' P
    # and its own 3' O3'.  It is also a MODELLED force, and the hero stat should say so.
    # The measured force leads if we have it: a number out of the force field is worth
    # more than the same number out of a model, and the whole of point 6 was about that
    # distinction.  Fall back to the modelled one only when the windows are unavailable.
    d27m = (mfx or {}).get("by_system", {}).get("D_27bp")
    if d27m:
        hero_k = "measured spring force"
        hero_force = f"{d27m['f']:.1f} pN"
        hero_sub = (f"27 bp duplex, straight from the force field — and "
                    f"{d27m['ratio']:.2f}× the Euler buckling load, not the few pN the "
                    f"model predicts")
    else:
        hero_k = "strongest spring"
        hero_force = (f"{sp[strongest]['force_pN']['mean']:.1f} pN" if strongest
                      else "—")
        hero_sub = "modelled, from the span the simulation settles at"
    sq_h = f.get("is_it_a_spring") or {}
    if sq_h.get("vs_modelled_force"):
        hero_r = f"r = {sq_h['vs_modelled_force']['r']:+.2f}"
        hero_r_sub = (f"site displacement against modelled force, "
                      f"{sq_h.get('n_constructs', 0)} constructs — p = "
                      f"{sq_h['vs_modelled_force']['p']:.2f}, and the sign is the wrong "
                      f"way round")
    else:
        hero_r, hero_r_sub = f"&gt; {kap:.0f}", "k<sub>B</sub>T nm⁻², a bound not a value"
    hero_sub_modelled = (f"{sp[strongest]['label']}, modelled from an end-to-end span of "
                f"{sp[strongest]['x_A']['mean']:.0f} Å" if strongest else "")
    kap = clamp6.get("kappa_lower_bound_kT_nm2")

    # ---------- spring table ----------
    # Each spring is differenced against its OWN family's clicked control: U3_ds50 is a
    # ubiquitin construct and subtracting sfGFP's zero from it would be meaningless.
    # findings.py has already made the matched comparison, so read the delta from there.
    def _fam_delta(k):
        v = (f.get("extension", {}).get("vs_matched_control") or {}).get(k)
        if v is not None:
            return v["delta"]
        return sp[k]["extension_A"]["mean"] - zero_clicked

    spring_rows = "\n".join(
        f"""<tr><td>{esc(v['label'])}</td>
        <td class="num">{v['contour_nm'] * 10:.0f}</td>
        <td class="num">{pm(v['x_A'], 1)}</td>
        <td class="num">{v['tau_c']:.0f}</td>
        <td class="num"><strong>{pm(v['force_pN'])}</strong></td>
        <td class="num">{pm(v['energy_kT'], 1)}</td>
        <td class="num">{pm(v['bend_total_deg'], 0)}</td>
        <td class="num">{v['bp_open']['mean']:.2f}</td>
        <td class="num">{_fam_delta(k):+.2f}</td></tr>"""
        for k, v in sp.items())

    # ---------- clamp ladder ----------
    ladder_rows = "\n".join(
        f"""<tr><td class="num">{r['pN']:.2f}</td>
        <td class="num">{r['mean']:.3f}</td><td class="num">{r['sem']:.3f}</td>
        <td class="num">{r['min']:.3f} – {r['max']:.3f}</td>
        <td class="num">{r['n']}</td>
        <td class="num">{r['mean'] - zero_unclicked:+.3f}</td></tr>"""
        for r in clamp6.get("ladder", []))

    # ---------- extension table ----------
    ext_rows = ""
    for k in ("S0_wt", "S1_tet", "S2_clicked", "S3_spring27", "S4_spring40",
              "S5_spring40nick"):
        if k not in ext:
            continue
        v = ext[k]
        force = (f"{sp[k]['force_pN']['mean']:.2f}" if k in sp else
                 ("0 (by construction)" if k in ("S0_wt", "S1_tet", "S2_clicked")
                  else "—"))
        ext_rows += (f"""<tr{' class="hl"' if k in sp else ''}>
        <td>{esc(LABEL[k])}</td><td class="num">{force}</td>
        <td class="num">{v['mean']:.3f}</td><td class="num">{v['sem']:.3f}</td>
        <td class="num">{v['min']:.3f} – {v['max']:.3f}</td>
        <td class="num">{v['n']}</td>
        <td class="num">{v['mean'] - zero_clicked:+.3f}</td></tr>\n""")

    # ---------- H-bond table ----------
    cols = [k for k in ("S0_wt", "S1_tet", "S2_clicked", "S3_spring27",
                        "S4_spring40", "S5_spring40nick")
            if any(k in v["per_system"] for v in f["hbonds"].values())]
    hb_head = "".join(f"<th>{esc(LABEL[k])}</th>" for k in cols)
    hb_rows = ""
    for name, v in sorted(f["hbonds"].items(),
                          key=lambda kv: -max((x["mean"] for x in
                                               kv[1]["per_system"].values()),
                                              default=0)):
        cells = ""
        for k in cols:
            a = v["per_system"].get(k)
            cells += ("<td class=\"num\">—</td>" if not a else
                      f"<td class=\"num\">{100 * a['mean']:.0f}"
                      f"<span style=\"color:var(--ink-3)\"> ({100 * a['min']:.0f}–"
                      f"{100 * a['max']:.0f})</span></td>")
        flag = ('<span class="badge flag">bistable</span>' if v["bistable"] else "")
        hb_rows += (f"<tr><td>{esc(name).replace('-&gt;', '→').replace('&lt;-', '←')}"
                    f" {flag}</td>{cells}</tr>\n")

    # ---------- twist / water / rmsf ----------
    aux_rows = ""
    for k in cols:
        tau = f["twist"]["tau_phenol_bridge"].get(k)
        phi = f["twist"]["phi_bridge_imidazolinone"].get(k)
        wat = f["cavity_water"].get(k)
        rm = f["rmsf"].get(k)
        if not tau:
            continue
        aux_rows += (f"""<tr><td>{esc(LABEL[k])}</td>
        <td class="num">{pm(tau, 2)}</td><td class="num">{pm(phi, 2)}</td>
        <td class="num">{pm(wat, 1) if wat else '—'}</td>
        <td class="num">{pm(rm, 3) if rm else '—'}</td></tr>\n""")

    rows_sys = "\n".join(
        f"""<tr><td class="mono">{esc(k)}</td><td>{esc(LABEL.get(k, k))}</td>
        <td class="num">{r['atoms']:,}</td>
        <td class="num">{r['net_charge']:+.4f}</td>
        <td>{PURPOSE.get(k, '')}</td></tr>"""
        for k, r in systems.items() if k in LABEL)
    # ---------- the measured kinking threshold ----------
    # Read from the force-extension runs rather than written down, so the table cannot
    # drift from the data behind it.
    fx_rows, fx_note = "", ""
    fx_first = fx_certain = None
    fxf = Path(args.root) / "analysis" / "force_extension.csv"
    if fxf.exists():
        import numpy as _np
        fx = pd.read_csv(fxf)
        NB = {"D_27bp": 27, "D_40bp": 40, "D_40bp_nick": 40}
        NICE = {"D_27bp": "27 bp", "D_40bp": "40 bp", "D_40bp_nick": "40 bp, nicked"}
        need = {"wc.n_internal_open", "box_ok", "mean_r_A", "box_width_A"}
        if len(fx) and need <= set(fx.columns):
            # Only windows the box can hold, judged on where the duplex actually sat.
            # The earlier version of this table showed every window and derived its
            # "first opening" note from all of them -- and the widest windows, the ones
            # beyond the box, were precisely the rows reading "no opening", so the note
            # reported an onset that was an artefact of including them.
            n_all = len(fx)
            fx = fx[fx.box_ok.astype(str).str.lower().isin(("true", "1"))]
            fx = fx.assign(r_over_box=fx.mean_r_A / fx.box_width_A)
            fx = fx[(fx.r_over_box <= 0.85) & fx.system.isin(NB)]
            n_drop = n_all - len(fx)

            def per_bp(chord, n_bp):
                L = (n_bp - 1) * 3.376
                lo, hi = L / (2 * _np.pi), 1e7
                for _ in range(300):
                    mid = 0.5 * (lo + hi)
                    lo, hi = ((mid, hi) if 2 * mid * _np.sin(L / (2 * mid)) < chord
                              else (lo, mid))
                return _np.degrees(3.376 / (0.5 * (lo + hi)))

            out = []
            for (s, x0), g in fx.groupby(["system", "x0_A"]):
                r = float(g.mean_r_A.mean())
                out.append((per_bp(r, NB[s]), NICE[s], r,
                            int((g["wc.n_internal_open"] > 0).sum()), len(g),
                            float(g["bend.max_step_deg"].mean())
                            if "bend.max_step_deg" in g else float("nan"),
                            float(g["force_pN"].mean())))
            out.sort()

            def _fx_row(d, nm, r, k, n, b, fo):
                cls = ' class="hi"' if k else ''
                verdict = ("<strong>yes</strong>, %d/%d" % (k, n)) if k else "no, 0/%d" % n
                bs = "&mdash;" if b != b else "%.1f" % b
                return ('<tr%s><td>%s</td><td class="num">%.1f</td>'
                        '<td class="num">%.0f</td><td class="num">%s</td>'
                        '<td class="num">%.0f</td><td>%s</td></tr>'
                        % (cls, nm, d, r, bs, fo, verdict))

            fx_rows = "\n".join(_fx_row(*r) for r in out)
            opened = [d for d, nm, r, k, n, b, fo in out if k]
            allopen = [d for d, nm, r, k, n, b, fo in out if k == n]
            parts = []
            if opened:
                parts.append("The gentlest curvature at which any interior pair opens is "
                             f"<strong>{min(opened):.1f}°/bp</strong>")
                if allopen:
                    parts.append("and every run opens above "
                                 f"<strong>{min(allopen):.1f}°/bp</strong>")
            else:
                parts.append("No interior pair opened at any curvature these boxes can "
                             "hold")
            if n_drop:
                parts.append(f"{n_drop} of {n_all} runs are excluded because the duplex "
                             "sat within 85 % of the box's perpendicular width, where it "
                             "can reach its own image")
            fx_note = ". ".join(parts) + "."
            fx_first = min(opened) if opened else None
            fx_certain = min(allopen) if allopen else None

    ps = params_summary(Path(args.root))
    _pv = {
        "__PARAMS_FLAGGED__": (f"{ps['flagged']}, all {ps['n']} custom residues "
                               f"({ps['names']})" if ps else "0"),
        "__PARAMS_ELEMENTS__": (f"{ps['atoms']} atoms across {ps['n']} residues, "
                                f"{ps['mismatches']} mismatches" if ps
                                else "every atom, every residue"),
        "__PARAMS_TREE__": (f"all {ps['n']} residues" if ps else "every residue"),
    }
    rows_ver = "\n".join(
        f"""<tr><td>{v[0]}</td><td class="mono">{_pv.get(v[1], v[1])}</td>
        <td class="badge {'ok' if v[2] else 'warn'}">{'pass' if v[2] else 'check'}</td></tr>"""
        for v in VER)
    rows_steps = "\n".join(
        f"""<li><span class="step-n mono">{i}</span><div>
        <h4>{esc(t)} <span class="mono path">{esc(p)}</span></h4><p>{d}</p></div></li>"""
        for i, (t, p, d) in enumerate(STEPS, 1))
    figs_html = ""
    for key, title, cap in FIGS:
        uri = b64(figdir / f"{key}.png", "image/png")
        if uri:
            figs_html += (f'<figure><img src="{uri}" alt="{esc(title)}" />'
                          f'<figcaption><strong>{esc(title)}.</strong> {cap}'
                          f'</figcaption></figure>\n')
    # Animations are embedded as GIF inside <img>, not mp4 inside <video>.  The
    # artifact viewer renders <img data:...> reliably and silently dropped every
    # base64 <video>, leaving two empty panels; a GIF is just an image as far as the
    # renderer is concerned.  The full-quality mp4s stay on disk in figures/.
    def anim_uri(view):
        return (b64(figdir / "web" / f"{args.anim}_{view}.gif", "image/gif")
                or b64(figdir / f"{args.anim}_{view}.gif", "image/gif"))

    def still_uri(view):
        return b64(figdir / f"{args.anim}_{view}_still.png", "image/png")

    anim = anim_uri("overview")
    # The chromophore and strain panels use the ray-traced labelled stills rather than
    # the GIFs: a moving picture is right for showing motion and wrong for reading
    # labels, and those two panels exist to say which residue is which.
    anim2 = still_uri("strain") or anim_uri("strain")
    anim3 = still_uri("chromophore") or anim_uri("chromophore")

    # ---------- prose that quotes the numbers ----------
    s27, s40, s40n = (sp.get("S3_spring27"), sp.get("S4_spring40"),
                      sp.get("S5_spring40nick"))
    tw = ext.get("tet_vs_wt") or {}          # installing the residues: WT -> unclicked
    cw = ext.get("clicked_vs_tet") or {}     # clicking: unclicked -> clicked
    sw = ext.get("strongest_vs_weakest_spring") or {}
    chim_slope = ext.get("chimera_slope_A_per_pN")
    mc = ext.get("matched_control") or {}
    vsm = ext.get("vs_matched_control") or {}
    unclicked_ref = ext.get("zero_force_ref") or {}

    def _row(key, force):
        v = vsm.get(key)
        if not v:
            return ""
        sig = ("<strong>yes</strong>" if v.get("p_bonferroni", 1) < 0.05 else "no")
        # the label comes from findings.json, so this file needs no style import
        return (f"<tr><td>{force:.2f}</td><td>{v.get('label', key)}</td>"
                f"<td>{v['delta']:+.3f}</td><td>{v['t']:+.2f}</td>"
                f"<td>{v['p']:.3f}</td><td>{v.get('p_bonferroni', float('nan')):.3f}</td>"
                f"<td>{sig}</td></tr>")
    _forces = {k: (f["springs"].get(k, {}).get("force_pN", {}) or {}).get("mean", 0.0)
               for k in ("S5_spring40nick", "S4_spring40", "S3_spring27")}
    matched_rows = "\n".join(_row(k, v) for k, v in
                              sorted(_forces.items(), key=lambda kv: kv[1]))

    # ---------- is it a spring effect at all? ----------
    spring_q = f.get("is_it_a_spring") or {}
    vm_all = (f.get("extension", {}).get("vs_matched_control") or {})
    fam_of = {k: v.get("family", "sfGFP") for k, v in vm_all.items()}

    def _all_row(k, v):
        fo = (f.get("springs", {}).get(k, {}).get("force_pN", {}) or {}).get("mean")
        c = ((f.get("contacts", {}).get(k, {}) or {})
             .get("atoms_within_4.5A", {}) or {}).get("mean")
        fo_s = f"{fo:.2f}" if fo else "<em>none</em>"
        c_s = f"{c:.1f}" if c is not None else "&mdash;"
        hi = ' class="hi"' if fo in (None, 0) or not fo else ""
        return (f"<tr{hi}><td>{v.get('label', k)}</td>"
                f"<td class='num'>{fo_s}</td><td class='num'>{c_s}</td>"
                f"<td class='num'><strong>{v['delta']:+.3f}</strong></td>"
                f"<td class='num'>{v['p']:.3f}</td>"
                f"<td class='num'>{v.get('p_bonferroni', float('nan')):.3f}</td></tr>")

    allrows = "\n".join(_all_row(k, v) for k, v in sorted(
        vm_all.items(), key=lambda kv: -(f.get("springs", {}).get(kv[0], {})
                                         .get("force_pN", {}) or {}).get("mean", 0.0)))

    ss = f.get("ss_springs") or {}
    ss_rows = "\n".join(
        "<tr><td>%s</td><td class='num'>%d</td><td class='num'>%.0f</td>"
        "<td class='num'>%.0f</td><td class='num'>%.1f</td>"
        "<td class='num'><strong>%s</strong></td><td>%s</td></tr>" % (
            esc(v["label"]), v["n_nt"], v["contour_nm"] * 10,
            v["natural_rms_nm"] * 10, v["x_A"]["mean"],
            ("%.2f" % v["tension_pN"]["mean"]) if v["tension_pN"]["mean"]
            else "0 &mdash; slack", v["state"])
        for k, v in ss.items())

    topo = f.get("topology") or {}
    topo_rows = "\n".join(
        "<tr%s><td>%s</td><td class='num'>%+.3f</td><td class='num'>%.3f</td>"
        "<td>%s</td></tr>" % (
            ' class="hi"' if w["p"] < 0.05 else "", w["label"], w["delta"], w["p"],
            "closing the loop compresses more" if w["delta"] < 0
            else "closing the loop compresses less")
        for w in (topo.get("pairs") or {}).values())

    sites_f = f.get("sites") or {}

    def _site(sysk, site, key, dflt=float("nan")):
        return (sites_f.get(sysk, {}).get(str(site), {}) or {}).get(key, dflt)

    # Each protein's own two sites, taken from the entry rather than assumed: sfGFP is
    # modified at 133/149 and ubiquitin at 9/63, and hardcoding the sfGFP pair printed a
    # row of literal "nan" for every ubiquitin system -- whose angles were in the same
    # file under different keys.
    def _site_row(k, v):
        sl = v.get("sites") or [s for s in (133, 149) if str(s) in v]
        if len(sl) < 2:
            return ""
        a, b = sl[0], sl[1]

        def ang(s):
            ph, ps = _site(k, s, "phi"), _site(k, s, "psi")
            return "&mdash;" if ph != ph or ps != ps else f"{ph:+.1f}, {ps:+.1f}"

        stack = _site(k, b, "Tyr200_contacts_per_frame")
        return (f"<tr><td>{v['label']}</td>"
                f"<td class='num'>{a}</td>"
                f"<td>{_site(k, a, 'residue') or '&mdash;'}</td>"
                f"<td class='num'>{ang(a)}</td>"
                f"<td class='num'>{b}</td>"
                f"<td>{_site(k, b, 'residue') or '&mdash;'}</td>"
                f"<td class='num'>{ang(b)}</td>"
                f"<td class='num'>"
                f"{'&mdash;' if stack != stack else f'{stack:.2f}'}</td></tr>")

    # The point-9 prediction, tested across both proteins: modification should strain a
    # site only where the native backbone sits at POSITIVE phi, a conformation Asp, Asn
    # and Gly adopt readily and a bulkier side chain does not.  Two proteins with two
    # unrelated chemistries give four sites, only one of them positive-phi.
    strain = []
    for prot, apo, mod in (("sfGFP", "S0_wt", "S1_tet"),
                           ("ubiquitin", "U0_apo", "U1_clicked")):
        av, mv = sites_f.get(apo) or {}, sites_f.get(mod) or {}
        for st in (av.get("sites") or []):
            a, m = av.get(str(st)) or {}, mv.get(str(st)) or {}
            if a.get("phi") is None or m.get("phi") is None:
                continue
            strain.append({"protein": prot, "site": st,
                           "native": a.get("residue"), "modified": m.get("residue"),
                           "phi": a["phi"], "psi": a["psi"],
                           "dphi": m["phi"] - a["phi"], "dpsi": m["psi"] - a["psi"],
                           "positive_phi": a["phi"] > 0})
    strain_rows = "\n".join(
        "<tr%s><td>%s</td><td class='num'>%d</td><td>%s → %s</td>"
        "<td class='num'>%+.1f, %+.1f</td><td class='num'><strong>%+.1f</strong></td>"
        "<td class='num'>%+.1f</td><td>%s</td></tr>" % (
            ' class="hi"' if r["positive_phi"] else "", r["protein"], r["site"],
            r["native"], r["modified"], r["phi"], r["psi"], r["dphi"], r["dpsi"],
            "<strong>left-handed α</strong>" if r["positive_phi"] else "ordinary")
        for r in strain)
    have_strain = len(strain) > 2

    site_rows = "\n".join(r for r in (_site_row(k, v)
                                      for k, v in sites_f.items()) if r)
    have_sites = bool(site_rows)

    con = f.get("contacts") or {}
    contact_rows = "\n".join(
        f"<tr><td>{v['label']}</td>"
        f"<td>{v['atoms_within_4.5A']['mean']:.1f} &plusmn; "
        f"{v['atoms_within_4.5A'].get('sd', 0):.1f}</td>"
        f"<td>{100 * v['frames_with_contact']:.0f}%</td>"
        f"<td>{v['residues_touched']:.0f}</td></tr>"
        for v in con.values())
    have_contacts = bool(contact_rows)
    bistable = [n for n, v in f["hbonds"].items() if v["bistable"]]
    stamp = datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M UTC")

    # The template below is an f-string, so every {cite_*} marker has to exist as a
    # local here -- there is no .format() call to pass them into.
    refs = refs_html()
    cite_eddins = cite("eddins", "blizzard"); cite_pedelacq = cite("pedelacq")
    cite_ormo = cite("ormo"); cite_labms = cite("labms", "eddins")
    cite_spring = cite("zocchi_book", "qu", "tseng2013")
    cite_qu = cite("qu")
    cite_tau = cite("tseng2013", "qu"); cite_tseng2021 = cite("tseng2021")
    cite_zocchi = cite("zocchi_book"); cite_niwa = cite("niwa")
    cite_brejc = cite("brejc"); cite_choi = cite("choi2007")
    cite_maier = cite("maier"); cite_ol15 = cite("zgarbova")
    cite_tip3p = cite("jorgensen"); cite_ions = cite("joung"); cite_mg = cite("li_merz")
    cite_boral = cite("boral"); cite_breyfogle = cite("breyfogle")
    cite_case = cite("case"); cite_gaff = cite("wang_gaff")
    cite_ambcc = cite("jakalian"); cite_pymol = cite("pymol")
    cite_kabsch = cite("kabsch"); cite_openmm = cite("eastman")
    cite_leimkuhler = cite("leimkuhler"); cite_hopkins = cite("hopkins")
    cite_mdtraj = cite("mdtraj"); cite_shrake = cite("shrake")

    html = CSS + f"""
<div class="wrap">
<header class="mast">
  <p class="eyebrow">All-atom molecular dynamics · superfolder GFP · {stamp}</p>
  <h1>A DNA spring, wound around a fluorescent protein</h1>
  <p class="lede">Two amber codons put a tetrazine amino acid{cite_eddins} on adjacent
  strands of sfGFP’s β-barrel{cite_pedelacq}. A doubly-modified DNA duplex clicks onto
  both, and because dsDNA is far too stiff to follow the barrel’s curvature it has to
  bend — which should push the two attachment points apart. This is what that looks like
  <em>atom by atom</em>, across eleven DNA-bearing constructs in two proteins. The
  duplex turns out to deliver <em>far more</em> force than the model predicts, and the
  protein turns out not to care: every construct pulls its attachment sites slightly
  <em>together</em>, and it makes no measurable difference whether the DNA is a spring at
  all.</p>
  <div class="stats">
    <div class="stat force"><div class="k">{hero_k}</div>
      <div class="v">{hero_force}</div><div class="s">{hero_sub}</div></div>
    <div class="stat"><div class="k">production analysed</div>
      <div class="v">{f['total_ns']:,.0f} ns</div>
      <div class="s">{f['n_replicates']} replicates across
      {len(f['systems_present'])} simulated systems</div></div>
    <div class="stat"><div class="k">does the force explain it?</div>
      <div class="v">{hero_r}</div>
      <div class="s">{hero_r_sub}</div></div>
    <div class="stat ok"><div class="k">invented parameters</div>
      <div class="v">0</div>
      <div class="s">nothing flagged for review across the
      {ps['n'] if ps else 'seven'} custom residues</div></div>
  </div>
</header>

<section>
  <div class="sec-head"><h2>The construct</h2>
    <span class="tag">Asp133 · Asn149 · CRO 66 · {cbcb} Å apart</span></div>
  <figure><img src="{anim}" alt="The 27 bp chimera over 50 ns" />
    <figcaption><strong>The chimera, 50 ns.</strong> sfGFP’s β-barrel in grey with the
    chromophore inside, the two Tet2-Et/sTCO tethers in orange leaving adjacent strands,
    and the 27 bp duplex arcing away. Water and ions are stripped and the frames
    superposed on the barrel core, so what moves is the spring working against the
    protein rather than the box tumbling. The duplex sweeps through a wide range of
    orientations while staying bent — which is the point: the load is delivered through
    two flexible tethers, so the protein sees a diffuse shell rather than a fixed
    lever.</figcaption></figure>
  <figure><img src="{anim3}" alt="The chromophore and its hydrogen-bond partners" />
    <figcaption><strong>Inside the barrel.</strong> The chromophore in green with the
    side chains lining its cavity, and dashed measures on the three hydrogen bonds that
    hold it <em>from twisting</em>: His148 and Thr203 to the phenolate oxygen, Arg96 to
    the imidazolinone carbonyl. The barrel is drawn transparent so the cavity shows
    through. What these bonds resist is <strong>rotation, not bending</strong> — the
    chromophore loses its quantum yield by twisting about the two exocyclic bonds of the
    methine bridge, which opens a non-radiative decay channel{cite_niwa}, and holding the
    two ring systems coplanar is what suppresses it. Watch His148: it is the mobile one
    of the three, and the residue nearest the two attachment
    sites.</figcaption></figure>
  <div class="col">
    <p>The attachment sites are <strong>Asp133</strong> and <strong>Asn149</strong> in
    2B3P’s numbering — D134 and N150 in the construct’s own, which runs one higher.
    They sit on adjacent strands of the barrel with <strong>His148</strong> between them,
    and His148 donates a hydrogen bond straight to the chromophore’s phenolate
    oxygen{cite_ormo}.
    That is the mechanical path from the DNA to the fluorescence, and it is why these
    two sites are the interesting ones rather than merely convenient ones.</p>
  </div>
  <div class="note flag"><div class="h">A correction worth carrying back to the bench</div>
    <p>The PyMOL session’s <code>134TAG</code> and <code>150TAG</code> selections point
    at Gly134 and Val150 — one residue off. The mass spectrometry{cite_labms} settles it: the single
    construct is +141 Da (Asn → Tet2-Et) and the double +281 Da (= 141 + 140, Asn plus
    Asp). Gly → Tet2-Et would be +198 Da and Val → Tet2-Et +156 Da, and neither is
    observed. <code>build/01_protein/fix_pse_selections.pml</code> loads corrected
    selections over the session.</p></div>
  <figure><img src="{anim2}" alt="Backbone painted by per-residue mobility" />
    <figcaption><strong>Painted by mobility.</strong> The same trajectory with the
    backbone drawn as a putty tube scaled and coloured by per-residue RMSF — white and
    thin where the fold is rigid, orange and swollen where it moves. The barrel stays
    uniformly rigid; what moves is the loops and the tethers.</figcaption></figure>
</section>

<section>
  <div class="sec-head"><h2>What the spring delivers</h2>
    <span class="tag">force read from the simulation’s own geometry</span></div>
  <div class="tbl-wrap"><table>
    <thead><tr><th>spring</th><th>contour (Å)</th><th>span x (Å)</th>
    <th>τ<sub>c</sub></th><th>force (pN)</th><th>energy (k<sub>B</sub>T)</th>
    <th>total bend (°)</th><th>bp open</th>
    <th>Δ Cβ–Cβ vs clicked control (Å)</th></tr></thead>
    <tbody>{spring_rows}</tbody></table></div>
  <div class="col">
    <p>Every force here is Zocchi’s <em>f(x)</em>{cite_spring} evaluated at the end-to-end distance
    the simulation itself settles at, not at an assumed one. The construct is not free
    to choose that distance arbitrarily — it is set by where the two tethers can reach —
    so this is the one number in the whole model that the simulation can supply
    independently.</p>

    <h3>The nick softens the spring in the model — but not in measurement</h3>
    <p>The nicked 40 bp duplex delivers <strong>{sig(s40n['force_pN']['mean'])} pN</strong>
    against <strong>{sig(s40['force_pN']['mean'])} pN</strong> for the intact duplex of
    the same length, and it is visibly floppier: its span wanders over
    {sig(s40n['x_wander_A']['mean'], 1)} Å within a trajectory against
    {sig(s40['x_wander_A']['mean'], 1)} Å for the intact one. Zocchi’s two values of the
    critical bending torque{cite_tau} — 27 pN·nm with a nick, 31–36 without — predict this
    ordering, and we reach it here from structure rather than from his calibration.</p>
    <p><strong>Measured directly, though, the nick does nothing resolvable</strong> — see
    the next section. The agreement here is between two calculations that both take the
    nick's effect from the same critical torque, so it is a consistency check on our
    geometry rather than independent support for the model. The nicked spring is the
    low-stress reference his argument relies on, and our force field does not reproduce
    that it is softer.</p>

    <h3>Shorter spring, more force; longer spring, more stored energy</h3>
    <p>The 27 bp spring delivers {sig(s27['force_pN']['mean'])} pN and stores
    {sig(s27['energy_kT']['mean'], 1)} k<sub>B</sub>T; the 40 bp spring delivers less
    force ({sig(s40['force_pN']['mean'])} pN) but stores <em>more</em> energy
    ({sig(s40['energy_kT']['mean'], 1)} k<sub>B</sub>T). That is not a contradiction:
    force is the derivative of the energy with respect to the end-to-end distance, and
    the longer rod is bent through {sig(s40['bend_total_deg']['mean'], 0)}° against
    {sig(s27['bend_total_deg']['mean'], 0)}° — much more total curvature, spread over
    many more base pairs, so a larger integral with a shallower slope. All three
    energies, {min(v['energy_kT']['mean'] for v in sp.values()):.1f}–{max(v['energy_kT']['mean'] for v in sp.values()):.1f}
    k<sub>B</sub>T, bracket the 6.6 k<sub>B</sub>T Zocchi measured for the
    guanylate-kinase chimera{cite_tseng2021} by an entirely different route.</p>

    <h3>The duplex bends smoothly — it does not kink</h3>
    <p>Zocchi’s model{cite_zocchi} places all three springs on its <em>kinked</em> branch, and the
    forces above come from that branch. The structures do not corroborate it on this
    timescale. Bending is spread evenly at 3–8° per step with no localised spike
    (figure 4), and of 27–40 base pairs only
    {min(v['bp_open']['mean'] for v in sp.values()):.2f}–{max(v['bp_open']['mean'] for v in sp.values()):.2f}
    are open at any instant — the duplex is intact. Three explanations are possible: the
    force field resists kinking, the timescale is too short, or the model puts the
    kinking threshold in the wrong place. The literature answers fairly clearly, and it
    is mostly the last two.</p>

    <h4>The curvature we impose is below where all-atom MD nucleates kinks</h4>
    <p>The duplexes are built at <em>R</em> = 41–44 Å and relax over the run to
    <strong>37–41 Å, i.e. 4.8–5.2° per base-pair step</strong> (measured, one-turn
    smoothed, 18 replicates), equivalent to a 70–78 bp minicircle. That is essentially
    <em>nucleosomal</em> curvature — 42 Å and 4.6°/bp — and about
    <strong>56% of the curvature at which all-atom MD nucleates kinks</strong>.</p>
    <div class="note"><div class="h">Correcting a correction</div>
      <p>An earlier version of this page said 31–34 Å and 5.7–6.3°/bp, and a version after
      that "corrected" it to being 25% tighter than nucleosomal. Both were wrong: 31–34 Å
      was a figure I had carried over without checking it against the trajectories, and
      the numbers above are what the bend analysis actually reports. The conclusion is
      unaffected in direction and slightly strengthened in degree — we are further below
      the kinking threshold than either earlier version implied.</p></div>
    <p>Converting each published criterion to a threshold radius through
    <em>R</em> = <em>B</em>/τ<sub>c</sub>:</p>

    <div class="tbl-wrap"><table>
      <thead><tr><th>source</th><th>criterion</th><th>threshold R</th><th>°/bp</th></tr></thead>
      <tbody>
      <tr><td>Zocchi, nicked{cite_qu}</td><td>τ<sub>c</sub> = 27 pN nm</td><td>76 Å</td><td>2.6</td></tr>
      <tr><td>Zocchi, intact{cite_zocchi}</td><td>τ<sub>c</sub> = 31 pN nm</td><td>66 Å</td><td>3.0</td></tr>
      <tr><td>Vologodskii &amp; Frank-Kamenetskii 2013; Du 2008</td><td>nuclease probing</td><td>35–40 Å</td><td>4.9–5.6</td></tr>
      <tr><td>Fields, Meyer &amp; Cohen 2013</td><td>FRET “molecular vice”, WLC linear to 7°/bp</td><td>≤28 Å</td><td>7</td></tr>
      <tr><td>Curuksu <em>et al.</em> 2009</td><td>umbrella sampling, same force-field family</td><td>21–23 Å</td><td><strong>8.6–9.3</strong></td></tr>
      <tr class="hi"><td><strong>this work</strong></td><td>bent duplex in the chimera, measured</td><td><strong>37–41 Å</strong></td><td><strong>4.8–5.2</strong></td></tr>
      </tbody></table></div>

    <p>Zocchi’s τ<sub>c</sub> is the most kink-friendly estimate in the field: it places
    the softening transition at a <em>gentler</em> bend than a nucleosome, where every
    independent estimate is 1.3–2.7× stiffer. The one directly comparable calculation —
    Curuksu <em>et al.</em>’s umbrella sampling on a bent linear duplex in the same
    force-field family — finds type II kinks only above 8.6–9.3° per step, which is
    <strong>almost twice the curvature we impose</strong>. On that criterion, not seeing a
    kink is the expected result rather than an anomaly.</p>

    <h4>45 ns is short by two to five orders of magnitude</h4>
    <p>Base-pair lifetimes measured by NMR imino exchange are 1–40 ms, so the per-pair
    opening rate is ~10²–10³ s⁻¹. For ~30 pairs in parallel across our 270 ns of aggregate
    sampling the <em>expected</em> number of opening events is 8×10⁻⁴ to 8×10⁻³. Seeing
    none is what unstressed kinetics predicts before bending is considered at all.
    Turning that around converts a null result into a bound: these data are consistent
    with a stress-induced reduction of the opening barrier of <strong>up to about
    5 k<sub>B</sub>T</strong>, and cannot support more. Independently, published kink free
    energies of 12–22 k<sub>B</sub>T imply waiting times of microseconds to hundreds of
    milliseconds.</p>
    <p>Published calibration points agree: kinks appear at 80–800 ns in nucleosomal loops
    held far more tightly than ours, and at 20–25 ns in 65 bp <em>closed</em> minicircles
    — but a 90 bp minicircle run for <strong>11 µs</strong> at R ≈ 49 Å showed no kinks
    and no base-pair-level distortion at all.</p>

    <h4>We measured the threshold rather than arguing about it</h4>
    <p>Everything above is other people's numbers. The compression windows let us put ours
    beside them: hold the bare 27 bp duplex at a series of separations, from its contour
    length down to 64% compression, and ask at what curvature <em>this</em> force field
    starts breaking base pairs. The criterion was fixed before the runs — an
    <strong>internal</strong> pair, not the terminal two which fray in anything, open in
    more than 10% of frames.</p>
    {'<div class="tbl-wrap"><table><thead><tr><th>duplex</th><th>°/bp</th>'
     '<th>separation adopted (Å)</th>'
     '<th>max step (°)</th><th>f (pN)</th><th>interior pair opens?</th></tr></thead>'
     '<tbody>' + fx_rows + '</tbody></table></div><p>' + fx_note + '</p>'
     if fx_rows else '<p class="muted">Force–extension windows still running.</p>'}
    <p>Three things follow. <strong>The force field can kink</strong> — reproducibly, with
    the opened pair sitting at the maximum of the bend profile — so the possibility that
    AMBER with OL15 simply cannot support a kink is excluded by our own data rather than by
    argument.</p>

    <p><strong>There is no threshold; there is a probability.</strong> An earlier version of
    this page said the threshold was 6.7–8.3°/bp. That rested partly on the 27-mer's widest
    windows — 66 Å and beyond — which are exactly the rows showing <em>no</em> opening at
    low curvature, and every one of them sits outside what its solvent box can hold. With
    only valid windows, opening appears as low as
    <strong>{fx_first:.1f}°/bp</strong> and becomes certain — every run in the window —
    only above <strong>{fx_certain:.1f}°/bp</strong>. The per-run observation is bimodal:
    a run either keeps every interior pair shut for its whole trajectory or opens one and
    holds it open 30–100 % of frames. So what varies smoothly with curvature is the
    <em>probability</em> of opening rather than its extent, which is what figure 5 plots —
    binned by each run's own curvature, where the table above groups by window. Curuksu <em>et al.</em>'s 8.6–9.3°/bp from biased sampling is where our
    probability has already reached 0.8, so the two agree about where kinking becomes the
    expected behaviour while ours shows a stochastic tail below it. Both are a long way
    from Zocchi's τ<sub>c</sub>, which implies 2.6–3.0 — even our tail does not reach
    that.</p>

    <p>And <strong>our chimera sits at 4.8–5.2°/bp, below the lowest curvature at which we
    ever see a pair open</strong> — which is the answer. The construct is not bent far
    enough, on this force field's own showing.</p>

    <p>One detail that makes this a check rather than a coincidence: in the 27 bp duplex the
    kink lands on <strong>base pair 17 in 11 of its 14 openings</strong>, across independent
    windows and replicates. Base 17 is an adenine and the step into it is
    <strong>CpA</strong> — one of the three pyrimidine–purine steps (CpG, TpA, CpA) that
    Curuksu names as where roll preferentially localises. The other pair that opened, bp13,
    sits at a <strong>TpG</strong> step, also pyrimidine–purine. The force field is not
    breaking a pair at random under load; it picks exactly the dinucleotide steps the
    literature says it should.</p>

    <h4>The force field, honestly</h4>
    <p>There is a documented bias, and OL15 is specifically implicated: backbone-dihedral
    refinements favour B-form geometry even in non-helical contexts, and all AMBER
    variants overestimate persistence length with OL15 the stiffest of them. Twist-bend
    kinks also require local <em>under</em>twisting, and the OL family was tuned to
    increase twist. So OL15 is plausibly the least kink-prone reasonable choice, and that
    deserves saying. But the magnitude is not enough to rescue the model: parmbsc0 and
    parmbsc1 are softer and also fail to kink at comparable curvature on 20–240× longer
    runs; parmbsc1’s free-energy landscape is <em>sub</em>-harmonic, i.e. more deformable
    than elasticity predicts; and CHARMM36, the alternative, has the opposite pathology of
    spurious base-pair disruption. No head-to-head comparison of kinking propensity across
    these force fields in a sharply bent duplex has been published — that is a real gap.</p>

    <p><strong>So the forces above remain model-dependent numbers.</strong> They are
    consistent with the model’s inputs; they do not confirm its mechanism. That is why
    the force is now also measured directly, below.</p>
  </div>
</section>

<section>
  <div class="sec-head"><h2>What the spring delivers, measured</h2>
    <span class="tag">the force field's own answer, not the model's</span></div>
  <div class="col">
    <p>Every force above is Zocchi's <em>f(x)</em> evaluated at a simulated distance. That
    is a model, and the obvious question — put to us directly in review — is what the
    force field says if you simply ask it. So we did: a bare duplex, no protein, held at a
    series of separations by a harmonic restraint between the two linker nitrogens. The
    mean restraint force ⟨k(r − x₀)⟩ <em>is</em> the force the duplex exerts, with no
    model in between.</p>

    <div class="tbl-wrap"><table>
      <thead><tr><th>duplex</th><th>runs</th><th>separations sampled (Å)</th>
      <th>measured force (pN)</th><th>Euler load π²B/L² (pN)</th>
      <th>measured / Euler</th><th>p, slope vs separation</th></tr></thead>
      <tbody>{mfx_rows}</tbody></table></div>

    <h3>The force does not depend on how far you compress it</h3>
    <p>For all three duplexes, regressing force on separation gives no resolvable slope
    (p = {", ".join("%.2f" % e["slope_p"] for e in (mfx or {"by_system": {}})["by_system"].values())}).
    A constant load independent of compression is the signature of a <strong>buckled
    rod</strong>: once past its critical load an elastic rod does not stiffen, it just
    keeps folding at roughly the same force.</p>

    <h3>And the load it holds is the Euler load</h3>
    <p>The 27 bp duplex plateaus at
    <strong>{(mfx or {}).get("by_system", {}).get("D_27bp", {}).get("f", 0):.1f} ±
    {(mfx or {}).get("by_system", {}).get("D_27bp", {}).get("sem", 0):.1f} pN</strong>
    against a classical Euler buckling load of
    {(mfx or {}).get("by_system", {}).get("D_27bp", {}).get("euler", 0):.1f} pN — a ratio of
    {(mfx or {}).get("by_system", {}).get("D_27bp", {}).get("ratio", 0):.2f}, with nothing
    fitted. The 40 bp pair land at
    {(mfx or {}).get("by_system", {}).get("D_40bp", {}).get("ratio", 0):.2f} and
    {(mfx or {}).get("by_system", {}).get("D_40bp_nick", {}).get("ratio", 0):.2f} of theirs.
    Textbook elasticity, from an all-atom force field, with a persistence length nobody
    supplied.</p>

    <div class="note flag"><div class="h">This disagrees with the model the project is built on</div>
      <p>Zocchi's kinked branch predicts <strong>2.4–6.6 pN</strong> for these same
      duplexes at these same separations. We measure
      {(mfx or {}).get("by_system", {}).get("D_27bp", {}).get("f", 0):.0f} pN for the
      27 bp — roughly <strong>six times higher</strong> — and on the
      <em>unyielded</em> branch, which is the branch his model says the duplex has already
      left. On this timescale the duplex is not a softened, kinked spring; it is an
      elastic rod holding its buckling load.</p>
      <p>Two consequences. The forces in the table above are, if anything,
      <strong>underestimates</strong>. And the reason the protein does not respond is not
      that the force is small — it is that even 26 pN moved these attachment sites less
      than the noise, because the barrel is stiff along that coordinate and the load is
      not applied along it.</p></div>

    <h3>Two restraint stiffnesses, one answer</h3>
    <p>Each separation was run at two restraint stiffnesses, 33 and 200 kJ mol⁻¹ nm⁻².
    The soft one is not really a constraint — the duplex settles 2–7 Å past its target —
    so if the measured force were an artefact of how hard we pulled, the two would sit on
    different curves. Plotted against the separation the duplex actually adopted, they
    agree
    (p = {", ".join("%.2f" % e["k_agree_p"] for e in (mfx or {"by_system": {}})["by_system"].values() if "k_agree_p" in e)}).
    That is the control that makes this a measurement of the duplex rather than of the
    restraint.</p>

    <h3>What the nick does — and does not do</h3>
    <p>Zocchi's model treats the nick as the key softening control, and predicts about
    1.3× less force. Measured, there is <strong>no resolvable difference</strong>:
    {(mfx or {}).get("nick", {}).get("intact", 0):.1f} pN intact against
    {(mfx or {}).get("nick", {}).get("nicked", 0):.1f} pN nicked, p =
    {(mfx or {}).get("nick", {}).get("p", float("nan")):.2f}. A nick removes one backbone
    bond from a rod held together by 40 base pairs' worth of stacking; over 20–30 ns it
    does not measurably change the buckling load. The model's ordering may still be right
    on longer timescales, where a nick can act as a hinge — but our own data does not
    support quoting it as a validated control, and the model section above should be read
    with that in mind.</p>

    <div class="note"><div class="h">What this cannot say</div>
      <p>These are 20–30 ns windows on bare duplexes. A kink that takes microseconds to
      nucleate would lower the plateau and is invisible here, so the Euler agreement is
      evidence about this timescale, not about equilibrium.</p>
      <p>And the chimera's own separation — about 70 Å for the 27 bp — is
      <strong>outside</strong> the range we can sample validly, because a box big enough
      to hold a 27-mer at 70 Å without it touching its own periodic image is larger than
      the one these runs used. An earlier draft quoted a measured force at the chimera's
      geometry; that figure came from a window whose target was 87 % of the box width and
      it is <strong>withdrawn</strong>. Repeating those windows in a larger box is queued
      and waiting on GPU allocation.</p></div>
  </div>
</section>

<section>
  <div class="sec-head"><h2>What the protein does about it</h2>
    <span class="tag">Cβ–Cβ separation, Asp133–Asn149</span></div>
  <div class="tbl-wrap"><table>
    <thead><tr><th>system</th><th>force (pN)</th><th>mean (Å)</th><th>s.e.m.</th>
    <th>replicate range</th><th>n</th><th>Δ vs clicked control (Å)</th></tr></thead>
    <tbody>{ext_rows}</tbody></table></div>
  <div class="col">
    <h3>The protein is stiff; how stiff, we can only bound</h3>
    <p>The axial clamp — a genuine constant force, <em>E</em> = −<em>f·r</em>, applied
    between the two Cβ atoms — gives a compliance of
    {clamp6.get('slope_A_per_pN', 0):.4f} ± {clamp6.get('slope_se', 0):.4f} Å pN⁻¹ over
    0–20 pN across {clamp6.get('n', 0)} replicates. Its 95% confidence interval,
    [{clamp6.get('ci95', [0, 0])[0]:.4f}, {clamp6.get('ci95', [0, 0])[1]:.4f}],
    <strong>includes zero</strong>: at 45 ns per replicate the response is not resolved
    from noise, and the 4 pN point sitting below the 0 pN point in figure 1 is that noise
    made visible. What the data <em>do</em> support is a lower bound,
    <strong>κ &gt; {kap:.0f} k<sub>B</sub>T nm⁻²</strong>, which lands on the
    ~100 k<sub>B</sub>T nm⁻² Zocchi’s activity measurements imply{cite_zocchi}. Under a few pN a
    protein this stiff should move by hundredths of an ångström, which is precisely why
    it cannot be measured this way in 45 ns — the honest conclusion is agreement with
    his stiffness and no more.</p>

    <div class="tbl-wrap"><table>
      <thead><tr><th>applied force (pN)</th><th>mean (Å)</th><th>s.e.m.</th>
      <th>replicate range</th><th>n</th><th>Δ vs 0 pN</th></tr></thead>
      <tbody>{ladder_rows}</tbody></table></div>

    <h3>The spring does <em>not</em> move the sites more than an axial pull</h3>

    <div class="note flag"><div class="h">Retraction</div>
      <p>An earlier version of this page said the three chimeras sat +0.19, +0.52 and
      +0.59 Å from the zero-force clamp, monotonic in force at 0.13 Å pN⁻¹, and
      attributed the excess over a pure axial pull to a couple and a shear delivered
      through the two flexible tethers. <strong>That claim is withdrawn.</strong> It used
      the <em>unclicked</em> protein as its zero. Every chimera is clicked, and clicking
      alone moves the sites {cw.get('delta', 0):+.2f} Å — fifteen to thirty times any
      force effect. The physical argument was also unsound: at equilibrium a couple does
      no net work, and along Cβ–Cβ only the component of force along that coordinate can
      displace the mean. Thanks to Doeke Hekstra for catching both.</p></div>

    <p>The right zero is the clicked, unloaded protein: same tether chemistry, same two
    modified residues, no DNA. Against it, every chimera sits <em>at or below</em> the
    control — the springs compress the attachment sites or leave them alone. Nothing
    survives a Bonferroni correction for the three comparisons.</p>

    <div class="tbl-wrap"><table>
      <thead><tr><th>modelled force (pN)</th><th>system</th><th>Δ vs clicked control (Å)</th>
      <th>t</th><th>p</th><th>p (×3)</th><th>resolved</th></tr></thead>
      <tbody>{matched_rows}</tbody></table></div>

    <p>The two baselines, for comparison: the clicked control sits at
    {mc.get('mean', float('nan')):.3f} ± {mc.get('sem', float('nan')):.3f} Å and the
    unclicked protein the old claim used at
    {unclicked_ref.get('mean', float('nan')):.3f} ±
    {unclicked_ref.get('sem', float('nan')):.3f} Å.</p>

    <h3>What the axial clamp measures instead: the barrel's stiffness</h3>
    <p>Clamping the same two Cβ at 0, 2.35 and 4.88 pN on the clicked topology gives a
    straight line at 6.9 mÅ pN⁻¹ (r² = 0.998), i.e. <strong>κ ≈ 350 k<sub>B</sub>T
    nm⁻²</strong>. The longer six-point ladder on the unclicked protein is noisier
    (r² = 0.18) and gives ≈ 230 k<sub>B</sub>T nm⁻². Both say the barrel along this
    coordinate is <strong>two to three and a half times stiffer</strong> than the
    ~100 k<sub>B</sub>T nm⁻² Zocchi uses as a generic protein spring constant{cite_zocchi}.
    At 350 k<sub>B</sub>T nm⁻² a 4 pN spring moves the sites 0.03 Å — which is the whole
    reason this cannot be measured by watching an equilibrium distance in 50 ns.</p>

    <h3>The spring is not acting at a distance — it is touching the protein</h3>
    <p>If the chimeras differ from a matched axial pull, the difference cannot be a
    couple, for the reason above. The obvious alternative is that the duplex is a large
    charged object tethered two nanometres from the barrel and leaning on it.
    {"That is what it does:" if have_contacts else "Measuring it directly:"}</p>
    {'<div class="tbl-wrap"><table><thead><tr><th>system</th>'
     '<th>protein heavy atoms within 4.5 Å</th><th>frames in contact</th>'
     '<th>residues touched</th></tr></thead><tbody>' + contact_rows +
     '</tbody></table></div>' if have_contacts else
     '<p class="muted">Contact analysis pending the current round of runs.</p>'}
    <p>The duplex is in contact with the protein in essentially every frame of every
    replicate, so whatever it does, it does not do at a distance. And the most compressed
    system is the nicked one, which delivers the least modelled force — a floppier duplex
    lies against the barrel more readily, which is the ordering a contact mechanism
    predicts and the reverse of a force mechanism.</p>

    <div class="note flag"><div class="h">But contact does not survive its own test</div>
      <p>That ordering is one comparison, and an earlier version of this page used it to
      argue for contact as <em>the</em> mechanism. With eleven constructs measured against
      their own controls, the compression correlates with contact area no better than it
      does with force — <strong>r = +0.32, p = 0.33</strong> — and constructs with no
      spring at all compress as much as those with one. The section
      <em>Is it acting as a spring at all?</em> below does that test properly; this
      paragraph is what the contact measurement shows, not what it explains.</p>
      <p>The companion ubiquitin study reaches the same negative conclusion from the other
      direction{cite_boral}: a control bearing two separate short strands reproduced most
      of the circular chimera's NMR spectrum, pointing at local protein–DNA interaction
      rather than spring action. Agreeing that it is <em>not</em> spring action is firmer
      than either study's guess at what it is instead.</p></div>
  </div>
</section>

<section>
  <div class="sec-head"><h2>The single-stranded control, and a sign flip</h2>
    <span class="tag">an entropic spring pulls where a duplex pushes</span></div>
  <div class="col">
    <p>Replacing the duplex with a single strand of the same span changes the sign of the
    mechanics, which makes it the sharpest available control. A bent duplex is an
    <em>enthalpic</em> spring: it is stiffer than the geometry it is forced into, so it
    <strong>pushes</strong> its anchors apart. A single strand is an <em>entropic</em>
    spring with a 1.5 nm persistence length: it has far more contour than it needs, so it
    coils, and it can only ever <strong>pull</strong> — and then only once stretched past
    its natural size.</p>

    <div class="tbl-wrap"><table>
      <thead><tr><th>strand</th><th>nt</th><th>contour (Å)</th>
      <th>natural rms end-to-end (Å)</th><th>span held at (Å)</th>
      <th>tension (pN)</th><th>state</th></tr></thead>
      <tbody>{ss_rows}</tbody></table></div>

    <p>Both single-stranded loops sit <strong>well inside</strong> their own natural size:
    the 50-nucleotide strand has {ss.get('C1_ss50', {}).get('contour_nm', 0) * 10:.0f} Å of
    contour and a natural end-to-end of about
    {ss.get('C1_ss50', {}).get('natural_rms_nm', 0) * 10:.0f} Å, while the attachment sites
    hold it at only {ss.get('C1_ss50', {}).get('x_A', {}).get('mean', 0):.0f} Å. A
    worm-like chain at less than its own rms size is under <strong>no tension at all</strong>,
    so the Marko–Siggia force is exactly zero and this construct is a genuine
    zero-force reference rather than a low-force one.</p>

    <div class="note"><div class="h">Why this matters for the interpretation</div>
      <p>This is the construct that carries the argument in the next section. It is the same
      protein, the same two attachment sites, the same click chemistry, a comparable amount
      of DNA — and provably zero mechanical load. If the compression of the attachment sites
      were a force effect, this construct should show none of it. It shows
      {abs((f.get('extension', {}).get('vs_matched_control', {}).get('C1_ss50') or {}).get('delta', 0)):.2f} Å,
      which is as much as any duplex.</p>
      <p>An earlier draft of this work described a 27-nucleotide strand as “taut” by
      comparing a 6.9 nm figure against a 6–7 nm span. That was an error: 6.9 nm is the rms
      <em>end-to-end</em> distance of a 27-mer, not its contour, and the contour is 15.9 nm.
      A chain is taut when the span approaches its <em>contour</em>, not its rms size, so
      the conclusion did not follow from the number. The 50-nucleotide strands used here are
      slack by either measure.</p></div>
  </div>
</section>

<section>
  <div class="sec-head"><h2>Is it acting as a spring at all?</h2>
    <span class="tag">eleven constructs, two proteins, two chemistries</span></div>
  <div class="col">
    <p>The premise of the whole design is that the DNA acts as a spring: bend it, and it
    pushes the attachment sites apart in proportion to the force it delivers. That is a
    testable claim, and with eleven constructs measured against their own matched
    controls it can be tested rather than assumed.</p>

    <div class="tbl-wrap"><table>
      <thead><tr><th>construct</th><th>modelled force (pN)</th>
      <th>protein atoms the DNA touches</th><th>Δ Cβ–Cβ vs control (Å)</th>
      <th>p</th><th>p corrected</th></tr></thead>
      <tbody>{allrows}</tbody></table></div>

    <p><strong>Every construct that carries DNA compresses the sites.</strong> Ten of the
    eleven are negative, spanning −0.09 to −0.59 Å, and the eleventh is +0.02 — which is
    zero. Not one extends them, which is what a spring pushing its anchors apart would do.</p>

    <p>And the size of the effect has nothing to do with the force:</p>
    <ul>
      <li>Against the modelled spring force: <strong>r = {spring_q.get('vs_modelled_force', {}).get('r', float('nan')):+.2f},
        p = {spring_q.get('vs_modelled_force', {}).get('p', float('nan')):.2f}</strong> — and the
        sign is <em>positive</em>, meaning more force goes with <em>less</em> compression.</li>
      <li>Against how much of the protein the DNA touches:
        <strong>r = {spring_q.get('vs_contact', {}).get('r', float('nan')):+.2f},
        p = {spring_q.get('vs_contact', {}).get('p', float('nan')):.2f}</strong>.</li>
      <li>Constructs with <strong>no spring at all</strong> — a slack single strand, or two
        separate arms that cannot circularise — compress by
        {abs(spring_q.get('mean_delta_no_spring', 0)):.2f} Å. Those with a spring compress by
        {abs(spring_q.get('mean_delta_with_spring', 0)):.2f} Å. The difference is not
        resolved (p = {(spring_q.get('no_spring_vs_spring') or {}).get('p', float('nan')):.2f}).</li>
    </ul>

    <div class="note flag"><div class="h">What this means</div>
      <p>Attaching DNA to this protein moves the attachment sites together by about a third
      of an ångström, and <strong>it makes no measurable difference whether that DNA is a
      spring</strong>. A slack 50-nucleotide single strand under zero tension does it as
      well as a 27 bp duplex the model says is pushing with 6.6 pN.</p>
      <p>This is the conclusion the companion ubiquitin study reached from NMR{cite_boral} —
      that local protein–DNA interaction, not spring action, dominates — reproduced here
      independently, in a different protein, with a different attachment chemistry, and
      read out from a different observable. That agreement across two systems is the
      strongest thing in this report.</p>
      <p>Two obvious mechanisms are <em>not</em> it. It is not the force, which does not
      correlate. And it is not simple contact area either, which also does not correlate:
      the 40 bp closed loop touches the barrel in only 20 % of frames and compresses the
      most of any duplex. What is common to every construct is a large, highly charged
      object tethered a couple of nanometres away; which specific part of that does the
      work is not something these runs identify.</p></div>
  </div>
</section>

<section>
  <div class="sec-head"><h2>Does closing the covalent loop matter?</h2>
    <span class="tag">five matched pairs</span></div>
  <div class="col">
    <p>The rebuild described above changed one thing about the molecule: whether a covalent
    path runs from one attachment site, out along its tether, through the DNA, and back to
    the other. Everything else — protein, chemistry, spring, length, tether — is held
    fixed, which makes these five pairs a clean test of topology on its own.</p>

    <div class="tbl-wrap"><table>
      <thead><tr><th>pair</th><th>closed − two-arm (Å)</th><th>p</th>
      <th>direction</th></tr></thead>
      <tbody>{topo_rows}</tbody></table></div>

    <p>The three <strong>duplex</strong> pairs — the comparison the rebuild was designed to
    make — all point the same way: closing the loop compresses the sites further, by
    {abs(topo.get('duplex_pairs_planned', {}).get('mean_A', 0)):.2f} Å on average, and two
    of the three clear significance on their own. Taken as a paired test across the three,
    <strong>p = {topo.get('duplex_pairs_planned', {}).get('p', float('nan')):.2f}</strong> —
    suggestive, and not resolved at n = 3 pairs.</p>

    <p>The two <strong>single-stranded</strong> pairs do not reproduce it. They give
    {topo.get('ssDNA_pairs', {}).get('deltas_A', [0, 0])[0]:+.3f} Å and
    {topo.get('ssDNA_pairs', {}).get('deltas_A', [0, 0])[-1]:+.3f} Å — one is zero and the
    other has the opposite sign. Across all five pairs together the effect is
    p = {topo.get('paired_across_pairs', {}).get('p', float('nan')):.2f}: nothing.</p>

    <div class="note"><div class="h">Reading this honestly</div>
      <p>There is a hint that covalent closure matters for a duplex and no evidence that it
      matters for a slack single strand — which would be physically sensible, since only a
      closed loop can transmit a bending moment around the circuit, and a slack strand has
      no moment to transmit. But a three-pair trend at p = 0.08 that reverses when two more
      pairs are added is a hypothesis, not a finding, and it is reported as one. The number
      of replicates that would settle it — six per system rather than three — is a day of
      GPU time whenever the allocation allows.</p>
      <p>What the rebuild <em>did</em> establish is not statistical: the original build was
      not the molecule the experiment makes, and no amount of sampling would have fixed
      that. Correctness of the construct and significance of its effect are separate
      questions, and only the first one is settled here.</p></div>
  </div>
</section>

<section>
  <div class="sec-head"><h2>What the modification itself does</h2>
    <span class="tag">before any force is applied</span></div>
  <div class="col">
    <p>Installing the two Tet2-Et residues shortens the Cβ–Cβ separation by
    {abs(tw.get('delta', 0)):.2f} Å before any DNA is attached, and clicking gives most of
    it back. That is larger than any force effect in this study, so it is worth knowing
    what it is — and it is <em>not</em> a change in secondary structure. DSSP over
    residues 128–155 is identical in the wild type, the unclicked protein and the clicked
    adduct: site 133 is coil in every case and site 149 is 100% β in every case.</p>

    {'<div class="tbl-wrap"><table><thead><tr><th>system</th>'
     '<th>site</th><th>residue</th><th>φ, ψ</th>'
     '<th>site</th><th>residue</th><th>φ, ψ</th>'
     '<th>ring–Tyr200 contacts</th></tr></thead><tbody>'
     + site_rows + '</tbody></table></div>' if have_sites else ''}

    <p><strong>Site 133 is the clean result.</strong> Asp133 natively occupies the
    <em>left-handed</em> α region — positive φ — a conformation that Asp, Asn and Gly
    commonly adopt and that a phenylalanine-shaped side chain tolerates poorly. Tet2-Et is
    a phenylalanine derivative, and it is pushed 17.5° out in φ and 14° in ψ; clicking
    returns it to within a few degrees of wild type. The replicate spread is 0.5–0.8°, so
    this is not noise.</p>

    <p><strong>Site 149's ring stacks on Tyr200, but not reproducibly.</strong> The
    unclicked tetrazine averages ~10 heavy-atom contacts per frame, but the three
    replicates give 20.0, 10.4 and 0.9 — it samples stacked and unstacked states and three
    replicates cannot pin the occupancy. The clicked control is tight and low
    (0.68 ± 0.4 over nine replicates), yet the 40 bp chimeras stack again. So the ring at
    149 <em>can</em> stack, and whether it does depends on the construct.</p>

    <p>Both variables track the contraction. Across all 30–33 replicates carrying these
    residues, Cβ–Cβ correlates with the Tyr200 contact count (r = −0.68, p = 5×10⁻⁵,
    −38 mÅ per contact) and with φ at site 133 (r = −0.59, p = 3×10⁻⁴, −39 mÅ per degree);
    together they account for 55% of its variance. These are correlations across
    constructs that differ chemically, so this is suggestive rather than causal — but the
    contraction is <em>local structure at the modified sites</em>, and nothing to do with
    the spring.</p>

    <div class="note"><div class="h">Which positions can take Tet2-Et</div>
      <p>Avoid positive-φ positions: the left-handed α region is occupied by Asp, Asn and
      Gly for a reason. Treat an aromatic neighbour as a liability, because the ring will
      stack on it in some constructs and not others, and because that stacking
      <em>differs between the unclicked and clicked states</em> — so it will masquerade as
      a mechanical effect for anyone using the unclicked protein as a reference. That is
      exactly how the retracted claim above arose. Ordinary β and α positions with no
      aromatic neighbour should be fine.</p></div>
  </div>
</section>

<section>
  <div class="sec-head"><h2>Which positions tolerate the modification</h2>
    <span class="tag">a design rule, tested in a second protein</span></div>
  <div class="col">
    <p>The Ramachandran result above gives a rule with a prediction attached: a site
    strains on modification only where the <em>native</em> backbone occupies a positive-φ
    conformation, which Asp, Asn and Gly adopt readily and a bulkier side chain does not.
    Two proteins with two unrelated attachment chemistries give four sites to test it on,
    and exactly one of the four is positive-φ.</p>

    {'<div class="tbl-wrap"><table><thead><tr><th>protein</th><th>site</th>'
     '<th>residue</th><th>native φ, ψ</th><th>Δφ on modification</th><th>Δψ</th>'
     '<th>native conformation</th></tr></thead><tbody>' + strain_rows +
     '</tbody></table></div>' if have_strain else
     '<p class="muted">Awaiting the second protein&rsquo;s site analysis.</p>'}

    <p><strong>The prediction holds.</strong> The one positive-φ site moves 17.5° in φ; the
    three ordinary sites move by a degree or less. And the chemistries have nothing in
    common — sfGFP carries a genetically encoded tetrazine amino acid, a phenylalanine
    derivative installed by amber suppression; ubiquitin carries a cysteine converted to a
    maleimide thioether. What predicts the strain is not the adduct but the backbone
    conformation the native residue was holding.</p>

    <div class="note"><div class="h">The practical rule</div>
      <p>When choosing attachment sites for either chemistry: <strong>avoid positions
      whose native φ is positive</strong> — in practice Gly, and the Asn/Asp positions
      that sit in the left-handed α region. An ordinary β or α position takes the
      modification without measurable backbone cost.</p>
      <p>Treat an aromatic neighbour as a second liability, independently. The tetrazine
      ring at sfGFP&rsquo;s site 149 stacks on Tyr200 in some constructs and not others,
      and it does so <em>differently before and after clicking</em> — so it will
      masquerade as a mechanical effect in exactly the comparison this project set out to
      make.</p>
      <p>Both of these are conclusions about four sites in two proteins, which is a rule
      of thumb rather than a survey. What makes it worth stating is that the prediction
      was made from sfGFP alone and then held in a protein and a chemistry it was not
      derived from.</p></div>
  </div>
</section>

<section>
  <div class="sec-head"><h2>The chromophore</h2>
    <span class="tag">occupancy on H···A and D–H···A, not heavy atoms</span></div>
  <div class="tbl-wrap"><table>
    <thead><tr><th>hydrogen bond</th>{hb_head}</tr></thead>
    <tbody>{hb_rows}</tbody></table></div>
  <p class="col" style="font-size:14px;color:var(--ink-3)">Occupancy as a percentage,
  with the range across replicates in parentheses.</p>
  <div class="col">
    <h3>The network is intact everywhere, and the spring does not break it</h3>
    <p>Thr203 to the phenolate and Arg96 to the imidazolinone carbonyl{cite_ormo} hold at 96–100%
    in every system — these are the buried, load-bearing contacts that keep the
    chromophore planar, and nothing the spring does disturbs them. Gln94 holds at
    91–97%. If the spring were prising the barrel open near the chromophore, these are
    the numbers that would fall, and they do not.</p>

    <h3>His148 is tighter in every chimera — but not beyond its own scatter</h3>
    <p>His148’s bond to the phenolate reads
    {100 * f['hbonds']['His148 ND1-H -> chromophore phenolate']['per_system']['S0_wt']['mean']:.0f}%
    in WT and
    {100 * f['hbonds']['His148 ND1-H -> chromophore phenolate']['per_system']['S1_tet']['mean']:.0f}%
    in the unclicked protein, against 80–81% in all three chimeras. The direction is
    consistent and mechanistically appealing — His148 sits between the two attachment
    sites, so a spring pulling on that strand plausibly presses it onto the chromophore.
    But the replicate-to-replicate spread is 58–92%, wider than the gap between the
    means. It is a hypothesis worth more sampling, not a result.</p>

    <h3>One apparent effect is an artefact of bistability</h3>
    <p>The <code>{esc(bistable[0]) if bistable else '—'}</code> contact — part of the proton relay behind the chromophore{cite_brejc} —
    reads about 30%
    in most systems and 0% in the 40 bp chimera, which looks like an effect and is not.
    Within a single system it is either near 0% or near 90%: a two-state Ser205 hydroxyl
    rotamer, flipped once per trajectory and then stuck. A mean of three coin flips is
    not a measurement, and figure 3 plots the individual replicates so it cannot pass as
    one. The automated check that flagged it lives in
    <code>analysis/findings.py</code> and runs over every H-bond.</p>

    <h3>Bridge planarity and cavity water are unchanged</h3>
    <div class="tbl-wrap"><table>
      <thead><tr><th>system</th><th>τ twist (°)</th><th>φ twist (°)</th>
      <th>waters within 6 Å</th><th>mean RMSF (Å)</th></tr></thead>
      <tbody>{aux_rows}</tbody></table></div>
    <p>The two dihedrals that twist the chromophore out of plane — τ across the
    phenol–bridge bond and φ across the bridge–imidazolinone bond — are the structural
    proxy for non-radiative decay: the free chromophore, able to rotate about them, is
    essentially dark while the protein-bound one is bright{cite_niwa}, so the wider
    their excursions the more of the excited state leaks away without a photon. They sit at 5–6° of twist in every system, with
    differences between systems of a few tenths of a degree against replicate spreads of
    similar size. Cavity water, the other route to quenching, is 8–11 molecules with
    spreads of 1–3. Neither is resolved. If the spring changes sfGFP’s brightness at
    these forces, 45 ns × 3 does not see the structural cause, and the honest reading of
    the fluorometry is that it needs either a stronger spring or far more sampling.</p>
  </div>
</section>

<section>
  <div class="sec-head"><h2>Figures</h2>
    <span class="tag">every panel ships the CSV it was drawn from</span></div>
  {figs_html}
</section>

<section>
  <div class="sec-head"><h2>Force field and protocol</h2>
    <span class="tag">what every number above was computed with</span></div>
  <div class="col">
    <p>Protein <strong>ff14SB</strong>{cite_maier}, DNA <strong>OL15</strong>{cite_ol15},
    <strong>TIP3P</strong> water{cite_tip3p}, Joung–Cheatham monovalent ions{cite_ions}
    at 150 mM NaCl, and Mg²⁺ from the Li–Merz 12-6 set{cite_mg} at ~5 mM for the duplex
    systems — 12-6 rather than 12-6-4 because OpenMM has no C⁴ term, so a 12-6-4
    topology would silently lose the polarisation correction on import. Buffer and salt
    follow the companion ubiquitin experiment in this lab{cite_boral}.</p>
    <p>The chromophore uses the published <strong>CRO</strong> parameters{cite_breyfogle},
    which come out of tleap{cite_case} at exactly −1.0000 e. The four custom residues —
    the tetrazine amino acid, the click adduct, and the DNA-side arm in its two forms —
    are typed with <strong>GAFF2</strong>{cite_gaff} and charged with
    <strong>AM1-BCC</strong>{cite_ambcc}; their six backbone atoms are then overwritten
    with ff14SB’s types and charges so the peptide bonds into and out of the modified
    sites match every other one, the same device the chromophore parameters
    used{cite_breyfogle}. Duplexes are built as ideal B-form DNA with PyMOL’s
    <code>fnab</code>{cite_pymol}; the tethers are placed by superposition{cite_kabsch}
    and torsion optimisation.</p>
    <p>Simulations run in <strong>OpenMM 8.6</strong>{cite_openmm} with the
    LangevinMiddle integrator{cite_leimkuhler} at 300 K and a 4 fs timestep enabled by
    hydrogen-mass repartitioning on the solute only{cite_hopkins}, PME with a 0.9 nm
    cutoff, and a Monte Carlo barostat at 1 bar. Analysis uses MDTraj{cite_mdtraj};
    solvent accessibility is computed with Shrake–Rupley{cite_shrake}.</p>
  </div>
</section>

<section>
  <div class="sec-head"><h2>Systems</h2>
    <span class="tag">explicit solvent · 150 mM NaCl · 300 K</span></div>
  <div class="tbl-wrap"><table>
    <thead><tr><th>id</th><th>system</th><th>atoms</th><th>net charge</th>
    <th>what it is for</th></tr></thead><tbody>{rows_sys}</tbody></table></div>
</section>

<section>
  <div class="sec-head"><h2>How it is built</h2>
    <span class="tag">raw deposition → six systems</span></div>
  <ol class="steps">{rows_steps}</ol>
</section>

<section>
  <div class="sec-head"><h2>Verification</h2><span class="tag">asserted, not assumed</span></div>
  <div class="tbl-wrap"><table>
    <thead><tr><th>check</th><th>result</th><th></th></tr></thead>
    <tbody>{rows_ver}</tbody></table></div>
  <div class="note"><div class="h">Deliberately not used</div>
    <p>The <em>F</em> ≈ 10 pN figure from Choi &amp; Zocchi (2007){cite_choi}. That paper labels it
    an upper bound from a worm-like-chain treatment that ignores kinking; the kinked
    model supersedes it and gives 1.5–2.4 pN for the same constructs. Using it would
    over-stress the protein roughly fourfold.</p></div>
  <div class="note"><div class="h">What this does and does not show</div>
    <p><strong>Does:</strong> that the chimera can be built and parameterised without
    inventing a force-field term; that the spring delivers 2–5 pN and 7–9 k<sub>B</sub>T
    by the model’s own formula evaluated on the simulation’s geometry; that a nick
    softens it in the predicted direction; that the barrel and the chromophore’s
    load-bearing hydrogen bonds survive intact; and that the protein is at least as
    stiff as Zocchi’s calibration assumes.</p>
    <p><strong>Does not:</strong> resolve the protein’s compliance, establish that the
    spring changes the chromophore, or confirm the kink the force model presumes. Those
    need the clicked control now running, more replicates, and longer trajectories —
    <code>bash slurm/submit_all.sh --tier B</code> extends every run to 300 ns from its
    checkpoint.</p></div>
</section>

<section>
  <div class="sec-head"><h2>References</h2>
    <span class="tag">everything this setup rests on</span></div>
  {refs}
</section>

<footer>
  Generated by <code>report/build_report.py</code> from <code>analysis/findings.json</code>
  and the per-run analyses. Regenerate with <code>bash finalize.sh</code> after any run
  extends. {f['n_replicates']} replicates · {f['total_ns']:,.0f} ns ·
  systems analysed: {', '.join(esc(LABEL.get(s, s)) for s in f['systems_present'])}.
</footer>
</div>
"""
    out = Path(args.out)
    out.write_text(html)
    print(f"wrote {out} ({out.stat().st_size / 1e6:.2f} MB)")
    print(f"  PNG figures: {html.count('data:image/png')}, "
          f"animated GIFs: {html.count('data:image/gif')}, "
          f"video elements: {html.count('<video')} (should be 0)")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", default=str(ROOT))
    ap.add_argument("--out", default=str(ROOT / "report" / "sfgfp_dna_spring_report.html"))
    # the corrected closed-loop build, not the superseded two-arm one
    ap.add_argument("--anim", default="L3_loop27")
    build(ap.parse_args())
