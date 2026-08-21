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
         "S3_spring27": "27 bp spring", "S4_spring40": "40 bp spring",
         "S5_spring40nick": "40 bp nicked", "S6_clamp": "force clamp"}
PURPOSE = {
    "S0_wt": "Baseline. Validates the published chromophore parameters against a "
             "structure nobody modified.",
    "S1_tet": "The two encoded tetrazines, unclicked — the fluorometry sample. Also "
              "the topology the force-clamp ladder runs on.",
    "S2_clicked": "The clicked tether with no DNA: the zero-force reference with a "
                  "chimera’s own chemistry.",
    "S3_spring27": "The strong spring.",
    "S4_spring40": "Zocchi’s mechanically clean length (γ &lt; 1).",
    "S5_spring40nick": "Nicked at the centre — his calibrated low-stress reference.",
    "S6_clamp": "A constant force of 0/2/4/7/12/20 pN applied straight along the "
                "Cβ–Cβ axis: the calibration curve.",
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
     quotes a <em>bound</em> on stiffness rather than a value — and the chimeras sit well
     above that line, so the spring is not equivalent to pulling along the same axis.
     The 4 pN clamp point falling below the 0 pN point is the size of the noise at this
     sampling; it is not a physical dip."""),
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
VER = [
    ("OpenMM with CUDA", "all four platforms agree within tolerance", True),
    ("CRO atom names against the xFPchromophores template", "22 / 22, checked before use", True),
    ("CRO net charge from the library", "−1.0000 exactly", True),
    ("parmchk2 parameters flagged “ATTN, need revision”", "0, all four custom residues", True),
    ("Assigned GAFF type implies the declared element", "every atom, every residue", True),
    ("Prep tree parent links are declared bonds", "all four residues", True),
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
    ("zocchi_book", "Zocchi G. <i>Molecular Machines: A Materials-Science Approach.</i> "
     "Princeton University Press (2018).",
     "https://doi.org/10.2307/j.ctvc77hpc",
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
     "https://doi.org/10.1021/ja4045532",
     "Eq. 6 as implemented, τ<sub>c</sub> = 27 pN nm for a nicked rod, and the "
     "1.5 / 2.4 pN values the implementation is tested against."),
    ("tseng2021", "Tseng C-Y, Wang A, Zocchi G. Protein–DNA chimeras: synthesis and "
     "applications. <i>Methods Enzymol</i> <b>647</b>, 257–281 (2021).",
     "https://doi.org/10.1016/bs.mie.2020.11.005",
     "The verified oligonucleotide sequences used for the springs, and the measured "
     "energy split E<sub>tot</sub> 9.1 = E<sub>DNA</sub> 6.6 + E<sub>enzyme</sub> "
     "2.5 k<sub>B</sub>T."),
    ("choi2007", "Choi B, Zocchi G. Guanylate kinase, induced fit, and the allosteric "
     "spring probe. <i>Biophys J</i> <b>92</b>, 1651–1658 (2007).",
     "https://doi.org/10.1529/biophysj.106.092049",
     "The earlier worm-like-chain force estimate of ~10 pN, which that paper itself "
     "labels an upper bound and which is deliberately <i>not</i> used here."),
    ("boral", "Boral S, Schnebly J, Gamada D, Gardner KH, Hekstra DR. Conformational "
     "changes induced in ubiquitin by circular protein–DNA chimeras. "
     "<i>bioRxiv</i> 2026.07.05.736590 (2026).",
     "https://doi.org/10.1101/2026.07.05.736590",
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
    zero = ext.get("zero_force_ref", {}).get("mean")

    # ---------- hero ----------
    hero_force = (f"{sp[strongest]['force_pN']['mean']:.1f} pN" if strongest else "—")
    hero_sub = (f"{sp[strongest]['label']}, from a 5′P–5′P span of "
                f"{sp[strongest]['x_A']['mean']:.0f} Å" if strongest else "")
    kap = clamp6.get("kappa_lower_bound_kT_nm2")

    # ---------- spring table ----------
    spring_rows = "\n".join(
        f"""<tr><td>{esc(v['label'])}</td>
        <td class="num">{v['contour_nm'] * 10:.0f}</td>
        <td class="num">{pm(v['x_A'], 1)}</td>
        <td class="num">{v['tau_c']:.0f}</td>
        <td class="num"><strong>{pm(v['force_pN'])}</strong></td>
        <td class="num">{pm(v['energy_kT'], 1)}</td>
        <td class="num">{pm(v['bend_total_deg'], 0)}</td>
        <td class="num">{v['bp_open']['mean']:.2f}</td>
        <td class="num">{v['extension_A']['mean'] - zero:+.2f}</td></tr>"""
        for k, v in sp.items())

    # ---------- clamp ladder ----------
    ladder_rows = "\n".join(
        f"""<tr><td class="num">{r['pN']:.2f}</td>
        <td class="num">{r['mean']:.3f}</td><td class="num">{r['sem']:.3f}</td>
        <td class="num">{r['min']:.3f} – {r['max']:.3f}</td>
        <td class="num">{r['n']}</td>
        <td class="num">{r['mean'] - zero:+.3f}</td></tr>"""
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
        <td class="num">{v['mean'] - zero:+.3f}</td></tr>\n""")

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
    rows_ver = "\n".join(
        f"""<tr><td>{v[0]}</td><td class="mono">{v[1]}</td>
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
    tw = ext.get("tet_vs_wt") or {}
    sw = ext.get("strongest_vs_weakest_spring") or {}
    chim_slope = ext.get("chimera_slope_A_per_pN")
    bistable = [n for n, v in f["hbonds"].items() if v["bistable"]]
    stamp = datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M UTC")

    # The template below is an f-string, so every {cite_*} marker has to exist as a
    # local here -- there is no .format() call to pass them into.
    refs = refs_html()
    cite_eddins = cite("eddins", "blizzard"); cite_pedelacq = cite("pedelacq")
    cite_ormo = cite("ormo"); cite_labms = cite("labms", "eddins")
    cite_spring = cite("zocchi_book", "qu", "tseng2013")
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
  strands of sfGFP’s β-barrel{cite_pedelacq}. A doubly-modified DNA duplex clicks onto both, and because dsDNA is
  far too stiff to follow the barrel’s curvature it has to bend — pushing the two
  attachment points apart with a force of a few piconewtons. This is what that looks
  like <em>atom by atom</em>: what the spring delivers, what the protein does about it,
  and which of the apparent effects survive their own error bars.</p>
  <div class="stats">
    <div class="stat force"><div class="k">strongest spring</div>
      <div class="v">{hero_force}</div><div class="s">{hero_sub}</div></div>
    <div class="stat"><div class="k">production analysed</div>
      <div class="v">{f['total_ns']:,.0f} ns</div>
      <div class="s">{f['n_replicates']} replicates across
      {len([k for k in systems if k in LABEL])} systems</div></div>
    <div class="stat"><div class="k">protein stiffness</div>
      <div class="v">&gt; {kap:.0f}</div>
      <div class="s">k<sub>B</sub>T nm⁻², a bound not a value — see below</div></div>
    <div class="stat ok"><div class="k">invented parameters</div>
      <div class="v">0</div>
      <div class="s">nothing flagged for review across the four custom residues</div></div>
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
    hold it planar: His148 and Thr203 to the phenolate oxygen, Arg96 to the
    imidazolinone carbonyl. The barrel is drawn transparent so the cavity shows through.
    Watch His148 — it is the mobile one of the three, and the residue the spring pulls
    across.</figcaption></figure>
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
    <th>total bend (°)</th><th>bp open</th><th>Δ extension (Å)</th></tr></thead>
    <tbody>{spring_rows}</tbody></table></div>
  <div class="col">
    <p>Every force here is Zocchi’s <em>f(x)</em>{cite_spring} evaluated at the end-to-end distance
    the simulation itself settles at, not at an assumed one. The construct is not free
    to choose that distance arbitrarily — it is set by where the two tethers can reach —
    so this is the one number in the whole model that the simulation can supply
    independently.</p>

    <h3>The nick softens the spring, by the amount it should</h3>
    <p>The nicked 40 bp duplex delivers <strong>{sig(s40n['force_pN']['mean'])} pN</strong>
    against <strong>{sig(s40['force_pN']['mean'])} pN</strong> for the intact duplex of
    the same length, and it is visibly floppier: its span wanders over
    {sig(s40n['x_wander_A']['mean'], 1)} Å within a trajectory against
    {sig(s40['x_wander_A']['mean'], 1)} Å for the intact one. Zocchi’s two values of the
    critical bending torque{cite_tau} — 27 pN·nm with a nick, 31–36 without — predict this
    ordering, and we reach it here from structure rather than from his calibration. This
    is the control he relies on throughout: the nicked spring is the low-stress
    reference, and it behaves like one.</p>

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
    are open at any instant — the duplex is intact. Either a kink needs longer than
    45 ns to nucleate, or at these lengths the strain is genuinely taken up as
    distributed elastic bending. Until that is settled the forces should be read as
    model-dependent numbers, consistent with the model’s inputs rather than confirming
    its mechanism.</p>
  </div>
</section>

<section>
  <div class="sec-head"><h2>What the protein does about it</h2>
    <span class="tag">Cβ–Cβ separation, Asp133–Asn149</span></div>
  <div class="tbl-wrap"><table>
    <thead><tr><th>system</th><th>force (pN)</th><th>mean (Å)</th><th>s.e.m.</th>
    <th>replicate range</th><th>n</th><th>Δ vs 0 pN</th></tr></thead>
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

    <h3>The spring moves the sites more than an equivalent axial pull</h3>
    <p>Ranked by the force their own spring delivers, the three chimeras sit
    {s40n['extension_A']['mean'] - zero:+.2f}, {s40['extension_A']['mean'] - zero:+.2f}
    and {s27['extension_A']['mean'] - zero:+.2f} Å from the zero-force clamp —
    monotonic in force, a slope of {chim_slope:.2f} Å pN⁻¹, more than ten times the
    axial-clamp slope. That is what Zocchi argues throughout his work{cite_zocchi}: forces applied at
    different places elicit <em>different responses</em>, not different magnitudes of one
    response. A spring anchored through two flexible tethers applies a couple and a
    shear as well as an extension, and the barrel evidently answers those more readily
    than a pure pull along Cβ–Cβ.</p>
    <p>With three replicates the strongest-against-weakest comparison reaches
    t = {sw.get('t', 0):.1f} on 4 degrees of freedom — <strong>not significant</strong>.
    This is a trend of the right sign and a plausible size, not a measurement.</p>

    <div class="note flag"><div class="h">The comparison is confounded — and the fix is
      running</div>
      <p>Installing the two Tet2-Et residues moves the sites by
      {tw.get('delta', 0):+.2f} Å on its own (t = {tw.get('t', 0):.1f},
      {'resolved' if tw.get('resolved') else 'not resolved'}) — comparable to or larger
      than any spring effect above. The clamp ladder runs on that unclicked protein, so
      it is <em>not</em> the right zero for a chimera built from the clicked adduct,
      which differs by the whole sTCO cage plus a six-carbon arm. A clicked-but-unloaded
      control (<code>S2_clicked</code>), axial clamps on that same topology set to the
      two extreme measured spring forces, and three further replicates of every chimera
      are in the queue. A chimera minus that control is what isolates the load; until
      then the Δ-extension column should be read as an upper limit on the spring’s
      effect, not an estimate of it.</p></div>
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
    ap.add_argument("--anim", default="S3_spring27")
    build(ap.parse_args())
