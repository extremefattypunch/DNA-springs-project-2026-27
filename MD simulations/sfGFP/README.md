# All-atom MD of the sfGFP–DNA spring chimera

An OpenMM pipeline that builds a Zocchi-style DNA-spring chimera on superfolder GFP
from the raw PDB deposition, runs it in explicit solvent, and measures what the spring
does to the barrel and the chromophore.

Two amber (TAG) codons put the tetrazine amino acid **Tet2-Et** at two solvent-exposed
sites; a doubly-**sTCO**-modified DNA duplex is clicked on, closing a covalent loop.
Because dsDNA's persistence length (~50 nm) dwarfs the barrel, the duplex must bend,
and it pushes the two attachment points apart with a force of a few pN.


## Quick start

```bash
source env/activate.sh          # sets MAMBA_*, $DNASPRING_PY, $DNASPRING_SCRATCH
./env/create_envs.sh            # one-time, ~1 h
./build/build_all.sh            # raw 2B3P -> five parameterised, solvated systems
bash slurm/submit_all.sh --tier A          # 11 configurations x 3 replicates x 50 ns
bash slurm/status.sh                       # queue + per-run progress
$DNASPRING_PY analysis/collect.py          # every analysis over every replicate
$DNASPRING_PY analysis/figures.py          # report figures + the CSVs behind them
$DNASPRING_PY analysis/animate.py --traj … --top … --name S3_spring27
bash slurm/submit_all.sh --tier B          # extend to 300 ns from the checkpoints
```

Trajectories go to `/n/netscratch/hekstra_lab/Lab/ian_poon/sfGFP-md` (symlinked as
`data/`, 90-day purge). Code, analysis CSVs and figures are versioned; trajectories
and checkpoints are not — this tree sits in a git-backed Obsidian vault with a public
GitHub remote.

## The systems

| ID | System | Atoms | What it is for |
|---|---|---|---|
| `S0_wt` | WT sfGFP | 44,016 | baseline; validates the chromophore parameters |
| `S1_tet` | sfGFP + 2× Tet2-Et (unclicked) | 44,011 | matches the 2-tet fluorometry sample; also the force-clamp topology |
| `S3_spring27` | 27 bp spring | 86,126 | the strong spring |
| `S4_spring40` | 40 bp spring | 121,404 | Zocchi's mechanically clean length (γ<1) |
| `S5_spring40nick` | 40 bp, nicked | 121,420 | his calibrated low-stress reference |
| `S6_clamp` | constant force 0/2/4/7/12/20 pN on `S1_tet` | 44,011 | the calibration curve |

The clamp ladder runs on `S1_tet`'s topology — the same two residues the spring pulls
on — so its force–response curve and the chimeras are directly comparable. The
chimeras' tension is read off their own geometry through Zocchi's formula; reading
their measured extension against the clamp curve is the independent second route to
the same number.

`--clamp-pN` applies a genuine constant force, `E = −f·r` between the two Cβ atoms, so
`F = −dE/dr = +f` pushes them apart — the sign a compressed leaf spring applies. A
moving harmonic restraint would instead impose a velocity and report a rate-dependent
force.

## Where the numbers come from

| Quantity | Value | Source |
|---|---|---|
| Structure | PDB **2B3P**, sfGFP, 1.40 Å | Pédelacq et al., *Nat Biotechnol* 24:79 (2006) |
| Chromophore | **CRO** (Thr65-Tyr66-Gly67), anionic | 2B3P; parameters from `leaprc.xFPchromophores` |
| Chromophore parameters | mixed GAFF/parm94, net charge −1.0000 | Breyfogle et al., *JPCB* 127:5772 (2023) |
| Attachment sites | **Asp133 / Asn149** in 2B3P numbering | see *Numbering*, below |
| Tet2-Et | 4-(6-ethyl-1,2,4,5-tetrazin-3-yl)-L-phenylalanine, C₁₃H₁₃N₅O, 255.27 Da | Eddins et al., *Bio-protocol* 14(16):e5048 (2024) |
| Click product | 4,5-dihydropyridazine fused to *trans*-bicyclo[6.1.0]nonane, N₂ lost | Eddins Fig. 1B; graphical overview p. 2 |
| Spring model | `f(x) = τ_c / [2R√(1−(x/2R)²)]` (kinked); `10B/L² − T/(L−x)` (bent) | Zocchi, *Molecular Machines* (2018), Eqs. 2.188–2.206, 3.51–3.52 |
| `B`, `τ_c`, rise | 200 pN·nm², 27 (nicked) / 31–36 (continuous) pN·nm, 0.33 nm/bp | Book pp. 73, 80; Tseng & Zocchi *JACS* 135:11879 (2013) p. 11884 |
| DNA sequence | centre of Zocchi's verified 2021 review 60mer | *Methods Enzymol* 647:257 (2021) p. 269 |
| Solvent | 150 mM NaCl, ~5 mM MgCl₂ (DNA systems), pH 7.4, 300 K | Boral et al. ubiquitin preprint pp. 16–18 |
| Force field | ff14SB + DNA.OL15 + TIP3P + Joung–Cheatham, Mg²⁺ Li/Merz 12-6 | — |

`analysis/test_spring_model.py` reproduces Zocchi's own published outputs from the
implementation: RLuc 60-mer 1.47 vs 1.5 pN, 40-mer 2.33 vs 2.4 pN, GK `E_DNA` 6.62 vs
6.6 k<sub>B</sub>T.

**Deliberately not used:** the *F* ≈ 10 pN / 25 k<sub>B</sub>T figure from Choi &
Zocchi, *Biophys J* 92:1651 (2007) p. 1657. That paper labels it "an upper bound" from
a worm-like-chain treatment that ignores kinking; the kinked model supersedes it and
gives 1.5–2.4 pN for the same constructs. Using it would over-stress the protein ~4×.

## Numbering — and a correction to the `.pse`

2B3P numbers the chain 1–246 with the Thr65-Tyr66-Gly67 chromophore collapsed into
`CRO 66`, so 65 and 67 are absent. The construct's own numbering runs **+1** relative
to this: construct D134/N150 are 2B3P **Asp133/Asn149**.

The session `pdb files/sfGFP 2b3p (150TAG, 134TAG).pse` bookmarks `134TAG` → Gly134
and `150TAG` → Val150, which is one residue off. The ESI-MS in the lab slides settles
it: Δ = **+141 Da** single (Asn→Tet2-Et) and **+281 Da** double (= 141 + 140, Asn +
Asp). Gly→Tet2-Et would be +198 and Val→Tet2-Et +156; neither is observed. Load
`build/01_protein/fix_pse_selections.pml` over the session for corrected selections.

Everything downstream reads `build/01_protein/attachment_sites.json`, so nothing
hard-codes a residue number or a distance.

## Design decisions, and their reasons

- **His148 is HID, not Amber's default HIE.** Its ND1 sits 2.98 Å from the chromophore
  phenolate (an acceptor, since CRO is modelled anionic) and its NE2 2.96 Å from the
  Arg168 backbone NH (a donor); only the delta tautomer satisfies both. The dynamics
  confirm it: NE2→Arg168 holds at 2.29 Å / 81% occupancy through the run. All ten
  histidines are assigned this way, with their evidence recorded.
- **The DNA-side linker chemistry is ours, not the literature's.** No source in the
  reading list says how sTCO reaches DNA — the only characterised reagent is
  sTCO-PEG5000, a gel-shift diagnostic. We use the standard route from sTCO's own
  precursor, activated as an NHS carbonate and reacted with a 5′-amino-modifier-C6
  oligo, giving a carbamate: `…C9(H)–CH₂–O–C(=O)–NH–(CH₂)₆–O–PO₂⁻–5′DNA`. Swap it by
  editing `build/02_params/residue_defs.py`; the force ladder is recomputed from the
  measured arm length.
- **The 4,5-dihydropyridazine tautomer**, as drawn in Eddins Fig. 1B. Aromatisation to
  pyridazine is blocked by the sp³ ring-fusion carbons of the bicyclononane.
- **One sTCO diastereomer**, (1R,8S,9S). The cycloaddition can occur on either
  tetrazine face and the literature does not resolve the product.
- **Terminal ACE/NME caps** at Ser2 and Gly232. Both positions are internal in the real
  construct (Met1 precedes; GSHHHHHH follows), so charged termini would be an invented
  pair of charges. Residue 1 and 233–246 are zero-occupancy stubs in the deposition.
- **Charges: AM1-BCC side chains, ff14SB backbone.** AM1-BCC is the self-consistent
  choice for GAFF2-typed atoms — GAFF targets AM1-BCC charges, themselves fitted to
  HF/6-31G(d) RESP. The six backbone atoms are then overwritten with ff14SB's types and
  charges and the residual spread over the side chain, so the peptide bonds into and
  out of the modified sites are described like every other. The correction is small
  (−0.003 e/atom for TET), which is itself reassuring.
- **A 27 bp spring, not 20.** With the sites 31.3 Å apart and the tether reaching
  22.3 Å, a 20 bp duplex's 6.6 nm contour is shorter than the span it must cross: it
  would be stretched, and the linkers rather than the DNA would carry the strain.
- **12-6 Mg²⁺, not 12-6-4.** OpenMM has no native C4 term, so a 12-6-4 prmtop would
  silently lose the polarisation correction on the way in.

`build/02_params/PARAMS_REPORT.json` records every residue's formula, mass, charge
method, backbone correction and frcmod audit. `build/02_params/junction.frcmod`
documents each bridging parameter and the Amber/GAFF term it was copied from.

## Verification

| Check | Status |
|---|---|
| `python -m openmm.testInstallation`, CUDA | passes; all four platforms agree |
| CRO atom names match the xFPchromophores template | 22/22, checked before use |
| CRO net charge from the library | −1.0000 exactly |
| `parmchk2` parameters flagged "ATTN, need revision" | **0** across TET, TDP, DNL |
| Every assigned GAFF type implies the declared element | 111/111 atoms |
| Every prep tree parent link is a declared bond | passes for all three residues |
| tleap errors on the assembled chimera | **0** |
| Net charge, every system | integral to <1e-4 |
| Both carbamate bonds present in the prmtop | verified per system |
| Linker bonds after minimisation | C–N 1.365/1.366 Å (ideal ~1.38); O–P 1.599/1.594 Å (ideal 1.61) |
| Spring model vs Zocchi's published forces | 1.47/1.5, 2.33/2.4 pN; 6.62/6.6 k<sub>B</sub>T |
| Bend transform is the identity at zero curvature | O3′–P 1.601 Å at R = 566 Å |
| Chromophore H-bond network, 50 ns WT | Thr203 97%, Arg96 100%, Gln94 95% |

## Layout

```
env/          pinned conda specs, activate.sh (works around the dead /n/hekstra_lab
              paths still in ~/.bashrc)
build/01_protein/   2B3P -> tleap-ready structure; sites, SASA, tautomers, force ladder
build/02_params/    TET / TDP / DNL: tables -> RDKit -> antechamber -> prepgen -> frcmod
build/03_dna/       ideal B-form duplex, bent to a target span
build/04_assemble/  tether posing, duplex docking, tleap driver
build/systems/      one prmtop/inpcrd/residue_map per system
runs/               the OpenMM engine (minimise / equilibrate / produce, restart-safe)
slurm/              one submit template, submit_all.sh, status.sh
analysis/           spring model + its tests, chromophore, mechanics, barrel,
                    collect (aggregate), figures, animate
figures/            report figures, the CSV behind each, and the animations
```

## Notes for whoever runs this next

- `~/.bashrc`'s mamba block still points at `/n/hekstra_lab/…`, retired when the share
  moved to `/n/lab_storage/hekstra_lab`. `env/activate.sh` works around it rather than
  editing your dotfiles; fixing the `.bashrc` is worth doing separately.
- **tleap cannot handle a path containing a space**, and this directory is literally
  `MD simulations`. Every driver copies its inputs into the working directory and
  refers to them by basename.
- `gpu_test` is entirely MIG `a100_3g.20gb` slices: fine for correctness, useless for
  timing. Production runs at ~427 ns/day on one such slice at 4 fs with HMR.
- Never write trajectories to `/n/lab_storage` — 34 TB filesystem, 307 GB free.
