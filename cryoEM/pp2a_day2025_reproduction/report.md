# Report: reproducing the key cryo-EM and MD results of Day et al. 2025

Paper: *Regulatory mechanisms of PP2A complex assembly driven by physicochemical differences in A-subunit
isoforms*, Structure 33:1688–1699. Plan: `../gpt-plans/pp2a-day2025-reproduction.md`. Numbers below are in
`results/*.csv`; figures in `figures/`; tool versions in `software_versions.json`; data paths in `DATA.md`.

**Scope.** No movies, micrographs or particle stacks were deposited (no EMPIAR entry), so the cryo-EM half is a
map- and model-level reproduction from EMD-48300 / PDB 9MIP (Aα-B56γ-Cα) and EMD-48224 / PDB 9MF5
(Aβ-B56γ-Cα): primary maps, half maps, deposited FSC curves and models. Everything particle-based (Fig S5B/D/F/H
2D classes and viewing distributions, Fig 4 and S7 3D classification and C-tail counts, Fig S6 particle counts,
picking) is not reproducible and is not attempted. The MD half (Fig 1A–C, Videos S1/S2) is Tier 1: one 100 ns
run per isoform, in OpenMM because NAMD is behind a registration wall on this cluster.

**Status (2026-09-26, 18:20).** Cryo-EM: complete. MD: systems built (Aα 163,423 particles, Aβ 225,676; 16 Å
padding, OPC, 150 mM NaCl); the Aα 1 ns benchmark ran at 53.8 ns/day (2 fs, one gpu_test A100 slice), so each
100 ns run needs ≈1.9 days; the Aα production chain started 18:19 and an orchestrator job feeds 12 h chunks one
after another (Aα, then Aβ). §A is filled in when they finish (≈2026-09-30).

## Headline comparison with the paper

| Quantity | Paper Aα | Ours Aα | Paper Aβ | Ours Aβ | Method |
|---|---|---|---|---|---|
| Map resolution, FSC 0.143 (Å) | 3.4 | 3.42 (deposited curve) · 3.50 (CryoSPARC v5) · 3.50 (mtriage) · 3.54 (RELION, loose mask) | 3.2 | 3.22 · 3.33 · 3.34 · 3.33 | §C |
| Model resolution, FSC 0.5 (Å) | 3.8 | 3.80 | 3.6 | 3.63 | phenix.mtriage, masked |
| MolProbity / clashscore | 1.25 / 4.88 | 1.25 / 4.88 | 1.15 / 3.59 | 1.15 / 3.59 | phenix.molprobity |
| Ramachandran favoured / allowed / outliers (%) | 98.08 / 1.92 / 0 | 98.08 / 1.92 / 0 | 98.22 / 1.78 / 0 | 98.22 / 1.78 / 0 | phenix.molprobity |
| RMS bonds (Å) / angles (°) | 0.010 / 1.291 | 0.0102 / 1.29 | 0.009 / 1.241 | 0.0089 / 1.24 | phenix.molprobity |
| HEAT1–HEAT15 CoM distance (Å) | 62.7 | 63.3 | 59.6 | 60.4 | Cα centroids, 39-residue HEAT grid |
| Salt-bridge Cα–Cα (Å) | 13.2 (R113–E191) | 12.7 | 12.4 (R125–E203) | 12.4 | deposited models |
| Largest per-HEAT CoM shift after HEAT 12–15 fit | HEAT 3–6, 2–3 Å | HEAT 1 (3.2), 4 (2.9), 6 (2.8), 5 (2.7); HEAT 13–14 < 0.4 | – | – | Fig 3B analogue |
| B56γ C-helix density at 6 Å, helix/B56γ-core ratio | visible | 0.60 (6.2 σ) | absent | 0.30 (3.3 σ) | Fig 2B analogue, controls 1.0 |
| Interface residues of Fig 3D/E in contact (≤ 4 Å) | 11/11 | 11/11 | 11/11 | 11/11 | deposited models |

## A. Do N-terminal RMSD/RMSF separate Aβ from Aα as in Fig 1?

*Pending.* Systems are built (`md/systems/build.json`): Aα box 129.3 × 144.3 × 77.9 Å, 163,423 particles,
24 neutralising Na⁺ (paper 23) + 105 NaCl pairs; Aβ box 128.7 × 136.6 × 112.4 Å, 225,676 particles (the
low-confidence 25-residue tail sticks out of the horseshoe), 25 Na⁺ (paper 25) + 147 NaCl pairs; padding
16.8–17.7 Å per face (paper 18–19.5 Å); 150 mM salt. Engine: OpenMM 8.6 port of the NAMD protocol (`scripts/03_openmm_pp2a_engine.py`):
restrained minimisation (10 kcal/mol/Å² on protein heavy atoms), 30 × 10,000-step heating to 300 K, five 50,000-step
restraint-release stages under a Monte Carlo barostat, then 100 ns NPT at 2 fs with Langevin friction 5 ps⁻¹.
Deviations from NAMD 3.0.1: barostat type (MC vs Langevin piston), Amber 10 Å hard cutoff with dispersion
correction instead of NAMD's 12/10 Å switching, timestep unstated in the paper (2 fs chosen). Starting models are
AlphaFold DB v6 (AF2) full-length predictions (mean pLDDT 94.9 / 92.8; N-terminal pLDDT 54.9 for Aα 1–13 and 43.5
for Aβ 1–25) — the AlphaFold Server (AF3) models the paper used were not supplied. Results (RMSD of residues 1–13 /
1–25 and of the rest, per-residue RMSF, Fig 1C overlay, Videos S1/S2) will be added here.

## B. Which engine and starting models, and how does a single replicate limit A?

OpenMM 8.6 + AmberTools 24.8 tleap (ff19SB, OPC), AlphaFold DB v6 models, one replicate per isoform (Tier 1, Ian's
choice). The paper's Fig S3 shows three AF3 replicates plus a cryo-EM-start run per isoform; with one replicate we
can show whether the qualitative separation appears, not its reproducibility across starting conformations of a
low-confidence tail (pLDDT 43–55). Any statement about *how much* more mobile Aβ's N-terminus is carries the
usual single-trajectory caveat.

## C. Are the deposited map resolutions and Table S1 statistics reproduced?

Yes for everything that does not depend on a mask choice, and the mask dependence is quantified.

* **Model validation is reproduced exactly.** `phenix.molprobity` on 9MIP / 9MF5 gives MolProbity 1.25 / 1.15,
  clashscore 4.88 / 3.59, Ramachandran favoured 98.08 / 98.22 %, no outliers, RMS bonds 0.0102 / 0.0089 Å, angles
  1.29 / 1.24° — every Table S1 validation row to the printed precision. `phenix.mtriage` gives model-vs-map
  FSC 0.5 at 3.80 / 3.63 Å (Table S1: 3.8 / 3.6).
* **Half-map resolution.** The deposited FSC curves cross 0.143 at 3.42 / 3.22 Å. Recomputing from the deposited
  half maps: CryoSPARC v5 Validation (FSC) with its auto mask 3.50 / 3.33 Å; Phenix mtriage masked 3.50 / 3.34 Å;
  RELION 3.1 postprocess (phase-randomisation corrected) 3.54 / 3.33 Å with a loose mask (15 Å-lowpassed
  half-map average, threshold 0.2× the contour-matched value, 3 px extension, 6 px soft edge) but 3.93 / 3.55 Å with a
  tight mask matched to the EMDB contour volume. A sweep of eight masks (`22b_relion_mask_sweep.sh`) plateaus at
  3.56 / 3.34 Å once the mask stops cutting into the molecule. All three independent tools therefore agree at
  ≈3.5 / 3.33 Å; the deposited 3.4 / 3.2 Å come from the refinement's own tighter FSC mask. Within the paper's
  precision the Aβ value is reproduced (±0.15 Å criterion met: 3.33 vs 3.2 is 0.13); Aα misses the criterion by
  0.1 Å (3.50 vs 3.4). `figures/figS5_fsc.png` overlays all curves.
* **Local resolution** (Fig S5A/E). RELION 3.1 `--locres` (the paper's tool) within 3 Å of the model spans
  3.48–4.47 Å (Aα) and 3.30–3.91 Å (Aβ) at the 1st–99th percentile; CryoSPARC v5 Local Resolution spans
  3.21–6.67 / 3.00–4.99 Å. The paper's Table S1 ranges (5.87–3.04 / 5.22–2.82 Å) match the CryoSPARC-style spread
  better than RELION's compressed one; both tools agree that Aβ is uniformly better resolved and that the
  B56γ periphery is the worst region (`figures/figS5_localres.png`, blue lobes).

## D. Is the B56γ C-terminal helix density present at 6 Å in Aα and absent in Aβ?

Yes, qualitatively and quantitatively, with the caveat that the ratios sit just outside the pre-registered cut-offs.

9MIP models B56γ to residue 449 and 9MF5 only to 403; residues 402–436 of 9MIP chain B form a helix (DSSP) whose Cα
atoms come within 6.9 Å of the Aα N-terminus (residues 7–20) — this is the "B56γ C-terminal helix" of the paper.
Its atoms were carried into the Aβ frame by superposing the B56γ cores (38–403 Cα, RMSD 1.37 Å). In the half-map
averages low-passed to 6 Å (`relion_image_handler`), the mean density in a 3 Å envelope around the helix, relative to
the B56γ core, is **0.60 for Aα and 0.30 for Aβ** (6.2 σ vs 3.3 σ above solvent; identical on the low-passed
deposited maps, 0.60 / 0.31). Control regions (HEAT 7 of the A subunit, the C-subunit core) give 1.07 / 1.00 and
1.06 / 1.06, so the difference is specific to the helix. The pre-registered thresholds were r_Aα ≥ 0.6 (got 0.595)
and r_Aβ ≤ 0.25 (got 0.30): the two-fold contrast the paper describes is there, but "absent" overstates it — a
third of the core density remains in the Aβ helix region, consistent with the paper's own description that the
helix appears only at low contour and was missing from every Aβ 3D class. `figures/fig2b_lowpass6.png` shows the
ribbon fully enclosed by 6 Å mesh in Aα and outside the density in Aβ.

## E. Are the HEAT CoM distances, per-repeat shifts and salt-bridge distances reproduced within 0.5 Å?

Mostly. With HEAT repeats defined on the canonical 39-residue grid (HEAT k = 8+39(k−1) … 8+39k−1 in Aα numbering,
anchored by R113 ∈ HEAT3, E191 ∈ HEAT5, W257 ∈ HEAT7, K416/R418 ∈ HEAT11, D572/D574 ∈ HEAT15; Aβ = +12):

* HEAT1–HEAT15 Cα-centroid distance 63.3 Å (Aα) vs 60.4 Å (Aβ); paper 62.7 vs 59.6. Both are 0.6–0.8 Å above the
  paper (outside the 0.5 Å criterion) but the Aα−Aβ difference, 2.9 Å, matches the paper's 3.1 Å. The offset is a
  boundary-definition effect: a DSSP-paired definition gives 68.9 Å for Aα and fails to pair 30 helices for Aβ,
  so the grid definition is reported.
* Per-repeat CoM shift after superposing on HEAT 12–15 Cα (fit RMSD 0.79 Å over 153 atoms): 3.2, 2.3, 2.5, 2.9, 2.7,
  2.8, 2.5, 2.2, 1.8, 1.2, 0.9, 0.7, 0.3, 0.2, 0.8 Å for HEAT 1–15 (`figures/fig3b_com_shift.png`). The paper's
  Fig 3B has the same shape — 2–3 Å across HEAT 1–7 decaying to ~0 at HEAT 13–14 — with its maximum at HEAT 3–6;
  ours peaks at HEAT 1 and 4–6.
* Salt bridge: Cα–Cα R113–E191 = 12.7 Å (paper 13.2) in Aα; R125–E203 = 12.4 Å (paper 12.4) in Aβ. Both bridges are
  formed (closest Arg N – Glu O 2.9 / 2.8 Å). The Aα value is 0.5 Å from the paper's; the paper's 0.8 Å Aα–Aβ
  difference is 0.3 Å here.

## F. Are the eleven interface residues in contact in both models?

Yes. Every residue the paper lists is within 4 Å of its partner in both deposited models: A–B56γ E100/E112,
E101/E113, R183/R195, W257/W269 and A–Cα K416/K428, R418/R430, Y456/Y468, Y495/Y507, R498/R510, D572/D584,
D574/D586. The full contact sets are 32 (Aα) vs 17 (Aβ) residues at the B56γ interface — the Aβ complex has ~60
fewer ordered B56γ residues (38–403 vs 21–449) — and 28 vs 30 at the Cα interface (`figures/fig3de_interfaces.png`).

## G. Does EMReady on the half-map average reproduce the deposited sharpened map?

Partly. EMReady v1.1 (the version the paper cites) applied to the half-map averages gives maps whose real-space
correlation with the deposited primary maps is 0.89 (Aα) and 0.80 (Aβ) over all voxels, 0.86 and 0.69 within the
top 1 % of deposited density, and the FSC between our map and the deposited one falls to 0.5 at 3.98 Å (Aα) and
4.51 Å (Aβ). That is EMReady-like sharpening (RELION's B-factor-sharpened map looks quite different) but not the
same map: the deposited maps were presumably made from the refinement's own unsharpened output, possibly with a
later EMReady release, and our resampling from EMReady's 1.0 Å grid back to the deposited grid adds interpolation
error. The pre-registered CC > 0.95 is not met for either isoform, and the Aβ agreement is clearly poorer, which
is consistent with its weaker B56γ density (EMReady's local enhancement is most sensitive where the input is weakest).

## H. What could not be reproduced, and why

* **Anything requiring particles**: Fig S5B/F 2D classes, Fig S5D/H viewing distributions, Fig 4 and S7 (10-class 3D
  classification of the final particles and the 9/10 vs 1/10 C-tail occupancy), Fig S6 particle counts, TOPAZ and
  blob picking, per-particle CTF. No EMPIAR deposition exists and the authors deposited only maps, half maps, FSC
  curves and models. A synthetic-particle test was deliberately not done (Ian's decision).
* **Model building**: the deposited models were validated, not rebuilt (Coot/ISOLDE/Phenix real-space refinement not
  rerun).
* **Tool substitutions**: NAMD → OpenMM 8.6 (same ff19SB/OPC Hamiltonian; barostat and cutoff differ); VMD →
  MDAnalysis and ChimeraX for RMSD, overlays and movies; AlphaFold3 → AlphaFold DB v6 (AF2) models; EMReady run
  on the half-map average rather than the authors' unsharpened refinement output. RELION 3.1, GROMACS 2023.3,
  ChimeraX 1.12, Phenix 2.2 and EMReady v1.1 are the paper's own tools (or newer point releases) and were installed
  on holylabs for this work.
* **Tier 2 MD** (Fig S3: three AF3 replicates and cryo-EM-start runs per isoform) was deferred.

## Deviations table (paper → here)

| Item | Paper | Here |
|---|---|---|
| MD engine / dt | NAMD 3.0.1, dt unstated | OpenMM 8.6, 2 fs, rigid bonds |
| Barostat | Langevin piston 1 atm (100 fs / 50 fs) | MonteCarloBarostat 1 atm, every 25 steps |
| Cutoff | NAMD default (12/10 Å switching implied) | 10 Å, PME, dispersion correction |
| Starting models | AlphaFold3 de novo | AlphaFold DB v6 (AF2) |
| Replicates | 3 AF3 + 1 EM-start per isoform | 1 AF per isoform |
| RMSD tool | VMD | MDAnalysis (same convention: fit and RMSD on the selection vs frame 1) |
| RMSF tool | GROMACS 2023.3 | GROMACS 2023.3 (`gmx rmsf -res -fit`) |
| Movies / figures | VMD / ChimeraX | ChimeraX 1.12 (offscreen) |
| Raw (unsharpened) map | refinement output | (half1 + half2)/2 |
| Sharpening | EMReady | deposited map used; EMReady v1.1 re-run for comparison |
| Local resolution | RELION 3.1 | RELION 3.1.4 and CryoSPARC v5.0.7 |
