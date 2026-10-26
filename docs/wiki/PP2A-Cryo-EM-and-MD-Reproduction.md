# Reproducing Day et al. 2025 — PP2A Aα vs Aβ holoenzymes (cryo-EM + MD)

A reproduction of the key computational results of **Day, Huang, Leonard, O'Connor, Narla & Taylor,
*Structure* 33:1688–1699 (2025)**, "Regulatory mechanisms of PP2A complex assembly driven by
physicochemical differences in A-subunit isoforms", on the FASRC Cannon cluster with the lab's CryoSPARC
instance and open-source tools installed for the purpose. Separate from the DNA-spring work on this wiki.

**Full write-up:** [`report.md`](https://github.com/extremefattypunch/DNA-springs-project-2026-27/blob/main/cryoEM/pp2a_day2025_reproduction/report.md)
 · [runbook / README](https://github.com/extremefattypunch/DNA-springs-project-2026-27/blob/main/cryoEM/pp2a_day2025_reproduction/README.md)
 · [plan](https://github.com/extremefattypunch/DNA-springs-project-2026-27/blob/main/cryoEM/gpt-plans/pp2a-day2025-reproduction.md)
 · [scripts](https://github.com/extremefattypunch/DNA-springs-project-2026-27/tree/main/cryoEM/pp2a_day2025_reproduction/scripts)
 · [result tables](https://github.com/extremefattypunch/DNA-springs-project-2026-27/tree/main/cryoEM/pp2a_day2025_reproduction/results)
 · [figures](https://github.com/extremefattypunch/DNA-springs-project-2026-27/tree/main/cryoEM/pp2a_day2025_reproduction/figures)

*Status 2026-09-26:* cryo-EM part **done**; MD part **running** (see [Pending](#pending)).

---

## What the paper claims, and what could be tested

The paper compares the two PP2A scaffold isoforms Aα and Aβ (Aβ carries a 12-residue N-terminal extension,
MAGASELGTGPG) by all-atom MD of the free A subunits and by cryo-EM of the Aα-B56γ-Cα (3.4 Å, EMD-48300 /
9MIP) and Aβ-B56γ-Cα (3.2 Å, EMD-48224 / 9MF5) holoenzymes. Its computational claims:

1. **MD:** the Aβ N-terminus is far more mobile than Aα's (RMSD/RMSF, Fig 1; Videos S1/S2).
2. **Cryo-EM:** the B56γ C-terminal helix is visible at 6 Å in the Aα map but absent in Aβ (Fig 2B);
   Aβ is more concave (HEAT1–15 distance 59.6 vs 62.7 Å; salt bridge 12.4 vs 13.2 Å; Fig 3); the
   A–B56γ and A–Cα interface residues are conserved (Fig 3D/E); the C-subunit tail occupies the A–B56 groove
   in 9/10 Aβ 3D classes vs 1/10 for Aα (Fig 4/S7).

**No raw cryo-EM data were deposited** — no EMPIAR entry, no movies, micrographs or particle stacks; only
maps, half maps, FSC curves and models. So every particle-level result (2D classes, viewing distributions,
the Fig 4/S7 3D classification and C-tail counts, the Fig S6 particle numbers) is **out of reach**, and the
cryo-EM reproduction is a map- and model-level one. A synthetic-particle test of the CryoSPARC pipeline
was considered and deliberately not done.

---

## Cryo-EM results (done)

### Table S1 and resolution

| quantity | paper Aα | ours Aα | paper Aβ | ours Aβ | how |
|---|---|---|---|---|---|
| map resolution, FSC 0.143 (Å) | 3.4 | 3.42 (deposited curve) · 3.50 (CryoSPARC v5) · 3.50 (Phenix) · 3.54 (RELION, loose mask) | 3.2 | 3.22 · 3.33 · 3.34 · 3.33 | half maps re-analysed with three tools |
| model resolution, FSC 0.5 (Å) | 3.8 | 3.80 | 3.6 | 3.63 | phenix.mtriage |
| MolProbity / clashscore | 1.25 / 4.88 | **1.25 / 4.88** | 1.15 / 3.59 | **1.15 / 3.59** | phenix.molprobity on the deposited models |
| Ramachandran favoured / allowed / outliers (%) | 98.08 / 1.92 / 0 | **98.08 / 1.92 / 0** | 98.22 / 1.78 / 0 | **98.22 / 1.78 / 0** | |
| RMS bonds (Å) / angles (°) | 0.010 / 1.291 | 0.0102 / 1.29 | 0.009 / 1.241 | 0.0089 / 1.24 | |

Every validation row of Table S1 comes back to the printed digit. The half-map resolution recomputed from
the deposited half maps lands at ≈3.50 / 3.33 Å with CryoSPARC v5, Phenix and RELION alike; the extra
0.1 Å in the paper comes from the refinement's own tighter FSC mask — a RELION sweep of eight masks moves
the number from 3.93 Å (tight, contour-matched mask) to a 3.56 Å plateau (loose mask).

![FSC curves](images/pp2a/figS5_fsc.png)

### Fig 2 — the maps, and the B56γ C-terminal helix

The deposited (EMReady-sharpened) maps, coloured by nearest subunit (A red/grey, B56γ blue, Cα orange):

![Fig 2A analogue](images/pp2a/fig2a_maps.png)

**Fig 2B reproduces.** 9MIP models B56γ to residue 449, 9MF5 only to 403; residues 402–436 of 9MIP chain B
form the helix that reaches the Aα N-terminus (6.9 Å away). Carrying that helix into the Aβ frame (B56γ
cores superposed, RMSD 1.4 Å) and low-passing the half-map averages to 6 Å, the density around the helix
relative to the B56γ core is **0.60 in Aα vs 0.30 in Aβ** (6.2 σ vs 3.3 σ above solvent), while control
regions on the A and C subunits give 1.0 in both maps. Visually the ribbon is enclosed by continuous mesh in
Aα and sits outside the density in Aβ. Pre-registered thresholds (≥ 0.6 and ≤ 0.25) were missed by a hair
on both sides: the two-fold contrast is real, "absent" overstates it.

![Fig 2B analogue](images/pp2a/fig2b_lowpass6.png)

### Fig 3 — Aβ is more concave; the interfaces are conserved

| quantity | paper Aα | ours Aα | paper Aβ | ours Aβ |
|---|---|---|---|---|
| HEAT1–HEAT15 centroid distance (Å) | 62.7 | 63.3 | 59.6 | 60.4 |
| HEAT3–HEAT5 salt bridge, Cα–Cα (Å) | 13.2 (R113–E191) | 12.7 | 12.4 (R125–E203) | 12.4 |
| Fig 3D/E interface residues within 4 Å of B56γ / Cα | 11 / 11 | 11 / 11 | 11 / 11 | 11 / 11 |

The Aα–Aβ difference in the HEAT1–15 distance is 2.9 Å here against the paper's 3.1 Å (both of our absolute
values are 0.6–0.8 Å higher, a HEAT-boundary-definition effect). The per-repeat centroid shift after
superposing on HEAT 12–15 has the paper's profile — 2–3 Å across HEAT 1–7, falling to ~0 at HEAT 13–14:

![Fig 3B analogue](images/pp2a/fig3b_com_shift.png)

![Fig 3C analogue](images/pp2a/fig3c_salt_bridge.png)

### Fig S5 — local resolution

RELION 3.1 (the paper's tool) and CryoSPARC v5 agree that Aβ is the more uniformly resolved map and that the
B56γ periphery is the worst region. Within 3 Å of the model: RELION 3.48–4.47 Å (Aα) and 3.30–3.91 Å (Aβ);
CryoSPARC 3.21–6.67 / 3.00–4.99 Å (paper's Table S1 ranges 5.87–3.04 / 5.22–2.82 Å).

![Fig S5 analogue](images/pp2a/figS5_localres.png)

### EMReady

Re-running EMReady v1.1 on the half-map averages gives maps that correlate 0.89 (Aα) and 0.80 (Aβ) with the
deposited sharpened maps — EMReady-like, not identical (the authors started from the refinement's own
unsharpened map, possibly with a later EMReady release).

---

## MD (Tier 1: one 100 ns run per isoform)

Set-up follows the paper: AlphaFold full-length models (AlphaFold DB v6, since the AlphaFold3 Server
models were not supplied), ff19SB + OPC via tleap, 16 Å padding, 150 mM NaCl; restrained minimisation,
30-stage heating to 300 K, five-stage restraint release, then 100 ns NPT. The engine is **OpenMM 8.6**
rather than NAMD (NAMD's download is registration-gated; the Hamiltonian is identical, the barostat and
cutoff conventions differ and are listed in the report).

| system | particles | box (Å) | ions | status |
|---|---|---|---|---|
| Aα (P30153, 589 aa) | 163,423 | 129 × 144 × 78 | 24 Na⁺ + 105 NaCl | benchmark 53.8 ns/day; 100 ns production started 2026-09-26 18:19 |
| Aβ (P30154, 601 aa) | 225,676 | 129 × 137 × 112 | 25 Na⁺ + 147 NaCl | queued behind Aα |

### Pending

- **MD production** — ≈1.9 days per isoform on one free `gpu_test` A100 slice; the two run back to back
  because the partition allows one GPU job at a time next to the Open OnDemand desktop. Expected to finish
  ≈ **2026-09-30**. An orchestrator job resubmits 12 h chunks automatically.
- **Fig 1A/1B** (RMSF; RMSD of the N-terminus and of the rest), **Fig 1C** (N-terminus at 0–100 ns in 20 ns
  steps, red→blue) and **Videos S1/S2** (smoothed trajectories, N-terminus blue) — scripts written
  (`10`–`15`), run as soon as the trajectories exist.
- Report §A (the MD answer) and §B will be filled in then. Not planned: the Fig S3 replicate set (three AF3
  replicates + cryo-EM-start runs per isoform).

---

## How it was done

- **Deposited data only**: EMD-48300 / 48224 maps, half maps and FSC XML; PDB 9MIP, 9MF5, 3FGA; UniProt
  P30153 / P30154; AlphaFold DB models.
- **CryoSPARC v5.0.7** (lab instance, free `Cannon-test` lane): Import Volumes → Local Resolution
  Estimation → Local Filtering → Validation (FSC), driven with cryosparc-tools and watched in the GUI.
- **Installed on holylabs for this work** (no SBGrid): RELION 3.1.4, GROMACS 2023.3, ChimeraX 1.12 (el9 build
  run in a CentOS Stream 9 Singularity image with an rpm sysroot; offscreen OSMesa rendering), EMReady v1.1.
  Phenix 2.2 from a labmate's install. Everything on the free `test`/`gpu_test` partitions; zero paid GPU hours.
- Every number quoted above is a row in `results/*.csv`; `software_versions.json` records the versions.

## What is not reproducible

2D classes, viewing-direction plots, the 10-class 3D classification and its 9/10 vs 1/10 C-tail occupancy
(Fig 4, S7), the Fig S6 particle counts and any picking/refinement — all need particles that were never
deposited. Model building was not redone (the deposited models were validated, not rebuilt).
