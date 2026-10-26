# All-atom MD of the PP2A Aα–ssDNA spring chimera

The sfGFP DNA-spring pipeline (`../sfGFP`) applied to the scaffold subunit of PP2A: the
Aα monomer (PR65, UniProt P30153, 589 residues, fifteen HEAT repeats) with **Tet2-Et at
Asp7 and Asn553** — construct positions **D6X / N552X**; the construct counts from Ala2 —
clicked to a single 50 nt DNA strand modified at both ends, closing a covalent loop.
No B or C subunit. Wiki page: `docs/wiki/PP2A-DNA-Spring-MD.md`.

Nothing chemical is new. The custom residues (TET, TDP, DNL, DNP, DNH), the DNA
builder, the assemblers, the tleap driver and the OpenMM engine are the sfGFP ones,
called by path or symlinked (`build/02_params`, `runs/openmm_engine.py`). What is new
is the protein preparation, the run orchestration for this cluster's constraints, and
the analyses a HEAT solenoid and a single strand need that a β-barrel and a duplex did
not: RMSD in Day et al.'s split, bend/twist/aperture of the solenoid, the strand's shape
and modelled tension, and a measured force–extension curve for the bare strand.

## Quick start

```bash
source env/activate.sh                 # holylabs env, netscratch data root; no lab_storage
bash build/build_pp2a.sh               # protein, spring, P1_tet, P2_clicked, P3_ss50, D_ss50
sbatch slurm/orchestrator.sbatch       # walks slurm/jobs.txt: free benchmarks, then production
bash slurm/status.sh                   # queue, per-run progress, fairshare
bash finalize.sh                       # collect -> findings -> figures -> animations -> wiki page
```

## The systems

| ID | system | atoms | what it is for |
|---|---|---|---|
| `P0_apo` | AlphaFold Aα alone | 163,423 | the cryo-EM reproduction's finished 100 ns run (ff19SB/OPC); baseline for RMSD/RMSF/φψ |
| `P1_tet` | Aα + 2× Tet2-Et, unclicked | 121,808 | the 2-tet sample; φ/ψ before clicking |
| `P2_clicked` | tethers clicked to aminohexanol, no DNA | 125,624 | **the zero reference** for every chimera number |
| `P3_ss50` | 50 nt ssDNA, closed loop | 259,465 | the chimera |
| `D_ss50` | the strand alone | 277,063 | restraint windows x₀ = 60/80/100 Å: the force field's own force–extension curve |

One 100 ns trajectory per system (Ian's choice, 2026-10-26), so comparisons use block
means (ten blocks per run) and say so.

## Where the numbers come from

| quantity | value | source |
|---|---|---|
| start model | AlphaFold DB v6, AF-P30153-F1 | `$PP2A_DATA/md/systems/Aa/protein.pdb` (cryoEM reproduction) |
| numbering | construct n = UniProt n + 1 | `build/01_protein/attachment_sites.json["numbering_note"]` |
| Cβ7–Cβ553 | 54.4 Å (AF model) | `attachment_sites.json` |
| anchor span | 90 Å (chosen) | `build/build_pp2a.sh` — 88 % tether extension, 11.7 Å lift, as sfGFP |
| ssDNA model | Marko–Siggia, 0.59 nm/nt, ξ = 1.5 nm, zero below √(2ξL) | `../sfGFP/analysis/spring_model.py:SsSpring` |
| HEAT grid | repeat k = 8 + 39(k−1) … 8 + 39k − 1 | `cryoEM/pp2a_day2025_reproduction/scripts/25_geometry.py` |
| literature end-to-end | Cα Asn30–Phe578 (= Yilmaz "N29–F577") | `analysis/common.py` |
| force field | ff14SB + OL15 + TIP3P + Joung–Cheatham + Li/Merz Mg²⁺ + GAFF2 | `../sfGFP/build/04_assemble/build_system.py` |

## Compute policy

Paid GPUs are allowed on this project provided fairshare is spent well. The
orchestrator therefore runs 1 ns benchmarks on the free `gpu_test` MIG slices, sizes
every `gpu_requeue` production chunk from the measured ns/day (×1.15 — the Job Defense
Shield flags jobs that use < 50 % of their limit), refuses `gpu` unless asked in words,
and keeps the restraint windows on `gpu_test`. `gpu_requeue` is half the price of `gpu`;
preemption costs at most 200 ps because the engine checkpoints that often.

## Layout

```
env/activate.sh       paths: holylabs env, netscratch data ($PP2A_MD_SCRATCH), apo run ($PP2A_DATA)
build/01_protein/     prepare_pp2a.py -> pp2a_prepped.pdb + attachment_sites.json
build/02_params ->    ../../sfGFP/build/02_params   (custom residues, unchanged)
build/03_dna/         spring_50nt.{pdb,json}; bare_strand.py (D_ss50)
build/04_assemble/    the posed solutes and their tleap bond lists (*.json)
build/systems/        prmtop/inpcrd (not in git), residue_map.json, clamp_atoms.txt, build_report.json
runs/openmm_engine.py -> ../../sfGFP/runs/openmm_engine.py
slurm/                md_chain.sbatch (one run, restart-safe), orchestrator.sbatch, jobs.txt, status.sh
analysis/             common, rmsd, heat, ssdna, tet_sites, collect, findings, figures, animate,
                      force_windows; summary.csv, findings.json, force_extension.csv
figures/              pp2a_*.png with their CSVs; animations
report/               build_wiki_pp2a.py -> docs/wiki/PP2A-DNA-Spring-MD.md
```

Trajectories: `/n/netscratch/hekstra_lab/Lab/ian_poon/pp2a-md/runs/<TAG>/rep1/0{1,2,3}_*`
(90-day purge). Logs: `/n/netscratch/hekstra_lab/Lab/ian_poon/logs/pp2amd-*`.
