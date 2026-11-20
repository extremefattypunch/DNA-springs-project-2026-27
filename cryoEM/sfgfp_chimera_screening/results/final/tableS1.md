# Table S1 analogue: cryo-EM data collection and processing for the two final maps

sfGFP + 50 bp DNA screening session (12-11-2026), CryoSPARC v5 project P6. Columns follow Day et al. 2025 *Structure* Table S1. Every cell's source file and field: `tableS1_sources.json`; built by `scripts/50_figS6_tableS1.py`.

| | W3 rod_no_dot (“free DNA”) control, J39 | W4 dot-on-rod (J56) |
|---|---|---|
| **Data collection and processing** | | |
| Magnification | 165,000× (nominal) | 165,000× (nominal) |
| Voltage (kV) | 200 | 200 |
| Electron exposure (e–/Å²) | 39.72 (vacuum calibration; 28.2 detected through sample); 50 fractions, 2.94 s | 39.72 (vacuum calibration; 28.2 detected through sample); 50 fractions, 2.94 s |
| Defocus range (μm, underfocus) | 0.74–2.11 (p5–p95 of 8,036 curated micrographs; median 1.41) | 0.74–2.11 (p5–p95 of 8,036 curated micrographs; median 1.41) |
| Pixel size (Å) | 0.71 physical, as used in processing (EPU xml calibration 0.7175) / 1.42 processed (motion-corrected, F-crop ½) / 2.84 refined (128 px box); if 0.7175 is correct, all Å values ×1.0106 (see pixel-size note) | 0.71 physical, as used in processing (EPU xml calibration 0.7175) / 1.42 processed (motion-corrected, F-crop ½) / 2.84 refined (160 px box); if 0.7175 is correct, all Å values ×1.0106 (see pixel-size note) |
| Symmetry imposed | C1 | C1 |
| Micrographs (collected / curated / contributing to map) | 9,019 / 8,036 / 7,866 | 9,019 / 8,036 / 7,831 |
| Initial particle images (no.) | 1,583,293 rod-stream blob picks (J18); 606,539 extracted after power filter | 1,583,293 rod-stream blob picks (J18); 224,068 re-centred and re-extracted (J44) |
| Final particle images (no.) | 64,013 | 54,158 |
| Map resolution (Å) | 7.8 (CryoSPARC GSFSC, tight mask) / 8.7 (RELION 3.1, loose mask, phase-randomised) | 7.9 (CryoSPARC GSFSC, tight mask) / 10.8 (RELION 3.1, loose mask, phase-randomised) |
| &nbsp;&nbsp;FSC threshold | 0.143 | 0.143 |
| Map resolution, other masks (Å) | no mask 9.7 / spherical 8.9 (CryoSPARC); RELION with CryoSPARC mask 7.9 | no mask 11.2 / spherical 9.8 (CryoSPARC); RELION with CryoSPARC mask 8.0 |
| Map resolution range (Å) | 8.2–10.3 (RELION locres, p5–p95 in tight mask; median 9.1); 7.4–14.1 (CryoSPARC J42; median 11.4) | 8.7–11.6 (RELION locres, p5–p95 in tight mask; median 10.3); no CryoSPARC local-resolution job |
| Orientation: cFAR / SCF* | 0.27 / 0.94 | 0.17 / 0.70 |
| **Refinement** | | |
| Initial model used (PDB code) | not built: ~8–9 Å map of a rod with a terminal lobe, no secondary structure resolved, so an atomic model would be unconstrained | not built: ~11 Å anisotropic average (cFAR 0.17) over heterogeneous sfGFP–DNA assemblies; nothing to refine a model against |
| Model resolution (Å) | not built (no model) | not built (no model) |
| Map sharpening B factor (Å²) | -1048 (RELION auto-B, CryoSPARC mask) / -926 (loose mask); Guinier fit from 10 Å only, r = 0.74 / 0.92 | -714 (RELION auto-B, CryoSPARC mask) / -823 (loose mask); Guinier fit from 10 Å only, r = 0.80 / 0.76 |
| Model composition: non-hydrogen atoms / residues / ligands | not built (no model) | not built (no model) |
| B factors (Å²) | not built (no model) | not built (no model) |
| R.m.s. deviations: bond lengths (Å) / bond angles (°) | not built (no model) | not built (no model) |
| **Validation** | | |
| MolProbity score | not built (no model) | not built (no model) |
| Clashscore | not built (no model) | not built (no model) |
| Poor rotamers (%) | not built (no model) | not built (no model) |
| Ramachandran plot: favored / allowed / disallowed (%) | not built (no model) | not built (no model) |

Notes. Tight mask = CryoSPARC Validation (FSC) auto-tightened mask (column `fsc_resmask`); RELION values are `relion_postprocess` phase-randomisation-corrected FSC (`scripts/32_relion_postprocess.sh`). The W4 loose-mask value is the honest whole-particle figure (report §13b). Local-resolution ranges are the 5th–95th percentiles of RELION local-resolution voxels inside the CryoSPARC validation (tight FSC) mask, the same definition as the Fig S5 summary boxes (other voxel sets give other ranges; they are listed, labelled as secondary, in `local_resolution_stats.json`). The W3 column additionally lists CryoSPARC Local Resolution (J42, same tight mask) for reference; no equivalent CryoSPARC job was run for W4. RELION auto-B values are fitted from 10 Å to Nyquist (5.68 Å) on ~8–11 Å maps, so they are poorly determined and are listed only because Table S1 asks for the sharpening B factor.

Pixel size: processed with 0.71 Å/px (user screenshot); EPU xml calibration 0.7175 Å/px (all 9,019 movies) — if 0.7175 is correct, all Å values scale by ×1.0106; to be confirmed with the facility.
