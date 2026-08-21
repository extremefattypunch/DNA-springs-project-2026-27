# DNA springs on sfGFP — all-atom MD

**[▶ Read the full report](https://extremefattypunch.github.io/DNA-springs-project-2026-27/)**
 · [no-setup mirror](https://htmlpreview.github.io/?https://raw.githubusercontent.com/extremefattypunch/DNA-springs-project-2026-27/main/docs/index.html)
 · [code and structures](https://github.com/extremefattypunch/DNA-springs-project-2026-27/tree/main/MD%20simulations/sfGFP)

Two amber codons put the tetrazine amino acid **Tet2-Et** at Asp133 and Asn149 on
adjacent strands of sfGFP's β-barrel; a doubly-sTCO-modified DNA duplex clicks onto
both, closing a covalent loop. dsDNA is far too stiff to follow the barrel's curvature,
so it bends — pushing the two attachment points apart with a few piconewtons.

![the chimera](https://raw.githubusercontent.com/extremefattypunch/DNA-springs-project-2026-27/main/docs/figures/S3_spring27_overview_still.png)

## What the springs deliver

Forces are Zocchi's *f(x)* evaluated at the end-to-end distance each simulation settles
at, not at an assumed one. 1,650 ns across six systems, three replicates each.

| spring | span x | force | stored energy | base pairs open |
|---|---|---|---|---|
| 27 bp | 53.8 Å | **4.88 ± 0.06 pN** | 7.4 k<sub>B</sub>T | 0.03 of 27 |
| 40 bp | 61.1 Å | **2.97 ± 0.02 pN** | 8.8 k<sub>B</sub>T | 0.05 of 40 |
| 40 bp, nicked | 62.1 Å | **2.35 ± 0.08 pN** | 6.9 k<sub>B</sub>T | 0.03 of 40 |

The nick softens the spring in the direction Zocchi's two critical-bending-torque
values imply, and makes it visibly floppier. The stored energies bracket the
6.6 k<sub>B</sub>T he measured for the guanylate-kinase chimera by a different route.

## What survived its error bars, and what did not

- **Protein stiffness: a bound, not a value.** The axial clamp's compliance is
  0.0106 ± 0.0057 Å pN⁻¹ over 0–20 pN — the 95% interval includes zero. That supports
  κ > 107 k<sub>B</sub>T nm⁻², right on the ~100 his measurements imply, and nothing
  tighter.
- **The duplex bends, it does not kink.** The force model puts all three springs on its
  kinked branch, but curvature is spread at 3–8° per step with no localised spike and
  the duplex stays intact. The forces are model-dependent.
- **One apparent chromophore effect is an artefact.** Ser205–Glu222 occupancy is either
  ~0% or ~90% per replicate — a two-state rotamer. A mean of three coin flips is not a
  measurement.

![force response](https://raw.githubusercontent.com/extremefattypunch/DNA-springs-project-2026-27/main/docs/figures/fig1_force_response.png)
![chromophore hydrogen bonds](https://raw.githubusercontent.com/extremefattypunch/DNA-springs-project-2026-27/main/docs/figures/fig3_chromophore_hbonds.png)
![the chromophore pocket](https://raw.githubusercontent.com/extremefattypunch/DNA-springs-project-2026-27/main/docs/figures/S3_spring27_chromophore_still.png)

## Structures you can open

[`MD simulations/sfGFP/pdb_exports`](https://github.com/extremefattypunch/DNA-springs-project-2026-27/tree/main/MD%20simulations/sfGFP/pdb_exports)
holds every construct as a PDB with hydrogens, crystallographic numbering and CONECT
records for the non-standard residues, plus a tarball for one download.

The report carries inline citations for all 29 sources the setup rests on.
