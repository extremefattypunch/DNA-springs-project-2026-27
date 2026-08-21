# Viewable structures

Every construct in the study, as coordinates you can open anywhere. Written by
`build/export_pdbs.py`.

## What is here

| file | what it is |
|---|---|
| `S0_wt_*.pdb` | wild-type sfGFP |
| `S1_tet_*.pdb` | two Tet2-Et residues at 133/149, unclicked — the fluorometry sample |
| `S2_clicked_*.pdb` | both tethers clicked, no DNA — the unloaded control |
| `S3_spring27_*.pdb` | the chimera with a 27 bp spring (~4.9 pN) |
| `S4_spring40_*.pdb` | the chimera with a 40 bp spring (~3.0 pN) |
| `S5_spring40nick_*.pdb` | the same, nicked at the centre (~2.4 pN) |
| `dna_spring_*.pdb` | the isolated duplexes, bent, before assembly |
| `residue_*_model.mol2` | the four custom residues as parameterised, with bond orders |

`_start` is the equilibrated structure at the beginning of production; `_final` is the
last frame of the longest run. Solvent is stripped — the solvated systems are 44k–121k
atoms and live in `build/systems/*/system.prmtop` with the trajectories on netscratch.

## Things worth knowing before you look

Residues use **2B3P numbering**, so the attachment sites are `Asp133` and `Asn149`
(D134/N150 in the construct's own numbering, which runs one higher). Chains: `A` the
protein, `C`/`D` the two DNA strands, `E`/`F` the capped tethers in `S2_clicked`.
The linker residues are `TDP` (the clicked Tet2-Et/sTCO adduct, part of the protein
chain) and `DNL` (the amino-C6 arm, leading each DNA strand); `TET` is the unclicked
tetrazine and `DNH` the tether capped as a free alcohol. `CRO` is the chromophore.

CONECT records are included for every bond involving those residues. Keep them — no
viewer has a template for the fused bicyclononane cage, and distance-based bond
guessing mis-draws it.

## Opening them

PyMOL:

    pymol S3_spring27_final.pdb
    # then, to see the construct the way the report shows it:
    hide everything
    show cartoon, polymer.protein
    show cartoon, polymer.nucleic
    set cartoon_ring_mode, 3
    show sticks, resn CRO+TDP+DNL
    color grey80, polymer.protein
    color skyblue, polymer.nucleic
    color orange, resn TDP+DNL
    color limegreen, resn CRO
    orient

ChimeraX:

    open S3_spring27_final.pdb
    cartoon; nucleotides stubs
    show :CRO,TDP,DNL atoms; style :CRO,TDP,DNL stick
    color /A grey; color /C,D cornflowerblue; color :TDP,DNL orange; color :CRO green

To measure the spring the way the analysis does — the distance between the two 5′
phosphates, which is what sets the force:

    # PyMOL
    distance span, chain C and resi 1 and name P, chain D and resi 1 and name P

And the deformation coordinate:

    distance sites, chain A and resi 133 and name CB, chain A and resi 149 and name CB
