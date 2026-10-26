#!/bin/bash
# The PP2A A-alpha DNA-spring set.  One protein, one spring, three systems plus the
# bare strand -- the 50 nt ssDNA chimera and the two controls the sfGFP review asked for.
#
#   P1_tet      A-alpha with Tet2-Et at Asp7 and Asn553, unclicked (the 2-tet sample)
#   P2_clicked  both tethers clicked to an inert aminohexanol, no DNA: the matched
#               zero-force reference (sfGFP's S2_clicked).  A chimera is compared to
#               THIS, never to the unmodified protein -- clicking alone moved sfGFP's
#               sites by 0.65 A, far more than any spring could.
#   P3_ss50     the chimera: one 50 nt strand clicked at both ends, a covalent circle
#   D_ss50      the strand alone, for the restraint-window force measurement
#
# Sites are UniProt Asp7/Asn553 (= construct D6X/N552X; see 01_protein/prepare_pp2a.py).
# The apo baseline is the cryo-EM reproduction's finished 100 ns run of the same AF
# model (ff19SB/OPC), so no P0 is built here.
#
# Everything chemical is the sfGFP machinery, called by path: parameters in
# ../sfGFP/build/02_params (symlinked as build/02_params), the DNA builder, the
# assemblers and the tleap driver.  Nothing in them is protein-specific.
set -euo pipefail
cd "$(dirname "$0")/.."
source ./env/activate.sh
P="$DNASPRING_PY"; SF="$SFGFP_REPO/build"; PARAMS=build/02_params
SITES="7 553"
# Asp7 CB to Asn553 CB is 54.4 A in the AF model and the tethers reach 22.4 / 26.1 A
# extended.  A 90 A anchor span puts both arms at 88% of reach with an 11.7 A lift off
# the surface -- the same fractions as the sfGFP (66 A) and ubiquitin (60 A) builds.
# The assembler's own limit is 97 A; 90 was chosen, not measured, and the strand is free
# to pull its anchors closer once released (sfGFP's relaxed from 66 to 48-65 A).
SPAN=90

echo "##### protein"
$P build/01_protein/prepare_pp2a.py

echo "##### spring  (50 nt, bent to $SPAN A, single strand, the sfGFP/ubiquitin 50-mer)"
$P "$SF/03_dna/build_dna.py" --n-bp 50 --span-A $SPAN --single-strand --outdir build/03_dna

frc=("$PARAMS/TET/TET.frcmod" "$PARAMS/TDP/TDP.frcmod" "$PARAMS/DNL/DNL.frcmod" \
     "$PARAMS/DNP/DNP.frcmod" "$PARAMS/DNH/DNH.frcmod" "$PARAMS/junction.frcmod")
prep=("$PARAMS/TET/TET.prepin" "$PARAMS/TDP/TDP.prepin" "$PARAMS/DNL/DNL.prepin" \
      "$PARAMS/DNP/DNP.prepin" "$PARAMS/DNH/DNH.prepin")
bonds_of () { $P -c "import json;print(' '.join(json.load(open('$1'))['tleap_bonds']))"; }

echo "##### P1_tet"
$P "$SF/04_assemble/mutate_residue.py" --protein build/01_protein/pp2a_prepped.pdb \
   --residue TET --sites $SITES --out build/04_assemble/pp2a_tet.pdb
$P "$SF/04_assemble/build_system.py" --name P1_tet \
   --solute build/04_assemble/pp2a_tet.pdb --outdir build/systems/P1_tet \
   --mg-M 0.005 --frcmod "${frc[@]}" --prep "${prep[@]}" \
   --expect-residue TET --clamp-sites $SITES

echo "##### P2_clicked"
$P "$SF/04_assemble/assemble_clicked.py" --protein build/01_protein/pp2a_prepped.pdb \
   --sites $SITES --arm TDP DNH --out build/04_assemble/pp2a_clicked.pdb
$P "$SF/04_assemble/build_system.py" --name P2_clicked \
   --solute build/04_assemble/pp2a_clicked.pdb --outdir build/systems/P2_clicked \
   --mg-M 0.005 --frcmod "${frc[@]}" --prep "${prep[@]}" \
   --bond $(bonds_of build/04_assemble/pp2a_clicked.json) \
   --expect-residue TDP DNH --clamp-sites $SITES

echo "##### P3_ss50"
$P "$SF/04_assemble/assemble_chimera.py" --protein build/01_protein/pp2a_prepped.pdb \
   --dna build/03_dna/spring_50nt.pdb --sites $SITES --out build/04_assemble/pp2a_ss50.pdb
$P "$SF/04_assemble/build_system.py" --name P3_ss50 \
   --solute build/04_assemble/pp2a_ss50.pdb --outdir build/systems/P3_ss50 \
   --mg-M 0.005 --frcmod "${frc[@]}" --prep "${prep[@]}" \
   --bond $(bonds_of build/04_assemble/pp2a_ss50.json) \
   --expect-residue TDP DNL DNP --expect-loop TDP --clamp-sites $SITES

echo "##### D_ss50  (bare strand for the force windows)"
$P build/03_dna/bare_strand.py --spring build/03_dna/spring_50nt.pdb --name D_ss50 \
   --pad "${BARE_PAD:-24}"

echo "##### done"
for d in build/systems/*/; do
  $P -c "
import json
r = json.load(open('$d/build_report.json')); l = r.get('covalent_loop', {})
tag = ('loop ' + ('CLOSED' if l.get('closed') else 'OPEN') +
       f\" (wanted {l.get('expected')})\") if l else 'no loop test'
box = r.get('box', [None])[0]
print(f\"  {r['system']:<11} {r['atoms']:>7} atoms  box {box:7.2f} A  charge {r['net_charge']:+.4f}  \"
      f\"{tag:<32} problems {len(r['problems'])}\")"
done
