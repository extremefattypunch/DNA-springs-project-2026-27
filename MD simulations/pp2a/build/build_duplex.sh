#!/bin/bash
# The duplex spring (see the note at the end of build_pp2a.sh): P4_ds50 and D_ds50.
set -euo pipefail
cd "$(dirname "$0")/.."
source ./env/activate.sh
P="$DNASPRING_PY"; SF="$SFGFP_REPO/build"; PARAMS=build/02_params
SITES="7 553"; SPAN=90
frc=("$PARAMS/TET/TET.frcmod" "$PARAMS/TDP/TDP.frcmod" "$PARAMS/DNL/DNL.frcmod" \
     "$PARAMS/DNP/DNP.frcmod" "$PARAMS/DNH/DNH.frcmod" "$PARAMS/junction.frcmod")
prep=("$PARAMS/TET/TET.prepin" "$PARAMS/TDP/TDP.prepin" "$PARAMS/DNL/DNL.prepin" \
      "$PARAMS/DNP/DNP.prepin" "$PARAMS/DNH/DNH.prepin")
bonds_of () { $P -c "import json;print(' '.join(json.load(open('$1'))['tleap_bonds']))"; }

echo "##### duplex spring (50 bp, bent to $SPAN A, both strands)"
$P "$SF/03_dna/build_dna.py" --n-bp 50 --span-A $SPAN --outdir build/03_dna

echo "##### P4_ds50"
$P "$SF/04_assemble/assemble_chimera.py" --protein build/01_protein/pp2a_prepped.pdb \
   --dna build/03_dna/spring_50bp.pdb --sites $SITES --out build/04_assemble/pp2a_ds50.pdb
$P "$SF/04_assemble/build_system.py" --name P4_ds50 \
   --solute build/04_assemble/pp2a_ds50.pdb --outdir build/systems/P4_ds50 \
   --mg-M 0.005 --frcmod "${frc[@]}" --prep "${prep[@]}" \
   --bond $(bonds_of build/04_assemble/pp2a_ds50.json) \
   --expect-residue TDP DNL DNP --expect-loop TDP --clamp-sites $SITES

echo "##### D_ds50  (bare duplex for the force windows)"
$P build/03_dna/bare_strand.py --spring build/03_dna/spring_50bp.pdb --name D_ds50 --pad "${BARE_PAD:-24}"

for d in build/systems/P4_ds50 build/systems/D_ds50; do
  $P -c "
import json
r = json.load(open('$d/build_report.json')); l = r.get('covalent_loop', {})
tag = ('loop ' + ('CLOSED' if l.get('closed') else 'OPEN') + f\" (wanted {l.get('expected')})\") if l else 'no loop test'
print(f\"  {r['system']:<11} {r['atoms']:>7} atoms  box {r['box'][0]:7.2f} A  charge {r['net_charge']:+.4f}  {tag:<32} problems {len(r['problems'])}\")"
done
