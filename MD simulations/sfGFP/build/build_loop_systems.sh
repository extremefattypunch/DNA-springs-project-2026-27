#!/bin/bash
# Build the corrected, covalently closed chimeras (Doeke's point 2).
#
# These replace S3-S5 as the model of the experiment.  S3-S5 are kept and re-read as
# what they actually are: the two-arm construct, one linker on each of two different
# strands, which is Zocchi's pre-ligation chimera and the closest analogue we have to
# the non-circular control in the ubiquitin preprint.  The difference between the two
# sets is one covalent bond's worth of topology and nothing else -- same protein, same
# sites, same chemistry, same spring length -- so they compare directly.
set -euo pipefail
cd "$(dirname "$0")/.."
source ./env/activate.sh
export PATH="$DNASPRING_ENV/bin:$PATH" AMBERHOME="$DNASPRING_ENV"
P="$DNASPRING_PY"
PARAMS=build/02_params
SITES="133 149"
# The span between the two linker anchors is a property of the protein and the two
# crosslinkers -- Zocchi's h -- so it is the same for every spring length and only the
# contour changes.  That is what makes the force ladder a controlled comparison.
SPAN=66

echo "##### 1. springs: one strand modified at both ends"
for n in 27 40; do
  $P build/03_dna/build_dna.py --n-bp $n --span-A $SPAN --outdir build/03_dna
done
$P build/03_dna/build_dna.py --n-bp 40 --span-A $SPAN --nick-after 20 --outdir build/03_dna

echo "##### 2. assemble and parameterise"
common_frcmod=("$PARAMS/TDP/TDP.frcmod" "$PARAMS/DNL/DNL.frcmod" \
               "$PARAMS/DNP/DNP.frcmod" "$PARAMS/junction.frcmod")
common_prep=("$PARAMS/TDP/TDP.prepin" "$PARAMS/DNL/DNL.prepin" "$PARAMS/DNP/DNP.prepin")

build_loop () {   # $1 = system name, $2 = spring pdb
  local name=$1 spring=$2 pdb=build/04_assemble/chimera_${1#L?_}.pdb
  $P build/04_assemble/assemble_chimera.py \
     --protein build/01_protein/sfgfp_prepped.pdb --dna "$spring" --out "$pdb"
  local bonds
  bonds=$($P -c "import json;print(' '.join(json.load(open('${pdb%.pdb}.json'))['tleap_bonds']))")
  $P build/04_assemble/build_system.py --name "$name" --solute "$pdb" \
     --outdir "build/systems/$name" --mg-M 0.005 \
     --frcmod "${common_frcmod[@]}" --prep "${common_prep[@]}" \
     --bond $bonds --expect-residue CRO TDP DNL DNP --expect-loop TDP \
     --clamp-sites $SITES
}
build_loop L3_loop27     build/03_dna/spring_27bp.pdb
build_loop L4_loop40     build/03_dna/spring_40bp.pdb
build_loop L5_loop40nick build/03_dna/spring_40bp_nick.pdb

echo "##### done"
for d in build/systems/L?_*/; do
  $P -c "
import json
r = json.load(open('$d/build_report.json'))
loop = r.get('covalent_loop', {})
print(f\"  {r['system']:<16} {r['atoms']:>7} atoms  charge {r['net_charge']:+.4f}  \"
      f\"loop {'CLOSED' if loop.get('closed') else 'OPEN'}  \"
      f\"problems {len(r['problems'])}\")"
done
