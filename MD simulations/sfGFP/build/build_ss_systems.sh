#!/bin/bash
# The ssDNA spring (Doeke's point 7) and the non-circular control (point 3).
#
# One 50-nt strand, clicked at both ends, is the ssDNA chimera.  Cut that same strand
# in the middle and it becomes two 25-nt arms, each attached at one end -- a chimera
# that cannot close.  Building the control this way holds the DNA content, the sequence,
# the sites and the chemistry exactly fixed and changes only the topology, which is a
# tighter comparison than the two 12-mers used in the ubiquitin preprint and is Zocchi's
# own recommended design: change the mechanics, hold every contact constant.
#
# 50 nt, not 27: at 0.255 nm/base a 27-mer has ~6.9 nm of contour against the 6.6 nm
# the sites plus linkers impose, so it would be taut and would report stretching rather
# than the slack, zero-stress reference the control is for.  50 nt is also the length
# used in the ubiquitin work, which makes the two studies comparable.
set -euo pipefail
cd "$(dirname "$0")/.."
source ./env/activate.sh
export PATH="$DNASPRING_ENV/bin:$PATH" AMBERHOME="$DNASPRING_ENV"
P="$DNASPRING_PY"; PARAMS=build/02_params; SITES="133 149"

$P build/03_dna/build_dna.py --n-bp 50 --span-A 66 --single-strand --outdir build/03_dna
$P build/03_dna/build_dna.py --n-bp 50 --span-A 66 --single-strand \
   --nick-after 25 --nick-modified-strand --outdir build/03_dna

frc=("$PARAMS/TDP/TDP.frcmod" "$PARAMS/DNL/DNL.frcmod" "$PARAMS/DNP/DNP.frcmod" \
     "$PARAMS/junction.frcmod")
prep=("$PARAMS/TDP/TDP.prepin" "$PARAMS/DNL/DNL.prepin" "$PARAMS/DNP/DNP.prepin")

build_one () {   # $1 = name, $2 = spring pdb, $3 = loop flag
  local name=$1 spring=$2 flag=$3 pdb=build/04_assemble/${1#C?_}.pdb
  $P build/04_assemble/assemble_chimera.py \
     --protein build/01_protein/sfgfp_prepped.pdb --dna "$spring" --out "$pdb"
  local bonds
  bonds=$($P -c "import json;print(' '.join(json.load(open('${pdb%.pdb}.json'))['tleap_bonds']))")
  $P build/04_assemble/build_system.py --name "$name" --solute "$pdb" \
     --outdir "build/systems/$name" --mg-M 0.005 \
     --frcmod "${frc[@]}" --prep "${prep[@]}" \
     --bond $bonds --expect-residue CRO TDP DNL DNP $flag TDP --clamp-sites $SITES
}
build_one C1_ss50    build/03_dna/spring_50nt.pdb         --expect-loop
build_one C2_ssarms  build/03_dna/spring_50nt_nickmod.pdb --expect-open-loop

for d in build/systems/C?_*/; do
  $P -c "
import json
r = json.load(open('$d/build_report.json')); l = r.get('covalent_loop', {})
print(f\"  {r['system']:<12} {r['atoms']:>7} atoms  charge {r['net_charge']:+.4f}  \"
      f\"loop {'CLOSED' if l.get('closed') else 'OPEN':6s} (wanted {l.get('expected')})  \"
      f\"problems {len(r['problems'])}\")"
done
