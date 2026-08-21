#!/bin/bash
# Build every simulated system from the raw 2B3P deposition, end to end.
#
# S6 (the force-clamp ladder) deliberately reuses S1's topology: clamping the same
# two Tet2-Et residues the spring pulls on is what makes the calibration curve and
# the chimera directly comparable.
set -euo pipefail
cd "$(dirname "$0")/.."
source ./env/activate.sh
export PATH="$DNASPRING_ENV/bin:$PATH" AMBERHOME="$DNASPRING_ENV"
P="$DNASPRING_PY"
PARAMS=build/02_params
SITES="133 149"

echo "##### 1. protein preparation"
$P build/01_protein/prepare_protein.py

echo "##### 2. custom residue parameters"
$P $PARAMS/build_residues.py

echo "##### 3. DNA springs"
for n in 27 40; do
  $P build/03_dna/build_dna.py --n-bp $n --span-A 62 --outdir build/03_dna
done
$P build/03_dna/build_dna.py --n-bp 40 --span-A 62 --nick-after 20 --outdir build/03_dna

echo "##### 4. assemble and parameterise"
common_frcmod=("$PARAMS/TDP/TDP.frcmod" "$PARAMS/DNL/DNL.frcmod" "$PARAMS/junction.frcmod")
common_prep=("$PARAMS/TDP/TDP.prepin" "$PARAMS/DNL/DNL.prepin")

# S0: wild type
$P build/04_assemble/build_system.py --name S0_wt \
   --solute build/01_protein/sfgfp_prepped.pdb --outdir build/systems/S0_wt \
   --clamp-sites $SITES

# S1: two unclicked Tet2-Et residues (matches the 2-tet fluorometry sample), and the
#     topology the force-clamp ladder runs on
$P build/04_assemble/mutate_residue.py --protein build/01_protein/sfgfp_prepped.pdb \
   --residue TET --sites $SITES --out build/04_assemble/sfgfp_tet.pdb
$P build/04_assemble/build_system.py --name S1_tet \
   --solute build/04_assemble/sfgfp_tet.pdb --outdir build/systems/S1_tet \
   --frcmod "$PARAMS/TET/TET.frcmod" --prep "$PARAMS/TET/TET.prepin" \
   --expect-residue CRO TET --clamp-sites $SITES
mkdir -p build/systems/S6_clamp
for f in system.prmtop system.inpcrd residue_map.json clamp_atoms.txt build_report.json; do
  cp build/systems/S1_tet/$f build/systems/S6_clamp/$f
done

# S3-S5: the chimeras
build_chimera () {   # $1 = system name, $2 = spring pdb
  local name=$1 spring=$2 pdb=build/04_assemble/chimera_${name}.pdb
  $P build/04_assemble/assemble_chimera.py \
     --protein build/01_protein/sfgfp_prepped.pdb --dna "$spring" --out "$pdb"
  local bonds
  bonds=$($P -c "import json,sys;print(' '.join(json.load(open('${pdb%.pdb}.json'))['tleap_bonds']))")
  $P build/04_assemble/build_system.py --name "$name" --solute "$pdb" \
     --outdir "build/systems/$name" --mg-M 0.005 \
     --frcmod "${common_frcmod[@]}" --prep "${common_prep[@]}" \
     --bond $bonds --expect-residue CRO TDP DNL --clamp-sites $SITES
}
build_chimera S3_spring27     build/03_dna/spring_27bp.pdb
build_chimera S4_spring40     build/03_dna/spring_40bp.pdb
build_chimera S5_spring40nick build/03_dna/spring_40bp_nick.pdb

echo "##### done"
for d in build/systems/*/; do
  $P -c "
import json,sys
r=json.load(open('$d/build_report.json'))
print(f\"  {r['system']:<16} {r['atoms']:>7} atoms  charge {r['net_charge']:+.4f}  \"
      f\"waters {r['waters']:>6}  problems {len(r['problems'])}\")"
done
