#!/bin/bash
# The ubiquitin chimera set (Doeke's points 3 and 7).
#
# Why ubiquitin as well as sfGFP: it is the better-controlled experiment.  The companion
# preprint from this lab measures NMR chemical shifts on this exact construct along a
# perturbation ladder -- apo, crosslinker only, two short strands, ssDNA, dsDNA -- and
# finds most of the effect has already appeared before any duplex exists.  If that is
# right, the DNA's contribution here is chemistry and contact, not spring action, and
# contacts are something 50 ns of MD can measure directly.
#
# What these runs cannot do, stated up front: the beta-5 register shift the preprint is
# ultimately about is a slow, discrete equilibrium, and the displacement a ~4 pN spring
# would produce along it is ~0.1 A -- far below thermal noise and far slower than MD.
# These runs test the delivered force, the local structure at the attachment sites, and
# the protein-DNA contacts.  They do not test the retraction equilibrium.
#
# Five systems, matching the preprint's ladder as closely as our chemistry allows:
#   U0_apo      the F45W/T9C/K63C protein alone
#   U1_clicked  both crosslinkers in place, clicked to a short aminohexanol instead of
#               an oligo.  This is NOT the preprint's DBCO-Ub, which has an unreacted
#               alkyne; it is the matched-chemistry, no-DNA control, the ubiquitin
#               analogue of S2_clicked.
#   U2_ss50     circular ssDNA-Ub
#   U3_ds50     circular dsDNA-Ub
#   U4_ssarms   the same 50 nt cut in half: two 25-nt arms, no circle
set -euo pipefail
cd "$(dirname "$0")/.."
source ./env/activate.sh
export PATH="$DNASPRING_ENV/bin:$PATH" AMBERHOME="$DNASPRING_ENV"
P="$DNASPRING_PY"; PARAMS=build/02_params
SITES="9 63"
# Cys9 CB to Cys63 CB is 24.8 A in 1UBQ and the DBCO arms reach 30-34 A extended, so a
# 60 A anchor span leaves both at 88% of reach -- the same fraction the sfGFP build uses.
SPAN=60
DNA=build/03_dna/ubq
mkdir -p "$DNA"

echo "##### protein"
$P build/01_protein/prepare_ubiquitin.py

echo "##### springs"
$P build/03_dna/build_dna.py --n-bp 50 --span-A $SPAN --outdir "$DNA"
$P build/03_dna/build_dna.py --n-bp 50 --span-A $SPAN --single-strand --outdir "$DNA"
$P build/03_dna/build_dna.py --n-bp 50 --span-A $SPAN --single-strand \
   --nick-after 25 --nick-modified-strand --outdir "$DNA"

frc=("$PARAMS/CMS/CMS.frcmod" "$PARAMS/CDB/CDB.frcmod" "$PARAMS/DNL/DNL.frcmod" \
     "$PARAMS/DNP/DNP.frcmod" "$PARAMS/DNH/DNH.frcmod" "$PARAMS/junction.frcmod")
prep=("$PARAMS/CMS/CMS.prepin" "$PARAMS/CDB/CDB.prepin" "$PARAMS/DNL/DNL.prepin" \
      "$PARAMS/DNP/DNP.prepin" "$PARAMS/DNH/DNH.prepin")

echo "##### U0_apo"
$P build/04_assemble/build_system.py --name U0_apo \
   --solute build/01_protein/ubq_prepped.pdb --outdir build/systems/U0_apo \
   --mg-M 0.005 --frcmod "${frc[@]}" --prep "${prep[@]}" \
   --expect-residue CYS TRP --clamp-sites $SITES

echo "##### U1_clicked"
$P build/04_assemble/assemble_clicked.py --protein build/01_protein/ubq_prepped.pdb \
   --sites $SITES --arm CMS CDB DNH --out build/04_assemble/ubq_clicked.pdb
b=$($P -c "import json;print(' '.join(json.load(open('build/04_assemble/ubq_clicked.json'))['tleap_bonds']))")
$P build/04_assemble/build_system.py --name U1_clicked \
   --solute build/04_assemble/ubq_clicked.pdb --outdir build/systems/U1_clicked \
   --mg-M 0.005 --frcmod "${frc[@]}" --prep "${prep[@]}" --bond $b \
   --expect-residue CMS CDB DNH --clamp-sites $SITES

build_ubq () {   # $1 = name, $2 = spring pdb, $3 = loop assertion flag
  local name=$1 spring=$2 flag=$3 pdb=build/04_assemble/ubq_${1#U?_}.pdb
  $P build/04_assemble/assemble_chimera.py \
     --protein build/01_protein/ubq_prepped.pdb --dna "$spring" --sites $SITES \
     --arm-5prime CMS CDB DNL --arm-3prime CMS CDB DNP --out "$pdb"
  local bonds
  bonds=$($P -c "import json;print(' '.join(json.load(open('${pdb%.pdb}.json'))['tleap_bonds']))")
  $P build/04_assemble/build_system.py --name "$name" --solute "$pdb" \
     --outdir "build/systems/$name" --mg-M 0.005 \
     --frcmod "${frc[@]}" --prep "${prep[@]}" --bond $bonds \
     --expect-residue CMS CDB DNL DNP $flag CMS --clamp-sites $SITES
}
build_ubq U2_ss50   "$DNA/spring_50nt.pdb"         --expect-loop
build_ubq U3_ds50   "$DNA/spring_50bp.pdb"         --expect-loop
build_ubq U4_ssarms "$DNA/spring_50nt_nickmod.pdb" --expect-open-loop

echo "##### done"
for d in build/systems/U?_*/; do
  $P -c "
import json
r = json.load(open('$d/build_report.json')); l = r.get('covalent_loop', {})
tag = ('loop ' + ('CLOSED' if l.get('closed') else 'OPEN') +
       f\" (wanted {l.get('expected')})\") if l else 'no loop test'
print(f\"  {r['system']:<12} {r['atoms']:>7} atoms  charge {r['net_charge']:+.4f}  \"
      f\"{tag:<32} problems {len(r['problems'])}\")"
done
