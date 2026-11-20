#!/bin/bash -l
#SBATCH -J sfgfp-probe
#SBATCH -p test
#SBATCH -c 1
#SBATCH --mem=1G
#SBATCH -t 0-00:30
#SBATCH -o /n/netscratch/hekstra_lab/Lab/ian_poon/logs/%x-%j.out
#SBATCH -e /n/netscratch/hekstra_lab/Lab/ian_poon/logs/%x-%j.err
# Reachability probe + recursive listing of the facility rsync module (no data copied).
set -uo pipefail
source /n/home05/ianpoon/DNA-springs-project-2026-27/cryoEM/sfgfp_chimera_screening/env/activate_sfgfp.sh
echo "node $(hostname)"
if timeout 10 bash -c '</dev/tcp/10.245.3.142/873'; then echo "TCP 873 reachable"; else echo "TCP 873 UNREACHABLE"; exit 3; fi
out=$SFGFP_DATA/raw_inventory/listing.txt
timeout 1500 rsync --list-only -r --password-file="$SFGFP_RSYNC_PASS" "$SFGFP_RSYNC_URL/" > "$out.tmp" && mv "$out.tmp" "$out"
echo "rsync exit $?"; wc -l "$out"
awk '{s+=$2} END {printf "total bytes %.0f (%.3f TB)\n", s, s/1e12}' <(sed 's/,//g' "$out")
