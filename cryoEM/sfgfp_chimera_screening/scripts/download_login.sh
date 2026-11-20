#!/bin/bash
# Pull the 12-11-2026 sfGFP chimera screening session from the facility rsync server.
#
# Must run on a LOGIN node: Cannon compute nodes have no route to 10.245.3.142:873
# (probe job 60306535, 2026-11-13). It is I/O-only, niced, resumable and retries.
# Skips EPU's per-movie single-frame preview .mrc (311 GB, unused by processing).
#
#   ssh holylogin05.rc.fas.harvard.edu
#   setsid nohup bash .../scripts/download_login.sh > /dev/null 2>&1 &
set -uo pipefail
source /n/home05/ianpoon/DNA-springs-project-2026-27/cryoEM/sfgfp_chimera_screening/env/activate_sfgfp.sh >/dev/null 2>&1 || true
DEST=${SFGFP_RAW:-/n/netscratch/hekstra_lab/Lab/ian_poon/sfgfp_chimera_screening_20261112/raw}
LOG=${SFGFP_LOGS:-/n/netscratch/hekstra_lab/Lab/ian_poon/logs}/sfgfp-download-$(hostname -s)-$(date +%Y%m%dT%H%M%S).log
URL=rsync://em00000013@10.245.3.142/em00000013/
PASS=$HOME/.config/rsync/em00000013.pass
mkdir -p "$DEST"
echo $$ > "$DEST/../download.pid"
FILTER=(--include='*_Fractions.mrc' --exclude='**/Data/*.mrc')
{
  echo "start $(date -Is) on $(hostname) -> $DEST"
  for try in $(seq 1 20); do
    nice -n 10 ionice -c2 -n7 rsync -rlt --chmod=D750,F640 --partial-dir=.rsync-partial \
      --timeout=600 --info=stats2,progress2 --no-inc-recursive \
      "${FILTER[@]}" --password-file="$PASS" "$URL" "$DEST/"
    rc=$?
    echo "try $try rc=$rc $(date -Is)"
    [ $rc -eq 0 ] && { echo "DONE $(date -Is)"; break; }
    sleep 60
  done
} >> "$LOG" 2>&1
rm -f "$DEST/../download.pid"
