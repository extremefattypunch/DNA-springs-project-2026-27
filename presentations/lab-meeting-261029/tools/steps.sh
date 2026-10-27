#!/usr/bin/env bash
# Contact sheet of every step of one gel: tools/steps.sh GEL_ID  ->  _render/GEL_ID/steps.png
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"; ID="$1"; OUT="$ROOT/_render/$ID"
python3 "$ROOT/tools/build.py" --test "$ID" >/dev/null
N=$(python3 -c "import json,sys;print(len(json.load(open('$ROOT/specs/$ID.json'))['steps']))")
files=()
for k in $(seq 0 "$N"); do
  google-chrome-stable --headless=new --no-sandbox --disable-gpu --hide-scrollbars --force-device-scale-factor=1 \
    --window-size=1920,1080 --virtual-time-budget=6000 --user-data-dir="$OUT/.chrome" \
    --screenshot="$OUT/step-$k.png" "file://$OUT/test.html?test=$ID&mode=sim&step=$k" >/dev/null 2>&1
  files+=("$OUT/step-$k.png")
done
magick montage "${files[@]}" -tile 2x -geometry 960x540+6+6 -background '#888' "$OUT/steps.png"
echo "wrote $OUT/steps.png ($((N+1)) steps)"
