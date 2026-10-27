#!/usr/bin/env bash
# Every step of one slide of the built deck: tools/deckshots.sh SLIDE_NUMBER -> _render/deck/sN.png
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"; N="$1"; OUT="$ROOT/_render/deck"; mkdir -p "$OUT"
K=$(python3 - "$ROOT" "$N" <<'PY'
import json, sys, pathlib
root, n = pathlib.Path(sys.argv[1]), int(sys.argv[2])
deck = json.loads((root / 'deck.json').read_text()); sl = deck['slides'][n - 1]
k = 1
for p in sl.get('panels', []):
    k += len(json.loads((root / 'specs' / (p['gel'] + '.json')).read_text())['steps'])
print(k)
PY
)
files=()
for k in $(seq 0 $((K-1))); do
  google-chrome-stable --headless=new --no-sandbox --disable-gpu --hide-scrollbars --force-device-scale-factor=1 \
    --window-size=1920,1080 --virtual-time-budget=6000 --user-data-dir="$OUT/.chrome-$N" \
    --screenshot="$OUT/s$N-$k.png" "file://$ROOT/index.html?shot=$N&step=$k" >/dev/null 2>&1
  files+=("$OUT/s$N-$k.png")
done
magick montage "${files[@]}" -tile 3x -geometry 640x360+4+4 -background '#777' "$OUT/s$N.png"
echo "slide $N: $K steps -> $OUT/s$N.png"
