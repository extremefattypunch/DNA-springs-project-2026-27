#!/usr/bin/env bash
# Render one gel spec headlessly and build a comparison sheet.
#   tools/render.sh GEL_ID [STEP]
# Writes into deck/_render/GEL_ID/:
#   sim.png      simulated gel, final simulated step (or STEP if given)
#   real.png     same layout with the real image swapped in
#   overlay.png  real image with the simulated bands drawn over it at ~60% opacity
#   compare.png  real / sim / overlay stacked, downscaled (read this one first)
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
ID="$1"; STEP="${2:-}"
OUT="$ROOT/_render/$ID"
python3 "$ROOT/tools/build.py" --test "$ID" >/dev/null
URL="file://$OUT/test.html?test=$ID"
[ -n "$STEP" ] && URL="$URL&step=$STEP"
shot() {
  google-chrome-stable --headless=new --no-sandbox --disable-gpu --hide-scrollbars \
    --force-device-scale-factor=1 --window-size=1920,1080 --virtual-time-budget=6000 \
    --user-data-dir="$OUT/.chrome" --screenshot="$2" "$1" >/dev/null 2>&1
}
shot "$URL&mode=sim" "$OUT/sim.png"
shot "$URL&mode=real" "$OUT/real.png"
shot "$URL&mode=overlay" "$OUT/overlay.png"
magick "$OUT/real.png" "$OUT/sim.png" "$OUT/overlay.png" -resize 1280x720 -append "$OUT/compare.png"
echo "wrote $OUT/{sim,real,overlay,compare}.png"
