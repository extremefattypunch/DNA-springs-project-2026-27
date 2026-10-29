#!/bin/bash
# Render one cryoEM explainer frame at deck size: tools/cryoshot.sh PAGE STEP [OUT]
D=$(cd "$(dirname "$0")/.." && pwd); P=$1; K=$2; O=${3:-$D/_render/cryo/${P%.html}-$K.png}
mkdir -p "$(dirname "$O")"
google-chrome-stable --headless=new --no-sandbox --disable-gpu --hide-scrollbars --window-size=1178,748 \
  --virtual-time-budget=4000 --screenshot="$O" "file://$D/interactive/cryo/$P?step=$K&static" 2>/dev/null >/dev/null
echo "$O"
