#!/usr/bin/env python3
"""Inline deck.json + specs/*.json into template.html.

  build.py                      -> deck/index.html (the artifact page)
  build.py --test GEL_ID        -> deck/_render/GEL_ID/test.html (one gel, for render.py)
"""
import argparse, json, pathlib, sys

ROOT = pathlib.Path(__file__).resolve().parent.parent
ap = argparse.ArgumentParser()
ap.add_argument("--test")
a = ap.parse_args()

deck = json.loads((ROOT / "deck.json").read_text())
gels = {}
for p in sorted((ROOT / "specs").glob("*.json")):
    try:
        spec = json.loads(p.read_text())
    except json.JSONDecodeError as e:
        sys.exit(f"{p.name}: invalid JSON: {e}")
    gels[spec["id"]] = spec

def check(spec):
    errs = []
    lanes = {l["id"]: l for l in spec.get("lanes", [])}
    for l in spec.get("lanes", []):
        for i, b in enumerate(l.get("bands", [])):
            if b.get("sp") not in spec.get("species", {}):
                errs.append(f"lane {l['id']} band {i}: species {b.get('sp')!r} not in species map")
    for ar in spec.get("arrows", []):
        for end in ("from", "to"):
            lid, bi = ar[end]
            if lid not in lanes or bi >= len(lanes[lid].get("bands", [])):
                errs.append(f"arrow {end} {ar[end]} does not name an existing lane/band")
    seen = set()
    for st in spec.get("steps", []):
        for lid in st.get("lanes", []):
            if lid not in lanes:
                errs.append(f"step lane {lid!r} unknown")
            seen.add(lid)
        for g in st.get("glow", []):
            if g[0] not in lanes or g[1] >= len(lanes[g[0]].get("bands", [])):
                errs.append(f"glow {g} does not name an existing band")
    for lid in lanes:
        if lid not in seen and lanes[lid].get("bands"):
            errs.append(f"lane {lid} has bands but is never shown in a step")
    if not spec.get("steps") or not spec["steps"][-1].get("real"):
        errs.append("last step must be {\"real\": true}")
    return errs

bad = False
for gid, spec in list(gels.items()):
    if a.test and gid != a.test:
        continue
    errs = check(spec)
    for e in errs:
        print(f"[{gid}] {e}", file=sys.stderr)
    if errs:
        if a.test:
            bad = True
        else:
            print(f"[{gid}] skipped in the full build", file=sys.stderr)
            del gels[gid]
if bad:
    sys.exit("spec errors (see above)")

charts = {}
for p in sorted((ROOT / "chartspecs").glob("*.json")):
    c = json.loads(p.read_text()); charts[c["id"]] = c

tpl = (ROOT / "template.html").read_text()
if a.test:
    if a.test not in gels:
        sys.exit(f"no spec with id {a.test}")
    payload = {"slides": deck["slides"], "gels": {a.test: gels[a.test]}}
    out = ROOT / "_render" / a.test / "test.html"
    out.parent.mkdir(parents=True, exist_ok=True)
    html = tpl.replace("/*__DECK__*/null", json.dumps(payload, separators=(",", ":")))
    html = f'<base href="{ROOT.as_uri()}/">\n' + html
else:
    payload = {"slides": deck["slides"], "gels": gels, "charts": charts}
    out = ROOT / "index.html"
    html = tpl.replace("/*__DECK__*/null", json.dumps(payload, separators=(",", ":")))
out.write_text(html)
print(f"wrote {out}")
