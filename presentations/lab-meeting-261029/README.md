# Lab Meeting 261029 deck

Source for the animated deck. `index.html` is generated, do not edit it.

- `template.html`: deck engine (stage, navigation, gel renderer, chimera sketches)
- `deck.json`: slide list. `kind: "image"` slides show the original PDF page
  (`slides/p-NN.jpg`); `kind: "gel"` slides draw one or two gel panels from specs.
- `specs/<id>.json`: one simulated gel each (format: `SPEC.md`, science + reveal
  rules: `CONTEXT.md`); `gels/<id>.jpg`: the real gel image swapped in at the end.
- `tools/build.py` -> `index.html`; `tools/render.sh ID`, `tools/steps.sh ID` render
  checks into `_render/ID/`; `tools/play.mjs` drives the built deck in headless Chrome.

cryoEM section (last 14 slides): `cryo/build_cryo.py` writes it into `deck.json`
(replacing PDF pages 25-28). Ian's Manim clips are copied into `cryo/`; a media point
with `seek: [a, b]` plays that beat of the clip once and holds. The optics and
reconstruction slides are iframes in `interactive/cryo/` driven by `{deck: 'cryo', step}`;
`tools/cryoshot.sh PAGE STEP` renders one frame.

Navigation: right arrow / space / page down / click = next step, left arrow /
page up = back, f = full screen, `#sN` in the link opens slide N.

To replace another original slide with a built one, change its entry in
`deck.json` from `kind: "image"` to a new kind handled in `template.html`.

Published as the Claude artifact https://claude.ai/artifact/3WwWoCorh1GgpnCWQCaTWn.
To update it, run `python3 tools/build.py`, then publish `index.html` to that URL
with `gels/*.jpg` and the `slides/p-NN.jpg` files that `deck.json` uses as
supporting files at the same relative paths. `index.html` and `_render/` are
generated and git-ignored. `slides/` holds 192 dpi renders of
"Lab Meeting Ian 261029.pdf" (`pdftoppm -r 192 -jpeg`).
