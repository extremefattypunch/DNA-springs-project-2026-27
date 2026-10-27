# Gel spec format (`specs/<id>.json`)

All coordinates are **pixels of the gel image file** (`gels/<id>.jpg`), so the
simulated bands land exactly on the real bands when the real image is swapped in.
The worked example is `specs/p05-sfgfp.json` (copy its style).

```jsonc
{
  "id": "p15-pp2a-wrong",            // = file name, = gels/<id>.jpg
  "image": "gels/p15-pp2a-wrong.jpg",
  "size": [728, 380],                // image width, height (identify / PIL)
  "style": "dark",                   // "dark": bright bands on dark gel (fluorescence)
                                     // "mid": bright bands on a mid-grey photo
                                     // "light": dark bands on a light photo (DNA stain)
  "alt": "one sentence describing the real gel",
  "gel": {
    "poly": [[x,y],[x,y],[x,y],[x,y]], // gel slab outline TL, TR, BR, BL (only the gel,
                                       // not labels/ladder panels/background)
    "wellY": 20,                      // y of the well line (top of lanes) at the gel's centre x
    "tilt": 0,                        // degrees, clockwise positive, if the photo is rotated
    "cy": 190                         // optional: y at which lane x values are measured
                                      // (default image mid-height); matters only with tilt
  },
  "species": {                       // one entry per band identity used below
    "protein": { "name": "pp2a only", "note": "unreacted protein", "cls": "protein",
                 "icon": { "k": "protein", "tet": 1 } },
    "p1d_ss":  { "name": "pp2a + ssDNA", "note": "linear 1 : 1", "cls": "ss",
                 "icon": { "k": "p1d" } },
    "circ_ds": { "name": "Circular pp2a + dsDNA", "note": "runs slower · product of interest",
                 "cls": "circ", "icon": { "k": "circ", "ds": true } },
    "imp":     { "cls": "impurity" }  // unlabelled faint bands
    // optional "side": "left" | "right" forces the label rail; normally automatic
  },
  "lanes": [
    { "id": "L1", "x": 267, "w": 32,          // lane centre x (at gel.cy) and band width
      "label": "ssDNA", "chip": "ss",         // header text; chip: "ss" | "ds" | "ctrl" | omit
      "sub": "TCO",                           // optional 2nd header line
      "hl": [x, y, w, h],                     // optional rect (image px) highlighted while this
                                              // lane is the newest one, e.g. a table column
      "bands": [
        { "y": 145, "h": 12, "i": 0.74, "sp": "p1d_ss" },
        // y: band centre; h: thickness (about 1.3-1.5 x the measured FWHM);
        // i: brightness 0-1 relative to the brightest band on this gel;
        // optional: "w" width override, "x" centre override,
        //           "blur": extra px of softness,
        //           "smear": { "up": px, "down": px, "k": 0-1 } streak above/below
        { "y": 207, "h": 14, "i": 0.36, "sp": "protein" }
      ] }
  ],
  "groups": [ { "label": "1 tet", "lanes": ["L1", "L2", "L3"] } ],  // bracket over lanes
  "arrows": [ { "from": ["L2", 0], "to": ["L3", 1] } ],            // [lane id, band index]; recorded, not drawn
  "pins": [ [x, y, w, h] ],          // optional: parts of the real image shown from the start
                                     // (e.g. a conditions table that is not part of the gel)
  "blank": "caption for the empty-gel step",
  "steps": [
    { "lanes": ["L1"], "say": "<b>pp2a alone</b>: ..." },
    { "lanes": ["L2", "L4"], "say": "..." },
    { "glow": [["L6", 0]], "say": "..." },
    { "real": true, "say": "Real gel, as imaged" }
  ]
}
```

## Icons (`species.*.icon`) - preprint-style chimera sketches
`k` is one of:
- `protein` - protein alone; `tet`: 0 (no handle drawn), 1 or 2 handles
- `dna` - free DNA strand
- `p1d` - linear chimera, DNA attached at one end (`tet: 2` also shows the unused site)
- `circ` - circular chimera, both DNA ends on the same protein
- `p2d1` - one DNA strand bridging two proteins
- `p1d2` - one protein carrying two DNA arms
- `p2d2` - chain protein-DNA-protein-DNA
- `p2d2circ` - two proteins + two strands closed into one ring
- `gstp` - protein still carrying its GST tag; `gst` - free GST; `thrombin`
- `pegp` - protein with PEG-capped tet site(s) (`tet` 1 or 2)
Add `"ds": true` for the hybridized (dsDNA) form: a violet complementary strand is drawn.

## Colour classes (`species.*.cls`)
`protein` grey-white, `ss` teal (ssDNA chimeras), `ds` violet (dsDNA chimeras),
`circ` amber (the band of interest that glows: circular dsDNA product, or a candidate circular ssDNA "new band"), `other` pink (GST,
thrombin, PEG...), `impurity` dim grey, no label.

## Tools (run from the deck folder)
- `python3 tools/grid.py gels/X.jpg _render/X-grid.png --box x y w h --scale 4 --step 10`
  zoomed crop with a labelled pixel grid in ORIGINAL image coordinates. Read it.
- `python3 tools/profile.py gels/X.jpg cols --y0 A --y1 B` column means (lane centres)
- `python3 tools/profile.py gels/X.jpg lanes --xs 120,150,180 --w 12 --y0 A --y1 B [--invert] [--tilt T --cy C]`
  band peaks per lane with relative intensities comparable across the listed lanes.
- `tools/render.sh X` builds `_render/X/test.html` and screenshots it:
  `_render/X/compare.png` (real / simulated / overlay). In `overlay.png` the simulated
  bands are drawn over the real image: every simulated band must sit on its real band.
- `tools/steps.sh X` -> `_render/X/steps.png`, one frame per step: check the reveal order,
  labels (no overlaps, nothing clipped), arrows, glow.
- `python3 tools/build.py --test X` validates the spec (unknown species, bad arrow
  indices, lanes never shown, last step not real).
Do NOT edit `template.html`, `deck.json` or other gels' specs. If the engine cannot
express something, say so in your result.
