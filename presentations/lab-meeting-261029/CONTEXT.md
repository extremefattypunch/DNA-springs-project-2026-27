# What these gels are about (read fully before touching a spec)

Ian Poon's lab meeting (2026-10-29) reports DNA-spring protein chimeras made with
**sfGFP** and the **PP2A A subunit (pp2a)**. The work follows Boral, Schnebly et al.
(bioRxiv 2026.07.05, "Conformational changes induced in ubiquitin by circular
protein-DNA chimeras"; reference figures are in `pages/preprint-*.png`).

## Chemistry
- The protein carries **tetrazine (tet)** unnatural amino acids: **1-tet** (one
  reactive site, the control, like the paper's single-cysteine K63C) or **2-tet**
  (two sites, like the paper's double-cysteine T9C/K63C).
- The ssDNA (50 nt or 80 nt) carries a reactive group at **both** ends: **TCO** or
  **BCN** (click partners of tetrazine). "2-TCO ssDNA" = TCO on both ends.
- The complementary strand ("cDNA", "+cDNA", "hybridized") turns the attached
  ssDNA into a stiff dsDNA "spring".

## Chimera variants = the MAIN bands (paper Fig 2 / S1 / S2)
Top of gel = slow = large. Typical order, slow to fast, in protein stain:
1. multi-protein / multi-DNA species: **2 protein + 2 DNA** (chain or dimeric ring),
   **1 protein + 2 DNA arms** (2-site protein only), **2 protein + 1 DNA** (one strand
   bridging two proteins; the "2sfGFP+DNA" / "2 pp2a" rows)
2. **linear 1 protein + 1 DNA** (one end attached) - "sfGFP+DNA", "pp2a+DNA", "1xpp2a+DNA"
3. **circular 1 protein + 1 DNA** (both ends on the same protein; 2-tet only). Same
   mass as linear but runs slightly FASTER than linear (paper Fig 2).
4. free DNA (faint in protein stain)
5. **protein only** (unreacted) - "sfGFP", "pp2a"; pp2a often shows as a doublet.

## Hybridization (paper Fig 4 / S4) - the arrows
- Adding cDNA shifts chimera bands. **Linear** dsDNA chimeras run **faster** (arrow
  points DOWN). **Circular** dsDNA chimeras run **slower** (arrow points UP, often a
  large shift). An up-shifted band in a 2-tet +cDNA lane = **circular dsDNA chimera**
  = the band of interest that glows.
- On the slides, **arrows go from the unhybridized (ssDNA) band to the hybridized
  (dsDNA) version of the same species**. Record exactly the arrows drawn on the
  image in `arrows` (same count, tail band, head band; no invented arrows): they fix
  which band becomes which. The USER DOES NOT WANT ARROWS IN THE ANIMATION, so the
  engine does not draw them; they appear only on the real gel image at the end.

## Other species that appear on some slides
- GST-tagged pp2a ("uncleaved pp2a", "pp2a+GST"), cleaved pp2a, free GST, free
  thrombin (GST-cleavage gels, thrombin series).
- PEG-capped tet groups ("2 tet pp2a+PEG", "1 tet pp2a+PEG") on the resin slide.
- DNA-stain gels (SYBR) show DNA-containing bands; on the light photo they are dark.

## What to draw
- **Only main bands**: the chimera variants above, protein only, and any species the
  slide explicitly labels. Ladders / marker lanes are NOT drawn.
- Gels crowded with impurity bands: also add a few (at most ~4 per lane) clearly
  fainter unlabelled bands (`"sp": "imp"`, `i` <= 0.2) where the real image has
  them, so the lane does not look cleaner than it is. Skip trivial specks.
- Relative brightness matters: `i` is relative to the brightest band on that gel.
  Saturated / bloated bands get high `i` and a larger `h`.
- Smudged / smeared bands (e.g. DTT optimization: less smear as DTT increases) use
  `smear` (upward/downward tail) and `blur`; the progression must be visible.

## Labels (the audience must always know what each band is)
- Every non-impurity species gets a label with a preprint-style chimera sketch
  (`icon`). Names <= 24 characters, notes <= 34 characters, plain wording:
  "pp2a only", "pp2a + ssDNA", "pp2a + dsDNA", "2 pp2a + 2 DNA", "Circular pp2a + dsDNA".
- Use the protein name the slide uses (sfGFP, pp2a). Keep the slide's own band
  names where it gives them, refined with ss/ds and linear/circular.
- Lane headers: short chip text (<= 7 chars, e.g. "ssDNA", "+cDNA", "pp2a", "0.50 U",
  "1:2", "16 mM") and optional `sub` (<= 8 chars). Group brackets mirror the boxes /
  headings drawn on the image ("1 tet", "2 tet BCN", "50nt", "control"...).

## Order of reveal (the user's rules)
1. Protein-only / no-reaction control lanes first (the reference row).
2. **1-tet unhybridized** lanes.
3. **2-tet unhybridized** lanes.
4. **All remaining hybridized** lanes together.
5. **Glow** the hybridized circular product band(s) (only if the gel has them).
6. `{"real": true}` swaps in the real gel image. Always the last step.
- **Concentration / ratio / time series**: one lane per step from lowest to highest
  value (e.g. thrombin 0 -> 0.01 -> 0.05 -> 0.10 -> 0.25 -> 0.50 U), each lane labelled
  with its value. Combine with the rules above when a gel has both.
- Each step's `say` caption: one sentence <= 95 characters, lead phrase in <b>...</b>.
