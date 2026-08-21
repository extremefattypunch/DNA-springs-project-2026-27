# All-atom MD of the sfGFP–DNA spring chimera

`index.html` is the full report — self-contained, with every figure and animation
embedded, and inline citations for everything the setup rests on.

**Read it here:** https://extremefattypunch.github.io/DNA-springs-project-2026-27/

If that link 404s, GitHub Pages has not been switched on yet. One click:
**Settings → Pages → Build and deployment → Source: “Deploy from a branch” →
Branch: `main`, folder `/docs` → Save.** It goes live in about a minute.

Until then the same file renders through a proxy, no setup needed:
https://htmlpreview.github.io/?https://raw.githubusercontent.com/extremefattypunch/DNA-springs-project-2026-27/main/docs/index.html

The code, structures and analysis live in [`MD simulations/sfGFP`](../MD%20simulations/sfGFP);
viewable PDBs of every construct are in
[`MD simulations/sfGFP/pdb_exports`](../MD%20simulations/sfGFP/pdb_exports).
`figures/` here holds the same figures as loose PNGs, for slides.

Regenerate everything with `bash "MD simulations/sfGFP/finalize.sh"`, then
`cp "MD simulations/sfGFP/report/sfgfp_dna_spring_report.html" docs/index.html`.
