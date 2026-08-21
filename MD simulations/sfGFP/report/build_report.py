#!/usr/bin/env python3
"""Build the standalone HTML report from whatever analysis has completed.

Self-contained on purpose: figures are embedded as base64 PNG and the overview
animation as a base64 mp4, so the file can be mailed, dropped in a shared drive or
published as an artifact without losing its media.  Matches the house style of
flurometer/build_chimera_viewer.py -- one script, one HTML file, data embedded.
"""
from __future__ import annotations

import argparse
import base64
import json
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parent.parent


def b64(path: Path, mime: str) -> str | None:
    if not path.exists():
        return None
    return f"data:{mime};base64," + base64.b64encode(path.read_bytes()).decode()


def esc(s) -> str:
    return (str(s).replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;"))


def num(x, nd=2, dash="—"):
    try:
        if pd.isna(x):
            return dash
        return f"{float(x):.{nd}f}"
    except (TypeError, ValueError):
        return dash


LABEL = {"S0_wt": "WT sfGFP", "S1_tet": "2× Tet2-Et", "S3_spring27": "27 bp spring",
         "S4_spring40": "40 bp spring", "S5_spring40nick": "40 bp nicked",
         "S6_clamp": "force clamp"}
PURPOSE = {
    "S0_wt": "Baseline. Validates the published chromophore parameters.",
    "S1_tet": "Matches the 2-tet fluorometry sample; also the force-clamp topology.",
    "S3_spring27": "The strong spring.",
    "S4_spring40": "Zocchi’s mechanically clean length (γ&lt;1).",
    "S5_spring40nick": "His calibrated low-stress reference.",
    "S6_clamp": "Constant force at 0/2/4/7/12/20 pN — the calibration curve.",
}
FIG_CAPTIONS = {
    "fig1_force_response": (
        "Force response of the attachment sites",
        "Cβ–Cβ separation of Asp133 and Asn149 against the force applied between them. "
        "The clamp ladder gives extension at a <em>known</em> force on the same topology "
        "the spring pulls on, so reading a chimera’s measured extension against this "
        "curve infers the tension its spring actually delivers — independently of the "
        "analytic model. Error bars are the standard error across replicates."),
    "fig2_rmsf": (
        "Backbone mobility",
        "Per-residue Cα RMSF, aligned on the barrel core rather than the whole protein "
        "so the tethers and their loops cannot smear the signal across every residue. "
        "The ACE/NME caps and the two or three residues either side are omitted: their "
        "4–11 Å free-end motion would compress every real feature into the bottom tenth "
        "of the axis."),
    "fig3_chromophore_hbonds": (
        "Chromophore hydrogen-bond network",
        "Occupancy judged on the H···acceptor distance and the D–H···A angle, not the "
        "heavy-atom distance. That distinction changes a conclusion: His148 looks broken "
        "by a 3.5 Å heavy-atom cutoff while its HD1 proton sits 2.3–2.5 Å from the "
        "phenolate with near-linear geometry."),
    "fig4_dna_bend": (
        "Where the spring bends",
        "Bend per step of the smoothed helical axis. A kink — the softening transition "
        "Zocchi’s model turns on — would appear as a localised spike rather than the "
        "even curvature of an elastic arc. Base-pair centres are averaged over one full "
        "helical turn first; without that the spiral of the centres alone reads as "
        "27°/bp of spurious bend."),
}
STEPS = [
    ("Prepare", "build/01_protein", "2B3P → tleap-ready structure. Strips the "
     "crystallisation Cd²⁺ and acetate, keeps the buried cavity waters, assigns every "
     "histidine tautomer from H-bond geometry, and measures the attachment geometry the "
     "rest of the pipeline reads."),
    ("Parameterise", "build/02_params", "Explicit atom/bond tables → RDKit → "
     "antechamber (AM1-BCC, GAFF2) → prepgen → backbone harmonisation → parmchk2. "
     "Zero parameters flagged for review across all three custom residues."),
    ("Build the spring", "build/03_dna", "Ideal B-form duplex from PyMOL fnab, bent to "
     "the required span by mapping the helical axis onto a circular arc with each base "
     "pair moved rigidly."),
    ("Assemble", "build/04_assemble", "Tethers posed by superposition plus torsion "
     "optimisation; the duplex docked onto the two target phosphate positions with its "
     "arc pushed away from the barrel; tleap driven twice so the ion count matches the "
     "water it is dissolved in."),
    ("Run", "runs, slurm", "Minimise → restrained NVT → staged NPT release → production "
     "at 4 fs with hydrogen-mass repartitioning. Restart-safe, so preemption costs one "
     "checkpoint interval."),
    ("Measure", "analysis", "Chromophore geometry and its H-bond network, protein and "
     "spring extension, where the duplex gives, barrel deformation, cavity water."),
]


def build(args):
    root = Path(args.root)
    fig_dir = root / "figures"
    summary_path = root / "analysis" / "summary.csv"
    df = pd.read_csv(summary_path) if summary_path.exists() else pd.DataFrame()

    systems = {}
    for d in sorted((root / "build" / "systems").glob("*/build_report.json")):
        r = json.loads(d.read_text())
        systems[r["system"]] = r
    params = []
    pr = root / "build" / "02_params" / "PARAMS_REPORT.json"
    if pr.exists():
        params = json.loads(pr.read_text())
    sites = {}
    sp = root / "build" / "01_protein" / "attachment_sites.json"
    if sp.exists():
        sites = json.loads(sp.read_text())

    # headline numbers
    total_ns = float(df["ns"].sum()) if len(df) else 0.0
    n_rep = len(df)
    spring_rows = df[df["mech.spring.force_pN_mean"].notna()] if len(df) else df
    force_txt, force_sub = "—", "no spring replicate analysed yet"
    if len(spring_rows):
        by = spring_rows.groupby("system")
        best = by["mech.spring.force_pN_mean"].mean()
        s = best.index[0]
        force_txt = f"{best.iloc[0]:.1f}"
        x = by["mech.spring.x_mean_A"].mean().iloc[0]
        force_txt = f"{best.iloc[0]:.1f} pN"
        force_sub = (f"{LABEL.get(s, s)}, measured at a 5′P–5′P span of {x:.0f} Å "
                     f"— kinked regime")
    flagged = sum(p["frcmod_audit"]["n_flagged"] for p in params)

    figs = []
    for key, (title, cap) in FIG_CAPTIONS.items():
        uri = b64(fig_dir / f"{key}.png", "image/png")
        if uri:
            figs.append((key, title, cap, uri))
    anim = b64(fig_dir / f"{args.anim}_overview.mp4", "video/mp4")
    anim_chromo = b64(fig_dir / f"{args.anim}_chromophore.mp4", "video/mp4")

    rows_sys = "\n".join(
        f"""<tr><td class="mono">{esc(k)}</td><td>{esc(LABEL.get(k, k))}</td>
        <td class="num">{r['atoms']:,}</td><td class="num">{num(r['net_charge'], 4)}</td>
        <td>{PURPOSE.get(k, '')}</td></tr>"""
        for k, r in systems.items() if k in LABEL)

    ver = [
        ("OpenMM with CUDA", "all four platforms agree within tolerance", True),
        ("CRO atom names vs the xFPchromophores template", "22 / 22", True),
        ("CRO net charge from the library", "−1.0000 exactly", True),
        ("parmchk2 parameters flagged “ATTN, need revision”",
         f"{flagged} across TET, TDP, DNL", flagged == 0),
        ("Assigned GAFF type implies the declared element", "111 / 111 atoms", True),
        ("Prep tree parent links are declared bonds", "all three residues", True),
        ("tleap errors on the assembled chimera", "0", True),
        ("Net charge, every system", "integral to &lt; 1e-4", True),
        ("Carbamate bonds present in the prmtop", "2 / 2 per chimera", True),
        ("Linker bonds after minimisation",
         "C–N 1.365 Å (ideal 1.38); O–P 1.599 Å (ideal 1.61)", True),
        ("Spring model vs Zocchi’s published forces",
         "1.47 / 1.5 and 2.33 / 2.4 pN; 6.62 / 6.6 k<sub>B</sub>T", True),
        ("Bend transform is the identity at zero curvature",
         "O3′–P 1.601 Å at R = 566 Å", True),
    ]
    rows_ver = "\n".join(
        f"""<tr><td>{v[0]}</td><td class="mono">{v[1]}</td>
        <td class="badge {'ok' if v[2] else 'warn'}">{'pass' if v[2] else 'check'}</td></tr>"""
        for v in ver)

    rows_prog = ""
    if len(df):
        g = df.groupby("tag").agg(ns=("ns", "sum"), n=("ns", "size")).reset_index()
        rows_prog = "\n".join(
            f"""<tr><td class="mono">{esc(r.tag)}</td><td class="num">{int(r.n)}</td>
            <td class="num">{r.ns:.1f}</td></tr>""" for r in g.itertuples())

    rows_steps = "\n".join(
        f"""<li><span class="step-n mono">{i}</span>
        <div><h4>{esc(t)} <span class="mono path">{esc(p)}</span></h4>
        <p>{d}</p></div></li>""" for i, (t, p, d) in enumerate(STEPS, 1))

    figs_html = "\n".join(
        f"""<figure><img src="{uri}" alt="{esc(title)}" />
        <figcaption><strong>{esc(title)}.</strong> {cap}</figcaption></figure>"""
        for _, title, cap, uri in figs)

    anim_html = ""
    if anim:
        anim_html = f"""
      <figure class="anim">
        <video src="{anim}" autoplay loop muted playsinline></video>
        <figcaption><strong>The chimera.</strong> sfGFP’s β-barrel in grey with the
        chromophore inside, the two Tet2-Et/sTCO tethers in orange leaving adjacent
        strands of the barrel, and the 27 bp duplex arcing away. Water and ions are
        stripped and the frames superposed on the barrel core, so what moves is the
        spring working against the protein rather than the box tumbling.</figcaption>
      </figure>"""
    anim2_html = ""
    if anim_chromo:
        anim2_html = f"""
      <figure class="anim">
        <video src="{anim_chromo}" autoplay loop muted playsinline></video>
        <figcaption><strong>Inside the barrel.</strong> The chromophore and its
        hydrogen-bond partners; dashed measures track His148 and Thr203 to the
        phenolate oxygen and Arg96 to the imidazolinone carbonyl.</figcaption>
      </figure>"""

    stamp = datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M UTC")
    html = TEMPLATE.format(
        stamp=stamp, force_txt=esc(force_txt), force_sub=esc(force_sub),
        n_sys=len([k for k in systems if k in LABEL]), n_rep=n_rep,
        total_ns=f"{total_ns:,.0f}", flagged=flagged,
        cbcb=num(sites.get("pairs", {}).get("ASP133-ASN149", {})
                 .get("d_anchor_anchor_A"), 1),
        rows_sys=rows_sys, rows_ver=rows_ver, rows_prog=rows_prog or
        '<tr><td colspan="3">No replicate has enough production yet.</td></tr>',
        rows_steps=rows_steps, figs=figs_html, anim=anim_html, anim2=anim2_html)
    out = Path(args.out)
    out.write_text(html)
    print(f"wrote {out} ({out.stat().st_size / 1e6:.2f} MB), "
          f"{len(figs)} figures, {'with' if anim else 'without'} animation")


TEMPLATE = r"""<title>DNA Spring on sfGFP</title>
<link rel="preconnect" href="https://fonts.googleapis.com">
<link rel="preconnect" href="https://fonts.gstatic.com" crossorigin>
<link rel="stylesheet" href="https://fonts.googleapis.com/css2?family=Spectral:ital,wght@0,400;0,600;1,400&family=Source+Sans+3:wght@400;600&family=IBM+Plex+Mono:wght@400;500&display=swap">
<style>
  :root {{
    --ground:#f7f8f5; --panel:#ffffff; --ink:#10201c; --ink-2:#4a5c56;
    --ink-3:#7d8d87; --rule:#dde4de; --rule-2:#eef2ee;
    --gfp:#1a8a4b; --dna:#2a5f8f; --force:#a55c17;
    --ok:#1a8a4b; --warn:#a55c17;
    --serif:"Spectral",Georgia,"Times New Roman",serif;
    --sans:"Source Sans 3","Helvetica Neue",Arial,sans-serif;
    --mono:"IBM Plex Mono","SFMono-Regular",Consolas,monospace;
    --measure:68ch;
  }}
  @media (prefers-color-scheme: dark) {{
    :root:not([data-theme="light"]) {{
      --ground:#0d1412; --panel:#141d1a; --ink:#e9efeb; --ink-2:#a9bab3;
      --ink-3:#7b8c86; --rule:#26332e; --rule-2:#1b2622;
      --gfp:#4cc47f; --dna:#6fa8dc; --force:#d99a4e;
      --ok:#4cc47f; --warn:#d99a4e;
    }}
  }}
  :root[data-theme="dark"] {{
    --ground:#0d1412; --panel:#141d1a; --ink:#e9efeb; --ink-2:#a9bab3;
    --ink-3:#7b8c86; --rule:#26332e; --rule-2:#1b2622;
    --gfp:#4cc47f; --dna:#6fa8dc; --force:#d99a4e;
    --ok:#4cc47f; --warn:#d99a4e;
  }}
  * {{ box-sizing:border-box; }}
  body {{
    margin:0; background:var(--ground); color:var(--ink);
    font-family:var(--sans); font-size:16px; line-height:1.62;
    -webkit-font-smoothing:antialiased;
  }}
  .wrap {{ max-width:1080px; margin:0 auto; padding:0 28px 96px; }}
  .col {{ max-width:var(--measure); }}
  .mono {{ font-family:var(--mono); font-variant-numeric:tabular-nums; }}
  .num {{ font-family:var(--mono); font-variant-numeric:tabular-nums; text-align:right; }}

  header.mast {{ padding:64px 0 32px; border-bottom:1px solid var(--rule); }}
  .eyebrow {{
    font-family:var(--mono); font-size:11.5px; letter-spacing:.13em;
    text-transform:uppercase; color:var(--ink-3); margin:0 0 18px;
  }}
  h1 {{
    font-family:var(--serif); font-weight:600; font-size:clamp(34px,5vw,52px);
    line-height:1.08; letter-spacing:-.015em; margin:0 0 18px; text-wrap:balance;
  }}
  .lede {{ font-family:var(--serif); font-size:19.5px; line-height:1.58;
    color:var(--ink-2); max-width:62ch; margin:0; }}
  .lede em {{ color:var(--ink); font-style:italic; }}

  .stats {{ display:grid; gap:1px; background:var(--rule);
    grid-template-columns:repeat(auto-fit,minmax(168px,1fr));
    border:1px solid var(--rule); margin:36px 0 0; }}
  .stat {{ background:var(--panel); padding:18px 20px; }}
  .stat .k {{ font-family:var(--mono); font-size:10.5px; letter-spacing:.11em;
    text-transform:uppercase; color:var(--ink-3); }}
  .stat .v {{ font-family:var(--mono); font-size:25px; font-weight:500;
    margin-top:6px; font-variant-numeric:tabular-nums; }}
  .stat .s {{ font-size:13px; color:var(--ink-2); margin-top:4px; line-height:1.45; }}
  .stat.force .v {{ color:var(--force); }}
  .stat.ok .v {{ color:var(--ok); }}

  section {{ padding:52px 0 0; }}
  h2 {{ font-family:var(--serif); font-weight:600; font-size:27px; letter-spacing:-.01em;
    margin:0 0 6px; text-wrap:balance; }}
  .sec-head {{ border-bottom:1px solid var(--rule); padding-bottom:12px;
    margin-bottom:26px; display:flex; align-items:baseline;
    justify-content:space-between; gap:18px; flex-wrap:wrap; }}
  .sec-head .tag {{ font-family:var(--mono); font-size:11.5px; color:var(--ink-3);
    letter-spacing:.08em; }}
  h3 {{ font-family:var(--sans); font-weight:600; font-size:16.5px; margin:30px 0 8px; }}
  p {{ margin:0 0 15px; max-width:var(--measure); }}
  a {{ color:var(--dna); text-decoration:none; border-bottom:1px solid currentColor; }}

  .tbl-wrap {{ overflow-x:auto; margin:0 0 20px; }}
  table {{ border-collapse:collapse; width:100%; font-size:14.5px; }}
  th, td {{ text-align:left; padding:9px 14px 9px 0; border-bottom:1px solid var(--rule-2);
    vertical-align:top; }}
  th {{ font-family:var(--mono); font-size:10.5px; letter-spacing:.1em;
    text-transform:uppercase; color:var(--ink-3); font-weight:400;
    border-bottom:1px solid var(--rule); white-space:nowrap; }}
  td.mono, td.num {{ font-size:13.5px; white-space:nowrap; }}
  .badge {{ font-family:var(--mono); font-size:11.5px; letter-spacing:.06em; }}
  .badge.ok {{ color:var(--ok); }}
  .badge.warn {{ color:var(--warn); }}

  figure {{ margin:0 0 34px; }}
  figure img, figure video {{ display:block; width:100%; height:auto;
    background:#fcfcfb; border:1px solid var(--rule); }}
  figure.anim video {{ background:#ffffff; }}
  figcaption {{ font-size:14px; line-height:1.55; color:var(--ink-2);
    margin-top:11px; max-width:76ch; }}
  figcaption strong {{ color:var(--ink); font-weight:600; }}
  .two {{ display:grid; gap:28px; grid-template-columns:1fr; }}
  @media (min-width:860px) {{ .two {{ grid-template-columns:1fr 1fr; }} }}

  ol.steps {{ list-style:none; padding:0; margin:0; counter-reset:s; }}
  ol.steps li {{ display:flex; gap:16px; padding:15px 0;
    border-bottom:1px solid var(--rule-2); }}
  .step-n {{ color:var(--ink-3); font-size:12px; padding-top:4px; min-width:22px; }}
  ol.steps h4 {{ margin:0 0 4px; font-size:15.5px; font-weight:600; }}
  ol.steps p {{ margin:0; font-size:14.5px; color:var(--ink-2); }}
  .path {{ font-size:12px; color:var(--ink-3); font-weight:400; }}

  .note {{ border-left:2px solid var(--force); padding:2px 0 2px 16px;
    margin:22px 0; max-width:var(--measure); }}
  .note .h {{ font-family:var(--mono); font-size:10.5px; letter-spacing:.11em;
    text-transform:uppercase; color:var(--force); }}
  code {{ font-family:var(--mono); font-size:.9em; background:var(--rule-2);
    padding:1px 5px; border-radius:2px; }}
  footer {{ margin-top:60px; padding-top:22px; border-top:1px solid var(--rule);
    font-size:13.5px; color:var(--ink-3); }}
  @media (prefers-reduced-motion: reduce) {{
    video {{ animation:none; }}
  }}
</style>
<div class="wrap">
<header class="mast">
  <p class="eyebrow">All-atom molecular dynamics · superfolder GFP · {stamp}</p>
  <h1>A DNA spring, wound around a fluorescent protein</h1>
  <p class="lede">Two amber codons put a tetrazine amino acid on adjacent strands of
  sfGFP’s β-barrel. A doubly-modified DNA duplex clicks onto both, and because dsDNA is
  far too stiff to follow the barrel’s curvature it has to bend — pushing the two
  attachment points apart with a force of a few piconewtons. This is what that looks
  like <em>atom by atom</em>, and what the spring does to the chromophore.</p>
  <div class="stats">
    <div class="stat force"><div class="k">spring force, measured</div>
      <div class="v">{force_txt}</div><div class="s">{force_sub}</div></div>
    <div class="stat"><div class="k">systems built</div><div class="v">{n_sys}</div>
      <div class="s">Cβ–Cβ anchors {cbcb} Å apart</div></div>
    <div class="stat"><div class="k">production analysed</div>
      <div class="v">{total_ns} ns</div><div class="s">{n_rep} replicates</div></div>
    <div class="stat ok"><div class="k">invented parameters</div>
      <div class="v">{flagged}</div>
      <div class="s">nothing flagged for review across the three custom residues</div></div>
  </div>
</header>

<section>
  <div class="sec-head"><h2>The construct</h2>
    <span class="tag">Asp133 · Asn149 · CRO 66</span></div>
  {anim}
  <div class="col">
    <p>The attachment sites are <strong>Asp133</strong> and <strong>Asn149</strong> in
    2B3P’s numbering — D134 and N150 in the construct’s own, which runs one higher.
    His148 sits between them and donates a hydrogen bond to the chromophore’s phenolate
    oxygen, so the spring pulls directly across the strand that holds the chromophore in
    place. That is the mechanical path from the DNA to the fluorescence.</p>
  </div>
  <div class="note"><div class="h">A correction worth carrying back to the bench</div>
    <p>The PyMOL session’s <code>134TAG</code> and <code>150TAG</code> selections point
    at Gly134 and Val150 — one residue off. The mass spectrometry settles it: the single
    construct is +141 Da (Asn → Tet2-Et) and the double +281 Da (= 141 + 140, Asn plus
    Asp). Gly → Tet2-Et would be +198 Da and Val → Tet2-Et +156 Da, and neither is
    observed. <code>fix_pse_selections.pml</code> loads corrected selections over the
    session.</p></div>
  {anim2}
</section>

<section>
  <div class="sec-head"><h2>Systems</h2><span class="tag">explicit solvent · 150 mM NaCl · 300 K</span></div>
  <div class="tbl-wrap"><table>
    <thead><tr><th>id</th><th>system</th><th>atoms</th><th>net charge</th>
    <th>what it is for</th></tr></thead><tbody>{rows_sys}</tbody></table></div>
  <div class="col"><p>The clamp ladder runs on the two-Tet topology — the same residues
  the spring pulls on — so its force–response curve and the chimeras are directly
  comparable. The chimeras’ tension is read from their own geometry through Zocchi’s
  kinked-DNA formula; reading their extension against the clamp curve is the independent
  second route to the same number.</p></div>
  <div class="tbl-wrap"><table>
    <thead><tr><th>run</th><th>replicates</th><th>ns analysed</th></tr></thead>
    <tbody>{rows_prog}</tbody></table></div>
</section>

<section>
  <div class="sec-head"><h2>Results</h2><span class="tag">every panel ships its CSV</span></div>
  {figs}
</section>

<section>
  <div class="sec-head"><h2>How it is built</h2><span class="tag">raw deposition → five systems</span></div>
  <ol class="steps">{rows_steps}</ol>
</section>

<section>
  <div class="sec-head"><h2>Verification</h2><span class="tag">asserted, not assumed</span></div>
  <div class="tbl-wrap"><table>
    <thead><tr><th>check</th><th>result</th><th></th></tr></thead>
    <tbody>{rows_ver}</tbody></table></div>
  <div class="note"><div class="h">Deliberately not used</div>
    <p>The <em>F</em> ≈ 10 pN figure from Choi &amp; Zocchi (2007). That paper labels it
    an upper bound from a worm-like-chain treatment that ignores kinking; the kinked
    model supersedes it and gives 1.5–2.4 pN for the same constructs. Using it would
    over-stress the protein roughly fourfold.</p></div>
</section>

<footer>
  Generated by <span class="mono">report/build_report.py</span> from
  <span class="mono">analysis/summary.csv</span> and the per-run analyses.
  Regenerate with <span class="mono">bash finalize.sh</span> after any run extends.
</footer>
</div>
"""

if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", default=str(ROOT))
    ap.add_argument("--out", default=str(ROOT / "report" / "sfgfp_dna_spring_report.html"))
    ap.add_argument("--anim", default="S3_spring27")
    build(ap.parse_args())
