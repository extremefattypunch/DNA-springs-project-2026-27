#!/usr/bin/env python3
"""Build a standalone, interactive overlay viewer for the sfGFP-DNA chimera run.

The 260812 folder is an Origin export: twenty ``Dflt{Em,Ex}<n>_Data.csv`` files
with no sample names in them.  The identity of each scan lives only in the
operator's notebook, so it is declared here (see ``ASSIGNMENT``) rather than
guessed from the filenames.

Emission files are two-column (Wavelength / S1).  Excitation files are
four-column (A,S1,B,R1): S1 is the emission-channel photon count and R1 is the
reference photodiode in microamps, which tracks the xenon lamp output across
the sweep.  Both are carried into the payload so the viewer can offer a
lamp-corrected excitation spectrum (S1/R1) alongside the raw counts.

Usage:
    python3 build_chimera_viewer.py                       # defaults below
    python3 build_chimera_viewer.py --root DIR --out FILE
"""

from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
DEFAULT_ROOT = HERE / "260812_sfGFPDNAchimeraALL_1-3"
DEFAULT_OUT = DEFAULT_ROOT / "chimera_spectra_viewer.html"

# Scan number -> (construct, replicate).  Declared by the operator; the file
# names carry no sample identity.  Replicate order follows the scan order
# within each construct, i.e. Em2/Em5/Em8 are 0-tet replicates 1/2/3.
ASSIGNMENT = {
    1: ("calibration", None),
    2: ("0-tet sfGFP", 1), 5: ("0-tet sfGFP", 2), 8: ("0-tet sfGFP", 3),
    3: ("1-tet sfGFP", 1), 6: ("1-tet sfGFP", 2), 9: ("1-tet sfGFP", 3),
    4: ("2-tet sfGFP", 1), 7: ("2-tet sfGFP", 2), 10: ("2-tet sfGFP", 3),
}

GROUP_ORDER = ["calibration", "0-tet sfGFP", "1-tet sfGFP", "2-tet sfGFP"]

# The emission monochromator sat here during the excitation sweeps, so every
# excitation scan carries a Rayleigh stray-light spike at this wavelength.
# Detected below and reported, not assumed.
SCATTER_GUESS = 510.0


# --------------------------------------------------------------------------
# reading
# --------------------------------------------------------------------------

def read_origin_csv(path: Path):
    """Return (long_names, units, columns) from an Origin ``*_Data.csv``.

    Layout is three header rows (short names / long names / units) followed by
    the numeric block.  Columns are returned as parallel lists of floats.
    """
    with path.open(newline="", encoding="utf-8-sig") as fh:
        rows = list(csv.reader(fh))
    if len(rows) < 4:
        raise ValueError(f"{path.name}: fewer than 4 rows")
    names, units = rows[1], rows[2]
    ncol = len(rows[0])
    cols = [[] for _ in range(ncol)]
    for ln, row in enumerate(rows[3:], start=4):
        if not row or not row[0].strip():
            continue
        if len(row) < ncol:
            raise ValueError(f"{path.name}:{ln}: expected {ncol} cells, got {len(row)}")
        for i in range(ncol):
            cell = row[i].strip()
            if cell == "":
                raise ValueError(f"{path.name}:{ln}: empty cell in column {i}")
            cols[i].append(float(cell))
    return names, units, cols


def uniform_axis(xs):
    """Return (x0, dx, n) if ``xs`` is a uniform ascending grid, else None."""
    if len(xs) < 2:
        return None
    dx = xs[1] - xs[0]
    if dx <= 0:
        return None
    for i in range(1, len(xs)):
        if abs((xs[i] - xs[i - 1]) - dx) > 1e-6:
            return None
    return xs[0], dx, len(xs)


def find_scatter_line(x, y):
    """Locate a 1-3 nm stray-light spike, or None.

    A Rayleigh line is a point that towers over its own neighbourhood; a real
    spectral band does not.  Compare the maximum against the median of the
    window 8-20 nm away from it, which is outside the spike but still inside
    the band it sits on.
    """
    if len(y) < 40:
        return None
    pk = max(range(len(y)), key=lambda i: y[i])
    near = [y[i] for i in range(len(y)) if 8 <= abs(x[i] - x[pk]) <= 20]
    if len(near) < 8:
        return None
    near.sort()
    base = near[len(near) // 2]
    if base > 0 and y[pk] > 8 * base:
        return x[pk], y[pk] / base
    return None


# --------------------------------------------------------------------------
# payload
# --------------------------------------------------------------------------

def build_payload(root: Path):
    traces, notes = [], []
    missing = []

    for mode, prefix in (("emission", "Em"), ("excitation", "Ex")):
        for num in sorted(ASSIGNMENT):
            path = root / f"Dflt{prefix}{num}_Data.csv"
            if not path.exists():
                missing.append(path.name)
                continue
            names, units, cols = read_origin_csv(path)
            group, rep = ASSIGNMENT[num]

            x, s1 = cols[0], cols[1]
            axis = uniform_axis(x)
            if axis is None:
                raise ValueError(f"{path.name}: wavelength axis is not a uniform grid")
            x0, dx, n = axis

            r1 = None
            if len(cols) >= 4:
                if any(abs(a - b) > 1e-6 for a, b in zip(cols[0], cols[2])):
                    raise ValueError(f"{path.name}: S1 and R1 wavelength columns disagree")
                r1 = cols[3]

            trace = {
                "id": f"{prefix}{num}",
                "file": path.name,
                "mode": mode,
                "group": group,
                "rep": rep,
                "x0": x0, "dx": dx, "n": n,
                "s1": s1,
                "r1": r1,
                "yUnit": units[1] if len(units) > 1 else "CPS",
                "rUnit": units[3] if len(units) > 3 else None,
            }
            traces.append(trace)

            if mode == "excitation" and group != "calibration":
                hit = find_scatter_line(x, s1)
                if hit:
                    trace["scatter"] = round(hit[0], 1)

    if missing:
        notes.append("missing files: " + ", ".join(missing))

    # Report the scatter line once, from the data, rather than per file.
    lines = sorted({t["scatter"] for t in traces if "scatter" in t})
    if lines:
        where = ", ".join(f"{v:g} nm" for v in lines)
        notes.append(
            f"Every excitation scan of a sample carries a stray-light spike at {where} "
            "— the emission monochromator was parked there, so the sweep sees the "
            "lamp directly when λ_ex = λ_em. It is 10–50× the real "
            "excitation band and is masked by default.")
    notes.append(
        "Excitation files carry a reference photodiode (R1, µA) that tracks lamp "
        "output across the sweep. Raw S1 is uncorrected for the xenon lamp profile; "
        "S1÷R1 is the lamp-corrected excitation spectrum.")
    notes.append(
        "Scan→sample assignment is declared, not inferred: "
        "1 = calibration; 2/5/8 = 0-tet; 3/6/9 = 1-tet; 4/7/10 = 2-tet "
        "(replicates 1/2/3 in that order).")

    return {
        "source": root.name,
        "groups": GROUP_ORDER,
        "traces": traces,
        "notes": notes,
        "scatterDefault": lines[0] if lines else SCATTER_GUESS,
    }


# --------------------------------------------------------------------------
# page
# --------------------------------------------------------------------------

HTML = r"""<!doctype html>
<html lang="en">
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>__TITLE__</title>
<style>
  /* Palette: three categorical slots (blue / orange / aqua) for the three
     constructs, validated all-pairs in both light and dark. Calibration is
     chart chrome, not a fourth series -- muted ink, thin, off by default, and
     living in a disjoint wavelength region. Emission is solid and excitation
     dashed, so the measurement axis never rides on colour. */
  :root {
    color-scheme: light;
    --surface:#fcfcfb; --plane:#f9f9f7; --ink:#0b0b0b; --ink2:#52514e;
    --muted:#898781; --grid:#e1e0d9; --axis:#c3c2b7; --hair:rgba(11,11,11,.10);
    --accent:#2a78d6; --warnbg:#fdf6e3; --warnink:#6b4f12; --warnline:#e6d18a;
    --c-cal:#898781; --c-g0:#2a78d6; --c-g1:#eb6834; --c-g2:#1baf7a;
  }
  @media (prefers-color-scheme: dark) {
    :root:not([data-theme="light"]) {
      color-scheme: dark;
      --surface:#1a1a19; --plane:#0d0d0d; --ink:#ffffff; --ink2:#c3c2b7;
      --muted:#898781; --grid:#2c2c2a; --axis:#383835; --hair:rgba(255,255,255,.10);
      --accent:#3987e5; --warnbg:#241f12; --warnink:#e8d9a8; --warnline:#4a3f22;
      --c-cal:#898781; --c-g0:#3987e5; --c-g1:#d95926; --c-g2:#199e70;
    }
  }
  :root[data-theme="dark"] {
    color-scheme: dark;
    --surface:#1a1a19; --plane:#0d0d0d; --ink:#ffffff; --ink2:#c3c2b7;
    --muted:#898781; --grid:#2c2c2a; --axis:#383835; --hair:rgba(255,255,255,.10);
    --accent:#3987e5; --warnbg:#241f12; --warnink:#e8d9a8; --warnline:#4a3f22;
    --c-cal:#898781; --c-g0:#3987e5; --c-g1:#d95926; --c-g2:#199e70;
  }

  * { box-sizing:border-box; }
  body {
    margin:0; background:var(--surface); color:var(--ink);
    font:13px/1.45 system-ui, -apple-system, "Segoe UI", sans-serif;
  }
  header {
    padding:13px 20px; border-bottom:1px solid var(--hair);
    display:flex; align-items:baseline; gap:14px; flex-wrap:wrap;
  }
  header h1 { margin:0; font-size:15.5px; font-weight:650; letter-spacing:-.01em; }
  header .sub { color:var(--muted); font-size:12px; }
  header .spacer { margin-left:auto; }
  .wrap { display:flex; align-items:stretch; min-height:calc(100vh - 50px); }
  aside {
    width:272px; flex:0 0 272px; padding:15px; background:var(--plane);
    border-right:1px solid var(--hair); overflow-y:auto; max-height:calc(100vh - 50px);
  }
  main { flex:1 1 auto; padding:15px 20px 26px; min-width:0; }
  .grp { margin-bottom:17px; }
  .grp > h2 {
    margin:0 0 6px; font-size:10.5px; font-weight:700; text-transform:uppercase;
    letter-spacing:.07em; color:var(--muted);
    display:flex; justify-content:space-between; align-items:center; gap:8px;
  }
  .grp > h2 button {
    font:inherit; font-size:10px; text-transform:none; letter-spacing:0;
    padding:1px 6px; border:1px solid var(--hair); background:var(--surface);
    border-radius:4px; color:var(--muted); cursor:pointer;
  }
  .grp > h2 button:hover { color:var(--accent); border-color:var(--accent); }
  label.chk {
    display:flex; align-items:center; gap:7px; padding:2.5px 0;
    cursor:pointer; user-select:none;
  }
  label.chk:hover { color:var(--accent); }
  label.chk input { margin:0; accent-color:var(--accent); flex:0 0 auto; }
  label.chk .sw { width:17px; height:0; flex:0 0 auto; }
  label.chk .cnt { margin-left:auto; color:var(--muted); font-size:11px;
                   font-variant-numeric:tabular-nums; }
  select, input[type=range], input[type=number] { font:inherit; }
  select, input[type=number] {
    padding:4px 6px; border:1px solid var(--hair); border-radius:5px;
    background:var(--surface); color:var(--ink); font-size:12.5px;
  }
  select { width:100%; }
  input[type=range] { width:100%; accent-color:var(--accent); }
  .row { display:flex; gap:8px; align-items:center; margin-top:6px; }
  .row > span { color:var(--muted); font-size:11.5px; white-space:nowrap; }
  .btn {
    padding:5px 10px; border:1px solid var(--hair); background:var(--surface);
    border-radius:5px; font:inherit; font-size:12px; cursor:pointer; color:var(--ink);
  }
  .btn:hover { border-color:var(--accent); color:var(--accent); }
  .btn[aria-pressed="true"] { border-color:var(--accent); color:var(--accent); }
  #plotbox { position:relative; }
  canvas { display:block; width:100%; touch-action:none; cursor:crosshair; }
  #tip {
    position:absolute; pointer-events:none; display:none; z-index:5;
    background:var(--surface); border:1px solid var(--hair); border-radius:6px;
    padding:7px 9px; font-size:11.5px; max-width:310px;
    box-shadow:0 4px 16px rgba(0,0,0,.16);
  }
  #tip b { display:block; margin-bottom:4px; font-variant-numeric:tabular-nums; }
  #tip .r { display:flex; gap:6px; align-items:center; white-space:nowrap;
            font-variant-numeric:tabular-nums; }
  #tip .r i { width:14px; height:0; flex:0 0 auto; }
  #tip .r em { font-style:normal; color:var(--ink2); overflow:hidden;
               text-overflow:ellipsis; }
  #tip .r span { margin-left:auto; font-weight:600; }
  #legend { margin-top:12px; border-top:1px solid var(--hair); padding-top:10px; }
  #legend .lgrp { display:flex; flex-wrap:wrap; gap:2px 14px; margin-bottom:2px; }
  #legend .it {
    display:flex; align-items:center; gap:7px; cursor:pointer;
    padding:2px 4px; border-radius:4px; font-size:12px;
  }
  #legend .it:hover { background:var(--plane); }
  #legend .it.off { opacity:.34; }
  #legend .it .sw { width:22px; height:0; flex:0 0 auto; }
  #legend .it.mean { font-weight:600; }
  #status { margin-top:9px; color:var(--muted); font-size:11.5px; }
  .flag {
    margin-top:11px; padding:8px 11px; background:var(--warnbg);
    border:1px solid var(--warnline); border-radius:6px; font-size:11.5px;
    color:var(--warnink); display:flex; gap:10px; align-items:center;
  }
  .flag button { margin-left:auto; flex:0 0 auto; }
  #notes {
    margin-top:13px; font-size:11.5px; color:var(--ink2);
    border-top:1px solid var(--hair); padding-top:10px;
  }
  #notes ul { margin:5px 0 0; padding-left:17px; }
  #notes li { margin-bottom:3px; }
  table.tbl { border-collapse:collapse; margin-top:12px; width:100%; font-size:11.5px; }
  table.tbl th, table.tbl td {
    text-align:right; padding:4px 8px; border-bottom:1px solid var(--hair);
    font-variant-numeric:tabular-nums; white-space:nowrap;
  }
  table.tbl th { color:var(--muted); font-weight:600; text-align:right; }
  table.tbl th:first-child, table.tbl td:first-child { text-align:left; }
  table.tbl td.name { display:flex; align-items:center; gap:7px; }
  table.tbl td.name .sw { width:18px; height:0; flex:0 0 auto; }
  @media (max-width:880px) {
    .wrap { flex-direction:column; }
    aside { width:auto; flex:none; max-height:none; border-right:0;
            border-bottom:1px solid var(--hair); }
  }
</style>

<header>
  <h1>__TITLE__</h1>
  <span class="sub" id="subtitle"></span>
  <span class="spacer"></span>
  <button class="btn" id="b_theme" title="light / dark / follow system">theme: auto</button>
</header>

<div class="wrap">
<aside>
  <div class="grp"><h2>Measurement <button data-all="mode">all</button></h2><div id="f_mode"></div></div>
  <div class="grp"><h2>Construct <button data-all="group">all</button></h2><div id="f_group"></div></div>
  <div class="grp"><h2>Replicate <button data-all="rep">all</button></h2><div id="f_rep"></div></div>
  <div class="grp"><h2>Trace <button data-all="kind">all</button></h2><div id="f_kind"></div></div>

  <div class="grp">
    <h2>Excitation handling</h2>
    <div class="row"><span>channel</span>
      <select id="o_exch">
        <option value="s1">S1 raw (CPS)</option>
        <option value="corr">S1 &divide; R1 lamp-corrected</option>
        <option value="r1">R1 reference (&micro;A)</option>
      </select>
    </div>
    <label class="chk"><input type="checkbox" id="o_mask" checked>
      mask &lambda;<sub>ex</sub>=&lambda;<sub>em</sub> scatter</label>
    <div class="row" id="maskrow">
      <span>centre</span>
      <input type="number" id="o_maskc" step="1" style="width:5.2em">
      <span>&plusmn;</span>
      <input type="number" id="o_maskw" step="1" min="0" max="30" value="3" style="width:3.6em">
      <span>nm</span>
    </div>
  </div>

  <div class="grp">
    <h2>Processing</h2>
    <label class="chk"><input type="checkbox" id="o_sd" checked>
      &plusmn;1 SD band on means</label>
    <label class="chk"><input type="checkbox" id="o_peaks"> mark &lambda;<sub>max</sub></label>
    <div class="row"><span>normalise</span>
      <select id="o_norm">
        <option value="none">raw signal</option>
        <option value="peak">peak = 1</option>
        <option value="area">area = 1</option>
      </select>
    </div>
    <div class="row"><span>y&#8209;axis</span>
      <select id="o_scale">
        <option value="linear">linear</option>
        <option value="log">log</option>
      </select>
    </div>
    <div class="row"><span>smooth</span>
      <input type="range" id="o_smooth" min="1" max="21" step="2" value="1">
      <span id="o_smooth_v" style="min-width:3.4em;text-align:right">off</span>
    </div>
  </div>

  <div class="grp">
    <h2>View</h2>
    <div class="row">
      <button class="btn" id="b_table" aria-pressed="false">table</button>
      <button class="btn" id="b_reset">reset view</button>
    </div>
    <div class="row">
      <button class="btn" id="b_png">PNG</button>
      <button class="btn" id="b_csv">CSV</button>
    </div>
  </div>
</aside>

<main>
  <div id="plotbox"><canvas id="cv"></canvas><div id="tip"></div></div>
  <div id="flags"></div>
  <div id="status"></div>
  <div id="legend"></div>
  <div id="tablebox"></div>
  <div id="notes"></div>
</main>
</div>

<script>
"use strict";
const PAYLOAD = __PAYLOAD__;
const RAW    = PAYLOAD.traces;
const GROUPS = PAYLOAD.groups;
const MODES  = ["emission", "excitation"];

/* Wavelength grids are uniform, so they ship as {x0,dx,n} and are expanded once. */
for (const t of RAW) {
  t.x = new Array(t.n);
  for (let i = 0; i < t.n; i++) t.x[i] = t.x0 + i * t.dx;
  t.xlo = t.x[0]; t.xhi = t.x[t.n - 1];
  t.ix = new Map(t.x.map((v, i) => [v, i]));
  t.repKey = t.rep == null ? null : String(t.rep);
  t.label = t.group === "calibration"
    ? `calibration · ${t.mode === "emission" ? "Em" : "Ex"}`
    : `${t.group} · rep ${t.rep}`;
}

const COLOR_VAR = { "calibration":"--c-cal", "0-tet sfGFP":"--c-g0",
                    "1-tet sfGFP":"--c-g1", "2-tet sfGFP":"--c-g2" };
let CSSV = {};
function readTheme() {
  const cs = getComputedStyle(document.documentElement);
  CSSV = {};
  for (const k of ["--ink","--ink2","--muted","--grid","--axis","--surface",
                   "--c-cal","--c-g0","--c-g1","--c-g2"]) {
    CSSV[k] = cs.getPropertyValue(k).trim();
  }
}
const colorOf = g => CSSV[COLOR_VAR[g]] || CSSV["--ink"];
/* Emission solid, excitation dashed: the measurement never rides on colour. */
const dashOf  = m => m === "emission" ? [] : [7, 4];

/* ---------- state ---------- */
const FACETS = [
  { key:"mode",  box:"f_mode",
    vals: MODES,
    label: v => v === "emission" ? "Emission (Em)" : "Excitation (Ex)" },
  { key:"group", box:"f_group",
    vals: GROUPS,
    label: v => v === "calibration" ? "calibration (off-scale)" : v },
  { key:"rep",   box:"f_rep",
    vals: ["1","2","3"],
    label: v => "replicate " + v },
  { key:"kind",  box:"f_kind",
    vals: ["individual","mean"],
    label: v => v === "individual" ? "individual scans" : "group mean" },
];

const state = {
  on: {
    mode:  new Set(MODES),
    /* Calibration is 20-40x the sample signal and sits in a disjoint band;
       starting it on would flatten every sample trace to the baseline. */
    group: new Set(GROUPS.filter(g => g !== "calibration")),
    rep:   new Set(["1","2","3"]),
    kind:  new Set(["individual","mean"]),
  },
  hidden: new Set(),
  exch: "s1",
  mask: true, maskC: PAYLOAD.scatterDefault, maskW: 3,
  sd: true, peaks: false,
  norm: "none", scale: "linear", smooth: 1,
  xdom: null, table: false, theme: "auto",
};

/* ---------- signal pipeline ----------
   channel -> scatter mask -> smooth -> normalise.  Masking before normalising
   matters: otherwise the stray-light spike becomes "peak = 1" and every real
   band collapses to a few percent of full scale. */
function channel(t) {
  if (t.mode === "emission" || !t.r1) return t.s1;
  if (state.exch === "r1") return t.r1;
  if (state.exch === "corr") {
    return t.s1.map((v, i) => (t.r1[i] > 1e-6 ? v / t.r1[i] : null));
  }
  return t.s1;
}

function movavg(y, w) {
  if (w <= 1) return y;
  const h = (w - 1) / 2, out = new Array(y.length);
  for (let i = 0; i < y.length; i++) {
    let sum = 0, k = 0;
    for (let j = Math.max(0, i - h); j <= Math.min(y.length - 1, i + h); j++) {
      const v = y[j];
      if (v == null || !isFinite(v)) continue;
      sum += v; k++;
    }
    out[i] = k ? sum / k : null;
  }
  return out;
}

function normalise(y, x) {
  if (state.norm === "none") return y;
  const idx = [];
  for (let i = 0; i < y.length; i++) if (y[i] != null && isFinite(y[i])) idx.push(i);
  if (!idx.length) return y;
  let d = 0;
  if (state.norm === "peak") {
    d = Math.max(...idx.map(i => Math.abs(y[i])));
  } else {
    for (let k = 1; k < idx.length; k++) {
      const a = idx[k - 1], b = idx[k];
      d += (y[a] + y[b]) / 2 * (x[b] - x[a]);
    }
    d = Math.abs(d);
  }
  return (d > 0 && isFinite(d)) ? y.map(v => v == null ? null : v / d) : y;
}

function procY(t) {
  let y = channel(t).slice();
  /* The calibration excitation scan's 365 nm spike is a real xenon lamp line,
     not stray light, so it is never masked. */
  if (t.mode === "excitation" && state.mask && t.group !== "calibration" && state.maskW >= 0) {
    y = y.map((v, i) => Math.abs(t.x[i] - state.maskC) <= state.maskW ? null : v);
  }
  if (state.smooth > 1) y = movavg(y, state.smooth);
  return normalise(y, t.x);
}

/* ---------- group means ----------
   Computed from the *processed* replicate curves, so the mean always equals
   the mean of what is on screen (normalise each replicate, then average, is
   also the right thing for comparing band shape).  It follows the Replicate
   checkboxes, not legend clicks: legend clicks only hide a drawn line. */
function meanOf(mode, group) {
  const reps = RAW.filter(t => t.mode === mode && t.group === group &&
                              t.repKey != null && state.on.rep.has(t.repKey));
  if (!reps.length) return null;
  const prep = reps.map(t => ({ t, y: procY(t) }));

  let common = null;
  for (const p of prep) {
    const ok = new Set();
    for (let i = 0; i < p.t.n; i++) {
      const v = p.y[i];
      if (v != null && isFinite(v)) ok.add(p.t.x[i]);
    }
    common = common === null ? ok : new Set([...common].filter(v => ok.has(v)));
  }
  const x = [...common].sort((a, b) => a - b);
  if (!x.length) return null;

  const y = new Array(x.length), sd = new Array(x.length);
  for (let k = 0; k < x.length; k++) {
    const vals = prep.map(p => p.y[p.t.ix.get(x[k])]);
    const m = vals.reduce((a, b) => a + b, 0) / vals.length;
    y[k] = m;
    sd[k] = vals.length > 1
      ? Math.sqrt(vals.reduce((a, b) => a + (b - m) * (b - m), 0) / (vals.length - 1))
      : 0;
  }
  /* Replicates in a group were not always swept over the same range (Ex3 stops
     at 500 nm while Ex6/Ex9 run to 600). The mean only exists where all of
     them do, so when that is narrower than the widest replicate, the label
     says so rather than letting a short mean look like a finished one. */
  const uLo = Math.min(...reps.map(t => t.xlo));
  const uHi = Math.max(...reps.map(t => t.xhi));
  const clipped = reps.length > 1 && (x[0] > uLo + 1e-6 || x[x.length - 1] < uHi - 1e-6);

  return {
    id: `mean:${mode}:${group}`, kind: "mean", mode, group, rep: null,
    nrep: reps.length, x, y, sd, n: x.length,
    xlo: x[0], xhi: x[x.length - 1], clipped,
    label: `${group} · mean (n=${reps.length}` +
           (clipped ? `, ${fmt(x[0])}–${fmt(x[x.length - 1])} nm` : "") + ")",
    file: reps.map(t => t.id).join("+"),
  };
}

/* ---------- selection ---------- */
function passes(t) {
  return state.on.mode.has(t.mode) && state.on.group.has(t.group) &&
         (t.repKey == null || state.on.rep.has(t.repKey));
}

/** Everything the current filters admit, hidden or not (the legend pool). */
function pool() {
  const out = [];
  if (state.on.kind.has("individual")) {
    for (const t of RAW) if (passes(t)) out.push({ ser: t, kind: "individual" });
  }
  if (state.on.kind.has("mean")) {
    for (const mode of MODES) {
      if (!state.on.mode.has(mode)) continue;
      for (const g of GROUPS) {
        if (g === "calibration" || !state.on.group.has(g)) continue;
        const m = meanOf(mode, g);
        if (m) out.push({ ser: m, kind: "mean" });
      }
    }
  }
  for (const r of out) {
    r.id = r.ser.id;
    r.color = colorOf(r.ser.group);
    r.dash = dashOf(r.ser.mode);
    r.width = r.kind === "mean" ? 2.6 : 1.3;
    r.alpha = r.kind === "mean" ? 1 : (r.ser.group === "calibration" ? .8 : .55);
  }
  return out;
}

/* ---------- UI ---------- */
function swatch(color, dash, w) {
  const s = document.createElement("span");
  s.className = "sw";
  s.style.borderTop = `${w || 2.4}px ${dash && dash.length ? "dashed" : "solid"} ${color}`;
  return s;
}

for (const f of FACETS) {
  const box = document.getElementById(f.box);
  for (const v of f.vals) {
    const count = f.key === "kind"
      ? (v === "mean" ? MODES.length * (GROUPS.length - 1) : RAW.length)
      : RAW.filter(t => (f.key === "rep" ? t.repKey : t[f.key]) === v).length;
    const l = document.createElement("label"); l.className = "chk";
    const i = document.createElement("input");
    i.type = "checkbox"; i.checked = state.on[f.key].has(v); i.dataset.facet = f.key;
    i.addEventListener("change", () => {
      i.checked ? state.on[f.key].add(v) : state.on[f.key].delete(v);
      draw();
    });
    l.append(i);
    if (f.key === "group") l.append(swatch(`var(${COLOR_VAR[v]})`, [], 3));
    if (f.key === "mode")  l.append(swatch("var(--ink2)", dashOf(v), 2.4));
    l.append(document.createTextNode(f.label(v)));
    const c = document.createElement("span");
    c.className = "cnt"; c.textContent = count;
    l.append(c);
    box.append(l);
  }
}

for (const b of document.querySelectorAll("[data-all]")) {
  b.addEventListener("click", () => {
    const k = b.dataset.all, f = FACETS.find(f => f.key === k);
    const every = state.on[k].size === f.vals.length;
    state.on[k] = new Set(every ? [] : f.vals);
    for (const i of document.querySelectorAll(`#${f.box} input`)) i.checked = !every;
    draw();
  });
}

const bind = (id, ev, fn) => document.getElementById(id).addEventListener(ev, fn);
bind("o_exch", "change", e => { state.exch = e.target.value; draw(); });
bind("o_mask", "change", e => { state.mask = e.target.checked; draw(); });
bind("o_maskc", "input", e => { state.maskC = +e.target.value; draw(); });
bind("o_maskw", "input", e => { state.maskW = +e.target.value; draw(); });
bind("o_sd", "change", e => { state.sd = e.target.checked; draw(); });
bind("o_peaks", "change", e => { state.peaks = e.target.checked; draw(); });
bind("o_norm", "change", e => { state.norm = e.target.value; draw(); });
bind("o_scale", "change", e => { state.scale = e.target.value; draw(); });
bind("o_smooth", "input", e => {
  state.smooth = +e.target.value;
  document.getElementById("o_smooth_v").textContent =
    state.smooth > 1 ? state.smooth + " nm" : "off";
  draw();
});
bind("b_reset", "click", () => { state.xdom = null; state.hidden.clear(); draw(); });
bind("b_table", "click", e => {
  state.table = !state.table;
  e.currentTarget.setAttribute("aria-pressed", String(state.table));
  draw();
});
bind("b_theme", "click", e => {
  state.theme = state.theme === "auto" ? "light" : state.theme === "light" ? "dark" : "auto";
  if (state.theme === "auto") document.documentElement.removeAttribute("data-theme");
  else document.documentElement.setAttribute("data-theme", state.theme);
  e.currentTarget.textContent = "theme: " + state.theme;
  draw();
});
document.getElementById("o_maskc").value = state.maskC;

/* ---------- canvas ---------- */
const cv = document.getElementById("cv"), ctx = cv.getContext("2d");
const tip = document.getElementById("tip");
const PAD = { l: 82, r: 18, t: 14, b: 48 };
let plot = null, drag = null;

function niceTicks(lo, hi, want) {
  if (!isFinite(lo) || !isFinite(hi) || !(hi > lo)) return [lo];
  const raw = (hi - lo) / want, mag = Math.pow(10, Math.floor(Math.log10(raw)));
  const step = ([1, 2, 2.5, 5, 10].find(m => m * mag >= raw) || 10) * mag;
  if (!isFinite(step) || step <= 0) return [lo, hi];
  const out = [];
  for (let v = Math.ceil(lo / step) * step; v <= hi + step * 1e-9; v += step) out.push(v);
  return out;
}
function logTicks(lo, hi) {
  const out = [];
  if (!isFinite(lo) || !isFinite(hi) || lo <= 0 || hi <= lo) return out;
  for (let e = Math.floor(Math.log10(lo)); e <= Math.ceil(Math.log10(hi)); e++) {
    for (const m of [1, 2, 5]) {
      const v = m * Math.pow(10, e);
      if (v >= lo * .999 && v <= hi * 1.001) out.push(v);
    }
  }
  return out;
}
const fmt = v => {
  if (v == null || !isFinite(v)) return "–";
  const a = Math.abs(v);
  if (a === 0) return "0";
  if (a >= 1e5 || a < 1e-3) return v.toExponential(1).replace("e+", "e");
  if (a >= 100) return v.toFixed(0);
  if (a >= 1) return v.toFixed(2).replace(/\.?0+$/, "");
  return v.toPrecision(3).replace(/0+$/, "").replace(/\.$/, "");
};

function yLabel() {
  const exOnly = state.on.mode.has("excitation") && !state.on.mode.has("emission");
  let base;
  if (state.norm === "peak")      base = "Normalised intensity (peak = 1)";
  else if (state.norm === "area") base = "Normalised intensity (area = 1)";
  else if (exOnly && state.exch === "r1")   base = "Reference detector R1 (µA)";
  else if (exOnly && state.exch === "corr") base = "Lamp-corrected intensity (S1 ÷ R1)";
  else base = "Fluorescence intensity (CPS)";
  return base;
}

function draw() {
  readTheme();
  const all = pool();
  const vis = all.filter(r => !state.hidden.has(r.id));
  const prep = vis.map(r => ({ r, x: r.ser.x, y: r.kind === "mean" ? r.ser.y : procY(r.ser),
                               sd: r.kind === "mean" ? r.ser.sd : null }));

  let x0 = Infinity, x1 = -Infinity;
  for (const p of prep) { x0 = Math.min(x0, p.r.ser.xlo); x1 = Math.max(x1, p.r.ser.xhi); }
  if (!isFinite(x0)) { x0 = 400; x1 = 700; }
  if (state.xdom) { x0 = state.xdom[0]; x1 = state.xdom[1]; }
  if (!(x1 > x0)) x1 = x0 + 1;

  const logMode = state.scale === "log";
  let y0 = Infinity, y1 = -Infinity;
  for (const p of prep) {
    const band = state.sd && p.sd;
    for (let i = 0; i < p.y.length; i++) {
      let v = p.y[i];
      if (v == null || !isFinite(v)) continue;
      if (p.x[i] < x0 || p.x[i] > x1) continue;
      const hi = band ? v + p.sd[i] : v, lo = band ? v - p.sd[i] : v;
      if (!(logMode && hi <= 0)) y1 = Math.max(y1, hi);
      if (!(logMode && lo <= 0)) y0 = Math.min(y0, lo);
    }
  }
  const empty = !isFinite(y0);
  if (empty) { y0 = 0; y1 = 1; }
  if (y0 === y1) y1 = y0 + (Math.abs(y0) || 1) * .1;
  if (logMode) {
    if (!(y1 > 0)) { y0 = 1; y1 = 10; }
    else { if (!(y0 > 0)) y0 = y1 / 1e4; y0 /= 1.6; y1 *= 1.6; }
  } else {
    const p = (y1 - y0) * .06;
    y0 = y0 >= 0 ? 0 : y0 - p;
    y1 += p;
  }

  const dpr = window.devicePixelRatio || 1;
  const W = cv.parentElement.clientWidth || 900;
  const H = Math.max(340, Math.min(640, Math.round(W * .53)));
  cv.width = Math.round(W * dpr); cv.height = Math.round(H * dpr);
  cv.style.height = H + "px";
  ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
  ctx.clearRect(0, 0, W, H);

  const iw = W - PAD.l - PAD.r, ih = H - PAD.t - PAD.b;
  const sx = v => PAD.l + (v - x0) / (x1 - x0) * iw;
  const ly0 = logMode ? Math.log10(y0) : y0, ly1 = logMode ? Math.log10(y1) : y1;
  const sy = v => PAD.t + ih - ((logMode ? Math.log10(v) : v) - ly0) / (ly1 - ly0) * ih;

  const xt = niceTicks(x0, x1, 9);
  const yt = logMode ? logTicks(y0, y1) : niceTicks(y0, y1, 7);
  ctx.lineWidth = 1; ctx.strokeStyle = CSSV["--grid"];
  ctx.beginPath();
  for (const t of xt) { const p = Math.round(sx(t)) + .5; ctx.moveTo(p, PAD.t); ctx.lineTo(p, PAD.t + ih); }
  for (const t of yt) { const p = Math.round(sy(t)) + .5; ctx.moveTo(PAD.l, p); ctx.lineTo(PAD.l + iw, p); }
  ctx.stroke();

  ctx.strokeStyle = CSSV["--axis"]; ctx.beginPath();
  ctx.moveTo(PAD.l + .5, PAD.t); ctx.lineTo(PAD.l + .5, PAD.t + ih + .5);
  ctx.lineTo(PAD.l + iw, PAD.t + ih + .5); ctx.stroke();

  ctx.fillStyle = CSSV["--muted"];
  ctx.font = "11px system-ui, -apple-system, Segoe UI, sans-serif";
  ctx.textAlign = "center"; ctx.textBaseline = "top";
  for (const t of xt) ctx.fillText(fmt(t), sx(t), PAD.t + ih + 7);
  ctx.textAlign = "right"; ctx.textBaseline = "middle";
  for (const t of yt) ctx.fillText(fmt(t), PAD.l - 8, sy(t));

  ctx.fillStyle = CSSV["--ink"];
  ctx.font = "12px system-ui, -apple-system, Segoe UI, sans-serif";
  ctx.textAlign = "center"; ctx.textBaseline = "bottom";
  ctx.fillText("Wavelength (nm)", PAD.l + iw / 2, H - 6);
  ctx.save();
  ctx.translate(15, PAD.t + ih / 2); ctx.rotate(-Math.PI / 2);
  ctx.textBaseline = "top"; ctx.fillText(yLabel(), 0, 0);
  ctx.restore();

  ctx.save();
  ctx.beginPath(); ctx.rect(PAD.l, PAD.t - 2, iw, ih + 4); ctx.clip();
  ctx.lineJoin = "round"; ctx.lineCap = "butt";

  // ±SD ribbons sit behind every line so they never obscure a mean.
  if (state.sd) {
    for (const p of prep) {
      if (!p.sd) continue;
      ctx.globalAlpha = .16; ctx.fillStyle = p.r.color;
      let run = [];
      const flush = () => {
        if (run.length < 2) { run = []; return; }
        ctx.beginPath();
        ctx.moveTo(sx(run[0].x), sy(run[0].hi));
        for (let k = 1; k < run.length; k++) ctx.lineTo(sx(run[k].x), sy(run[k].hi));
        for (let k = run.length - 1; k >= 0; k--) ctx.lineTo(sx(run[k].x), sy(run[k].lo));
        ctx.closePath(); ctx.fill();
        run = [];
      };
      for (let i = 0; i < p.y.length; i++) {
        const v = p.y[i];
        if (v == null || !isFinite(v)) { flush(); continue; }
        let hi = v + p.sd[i], lo = v - p.sd[i];
        if (logMode) { if (hi <= 0) { flush(); continue; } if (lo <= 0) lo = y0; }
        run.push({ x: p.x[i], hi, lo });
      }
      flush();
    }
    ctx.globalAlpha = 1;
  }

  for (const pass of ["individual", "mean"]) {
    for (const p of prep) {
      if (p.r.kind !== pass) continue;
      ctx.globalAlpha = p.r.alpha;
      ctx.strokeStyle = p.r.color; ctx.lineWidth = p.r.width;
      ctx.setLineDash(p.r.dash);
      ctx.beginPath();
      let pen = false;
      for (let i = 0; i < p.y.length; i++) {
        const v = p.y[i];
        // A masked point or a non-positive value on a log axis breaks the
        // polyline rather than being bridged over -- a bridge would draw a
        // straight line through data that was deliberately removed.
        if (v == null || !isFinite(v) || (logMode && v <= 0)) { pen = false; continue; }
        const px = sx(p.x[i]), py = sy(v);
        if (!pen) { ctx.moveTo(px, py); pen = true; } else ctx.lineTo(px, py);
      }
      ctx.stroke();
    }
  }
  ctx.setLineDash([]); ctx.globalAlpha = 1;

  if (state.peaks) {
    ctx.font = "10.5px system-ui, -apple-system, Segoe UI, sans-serif";
    for (const p of prep) {
      const k = peakIndex(p, x0, x1, logMode);
      if (k < 0) continue;
      const px = sx(p.x[k]), py = sy(p.y[k]);
      ctx.fillStyle = p.r.color;
      ctx.beginPath(); ctx.arc(px, py, 3.2, 0, 6.2832); ctx.fill();
      ctx.strokeStyle = CSSV["--surface"]; ctx.lineWidth = 2;
      ctx.beginPath(); ctx.arc(px, py, 3.2, 0, 6.2832); ctx.stroke();
      ctx.textAlign = "center"; ctx.textBaseline = "bottom";
      ctx.fillText(fmt(p.x[k]), px, py - 6);
    }
  }
  ctx.restore();

  if (empty) {
    ctx.fillStyle = CSSV["--muted"];
    ctx.font = "13px system-ui, -apple-system, Segoe UI, sans-serif";
    ctx.textAlign = "center"; ctx.textBaseline = "middle";
    ctx.fillText("No traces match the current filters.", PAD.l + iw / 2, PAD.t + ih / 2);
  }
  if (drag) {
    ctx.fillStyle = CSSV["--c-g0"]; ctx.globalAlpha = .12;
    ctx.fillRect(Math.min(drag.a, drag.b), PAD.t, Math.abs(drag.b - drag.a), ih);
    ctx.globalAlpha = 1;
  }

  plot = { prep, x0, x1, sx, sy, iw, ih, W, H, logMode };
  renderLegend(all);
  renderFlags(vis);
  renderStatus(vis, prep);
  renderTable(prep, x0, x1, logMode);
}

function peakIndex(p, x0, x1, logMode) {
  let bi = -1, bv = -Infinity;
  for (let i = 0; i < p.y.length; i++) {
    const v = p.y[i];
    if (v == null || !isFinite(v)) continue;
    if (p.x[i] < x0 || p.x[i] > x1) continue;
    if (logMode && v <= 0) continue;
    if (v > bv) { bv = v; bi = i; }
  }
  return bi;
}

/* ---------- legend ---------- */
function renderLegend(all) {
  const box = document.getElementById("legend");
  box.textContent = "";
  if (!all.length) { box.textContent = "No traces match the current filters."; return; }
  for (const mode of MODES) {
    const items = all.filter(r => r.ser.mode === mode);
    if (!items.length) continue;
    const row = document.createElement("div"); row.className = "lgrp";
    const h = document.createElement("div");
    h.className = "it"; h.style.cursor = "default"; h.style.color = "var(--muted)";
    h.style.fontSize = "10.5px"; h.style.textTransform = "uppercase";
    h.style.letterSpacing = ".07em"; h.style.fontWeight = "700";
    h.textContent = mode === "emission" ? "Emission" : "Excitation";
    row.append(h);
    items.sort((a, b) => GROUPS.indexOf(a.ser.group) - GROUPS.indexOf(b.ser.group) ||
                         (a.kind === "mean" ? -1 : 1) - (b.kind === "mean" ? -1 : 1) ||
                         String(a.ser.rep).localeCompare(String(b.ser.rep)));
    for (const r of items) {
      const it = document.createElement("div");
      it.className = "it" + (state.hidden.has(r.id) ? " off" : "") +
                     (r.kind === "mean" ? " mean" : "");
      it.title = r.kind === "mean"
        ? `mean of ${r.ser.file} over ${fmt(r.ser.xlo)}–${fmt(r.ser.xhi)} nm`
        : `${r.ser.file}  (${fmt(r.ser.xlo)}–${fmt(r.ser.xhi)} nm, n=${r.ser.n})`;
      it.append(swatch(r.color, r.dash, r.kind === "mean" ? 3 : 2),
                document.createTextNode(r.ser.label));
      it.addEventListener("click", ev => {
        if (ev.shiftKey) {
          const others = all.filter(o => o.id !== r.id);
          const soloed = others.every(o => state.hidden.has(o.id));
          state.hidden = new Set(soloed ? [] : others.map(o => o.id));
        } else {
          state.hidden.has(r.id) ? state.hidden.delete(r.id) : state.hidden.add(r.id);
        }
        draw();
      });
      row.append(it);
    }
    box.append(row);
  }
}

/* ---------- flags ----------
   One y-axis only.  R1 (microamps) and S1/R1 are not the same quantity as
   emission counts, so if both are on screen unnormalised, say so and offer
   the one-click fix rather than silently plotting two units on one scale. */
function renderFlags(vis) {
  const box = document.getElementById("flags");
  box.textContent = "";
  const mixed = state.exch !== "s1" && state.norm === "none" &&
                vis.some(r => r.ser.mode === "emission") &&
                vis.some(r => r.ser.mode === "excitation");
  if (!mixed) return;
  const d = document.createElement("div"); d.className = "flag";
  const unit = state.exch === "r1" ? "R1 (µA)" : "S1 ÷ R1";
  d.append(document.createTextNode(
    `⚠ Two units on one axis: emission is CPS, excitation is ${unit}. ` +
    `Their heights are not comparable.`));
  const b = document.createElement("button");
  b.className = "btn"; b.textContent = "show excitation only";
  b.addEventListener("click", () => {
    state.on.mode = new Set(["excitation"]);
    for (const i of document.querySelectorAll('#f_mode input')) {
      i.checked = i.parentElement.textContent.includes("Excitation");
    }
    draw();
  });
  d.append(b);
  box.append(d);
}

function renderStatus(vis, prep) {
  const bits = [];
  const nm = vis.filter(r => r.kind === "mean").length;
  bits.push(`${vis.length} traces (${vis.length - nm} scans, ${nm} means)`);
  if (state.mask && state.on.mode.has("excitation")) {
    bits.push(`scatter masked ${fmt(state.maskC - state.maskW)}–${fmt(state.maskC + state.maskW)} nm`);
  }
  const clipped = vis.filter(r => r.kind === "mean" && r.ser.clipped);
  if (clipped.length) {
    bits.push(`⚠ ${clipped.length} mean${clipped.length > 1 ? "s" : ""} limited ` +
              `to the range all replicates share`);
  }
  if (state.norm !== "none") {
    bits.push("each scan normalised before averaging");
  }
  bits.push(state.xdom
    ? `zoomed ${fmt(state.xdom[0])}–${fmt(state.xdom[1])} nm · double-click to reset`
    : "drag to zoom · shift-click a legend entry to solo");
  document.getElementById("status").textContent = bits.join("  —  ");
}

/* ---------- table view ----------
   Also discharges the contrast relief rule: every series is named in text
   with its numbers, so nothing depends on telling two hues apart. */
function renderTable(prep, x0, x1, logMode) {
  const box = document.getElementById("tablebox");
  box.textContent = "";
  if (!state.table) return;
  const tb = document.createElement("table"); tb.className = "tbl";
  const head = ["Trace", "Measurement", "Range (nm)", "pts",
                "λₘₐₓ (nm)", "peak", "area", "source"];
  const thead = document.createElement("thead"), hr = document.createElement("tr");
  for (const h of head) { const th = document.createElement("th"); th.textContent = h; hr.append(th); }
  thead.append(hr); tb.append(thead);
  const tbody = document.createElement("tbody");
  for (const p of prep) {
    const k = peakIndex(p, x0, x1, logMode);
    let area = 0, prev = -1;
    for (let i = 0; i < p.y.length; i++) {
      const v = p.y[i];
      if (v == null || !isFinite(v)) { prev = -1; continue; }
      if (prev >= 0) area += (p.y[prev] + v) / 2 * (p.x[i] - p.x[prev]);
      prev = i;
    }
    const tr = document.createElement("tr");
    const c0 = document.createElement("td"); c0.className = "name";
    c0.append(swatch(p.r.color, p.r.dash, p.r.kind === "mean" ? 3 : 2),
              document.createTextNode(p.r.ser.label));
    tr.append(c0);
    for (const v of [p.r.ser.mode,
                     `${fmt(p.r.ser.xlo)}–${fmt(p.r.ser.xhi)}`,
                     String(p.r.ser.n),
                     k < 0 ? "–" : fmt(p.x[k]),
                     k < 0 ? "–" : fmt(p.y[k]),
                     fmt(area),
                     p.r.ser.file]) {
      const td = document.createElement("td"); td.textContent = v; tr.append(td);
    }
    tbody.append(tr);
  }
  tb.append(tbody);
  box.append(tb);
}

/* ---------- hover ---------- */
cv.addEventListener("mousemove", ev => {
  if (!plot) return;
  const r = cv.getBoundingClientRect();
  const mx = ev.clientX - r.left, my = ev.clientY - r.top;
  if (drag) { drag.b = Math.max(PAD.l, Math.min(plot.W - PAD.r, mx)); draw(); return; }
  if (mx < PAD.l || mx > plot.W - PAD.r || my < PAD.t || my > PAD.t + plot.ih) {
    tip.style.display = "none"; return;
  }
  const xv = plot.x0 + (mx - PAD.l) / plot.iw * (plot.x1 - plot.x0);
  const rows = [];
  for (const p of plot.prep) {
    if (xv < p.r.ser.xlo - 1 || xv > p.r.ser.xhi + 1) continue;
    let bi = -1, bd = Infinity;
    for (let i = 0; i < p.x.length; i++) {
      const d = Math.abs(p.x[i] - xv);
      if (d < bd) { bd = d; bi = i; }
    }
    if (bi < 0 || bd > 2) continue;
    const v = p.y[bi];
    if (v == null || !isFinite(v)) continue;
    rows.push({ p, xv: p.x[bi], v, sd: p.sd ? p.sd[bi] : null,
                dy: Math.abs(plot.sy(v) - my) });
  }
  if (!rows.length) { tip.style.display = "none"; return; }
  rows.sort((a, b) => b.v - a.v);
  const near = rows.reduce((a, b) => (b.dy < a.dy ? b : a));

  tip.textContent = "";
  const h = document.createElement("b");
  h.textContent = `${fmt(near.xv)} nm`;
  tip.append(h);
  for (const q of rows.slice(0, 12)) {
    const d = document.createElement("div");
    d.className = "r"; d.style.fontWeight = q === near ? "700" : "400";
    d.append(swatch(q.p.r.color, q.p.r.dash, q.p.r.kind === "mean" ? 3 : 2));
    const e = document.createElement("em"); e.textContent = q.p.r.ser.label;
    const s = document.createElement("span");
    s.textContent = fmt(q.v) + (state.sd && q.sd ? ` ± ${fmt(q.sd)}` : "");
    d.append(e, s); tip.append(d);
  }
  if (rows.length > 12) {
    const d = document.createElement("div");
    d.className = "r"; d.style.color = "var(--muted)";
    d.textContent = `+${rows.length - 12} more`;
    tip.append(d);
  }
  tip.style.display = "block";
  const tw = tip.offsetWidth, th = tip.offsetHeight;
  tip.style.left = Math.max(4, Math.min(mx + 14, plot.W - tw - 4)) + "px";
  tip.style.top  = Math.max(4, Math.min(my - th / 2, plot.H - th - 4)) + "px";
});
cv.addEventListener("mouseleave", () => { tip.style.display = "none"; });

cv.addEventListener("mousedown", ev => {
  const r = cv.getBoundingClientRect(), mx = ev.clientX - r.left;
  if (!plot || mx < PAD.l || mx > plot.W - PAD.r) return;
  drag = { a: mx, b: mx }; tip.style.display = "none";
});
window.addEventListener("mouseup", () => {
  if (!drag || !plot) { drag = null; return; }
  const { a, b } = drag; drag = null;
  if (Math.abs(b - a) > 6) {
    const toData = px => plot.x0 + (px - PAD.l) / plot.iw * (plot.x1 - plot.x0);
    state.xdom = [toData(Math.min(a, b)), toData(Math.max(a, b))];
  }
  draw();
});
cv.addEventListener("dblclick", () => { state.xdom = null; draw(); });

/* ---------- export ---------- */
function download(name, blob) {
  const a = document.createElement("a");
  a.href = URL.createObjectURL(blob); a.download = name;
  a.click(); setTimeout(() => URL.revokeObjectURL(a.href), 2000);
}
bind("b_png", "click", () => {
  // The live canvas is transparent; most viewers render that black.
  const out = document.createElement("canvas");
  out.width = cv.width; out.height = cv.height;
  const o = out.getContext("2d");
  o.fillStyle = CSSV["--surface"] || "#fff";
  o.fillRect(0, 0, out.width, out.height);
  o.drawImage(cv, 0, 0);
  out.toBlob(b => download("chimera_spectra.png", b));
});
bind("b_csv", "click", () => {
  if (!plot || !plot.prep.length) return;
  const cols = [];
  for (const p of plot.prep) {
    cols.push({ head: p.r.ser.label.replace(/[",\n]/g, " "), p, sd: false });
    if (p.sd && state.sd) cols.push({ head: p.r.ser.label.replace(/[",\n]/g, " ") + " SD", p, sd: true });
  }
  const grid = [...new Set(plot.prep.flatMap(p => p.x))].sort((a, b) => a - b);
  const lines = [["Wavelength_nm", ...cols.map(c => c.head)].join(",")];
  for (const x of grid) {
    const row = [x];
    for (const c of cols) {
      const i = c.p.x.indexOf(x);
      const v = i < 0 ? null : (c.sd ? c.p.sd[i] : c.p.y[i]);
      row.push(v == null || !isFinite(v) ? "" : v);
    }
    lines.push(row.join(","));
  }
  download("chimera_spectra_visible.csv",
           new Blob([lines.join("\n")], { type: "text/csv" }));
});

/* ---------- notes ---------- */
{
  const n = document.getElementById("notes");
  const esc = t => t.replace(/[&<>]/g, c => ({ "&":"&amp;", "<":"&lt;", ">":"&gt;" }[c]));
  n.innerHTML = "<b>Reading these data</b><ul>" +
    PAYLOAD.notes.map(t => `<li>${esc(t)}</li>`).join("") + "</ul>";
}
document.getElementById("subtitle").textContent =
  `${RAW.length} scans · ${GROUPS.length - 1} constructs × 3 replicates · ${PAYLOAD.source}`;

let rt;
window.addEventListener("resize", () => { clearTimeout(rt); rt = setTimeout(draw, 90); });
matchMedia("(prefers-color-scheme: dark)").addEventListener("change", () => draw());
draw();
</script>
</html>
"""


def render(payload, title):
    blob = json.dumps(payload, separators=(",", ":"), ensure_ascii=True)
    blob = blob.replace("</", "<\\/")          # never close the host <script>
    return (HTML
            .replace("__PAYLOAD__", blob)
            .replace("__TITLE__", title))


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--root", type=Path, default=DEFAULT_ROOT,
                    help="folder holding the Dflt*_Data.csv files")
    ap.add_argument("--out", type=Path, default=DEFAULT_OUT,
                    help="output HTML file")
    ap.add_argument("--title", default="sfGFP–DNA chimera · excitation & emission overlay")
    args = ap.parse_args(argv)

    root = args.root.resolve()
    if not root.is_dir():
        sys.exit(f"not a directory: {root}")

    payload = build_payload(root)
    if not payload["traces"]:
        sys.exit(f"no Dflt*_Data.csv files found under {root}")

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(render(payload, args.title), encoding="utf-8")

    em = sum(1 for t in payload["traces"] if t["mode"] == "emission")
    ex = len(payload["traces"]) - em
    print(f"{args.out}  ←  {em} emission + {ex} excitation scans")
    for n in payload["notes"]:
        print("  · " + n)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
