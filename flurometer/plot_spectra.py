#!/usr/bin/env python3
"""
Superimpose every fluorometer spectrum in this folder onto one interactive graph.

Walks the fluorometer directory tree, parses each exported CSV plus the metadata
encoded in its filename, and writes ONE self-contained HTML file with the data
embedded.  No dependencies -- not even numpy -- and no network access: the page
is pure vanilla JS and opens straight from disk.

Typical use
-----------
    python3 plot_spectra.py                 # write spectra_viewer.html next to this script
    python3 plot_spectra.py --open          # ...and open it in a browser
    python3 plot_spectra.py --list          # just print what was discovered
    python3 plot_spectra.py --root ../other-runs --out /tmp/v.html

Adding data
-----------
Drop a new run folder alongside the existing ones and re-run.  Nothing in this
script is hard-coded to the current experiments; the only tunable is
CONSTRUCT_HINTS below, which maps a substring of a run-folder name to a pretty
construct label.

Filename convention understood (all parts after the sample are optional):

    <sample>_Ex<excitation>_Em<lo>-<hi>[_<extra>].csv

      0_4_Ex485_Em500-600.csv           -> sample 0.4,  Ex 485,   Em 500-600
      blank_Ex485_Em500-600.csv         -> blank
      calibration_Ex350_5_Em365-450_5.csv -> calibration, Ex 350.5
      4_0_Ex485_Em500-600_2tetGFP.csv   -> sample 4.0, replicate 2, tag "tetGFP"

Underscores inside numbers are decimal points (`0_05` -> 0.05, `4_0` -> 4.0),
matching how the instrument software sanitises filenames.  Files whose names do
not match are still plotted -- they just get the stem as their label.
"""

import argparse
import csv
import hashlib
import json
import math
import re
import sys
import webbrowser
from pathlib import Path

# Substring of a run-folder (or filename tag) -> human-readable construct label.
# First match wins, so put the more specific patterns first.
CONSTRUCT_HINTS = [
    ("ntet", "nTET sfGFP"),
    ("tetgfp", "nTET sfGFP"),
    ("wtsfgfp", "WT sfGFP"),
    ("wt", "WT sfGFP"),
]

# Sample tokens that mean "this is not a dilution of the analyte".
ROLE_TOKENS = {
    "blank": "blank",
    "buffer": "blank",
    "water": "blank",
    "calibration": "calibration",
    "calib": "calibration",
    "raman": "calibration",
    "standard": "calibration",
}

FNAME_RE = re.compile(
    r"^(?P<sample>.+?)"
    r"_Ex(?P<ex>\d+(?:_\d+)?)"
    r"_Em(?P<lo>\d+)-(?P<hi>\d+(?:_\d+)?)"
    r"(?:_(?P<extra>.+))?$",
    re.IGNORECASE,
)


# --------------------------------------------------------------------------
# parsing
# --------------------------------------------------------------------------

def _num(token):
    """'0_05' -> 0.05, '4_0' -> 4.0, '485' -> 485.0.  None if not numeric."""
    try:
        v = float(token.replace("_", "."))
    except ValueError:
        return None
    # float() happily accepts 'nan'/'inf', which would later blow up
    # json.dumps(allow_nan=False).  Treat them as missing.
    return v if math.isfinite(v) else None


def read_csv(path):
    """Return (channels, xs, {channel: ys}) from an instrument CSV export.

    Handles the two-row header the fluorometer writes (names, then units) and
    any number of signal columns (S1, S2, R1, ...).  Rows that are not fully
    numeric are skipped rather than fatal, so a stray footer will not kill a run.
    """
    with open(path, newline="", encoding="utf-8-sig", errors="replace") as fh:
        rows = [r for r in csv.reader(fh) if any(c.strip() for c in r)]
    if len(rows) < 2:
        raise ValueError("fewer than two rows")

    header = [c.strip() for c in rows[0]]
    body = rows[1:]
    units = None
    # A second header row is present iff its first cell is not a number.
    if body and _num(body[0][0].strip()) is None:
        units = [c.strip() for c in body[0]]
        body = body[1:]

    channels = header[1:] or ["S1"]
    xs, ys = [], {c: [] for c in channels}
    for row in body:
        cells = [c.strip() for c in row]
        x = _num(cells[0]) if cells else None
        if x is None:
            continue
        vals = []
        for i, _ in enumerate(channels, start=1):
            vals.append(_num(cells[i]) if i < len(cells) else None)
        if all(v is None for v in vals):
            continue
        xs.append(x)
        for c, v in zip(channels, vals):
            ys[c].append(v)
    if not xs:
        raise ValueError("no numeric data rows")
    y_unit = units[1] if units and len(units) > 1 else "CPS"
    x_unit = units[0] if units else "nm"
    return channels, xs, ys, x_unit, y_unit


def parse_name(stem):
    """Pull sample / excitation / emission / replicate / tag out of a filename."""
    meta = {
        "sample": stem, "sample_value": None, "role": "sample",
        "ex": None, "em_lo": None, "em_hi": None,
        "replicate": None, "tag": None, "parsed": False,
    }
    m = FNAME_RE.match(stem)
    if not m:
        return meta
    meta["parsed"] = True

    sample = m.group("sample")
    meta["sample"] = sample
    meta["role"] = ROLE_TOKENS.get(sample.lower().strip("_ "), "sample")
    meta["sample_value"] = _num(sample)
    meta["ex"] = _num(m.group("ex"))
    meta["em_lo"] = _num(m.group("lo"))
    meta["em_hi"] = _num(m.group("hi"))

    extra = m.group("extra")
    if extra:
        # "1tetGFP" -> replicate 1, tag "tetGFP";  "rep2" -> replicate 2
        em = re.match(r"^(?:rep|r|n)?(\d+)(.*)$", extra, re.IGNORECASE)
        if em:
            meta["replicate"] = int(em.group(1))
            meta["tag"] = em.group(2).strip("_-") or None
        else:
            meta["tag"] = extra.strip("_-") or None
    return meta


def guess_construct(*texts):
    for text in texts:
        if not text:
            continue
        low = text.lower()
        for needle, label in CONSTRUCT_HINTS:
            if needle in low:
                return label
    return None


def pretty_run(folder):
    """'260805_nTETsfGFPdiluation' -> ('260805_nTETsfGFPdiluation', '2026-08-05')."""
    m = re.match(r"^(\d{2})(\d{2})(\d{2})[_-]?", folder)
    date = f"20{m.group(1)}-{m.group(2)}-{m.group(3)}" if m else None
    return folder, date


def uniform_axis(xs):
    """Return (x0, dx, n) if xs is an evenly spaced grid, else None."""
    if len(xs) < 2:
        return None
    dx = xs[1] - xs[0]
    if dx == 0:
        return None
    for i in range(1, len(xs)):
        if abs((xs[i] - xs[i - 1]) - dx) > 1e-9:
            return None
    return xs[0], dx, len(xs)


# --------------------------------------------------------------------------
# discovery
# --------------------------------------------------------------------------

def discover(root, keep_duplicates=False, verbose=True):
    """Walk `root` and return (series, notes).  One series per signal column."""
    root = Path(root).resolve()
    series, notes = [], []
    by_content = {}

    for path in sorted(root.rglob("*.csv")):
        rel = path.relative_to(root)
        try:
            channels, xs, ys, x_unit, y_unit = read_csv(path)
        except Exception as exc:              # a bad file must not abort the run
            notes.append(f"skipped {rel}: {exc}")
            continue

        meta = parse_name(path.stem)
        parts = rel.parts
        run_folder = parts[0] if len(parts) > 1 else root.name
        run, run_date = pretty_run(run_folder)
        subdir = "/".join(parts[1:-1])
        construct = guess_construct(meta["tag"], run_folder) or "unassigned"

        for chan in channels:
            yv = ys[chan]
            # Content fingerprint: identical numbers => identical series, no
            # matter what the file is called.  Catches the instrument's
            # default-named exports (s.csv, massexptest/Data.csv).
            digest = hashlib.md5(
                repr([(round(a, 6), None if b is None else round(b, 6))
                      for a, b in zip(xs, yv)]).encode()
            ).hexdigest()

            # Sweep group keyed on the OBSERVED range, not the filename's claim.
            lo, hi = min(xs), max(xs)
            ex_txt = "?" if meta["ex"] is None else f"{meta['ex']:g}"
            sweep = f"Ex {ex_txt} / Em {lo:g}–{hi:g}"

            # '0_05' is the instrument's filename-safe spelling of 0.05.  Undo
            # the substitution rather than reformatting, so '4_0' stays '4.0'
            # (the notebook's notation) instead of collapsing to '4'.
            sample_label = (meta["sample"].replace("_", ".")
                            if meta["sample_value"] is not None else meta["sample"])

            label_bits = [run_date or run, sample_label]
            if meta["replicate"] is not None:
                label_bits.append(f"rep {meta['replicate']}")
            if len(channels) > 1:
                label_bits.append(chan)

            rec = {
                "id": f"{rel.as_posix()}::{chan}",
                "label": " · ".join(str(b) for b in label_bits),
                "path": rel.as_posix(),
                "run": run,
                "run_date": run_date,
                "subdir": subdir,
                "construct": construct,
                "sample": sample_label,
                "sample_raw": meta["sample"],
                "sample_value": meta["sample_value"],
                "role": meta["role"],
                "replicate": meta["replicate"],
                "channel": chan,
                "ex": meta["ex"],
                "declared_em": (
                    None if meta["em_lo"] is None
                    else f"{meta['em_lo']:g}–{meta['em_hi']:g}"
                ),
                "sweep": sweep,
                "x_unit": x_unit,
                "y_unit": y_unit,
                "parsed": meta["parsed"],
                "_digest": digest,
                "_xs": xs,
                "_ys": yv,
            }

            if not keep_duplicates and digest in by_content:
                keeper = by_content[digest]
                # Prefer the better-named copy: parsed filename, then longer
                # name, then shallower path.
                challenger_score = (rec["parsed"], len(path.stem), -len(parts))
                keeper_score = (keeper["parsed"], len(Path(keeper["path"]).stem),
                                -len(Path(keeper["path"]).parts))
                if challenger_score > keeper_score:
                    notes.append(
                        f"duplicate: {rec['path']} == {keeper['path']} "
                        f"-> keeping {rec['path']}")
                    series.remove(keeper)
                    by_content[digest] = rec
                    series.append(rec)
                else:
                    notes.append(
                        f"duplicate: {rec['path']} == {keeper['path']} "
                        f"-> dropped {rec['path']}")
                continue

            by_content[digest] = rec
            series.append(rec)

    # Stable, meaningful ordering: run, construct, sweep, concentration.
    series.sort(key=lambda s: (
        s["run"], s["construct"], s["sweep"],
        {"blank": -2, "calibration": -1}.get(s["role"], 0),
        s["sample_value"] if s["sample_value"] is not None else math.inf,
        s["replicate"] or 0, s["channel"],
    ))
    return series, notes


def to_payload(series):
    """Strip private fields and compress the wavelength axis where possible."""
    out = []
    for s in series:
        rec = {k: v for k, v in s.items() if not k.startswith("_")}
        xs = s["_xs"]
        axis = uniform_axis(xs)
        if axis:
            rec["x0"], rec["dx"], rec["n"] = axis[0], axis[1], axis[2]
        else:
            rec["x"] = xs
            rec["n"] = len(xs)
        rec["y"] = s["_ys"]
        out.append(rec)
    return out


# --------------------------------------------------------------------------
# HTML emitter
# --------------------------------------------------------------------------

HTML = r"""<!doctype html>
<meta charset="utf-8">
<title>__TITLE__</title>
<style>
  :root {
    --bg:#ffffff; --panel:#f7f8fa; --line:#dfe3e8; --ink:#1c2530;
    --muted:#63707f; --accent:#2563eb;
  }
  * { box-sizing:border-box; }
  body {
    margin:0; background:var(--bg); color:var(--ink);
    font:13px/1.45 -apple-system, "Segoe UI", Roboto, "Helvetica Neue", Arial, sans-serif;
  }
  header {
    padding:14px 20px; border-bottom:1px solid var(--line);
    display:flex; align-items:baseline; gap:14px; flex-wrap:wrap;
  }
  header h1 { margin:0; font-size:16px; font-weight:650; letter-spacing:-0.01em; }
  header .sub { color:var(--muted); font-size:12px; }
  .wrap { display:flex; align-items:stretch; min-height:calc(100vh - 52px); }
  aside {
    width:280px; flex:0 0 280px; padding:16px; background:var(--panel);
    border-right:1px solid var(--line); overflow-y:auto;
    max-height:calc(100vh - 52px);
  }
  main { flex:1 1 auto; padding:16px 20px 24px; min-width:0; }
  .grp { margin-bottom:18px; }
  .grp > h2 {
    margin:0 0 7px; font-size:10.5px; font-weight:700; text-transform:uppercase;
    letter-spacing:.07em; color:var(--muted);
    display:flex; justify-content:space-between; align-items:center;
  }
  .grp > h2 button {
    font-size:10px; text-transform:none; letter-spacing:0; padding:1px 6px;
    border:1px solid var(--line); background:#fff; border-radius:4px;
    color:var(--muted); cursor:pointer;
  }
  .grp > h2 button:hover { color:var(--accent); border-color:var(--accent); }
  label.chk {
    display:flex; align-items:center; gap:7px; padding:2.5px 0;
    cursor:pointer; user-select:none;
  }
  label.chk:hover { color:var(--accent); }
  label.chk input { margin:0; accent-color:var(--accent); flex:0 0 auto; }
  label.chk .cnt { margin-left:auto; color:var(--muted); font-size:11px; font-variant-numeric:tabular-nums; }
  select, input[type=range] { width:100%; }
  select {
    padding:5px 7px; border:1px solid var(--line); border-radius:5px;
    background:#fff; color:var(--ink); font:inherit; font-size:12.5px;
  }
  .row { display:flex; gap:8px; align-items:center; margin-top:6px; }
  .row > span { color:var(--muted); font-size:11.5px; white-space:nowrap; }
  .btn {
    padding:5px 10px; border:1px solid var(--line); background:#fff;
    border-radius:5px; font:inherit; font-size:12px; cursor:pointer; color:var(--ink);
  }
  .btn:hover { border-color:var(--accent); color:var(--accent); }
  #plotbox { position:relative; }
  canvas { display:block; width:100%; touch-action:none; cursor:crosshair; }
  #tip {
    position:absolute; pointer-events:none; display:none; z-index:5;
    background:rgba(255,255,255,.97); border:1px solid var(--line);
    border-radius:6px; padding:7px 9px; font-size:11.5px; max-width:290px;
    box-shadow:0 4px 14px rgba(0,0,0,.10);
  }
  #tip b { display:block; margin-bottom:4px; font-variant-numeric:tabular-nums; }
  #tip .r { display:flex; gap:6px; align-items:center; white-space:nowrap;
            font-variant-numeric:tabular-nums; }
  #tip .r i { width:9px; height:2.5px; border-radius:2px; flex:0 0 auto; }
  #tip .r em { font-style:normal; color:var(--muted); overflow:hidden;
               text-overflow:ellipsis; }
  #tip .r span { margin-left:auto; font-weight:600; }
  #legend {
    margin-top:12px; display:flex; flex-wrap:wrap; gap:3px 16px;
    border-top:1px solid var(--line); padding-top:11px;
  }
  #legend .it {
    display:flex; align-items:center; gap:7px; cursor:pointer;
    padding:2px 4px; border-radius:4px; font-size:12px;
  }
  #legend .it:hover { background:var(--panel); }
  #legend .it.off { opacity:.34; }
  #legend .it .sw { width:20px; height:0; border-top-width:2.5px; flex:0 0 auto; }
  #status { margin-top:10px; color:var(--muted); font-size:11.5px; }
  #notes {
    margin-top:14px; padding:9px 11px; background:#fffbeb;
    border:1px solid #fde68a; border-radius:6px; font-size:11.5px; color:#78350f;
  }
  #notes ul { margin:5px 0 0; padding-left:17px; }
  @media (max-width:860px) {
    .wrap { flex-direction:column; }
    aside { width:auto; flex:none; max-height:none; border-right:0;
            border-bottom:1px solid var(--line); }
  }
</style>

<header>
  <h1>__TITLE__</h1>
  <span class="sub" id="subtitle"></span>
</header>

<div class="wrap">
<aside>
  <div class="grp">
    <h2>Emission sweep</h2>
    <div id="f_sweep"></div>
  </div>

  <div class="grp">
    <h2>Processing</h2>
    <label class="chk"><input type="checkbox" id="o_blank"> subtract blank</label>
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
      <span id="o_smooth_v" style="min-width:3.2em;text-align:right">off</span>
    </div>
  </div>

  <div class="grp"><h2>Run <button data-all="run">all</button></h2><div id="f_run"></div></div>
  <div class="grp"><h2>Construct <button data-all="construct">all</button></h2><div id="f_construct"></div></div>
  <div class="grp"><h2>Sample <button data-all="sample">all</button></h2><div id="f_sample"></div></div>
  <div class="grp" id="g_replicate"><h2>Replicate <button data-all="replicate">all</button></h2><div id="f_replicate"></div></div>
  <div class="grp" id="g_channel"><h2>Channel <button data-all="channel">all</button></h2><div id="f_channel"></div></div>

  <div class="grp">
    <h2>Export</h2>
    <div class="row">
      <button class="btn" id="b_png">PNG</button>
      <button class="btn" id="b_csv">CSV</button>
      <button class="btn" id="b_reset">reset view</button>
    </div>
  </div>
</aside>

<main>
  <div id="plotbox"><canvas id="cv"></canvas><div id="tip"></div></div>
  <div id="legend"></div>
  <div id="status"></div>
  <div id="notes" hidden></div>
</main>
</div>

<script>
"use strict";
const PAYLOAD = __PAYLOAD__;
const SERIES = PAYLOAD.series;
const NOTES  = PAYLOAD.notes;

/* ---------- wavelength axis: stored as {x0,dx,n} when uniform ---------- */
for (const s of SERIES) {
  if (!s.x) { s.x = new Array(s.n); for (let i = 0; i < s.n; i++) s.x[i] = s.x0 + i * s.dx; }
  s.xlo = s.x[0]; s.xhi = s.x[s.n - 1];
}

/* ---------- colour: hue by construct, lightness by concentration -------
   Encoding the dilution series as a lightness ramp inside one hue keeps a
   titration readable at a glance, and keeps constructs separable even when
   printed greyscale (dash pattern carries the construct too).            */
const HUES = [214, 25, 145, 288, 45, 190, 330, 110];
const constructs = [...new Set(SERIES.map(s => s.construct))].sort();
const hueOf = Object.fromEntries(constructs.map((c, i) => [c, HUES[i % HUES.length]]));
const dashOf = Object.fromEntries(constructs.map((c, i) => [c, [[], [7, 4], [2, 3], [10, 3, 2, 3]][i % 4]]));

for (const c of constructs) {
  const fam = SERIES.filter(s => s.construct === c && s.role === "sample" && s.sample_value != null);
  const vals = [...new Set(fam.map(s => s.sample_value))].sort((a, b) => a - b);
  for (const s of SERIES.filter(s => s.construct === c)) {
    if (s.role !== "sample") { s.color = s.role === "blank" ? "#9aa5b1" : "#5f6b78"; }
    else if (s.sample_value == null || vals.length < 2) { s.color = `hsl(${hueOf[c]} 70% 45%)`; }
    else {
      const t = vals.indexOf(s.sample_value) / (vals.length - 1);
      s.color = `hsl(${hueOf[c]} ${58 + 24 * t}% ${72 - 40 * t}%)`;
    }
    // Reference traces get neutral greys and their own dash so they never read
    // as a member of the dilution series they sit next to.
    s.dash = s.role === "sample" ? dashOf[c]
           : s.role === "blank" ? [2, 3] : [6, 3, 2, 3];
    s.width = s.replicate && s.replicate > 1 ? 1.4 : 2.1;
  }
}

/* ---------- state ---------- */
const facets = ["run", "construct", "sample", "replicate", "channel"];
const facetVals = {};
for (const f of facets) {
  const seen = new Map();
  for (const s of SERIES) {
    const k = keyOf(s, f);
    seen.set(k, (seen.get(k) || 0) + 1);
  }
  facetVals[f] = [...seen.entries()].sort((a, b) => cmpKey(f, a[0], b[0]));
}
function keyOf(s, f) {
  if (f === "replicate") return s.replicate == null ? "–" : String(s.replicate);
  return String(s[f]);
}
function cmpKey(f, a, b) {
  if (f === "sample") {
    const va = sampleOrder(a), vb = sampleOrder(b);
    if (va !== vb) return va - vb;
  }
  return a.localeCompare(b, undefined, { numeric: true });
}
function sampleOrder(name) {
  const s = SERIES.find(s => String(s.sample) === name);
  if (!s) return 0;
  if (s.role === "blank") return -2;
  if (s.role === "calibration") return -1;
  return s.sample_value == null ? 1e9 : s.sample_value;
}

const sweeps = [...new Set(SERIES.map(s => s.sweep))].sort();
const state = {
  sweep: "__ALL__",
  on: Object.fromEntries(facets.map(f => [f, new Set(facetVals[f].map(e => e[0]))])),
  hidden: new Set(),
  blank: false, peaks: false, norm: "none", scale: "linear", smooth: 1,
  xdom: null,
};

/* ---------- blank lookup: same sweep, prefer same run ---------- */
function blankFor(s) {
  const c = SERIES.filter(b => b.role === "blank" && b.sweep === s.sweep && b.channel === s.channel);
  return c.find(b => b.run === s.run) || c[0] || null;
}

/* ---------- transforms ---------- */
function movavg(y, w) {
  if (w <= 1) return y;
  const h = (w - 1) / 2, out = new Array(y.length);
  for (let i = 0; i < y.length; i++) {
    let sum = 0, k = 0;
    for (let j = Math.max(0, i - h); j <= Math.min(y.length - 1, i + h); j++) {
      if (y[j] == null || !isFinite(y[j])) continue;
      sum += y[j]; k++;
    }
    out[i] = k ? sum / k : null;
  }
  return out;
}

function transform(s) {
  let y = s.y.map(v => (v == null || !isFinite(v)) ? null : v);

  // Subtract from *every* trace including the blanks themselves (which then
  // sit at zero).  Leaving blanks uncorrected would make them the tallest
  // curve on a blank-subtracted plot, which reads as a bug.
  if (state.blank) {
    const b = blankFor(s);
    if (b) {
      // Interpolate the blank onto this series' grid so a 500-600 blank can
      // still correct the overlapping part of a 380-700 sweep.
      y = y.map((v, i) => {
        if (v == null) return null;
        const bv = interp(b.x, b.y, s.x[i]);
        return bv == null ? null : v - bv;
      });
    }
  }
  if (state.smooth > 1) y = movavg(y, state.smooth);

  if (state.norm !== "none") {
    const idx = [];
    for (let i = 0; i < y.length; i++) if (y[i] != null) idx.push(i);
    let d = 0;
    if (state.norm === "peak") {
      d = Math.max(...idx.map(i => Math.abs(y[i])));
    } else {                       // trapezoidal area over the series' own span
      for (let k = 1; k < idx.length; k++) {
        const i0 = idx[k - 1], i1 = idx[k];
        d += (y[i0] + y[i1]) / 2 * (s.x[i1] - s.x[i0]);
      }
      d = Math.abs(d);
    }
    if (d > 0 && isFinite(d)) y = y.map(v => v == null ? null : v / d);
  }
  return y;
}

function interp(xs, ys, x) {
  if (x < xs[0] || x > xs[xs.length - 1]) return null;
  let lo = 0, hi = xs.length - 1;
  while (hi - lo > 1) { const m = (lo + hi) >> 1; if (xs[m] <= x) lo = m; else hi = m; }
  const a = ys[lo], b = ys[hi];
  if (a == null || b == null) return a == null ? b : a;
  const t = xs[hi] === xs[lo] ? 0 : (x - xs[lo]) / (xs[hi] - xs[lo]);
  return a + (b - a) * t;
}

function visible() {
  return SERIES.filter(s =>
    (state.sweep === "__ALL__" || s.sweep === state.sweep) &&
    facets.every(f => state.on[f].has(keyOf(s, f))) &&
    !state.hidden.has(s.id));
}

/* ---------- UI construction ---------- */
function chk(label, checked, count, onChange) {
  const l = document.createElement("label"); l.className = "chk";
  const i = document.createElement("input"); i.type = "checkbox"; i.checked = checked;
  i.addEventListener("change", () => onChange(i.checked));
  l.append(i, document.createTextNode(label));
  if (count != null) {
    const c = document.createElement("span"); c.className = "cnt"; c.textContent = count;
    l.append(c);
  }
  return l;
}

{ // sweep selector -- "all" superimposes every sweep on one axis
  const box = document.getElementById("f_sweep");
  const sel = document.createElement("select");
  const all = document.createElement("option");
  all.value = "__ALL__"; all.textContent = `all sweeps (${sweeps.length})`;
  sel.append(all);
  for (const sw of sweeps) {
    const o = document.createElement("option");
    o.value = sw; o.textContent = `${sw} nm  (${SERIES.filter(s => s.sweep === sw).length})`;
    sel.append(o);
  }
  sel.addEventListener("change", () => { state.sweep = sel.value; state.xdom = null; draw(); });
  box.append(sel);
}

for (const f of facets) {
  const box = document.getElementById("f_" + f);
  const entries = facetVals[f];
  if (entries.length < 2) {                       // a one-value facet is noise
    const g = document.getElementById("g_" + f);
    if (g) { g.hidden = true; continue; }
  }
  for (const [val, n] of entries) {
    box.append(chk(val, true, n, on => {
      on ? state.on[f].add(val) : state.on[f].delete(val);
      draw();
    }));
  }
}

for (const btn of document.querySelectorAll("[data-all]")) {
  btn.addEventListener("click", () => {
    const f = btn.dataset.all;
    const every = state.on[f].size === facetVals[f].length;
    state.on[f] = new Set(every ? [] : facetVals[f].map(e => e[0]));
    for (const i of document.querySelectorAll(`#f_${f} input`)) i.checked = !every;
    draw();
  });
}

const bind = (id, ev, fn) => document.getElementById(id).addEventListener(ev, fn);
bind("o_blank", "change", e => { state.blank = e.target.checked; draw(); });
bind("o_peaks", "change", e => { state.peaks = e.target.checked; draw(); });
bind("o_norm", "change", e => { state.norm = e.target.value; draw(); });
bind("o_scale", "change", e => { state.scale = e.target.value; draw(); });
bind("o_smooth", "input", e => {
  state.smooth = +e.target.value;
  document.getElementById("o_smooth_v").textContent =
    state.smooth > 1 ? state.smooth + " nm" : "off";
  draw();
});
bind("b_reset", "click", () => {
  state.xdom = null; state.hidden.clear(); draw();
});

/* ---------- canvas ---------- */
const cv = document.getElementById("cv"), ctx = cv.getContext("2d");
const tip = document.getElementById("tip");
const PAD = { l: 78, r: 18, t: 14, b: 48 };
let plot = null;                                  // last-drawn frame, for hover

function niceTicks(lo, hi, want) {
  if (!isFinite(lo) || !isFinite(hi) || !(hi > lo)) return [lo];
  const raw = (hi - lo) / want, mag = Math.pow(10, Math.floor(Math.log10(raw)));
  const step = [1, 2, 2.5, 5, 10].find(m => m * mag >= raw) * mag;
  if (!isFinite(step) || step <= 0) return [lo, hi];
  const out = [];
  for (let v = Math.ceil(lo / step) * step; v <= hi + step * 1e-9; v += step) out.push(v);
  return out;
}
function logTicks(lo, hi) {
  const out = [];
  // A non-positive bound makes Math.log10 return -Infinity, which turns the
  // loop below into an infinite push().  Callers should not do that, but a
  // degenerate domain must never hang the page.
  if (!isFinite(lo) || !isFinite(hi) || lo <= 0 || hi <= lo) return out;
  for (let e = Math.floor(Math.log10(lo)); e <= Math.ceil(Math.log10(hi)); e++) {
    for (const m of [1, 2, 5]) {
      const v = m * Math.pow(10, e);
      if (v >= lo * 0.999 && v <= hi * 1.001) out.push(v);
    }
  }
  return out;
}
const fmt = v => {
  const a = Math.abs(v);
  if (a === 0) return "0";
  if (a >= 1e5 || a < 1e-3) return v.toExponential(1).replace("e+", "e");
  if (a >= 100) return v.toFixed(0);
  if (a >= 1) return v.toFixed(2).replace(/\.?0+$/, "");
  return v.toPrecision(3).replace(/\.?0+$/, "");
};

function draw() {
  const vis = visible();
  const prep = vis.map(s => ({ s, y: transform(s) }));

  // x domain
  let x0 = Infinity, x1 = -Infinity;
  for (const { s } of prep) { x0 = Math.min(x0, s.xlo); x1 = Math.max(x1, s.xhi); }
  if (!isFinite(x0)) { x0 = 0; x1 = 1; }
  if (state.xdom) { x0 = state.xdom[0]; x1 = state.xdom[1]; }

  // y domain over the x window only, so zooming rescales usefully
  const logMode = state.scale === "log";
  let y0 = Infinity, y1 = -Infinity;
  for (const { s, y } of prep) {
    for (let i = 0; i < y.length; i++) {
      const v = y[i];
      if (v == null || !isFinite(v)) continue;
      if (s.x[i] < x0 || s.x[i] > x1) continue;
      if (logMode && v <= 0) continue;
      if (v < y0) y0 = v;
      if (v > y1) y1 = v;
    }
  }
  const empty = !isFinite(y0);
  if (empty) { y0 = 0; y1 = 1; }
  if (y0 === y1) { y1 = y0 + (Math.abs(y0) || 1) * 0.1; }
  if (logMode) {
    // Nothing positive left to show -- e.g. every trace filtered out, or a
    // blank subtracted from itself.  Pick an arbitrary decade rather than
    // handing log10() a zero.
    if (!(y1 > 0)) { y0 = 1; y1 = 10; }
    else { if (!(y0 > 0)) y0 = y1 / 1e3; y0 /= 1.6; y1 *= 1.6; }
  } else {
    const p = (y1 - y0) * 0.06;
    y0 = y0 >= 0 ? 0 : y0 - p;
    y1 += p;
  }

  // sizing (hi-DPI)
  const dpr = window.devicePixelRatio || 1;
  const W = cv.parentElement.clientWidth || 900;
  const H = Math.max(340, Math.min(620, Math.round(W * 0.52)));
  cv.width = W * dpr; cv.height = H * dpr; cv.style.height = H + "px";
  ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
  ctx.clearRect(0, 0, W, H);

  const iw = W - PAD.l - PAD.r, ih = H - PAD.t - PAD.b;
  const sx = v => PAD.l + (v - x0) / (x1 - x0) * iw;
  const ly0 = logMode ? Math.log10(y0) : y0, ly1 = logMode ? Math.log10(y1) : y1;
  const sy = v => PAD.t + ih - ((logMode ? Math.log10(v) : v) - ly0) / (ly1 - ly0) * ih;

  // grid + axes
  const xt = niceTicks(x0, x1, 9);
  const yt = logMode ? logTicks(y0, y1) : niceTicks(y0, y1, 7);
  ctx.lineWidth = 1; ctx.strokeStyle = "#eceff3";
  ctx.beginPath();
  for (const t of xt) { const p = Math.round(sx(t)) + .5; ctx.moveTo(p, PAD.t); ctx.lineTo(p, PAD.t + ih); }
  for (const t of yt) { const p = Math.round(sy(t)) + .5; ctx.moveTo(PAD.l, p); ctx.lineTo(PAD.l + iw, p); }
  ctx.stroke();

  ctx.strokeStyle = "#b9c2cc"; ctx.beginPath();
  ctx.moveTo(PAD.l + .5, PAD.t); ctx.lineTo(PAD.l + .5, PAD.t + ih + .5);
  ctx.lineTo(PAD.l + iw, PAD.t + ih + .5); ctx.stroke();

  ctx.fillStyle = "#63707f";
  ctx.font = "11px -apple-system, Segoe UI, Roboto, sans-serif";
  ctx.textAlign = "center"; ctx.textBaseline = "top";
  for (const t of xt) ctx.fillText(fmt(t), sx(t), PAD.t + ih + 7);
  ctx.textAlign = "right"; ctx.textBaseline = "middle";
  for (const t of yt) ctx.fillText(fmt(t), PAD.l - 8, sy(t));

  ctx.fillStyle = "#1c2530";
  ctx.font = "12px -apple-system, Segoe UI, Roboto, sans-serif";
  ctx.textAlign = "center"; ctx.textBaseline = "bottom";
  ctx.fillText("Emission wavelength (nm)", PAD.l + iw / 2, H - 6);
  ctx.save();
  ctx.translate(14, PAD.t + ih / 2); ctx.rotate(-Math.PI / 2);
  ctx.textBaseline = "top";
  ctx.fillText(yLabel(), 0, 0);
  ctx.restore();

  // series
  ctx.save();
  ctx.beginPath(); ctx.rect(PAD.l, PAD.t - 2, iw, ih + 4); ctx.clip();
  ctx.lineJoin = "round"; ctx.lineCap = "round";
  for (const { s, y } of prep) {
    ctx.strokeStyle = s.color; ctx.lineWidth = s.width; ctx.setLineDash(s.dash);
    ctx.beginPath();
    let pen = false;
    for (let i = 0; i < y.length; i++) {
      const v = y[i];
      // A log axis cannot show <=0, which blank subtraction routinely produces.
      // Break the polyline there instead of silently clamping.
      if (v == null || !isFinite(v) || (logMode && v <= 0)) { pen = false; continue; }
      const px = sx(s.x[i]), py = sy(v);
      if (!pen) { ctx.moveTo(px, py); pen = true; } else ctx.lineTo(px, py);
    }
    ctx.stroke();
  }
  ctx.setLineDash([]);

  if (state.peaks) {
    for (const { s, y } of prep) {
      let bi = -1, bv = -Infinity;
      for (let i = 0; i < y.length; i++) {
        if (y[i] == null || !isFinite(y[i])) continue;
        if (s.x[i] < x0 || s.x[i] > x1) continue;
        if (logMode && y[i] <= 0) continue;
        if (y[i] > bv) { bv = y[i]; bi = i; }
      }
      if (bi < 0) continue;
      const px = sx(s.x[bi]), py = sy(bv);
      ctx.fillStyle = s.color;
      ctx.beginPath(); ctx.arc(px, py, 3, 0, 6.2832); ctx.fill();
      ctx.font = "10.5px -apple-system, Segoe UI, Roboto, sans-serif";
      ctx.textAlign = "center"; ctx.textBaseline = "bottom";
      ctx.fillText(fmt(s.x[bi]), px, py - 5);
    }
  }
  ctx.restore();

  if (empty) {
    ctx.fillStyle = "#8d97a2";
    ctx.font = "13px -apple-system, Segoe UI, Roboto, sans-serif";
    ctx.textAlign = "center"; ctx.textBaseline = "middle";
    ctx.fillText(SERIES.length ? "No series match the current filters."
                               : "No data found.", PAD.l + iw / 2, PAD.t + ih / 2);
  }

  if (drag) {
    ctx.fillStyle = "rgba(37,99,235,.10)";
    ctx.fillRect(Math.min(drag.a, drag.b), PAD.t, Math.abs(drag.b - drag.a), ih);
  }

  plot = { prep, x0, x1, sx, sy, iw, ih, W, H, logMode };
  renderLegend(vis);
  renderStatus(vis);
}

function yLabel() {
  const u = SERIES.length ? SERIES[0].y_unit : "CPS";
  let base = state.norm === "peak" ? "Normalised intensity (peak = 1)"
           : state.norm === "area" ? `Normalised intensity (area = 1, ${u}·nm⁻¹)`
           : `Fluorescence intensity (${u})`;
  if (state.blank) base += ", blank-subtracted";
  return base;
}

function renderLegend(vis) {
  const box = document.getElementById("legend");
  box.textContent = "";
  const pool = SERIES.filter(s =>
    (state.sweep === "__ALL__" || s.sweep === state.sweep) &&
    facets.every(f => state.on[f].has(keyOf(s, f))));
  for (const s of pool) {
    const it = document.createElement("div");
    it.className = "it" + (state.hidden.has(s.id) ? " off" : "");
    it.title = s.path + (s.declared_em ? `  (declared Em ${s.declared_em})` : "");
    const sw = document.createElement("span");
    sw.className = "sw";
    sw.style.borderTop = `2.5px ${s.dash.length ? "dashed" : "solid"} ${s.color}`;
    it.append(sw, document.createTextNode(
      s.label + (state.sweep === "__ALL__" ? `  [${s.sweep}]` : "")));
    it.addEventListener("click", ev => {
      if (ev.shiftKey) {                          // shift-click = solo
        const others = pool.filter(o => o.id !== s.id);
        const soloed = others.every(o => state.hidden.has(o.id));
        state.hidden = new Set(soloed ? [] : others.map(o => o.id));
      } else {
        state.hidden.has(s.id) ? state.hidden.delete(s.id) : state.hidden.add(s.id);
      }
      draw();
    });
    box.append(it);
  }
  if (!pool.length) box.textContent = "No series match the current filters.";
}

function renderStatus(vis) {
  const bits = [`${vis.length} of ${SERIES.length} series shown`];
  if (state.blank) {
    const missing = vis.filter(s => !blankFor(s));
    if (missing.length) bits.push(`⚠ no matching blank for ${missing.length} series (left uncorrected)`);
  }
  if (state.xdom) bits.push(`zoomed ${fmt(state.xdom[0])}–${fmt(state.xdom[1])} nm · double-click to reset`);
  else bits.push("drag to zoom · shift-click a legend entry to solo it");
  document.getElementById("status").textContent = bits.join("  —  ");
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
  for (const { s, y } of plot.prep) {
    if (xv < s.xlo || xv > s.xhi) continue;
    let bi = 0, bd = Infinity;
    for (let i = 0; i < s.n; i++) { const d = Math.abs(s.x[i] - xv); if (d < bd) { bd = d; bi = i; } }
    if (y[bi] == null || !isFinite(y[bi])) continue;
    rows.push({ s, xv: s.x[bi], v: y[bi], dy: Math.abs(plot.sy(y[bi]) - my) });
  }
  if (!rows.length) { tip.style.display = "none"; return; }
  rows.sort((a, b) => b.v - a.v);
  const near = rows.reduce((a, b) => (b.dy < a.dy ? b : a));

  tip.innerHTML = "";
  const h = document.createElement("b");
  h.textContent = `${fmt(near.xv)} nm`;
  tip.append(h);
  for (const r of rows.slice(0, 10)) {
    const d = document.createElement("div");
    d.className = "r";
    d.style.fontWeight = r === near ? "600" : "400";
    const i = document.createElement("i"); i.style.background = r.s.color;
    const e = document.createElement("em"); e.textContent = r.s.label;
    const v = document.createElement("span"); v.textContent = fmt(r.v);
    d.append(i, e, v); tip.append(d);
  }
  if (rows.length > 10) {
    const d = document.createElement("div");
    d.className = "r"; d.style.color = "#63707f";
    d.textContent = `+${rows.length - 10} more`;
    tip.append(d);
  }
  tip.style.display = "block";
  const tw = tip.offsetWidth, th = tip.offsetHeight;
  tip.style.left = Math.min(mx + 14, plot.W - tw - 4) + "px";
  tip.style.top = Math.max(4, Math.min(my - th / 2, plot.H - th - 4)) + "px";
});
cv.addEventListener("mouseleave", () => { tip.style.display = "none"; });

/* ---------- drag-to-zoom ---------- */
let drag = null;
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
document.getElementById("b_png").addEventListener("click", () => {
  // Re-composite onto an opaque canvas: the live one is transparent, which
  // renders as black in most image viewers.
  const out = document.createElement("canvas");
  out.width = cv.width; out.height = cv.height;
  const o = out.getContext("2d");
  o.fillStyle = "#fff"; o.fillRect(0, 0, out.width, out.height);
  o.drawImage(cv, 0, 0);
  out.toBlob(b => download("spectra.png", b));
});
document.getElementById("b_csv").addEventListener("click", () => {
  const vis = visible();
  if (!vis.length) return;
  const prep = vis.map(s => ({ s, y: transform(s) }));
  const grid = [...new Set(prep.flatMap(p => p.s.x))].sort((a, b) => a - b);
  const head = ["Wavelength_nm", ...prep.map(p => p.s.label.replace(/[,\n]/g, " "))];
  const lines = [head.join(",")];
  for (const x of grid) {
    const row = [x];
    for (const p of prep) {
      const i = p.s.x.indexOf(x);
      row.push(i < 0 || p.y[i] == null || !isFinite(p.y[i]) ? "" : p.y[i]);
    }
    lines.push(row.join(","));
  }
  download("spectra_visible.csv", new Blob([lines.join("\n")], { type: "text/csv" }));
});

/* ---------- notes ---------- */
if (NOTES.length) {
  const n = document.getElementById("notes");
  n.hidden = false;
  n.innerHTML = "<b>Import notes</b><ul>" +
    NOTES.map(t => `<li>${t.replace(/[&<>]/g, c => ({ "&": "&amp;", "<": "&lt;", ">": "&gt;" }[c]))}</li>`).join("") +
    "</ul>";
}
document.getElementById("subtitle").textContent =
  `${SERIES.length} series · ${sweeps.length} sweep${sweeps.length === 1 ? "" : "s"} · ${PAYLOAD.generated_from}`;

let rt;
window.addEventListener("resize", () => { clearTimeout(rt); rt = setTimeout(draw, 90); });
draw();
</script>
"""


def render(series, notes, title, source):
    payload = {
        "series": to_payload(series),
        "notes": notes,
        "generated_from": source,
    }
    # ensure_ascii=True (the default) already escapes U+2028/U+2029, the two
    # characters JSON allows raw but JS treats as line terminators.  All that
    # is left is '<', so a label containing "</script>" cannot break out.
    blob = json.dumps(payload, allow_nan=False, separators=(",", ":"))
    blob = blob.replace("<", "\\u003c")
    return (HTML
            .replace("__PAYLOAD__", blob)
            .replace("__TITLE__", title))


# --------------------------------------------------------------------------
# cli
# --------------------------------------------------------------------------

def main(argv=None):
    here = Path(__file__).resolve().parent
    ap = argparse.ArgumentParser(
        description=__doc__.strip().splitlines()[0],
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--root", default=str(here), type=Path,
                    help="directory to scan recursively (default: this script's folder)")
    ap.add_argument("--out", default=None, type=Path,
                    help="output HTML path (default: <root>/spectra_viewer.html)")
    ap.add_argument("--title", default="Fluorometer emission spectra")
    ap.add_argument("--open", action="store_true", help="open the result in a browser")
    ap.add_argument("--list", action="store_true", dest="just_list",
                    help="print the discovered series and exit")
    ap.add_argument("--keep-duplicates", action="store_true",
                    help="keep files whose numeric content is identical to another")
    args = ap.parse_args(argv)

    root = args.root.resolve()
    if not root.is_dir():
        ap.error(f"--root is not a directory: {root}")

    series, notes = discover(root, keep_duplicates=args.keep_duplicates)
    if not series:
        print(f"No readable CSV files under {root}", file=sys.stderr)
        for n in notes:
            print("  " + n, file=sys.stderr)
        return 1

    if args.just_list:
        w = max(len(s["path"]) for s in series)
        print(f"{len(series)} series under {root}\n")
        print(f"{'file':<{w}}  {'sweep':<24} {'construct':<12} {'sample':<12} rep")
        print("-" * (w + 56))
        for s in series:
            print(f"{s['path']:<{w}}  {s['sweep']:<24} {s['construct']:<12} "
                  f"{s['sample']:<12} {s['replicate'] if s['replicate'] is not None else '-'}")
        if notes:
            print("\nnotes:")
            for n in notes:
                print("  " + n)
        return 0

    out = (args.out or root / "spectra_viewer.html").resolve()
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(render(series, notes, args.title, root.name), encoding="utf-8")

    sweeps = sorted({s["sweep"] for s in series})
    print(f"wrote {out}  ({len(series)} series, {len(sweeps)} sweeps, "
          f"{out.stat().st_size / 1024:.0f} kB)")
    for sw in sweeps:
        print(f"    {sw} nm  — {sum(1 for s in series if s['sweep'] == sw)} series")
    for n in notes:
        print("  note: " + n)
    if args.open:
        webbrowser.open(out.as_uri())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
