#!/usr/bin/env python3
"""Build ../REPORT.pdf from ../REPORT.md with Ian's report LaTeX template.

    python3 build.py

Steps: trim the markdown front matter the cover page replaces, record image
pixel sizes for the filter (which also sizes figures and table columns), run
pandoc (template.tex + filter.lua) to report.tex, then latexmk.
REPORT.md stays the single source: edit it, then rebuild.
"""
import json, pathlib, re, shutil, subprocess
from PIL import Image

HERE = pathlib.Path(__file__).resolve().parent
ROOT = HERE.parent
BUILD = HERE / 'build'
BUILD.mkdir(exist_ok=True)

md = (ROOT / 'REPORT.md').read_text()

# 1. the cover page carries title, author and dates; the PDF has its own contents page
md = md[md.index('**How this report was made.**'):]
md = re.sub(r'\n## Contents\n.*?\n---\n', '\n', md, count=1, flags=re.S)

(BUILD / '_report.md').write_text(md)

# 2. pixel sizes of every embedded image (the filter sizes figures from them)
sizes = {}
for src in re.findall(r'!\[[^\]]*\]\(([^)]+)\)', md):
    p = ROOT / src
    if p.suffix.lower() in ('.png', '.jpg', '.jpeg'):
        with Image.open(p) as im: sizes[src] = list(im.size)
(BUILD / 'imgsizes.json').write_text(json.dumps(sizes))

# 3. pandoc -> report.tex -> latexmk
subprocess.run(['pandoc', str(BUILD / '_report.md'), '-f', 'gfm+implicit_figures-emoji', '-t', 'latex',
                '--template', 'template.tex', '--lua-filter', 'filter.lua',
                '--shift-heading-level-by=-1', '-o', 'report.tex'], cwd=HERE, check=True)
subprocess.run(['latexmk', '-pdf', '-interaction=nonstopmode', '-halt-on-error', '-outdir=build', 'report.tex'],
               cwd=HERE, check=True, stdout=subprocess.DEVNULL)
shutil.copy(BUILD / 'report.pdf', ROOT / 'REPORT.pdf')
print('wrote', ROOT / 'REPORT.pdf')
