#!/usr/bin/env python3
"""Count words a marker would count if they pasted the PDF into Word:
everything from 'Context' up to 'References', minus Figure 1 and its caption,
running headers and page numbers. Citations, headings, tables and table
captions ARE counted (the strictest convention)."""
import re, subprocess, sys
pdf = sys.argv[1] if len(sys.argv) > 1 else 'build/report.pdf'
t = subprocess.run(['pdftotext', '-layout', pdf, '-'], capture_output=True, text=True).stdout
body = t[t.index('Context'):t.rindex('References')]
a = body.find('1 Starting point'); b = body.find('A1 and A3 address.')
if a >= 0 and b > a:
    body = body[:a] + body[b + len('A1 and A3 address.'):]
lines = [l for l in body.splitlines()
         if not re.match(r'\s*(Year in Industry · Quarterly Report 1\s+Ian Poon|\d+)\s*$', l)]
words = [w for w in ' '.join(lines).split() if re.search(r'[A-Za-z0-9]', w)]
print(len(words))
