#!/usr/bin/env bash
# Build ../Ian_Poon_YiI_Quarterly_Report_1.pdf.
# The word count printed in the title block uses the STRICTEST convention:
# every word of the rendered PDF from "Context" to "References" (text, headings,
# tables, table captions and in-text citations), excluding Figure 1 and its
# caption, running headers and page numbers (count_rendered.py). texcount's
# figure, which skips citations, is printed for comparison. Limit: 500-750.
set -euo pipefail
cd "$(dirname "$0")"
[ -f wordcount.tex ] || printf '\\newcommand{\\wordcount}{0}\n' > wordcount.tex
run() { latexmk -pdf -interaction=nonstopmode -halt-on-error -outdir=build report.tex > build.log 2>&1 || { tail -40 build.log; exit 1; }; }
run
n=$(python3 count_rendered.py build/report.pdf)
if ! grep -q "{$n}" wordcount.tex; then printf '\\newcommand{\\wordcount}{%s}\n' "$n" > wordcount.tex; run; fi
tc=$(texcount -1 -sum -q report.tex 2>/dev/null | grep -o '^[0-9]*')
cp build/report.pdf ../Ian_Poon_YiI_Quarterly_Report_1.pdf
echo "words: $n rendered (texcount without citations: $tc)  ->  ../Ian_Poon_YiI_Quarterly_Report_1.pdf ($(pdfinfo build/report.pdf | awk '/Pages/{print $2}') pages)"
[ "$n" -le 750 ] || { echo "OVER THE 750-WORD LIMIT"; exit 1; }
