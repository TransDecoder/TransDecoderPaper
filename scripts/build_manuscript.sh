#!/bin/bash
# Build the manuscript PDF locally with Tectonic (a self-contained LaTeX engine; runs BibTeX itself).
# Tectonic env: micromamba create -p $P/env/latex -c conda-forge tectonic poppler
# usage: scripts/build_manuscript.sh [outdir]      (default: $P/work/manuscript_build)
# Prints LaTeX/BibTeX problems (undefined references or citations, BibTeX warnings) and the page count.
set -euo pipefail
P=/broad/hptmp/bhaas/TRANSDECODER_PAPER
T=$(cd "$(dirname "$0")/.." && pwd)
OUT=${1:-$P/work/manuscript_build}
mkdir -p "$OUT"
cd "$T/manuscript"
"$P/env/latex/bin/tectonic" -X compile --keep-logs --outdir "$OUT" main.tex > "$OUT/build.out" 2>&1 \
  || { tail -30 "$OUT/build.out"; exit 1; }
grep -i -E "undefined|multiply defined" "$OUT/main.log" || true
grep -i "warning" "$OUT/main.blg" || true
grep -c "Overfull" "$OUT/build.out" | sed 's/^/overfull hbox warnings: /'
"$P/env/latex/bin/pdfinfo" "$OUT/main.pdf" | grep Pages
echo "PDF: $OUT/main.pdf"
