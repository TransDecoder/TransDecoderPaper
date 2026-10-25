#!/bin/bash
# Knit the paper notebooks to HTML (notebooks/*.html) and write their figures (figures/*.pdf,png).
# Environment: environment/rnb.yml  (micromamba create -p $P/env/rnb -f environment/rnb.yml)
# usage: scripts/knit_notebooks.sh [notebook.Rmd ...]     (default: all notebooks/*.Rmd)
set -euo pipefail
T=$(cd "$(dirname "$0")/.." && pwd)
RSCRIPT=${RSCRIPT:-/broad/hptmp/bhaas/TRANSDECODER_PAPER/env/rnb/bin/Rscript}
export RSTUDIO_PANDOC=$(dirname "$RSCRIPT")
cd "$T/notebooks"
nbs=("$@"); [ ${#nbs[@]} -gt 0 ] || nbs=([0-9S]*.Rmd)
for nb in "${nbs[@]}"; do
  nb=$(basename "$nb")
  echo "knitting $nb"
  "$RSCRIPT" -e "rmarkdown::render('$nb', quiet = TRUE, envir = new.env())"
done
