#!/bin/bash
# Submit every (organism, input variant, tool) benchmark run to UGER.
# usage: submit_all.sh <inputs_dir> <runs_dir> [tool ...]
# Skips runs whose preds.tsv already exists. VARIANTS (space-separated input basenames
# without .fa, default all) restricts which inputs are submitted.
set -euo pipefail
here=$(cd "$(dirname "$0")" && pwd)
inputs=$(readlink -f "$1"); runs=$(mkdir -p "$2" && readlink -f "$2"); shift 2
tools=("$@"); [ ${#tools[@]} -gt 0 ] || tools=(pytd pytd_norefine perl_td perl_td_v600 perl_td_v571 prodigal codan)
mkdir -p "$runs/logs"
n=0
for orgdir in "$inputs"/*/; do
  org=$(basename "$orgdir")
  for fa in "$orgdir"/*.fa; do
    v=$(basename "$fa" .fa)
    if [ -n "${VARIANTS:-}" ] && [[ " $VARIANTS " != *" $v "* ]]; then continue; fi
    for t in "${tools[@]}"; do
      od=$runs/$org/$v/$t
      [ -s "$od/preds.tsv" ] && continue
      res=(-l h_vmem=12G); [ "$t" = pytd_homology ] && res=(-pe smp 4 -binding linear:4 -l h_vmem=4G)
      qsub -cwd -V -N "b_${org}_${t}" -l h_rt=12:00:00 "${res[@]}" -o "$runs/logs/${org}.${v}.${t}.log" -j y -b y \
        "$here/run_tool.sh" "$t" "$fa" "$od" "$org" > /dev/null
      n=$((n+1))
    done
  done
done
echo "submitted $n jobs"
