#!/bin/bash
# Run every (organism, input variant, tool) benchmark run locally, JOBS at a time.
# usage: run_all.sh <inputs_dir> <runs_dir> [tool ...]
# Skips runs whose preds.tsv already exists. VARIANTS (space-separated input basenames
# without .fa, default all) restricts which inputs are run. JOBS (default 16) sets how many
# runs execute concurrently; pytd_homology runs use THREADS (default 4) cores each.
set -euo pipefail
here=$(cd "$(dirname "$0")" && pwd)
inputs=$(readlink -f "$1"); runs=$(mkdir -p "$2" && readlink -f "$2"); shift 2
tools=("$@"); [ ${#tools[@]} -gt 0 ] || tools=(pytd pytd_norefine perl_td perl_td_v600 perl_td_v571 prodigal codan)
mkdir -p "$runs/logs"
for orgdir in "$inputs"/*/; do
  org=$(basename "$orgdir")
  for fa in "$orgdir"/*.fa; do
    v=$(basename "$fa" .fa)
    if [ -n "${VARIANTS:-}" ] && [[ " $VARIANTS " != *" $v "* ]]; then continue; fi
    for t in "${tools[@]}"; do
      [ -s "$runs/$org/$v/$t/preds.tsv" ] && continue
      printf '%s\t%s\t%s\t%s\n' "$t" "$fa" "$runs/$org/$v/$t" "$org"
    done
  done
done | xargs -P "${JOBS:-16}" -d '\n' -n 1 bash -c '
  IFS=$'"'"'\t'"'"' read -r t fa od org <<< "$1"
  log=$0/logs/$org.$(basename "$fa" .fa).$t.log
  if "'"$here"'/run_tool.sh" "$t" "$fa" "$od" "$org" > "$log" 2>&1; then echo "done $org $t $fa"; else echo "FAILED $org $t $fa (see $log)"; fi
' "$runs"
