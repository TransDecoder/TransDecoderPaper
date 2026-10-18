#!/bin/bash
# Build the trimmed-UTR (250 nt) and full-UTR truth sets for every organism in
# truth_sources.tsv, locally, JOBS (default 8) organisms at a time.
# usage: build_all_truth.sh <annot_dir> <out_dir>
set -euo pipefail
here=$(cd "$(dirname "$0")" && pwd)
annot=$(readlink -f "$1"); out=$(mkdir -p "$2" && readlink -f "$2")
PYTHON=${PYTHON:-python3}
mkdir -p "$out/truth" "$out/truth_fullutr"
tail -n +2 "$here/../truth_sources.tsv" | xargs -P "${JOBS:-8}" -d '\n' -n 1 bash -c '
  IFS=$'"'"'\t'"'"' read -r org g a extra <<< "$1"
  for m in 250:truth -1:truth_fullutr; do
    "'"$PYTHON"'" "'"$here"'/build_truth_from_annotation.py" --genome "'"$annot"'/$g" --annotation "'"$annot"'/$a" \
      --add-utr --max-utr ${m%%:*} --out "'"$out"'/${m#*:}/$org" $extra 2> /dev/null
  done && echo "done $org" || echo "FAILED $org"
' _
