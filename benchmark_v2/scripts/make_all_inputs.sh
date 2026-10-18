#!/bin/bash
# Derive every benchmark input set from the truth sets built by build_all_truth.sh.
# usage: make_all_inputs.sh <bench_dir>   (reads <bench_dir>/truth{,_fullutr}, writes <bench_dir>/inputs/<org>)
# complete_fullutr[.unstranded] come from the untrimmed-UTR truth; all others from the 250-nt one.
set -euo pipefail
here=$(cd "$(dirname "$0")" && pwd)
b=$(readlink -f "$1"); PYTHON=${PYTHON:-python3}
for t in "$b"/truth/*.truth.tsv; do
  org=$(basename "$t" .truth.tsv); od=$b/inputs/$org
  "$PYTHON" "$here/make_benchmark_inputs.py" --fasta "$b/truth/$org.fasta" --truth "$t" --outdir "$od" --seed 1
  tmp=$(mktemp -d "$b/inputs/.fullutr_$org.XXXX")
  "$PYTHON" "$here/make_benchmark_inputs.py" --fasta "$b/truth_fullutr/$org.fasta" \
    --truth "$b/truth_fullutr/$org.truth.tsv" --outdir "$tmp" --seed 1
  for v in complete complete.unstranded; do
    for ext in fa truth.tsv; do mv "$tmp/$v.$ext" "$od/${v/complete/complete_fullutr}.$ext"; done
  done
  rm -rf "$tmp"
  echo "done $org"
done
