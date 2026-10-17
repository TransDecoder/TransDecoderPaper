#!/bin/bash
# Run the start-refinement strategy grid on one benchmark dataset (one UGER job).
# usage: run_strategy_grid.sh <org> <variant> <out.tsv>
# Needs the benchmark's pytd_norefine run for <org>/<variant>.
set -euo pipefail
here=$(cd "$(dirname "$0")" && pwd)
ROOT=${BENCH_ROOT:-/broad/hptmp/bhaas/TRANSDECODER_PAPER}
PY=${PYTHON:-$ROOT/env/td/bin/python}
org=$1; v=$2; out=$3
run=$ROOT/work/bench_v2/runs/$org/$v/pytd_norefine
truth=$ROOT/data/bench_v2/inputs/$org/$v.truth.tsv
gcode=$(awk -F'\t' -v o="$org" '$1==o{print $4}' "$here/../../benchmark_v2/organisms.tsv")
E="$PY $here/start_strategy_experiment.py --run-dir $run --truth $truth --pytd $ROOT/work/trees/PY_dev --genetic-code $gcode"
{
  $E --strategy none --label "$org/$v"
  for s in 1 2 3 4 5; do $E --strategy legacy --seed $s --label "$org/$v" | tail -1; done
  for s in 1 2 3 4 5; do for o in f1_balanced f1_raw youden; do for r in best first; do
    $E --strategy kfold --seed $s --thr-objective $o --rule $r --label "$org/$v" | tail -1; done; done; done
  for s in 1 2 3; do for r in best first; do
    $E --strategy fixed --seed $s --rule $r --label "$org/$v" | tail -1; done; done
} > "$out.tmp" && mv "$out.tmp" "$out"
