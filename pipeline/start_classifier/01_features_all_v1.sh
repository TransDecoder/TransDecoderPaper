#!/bin/bash
# E33 step 1: candidate features for every PyTD release run (3a793ff) used in E33. Serial, local.
set -euo pipefail
P=/broad/hptmp/bhaas/TRANSDECODER_PAPER; T=$P/TransDecoderPaper
PY=$P/env/startclf/bin/python; export PYTHONPATH=$P/work/trees/PY_v7
S=$T/pipeline/start_classifier; O=$T/data/start_classifier/features; mkdir -p $O
log=$O/extract.log; : > $log
one() { # run_dir gcode truth dataset
  [ "$(cat $1/pytd_version.txt)" = 3a793ff ] || { echo "bad version $1" >&2; exit 1; }
  $PY $S/01_features.py "$1" "$2" "$3" "$4" $O/$4.tsv.gz 2>> $log
}
while IFS=$'\t' read -r org _ _ gcode _; do
  [ "$org" = org ] && continue
  for v in complete complete_fullutr fragments fragments_mixed_t1 fragments_mixed_t2 fragments_mixed_t3; do
    one $P/work/bench_v3/runs/$org/$v/pytd_v7 $gcode $P/data/bench_v3/inputs/$org/$v.truth.tsv bench.$org.$v
  done
done < $T/benchmark_v2/organisms.tsv
one $P/work/simtrinity/runs/pytd_v7 Universal $P/data/simtrinity/inputs/fragments_simtrinity.truth.tsv sim.human.d1
for d in d2 d4 d8 d16 d32; do
  one $P/work/simtrinity/depth/$d/runs/pytd_v7 Universal $P/data/simtrinity/depth/$d/inputs/fragments_simtrinity.truth.tsv sim.human.$d
done
one $P/work/simtrinity_athal/depth/d1/runs/pytd_v7 Universal $P/data/simtrinity_athal/depth/d1/inputs/fragments_simtrinity.truth.tsv sim.athal.d1
echo DONE >> $log
