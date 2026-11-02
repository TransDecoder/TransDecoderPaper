#!/usr/bin/env bash
# Score the held-out preset test at the pinned release code (E40; design E21). Output data/benchmark/presets_release/.
set -euo pipefail
P=/broad/hptmp/bhaas/TRANSDECODER_PAPER; T=$P/TransDecoderPaper; PY=$P/env/td/bin/python; W=$P/work/release_pin/presets
OUT=$T/data/benchmark/presets_release; mkdir -p $OUT/preset_meta
TESTS=$P/work/bench_v3/presets/tests
{ first=1; for o in human mouse athal dmel; do
    other=human; [ $o = human ] && other=mouse
    vars=$(ls $TESTS/$o/*.fa | xargs -n1 basename | sed 's/\.fa$//' | paste -sd,)
    $PY $T/benchmark_v2/scripts/score_negatives.py $TESTS "$vars" \
      self_decoy=$W/runs/{org}/{variant}/self/preds.tsv \
      preset_same=$W/runs/{org}/{variant}/preset_$o/preds.tsv \
      preset_samemix=$W/runs/{org}/{variant}/preset_${o}mix/preds.tsv \
      preset_other=$W/runs/{org}/{variant}/preset_$other/preds.tsv \
      | { if [ $first = 1 ]; then cat; else tail -n +2; fi; } | awk -F'\t' -v o=$o 'NR==1||$1==o||$1=="org"'
    first=0
  done; } > $OUT/preset_scores.tsv
for n in humanmix; do mkdir -p $OUT/preset_meta/$n; cp $W/presets/$n/preset.json $OUT/preset_meta/$n/; done
wc -l $OUT/preset_scores.tsv
