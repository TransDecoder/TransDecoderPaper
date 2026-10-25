#!/bin/bash
P=/broad/hptmp/bhaas/TRANSDECODER_PAPER
r=$P/work/bench_v3/frames/runs/$1/$2/all
$P/env/td/bin/python $P/TransDecoderPaper/benchmark_v2/scripts/score_predictions.py --truth $P/data/bench_v3/inputs/$1/$2.truth.tsv --preds $r/preds.tsv | tail -1 | cut -f2- | sed "s/^/$1\t$2\tall\t/" | sed "s/\$/\t$(cat $r/runtime_sec.txt)/"
