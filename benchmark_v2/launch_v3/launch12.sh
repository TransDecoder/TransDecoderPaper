#!/bin/bash
# Short-ORF frame-contrast guard: LongOrfs -m 50 with Predict --short-orf-min-contrast 150/200/250
cd /broad/hptmp/bhaas/TRANSDECODER_PAPER
VARIANTS="complete.wRand fragments.wRand complete.wLnc fragments.wLnc" JOBS=8 \
  TransDecoderPaper/benchmark_v2/scripts/run_all.sh data/bench_v3/inputs work/bench_v3/runs \
  pytd_dev_m50_sc200 pytd_dev_m50_sc150 pytd_dev_m50_sc250
echo ALL_DONE12
