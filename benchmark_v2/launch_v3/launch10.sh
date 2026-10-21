#!/bin/bash
# Minimum ORF length sweep (PyTD dev defaults, LongOrfs -m 90..50) on inputs with random negatives
cd /broad/hptmp/bhaas/TRANSDECODER_PAPER
VARIANTS="complete.wRand fragments.wRand" JOBS=8 TransDecoderPaper/benchmark_v2/scripts/run_all.sh \
  data/bench_v3/inputs work/bench_v3/runs pytd_dev_m90 pytd_dev_m80 pytd_dev_m70 pytd_dev_m60 pytd_dev_m50
echo ALL_DONE10
