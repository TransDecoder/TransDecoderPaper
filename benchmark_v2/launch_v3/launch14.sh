#!/bin/bash
# Comparators on the GENCODE lncRNA-negative inputs (human, mouse)
cd /broad/hptmp/bhaas/TRANSDECODER_PAPER
VARIANTS="complete.wLnc fragments.wLnc" JOBS=4 TransDecoderPaper/benchmark_v2/scripts/run_all.sh \
  data/bench_v3/inputs work/bench_v3/runs cpat cpc2 borf rnasamba gmst prodigal
echo ALL_DONE14
