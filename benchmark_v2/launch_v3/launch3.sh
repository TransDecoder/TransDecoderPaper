#!/bin/bash
# Third runner: the additional comparators (CPAT, CPC2, Borf, RNAsamba).
cd /broad/hptmp/bhaas/TRANSDECODER_PAPER
R=TransDecoderPaper/benchmark_v2/scripts/run_all.sh
JOBS=6 $R data/bench_v3/inputs work/bench_v3/runs cpat cpc2 borf
VARIANTS="complete complete.wRand complete_fullutr fragments fragments.wRand" \
  JOBS=6 $R data/bench_v3/inputs work/bench_v3/runs rnasamba
echo ALL_DONE3
