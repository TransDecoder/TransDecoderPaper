#!/bin/bash
cd /broad/hptmp/bhaas/TRANSDECODER_PAPER
R=TransDecoderPaper/benchmark_v2/scripts/run_all.sh
JOBS=16 $R data/bench_v3/inputs work/bench_v3/runs pytd pytd_norefine perl_td perl_td_v600 perl_td_v571 prodigal codan gmst
VARIANTS="complete complete.unstranded complete.wRand fragments fragments.unstranded complete_fullutr complete_fullutr.unstranded" \
  JOBS=6 $R data/bench_v3/inputs work/bench_v3/runs pytd_homology
echo ALL_DONE
