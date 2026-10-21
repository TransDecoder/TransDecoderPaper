#!/bin/bash
# Second runner for everything except CodAn, so those runs don't wait behind CodAn in launch.sh.
cd /broad/hptmp/bhaas/TRANSDECODER_PAPER
while ps -eo args | grep -v grep | grep -q 'run_tool.sh perl_td_v571 .*creinhardtii//fragments.fa'; do sleep 20; done
R=TransDecoderPaper/benchmark_v2/scripts/run_all.sh
JOBS=10 $R data/bench_v3/inputs work/bench_v3/runs pytd pytd_norefine perl_td perl_td_v600 perl_td_v571 prodigal gmst
VARIANTS="complete complete.unstranded complete.wRand fragments fragments.unstranded complete_fullutr complete_fullutr.unstranded" \
  JOBS=4 $R data/bench_v3/inputs work/bench_v3/runs pytd_homology
echo ALL_DONE2
