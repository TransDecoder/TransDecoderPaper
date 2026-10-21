#!/bin/bash
# Final defaults after unifying the GC estimate and the exact-tie order: PyTD and Perl, all inputs
cd /broad/hptmp/bhaas/TRANSDECODER_PAPER
JOBS=6 TransDecoderPaper/benchmark_v2/scripts/run_all.sh data/bench_v3/inputs work/bench_v3/runs pytd_dev_v2 perl_td_dev_v2
echo ALL_DONE16
