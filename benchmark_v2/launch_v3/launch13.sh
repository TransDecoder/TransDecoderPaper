#!/bin/bash
# Final defaults (additive Kozak start rule, -m 50 + short-ORF test): PyTD and Perl, all inputs
cd /broad/hptmp/bhaas/TRANSDECODER_PAPER
JOBS=10 TransDecoderPaper/benchmark_v2/scripts/run_all.sh data/bench_v3/inputs work/bench_v3/runs pytd_dev_v2 perl_td_dev_v2
echo ALL_DONE13
