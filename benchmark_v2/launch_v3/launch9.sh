#!/bin/bash
# PyTD dev: additive Kozak term in the start keep-or-trim rule
cd /broad/hptmp/bhaas/TRANSDECODER_PAPER
JOBS=8 TransDecoderPaper/benchmark_v2/scripts/run_all.sh data/bench_v3/inputs work/bench_v3/runs pytd_dev_kozak_add
echo ALL_DONE9
