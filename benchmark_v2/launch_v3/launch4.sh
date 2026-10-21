#!/bin/bash
# Fourth runner: ANGEL (self-training; THREADS cores per run).
cd /broad/hptmp/bhaas/TRANSDECODER_PAPER
THREADS=4 JOBS=3 TransDecoderPaper/benchmark_v2/scripts/run_all.sh data/bench_v3/inputs work/bench_v3/runs angel
echo ALL_DONE4
