#!/bin/bash
cd /broad/hptmp/bhaas/TRANSDECODER_PAPER/work/bench_v3/history
xargs -P 6 -L 1 bash -c '/broad/hptmp/bhaas/TRANSDECODER_PAPER/TransDecoderPaper/benchmark_v2/scripts/hist_run.sh "$0" "$1" "$2" "$3" > /dev/null 2>&1 && echo "done $2" || echo "FAILED $2"' < jobs.txt
echo ALL_DONE
