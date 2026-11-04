#!/bin/bash
# E30: run the ORF predictors on the simulated GM12878 Trinity assembly (stranded; human settings).
# 4 jobs at a time (shared host). run_tool.sh is called, not modified.
P=/broad/hptmp/bhaas/TRANSDECODER_PAPER
R=$P/work/simtrinity/runs
IN=$P/data/simtrinity/inputs/fragments_simtrinity.fa
export P R IN
printf "%s\n" pytd_v7 perl_td_v600 gmst prodigal cpat cpc2 rnasamba borf | \
  xargs -P 4 -I{} bash -c 'echo "[$(date)] start {}"; THREADS=1 bash $P/TransDecoderPaper/benchmark_v2/scripts/run_tool.sh {} $IN $R/{} human > $R/logs/{}.log 2>&1; echo "[$(date)] end {} exit $?"'
echo SIMTRINITY_RUNS_DONE
