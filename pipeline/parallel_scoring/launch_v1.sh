#!/bin/bash
# Parity + timing series for the parallel-scoring branch, run one job at a time (timing on a shared host).
# Benchmark: one input per organism at --cpu 1 and 8; human complete also at 2 and 4.
# GM12878 Trinity (largest real input) at --cpu 1, 2, 4, 8.
P=/broad/hptmp/bhaas/TRANSDECODER_PAPER; R=$P/TransDecoderPaper/pipeline/parallel_scoring/run_one_v1.sh
L=$P/work/parallel_scoring/launch_v1.log; mkdir -p $P/work/parallel_scoring/logs
jobs=(
"human complete 1" "human complete 8" "human complete 2" "human complete 4"
"human complete.wLnc 1" "human complete.wLnc 8"
"mouse fragments.unstranded 1" "mouse fragments.unstranded 8"
"dmel complete_fullutr 1" "dmel complete_fullutr 8"
"athal fragments_mixed_t1 1" "athal fragments_mixed_t1 8"
"spombe fragments 1" "spombe fragments 8"
"tthermophila complete.unstranded.wRand 1" "tthermophila complete.unstranded.wRand 8"
"pfalciparum fragments.wRand 1" "pfalciparum fragments.wRand 8"
"creinhardtii complete 1" "creinhardtii complete 8"
"gm trinity 8" "gm trinity 4" "gm trinity 2" "gm trinity 1"
)
echo "=== start $(date) load: $(cat /proc/loadavg)" >> $L
for j in "${jobs[@]}"; do
  set -- $j
  if bash $R $1 $2 $3 > $P/work/parallel_scoring/logs/$1.$2.cpu$3.log 2>&1; then s=ok; else s=FAIL; fi
  echo "$s $j $(date +%T) load: $(cut -d' ' -f1-3 /proc/loadavg)" >> $L
done
echo "=== end $(date)" >> $L
