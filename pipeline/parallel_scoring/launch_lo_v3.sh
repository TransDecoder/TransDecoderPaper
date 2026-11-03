#!/bin/bash
# LongOrfs --cpu series: (A) one job at a time: --cpu 8 on all 108 benchmark inputs and GM12878
# LRAA/StringTie/Trinity, GM12878 Trinity --cpu 1, pyTransdecoder --cpu 8 on four inputs;
# (B) --cpu 1 on the 108 benchmark inputs and GM12878 LRAA/StringTie, three at a time (parity only).
P=/broad/hptmp/bhaas/TRANSDECODER_PAPER; R=$P/TransDecoderPaper/pipeline/parallel_scoring/run_lo_v3.sh
W=$P/work/longorfs_cpu; L=$W/launch_lo_v3.log; mkdir -p $W/logs
A=(); B=()
for d in $P/work/release_pin/bench/*/*/rel; do
  i=$(basename $(dirname $d)); o=$(basename $(dirname $(dirname $d))); A+=("bench $o $i 8"); B+=("bench $o $i 1")
done
for a in lraa stringtie trinity; do A+=("gm $a 8"); done
A+=("gm trinity 1"); B+=("gm lraa 1" "gm stringtie 1")
for x in "human complete" "tthermophila complete.unstranded" "pfalciparum fragments.wRand" "spombe fragments"; do A+=("pipe $x 8"); done
run() { if bash $R $1 > $W/logs/$(echo $1 | tr ' ' '.').log 2>&1; then s=ok; else s=FAIL; fi
        echo "$s $1 $(date +%T) load: $(cut -d' ' -f1-3 /proc/loadavg)" >> $L; }
export -f run; export R W L
echo "=== start $(date) A=${#A[@]} B=${#B[@]}" >> $L
for j in "${A[@]}"; do run "$j"; done
printf '%s\n' "${B[@]}" | xargs -P 3 -I{} bash -c 'run "{}"'
echo "=== end $(date)" >> $L
