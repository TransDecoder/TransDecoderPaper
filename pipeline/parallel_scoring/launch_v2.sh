#!/bin/bash
# Parity + timing of dev/faster-serial (PY_fast), one job at a time:
#  - timing at --cpu 1 on the ten v1 inputs and GM12878 Trinity (for comparison with v1 / release);
#  - parity at --cpu 8 on all 108 benchmark inputs (release defaults), 16 homology runs
#    (retention-only and start assist), GM12878 LRAA/StringTie/Trinity, one preset build, two preset runs.
P=/broad/hptmp/bhaas/TRANSDECODER_PAPER; R=$P/TransDecoderPaper/pipeline/parallel_scoring/run_one_v2.sh
W=$P/work/fast_scoring; L=$W/launch_v2.log; mkdir -p $W/logs
jobs=()
for x in "human complete" "human complete.wLnc" "mouse fragments.unstranded" "dmel complete_fullutr" \
         "athal fragments_mixed_t1" "spombe fragments" "tthermophila complete.unstranded.wRand" \
         "pfalciparum fragments.wRand" "creinhardtii complete"; do jobs+=("bench $x rel 1"); done
jobs+=("gm trinity 1" "gm trinity 8" "gm lraa 8" "gm stringtie 8")
for d in $P/work/release_pin/bench/*/*/rel; do
  i=$(basename $(dirname $d)); o=$(basename $(dirname $(dirname $d))); jobs+=("bench $o $i rel 8")
done
for o in human mouse dmel athal spombe tthermophila pfalciparum creinhardtii; do
  for i in fragments complete; do jobs+=("bench $o $i retonly 8" "bench $o $i hassist 8"); done
done
jobs+=("presetbuild human A" "presettest human complete.wRand.B200.r1 human 8" "presettest human fragments.codingonly.C50.r1 humanmix 8")
echo "=== start $(date) ${#jobs[@]} jobs load: $(cat /proc/loadavg)" >> $L
for j in "${jobs[@]}"; do
  if bash $R $j > $W/logs/$(echo $j | tr ' ' '.').log 2>&1; then s=ok; else s=FAIL; fi
  echo "$s $j $(date +%T) load: $(cut -d' ' -f1-3 /proc/loadavg)" >> $L
done
echo "=== end $(date)" >> $L
