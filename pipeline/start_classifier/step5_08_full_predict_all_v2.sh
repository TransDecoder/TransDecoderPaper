#!/bin/bash
# E33 step 5 (v2, post hoc): all hybrid (Swiss-Prot + UniRef50 fallback) full-mode Predict runs, 4 at a time (single-threaded Python each), logs in step5_full/logs/.
set -uo pipefail
P=/broad/hptmp/bhaas/TRANSDECODER_PAPER; D=$P/TransDecoderPaper/pipeline/start_classifier
L=$P/work/start_classifier/step5_full/logs; mkdir -p $L
ds=(sim.human.d8 sim.athal.d1 sim.athal.d16)
for o in human athal dmel tthermophila; do
  for v in complete fragments complete.wRand fragments.wRand; do ds+=(bench.$o.$v); done
done
ds+=(bench.human.complete.wLnc bench.human.fragments.wLnc)
printf "%s\n" "${ds[@]}" | xargs -P 4 -I{} bash -c "bash $D/step5_08_full_predict_v2.sh {} hybrid > $L/hybrid.{}.log 2>&1; echo \"\$(date +%T) {} rc=\$?\""
echo "$(date +%T) ALLDONE"
