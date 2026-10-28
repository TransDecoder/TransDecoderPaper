#!/bin/bash
# E33 step 2: self-trained decisions for every input in the manifest (serial, <= 3 threads).
set -euo pipefail
P=/broad/hptmp/bhaas/TRANSDECODER_PAPER; T=$P/TransDecoderPaper; D=$T/data/start_classifier
export OMP_NUM_THREADS=3 OPENBLAS_NUM_THREADS=3 MKL_NUM_THREADS=3 PYTHONPATH=$P/work/trees/PY_v7
O=$D/selftrain; mkdir -p $O; log=$O/run.log; : > $log
tail -n +2 $D/step2_manifest.tsv | while IFS=$'\t' read -r ds run truth excl gcode; do
  $P/env/startclf/bin/python $T/pipeline/start_classifier/04_selftrain.py $run $gcode $D/features/$ds.tsv.gz $ds $O/$ds 2>> $log
done
$P/env/startclf/bin/python $T/pipeline/start_classifier/05_evaluate.py $D/features $O $D/step2_manifest.tsv $D/step2 2>> $log
echo ALLDONE >> $log
