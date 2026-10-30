#!/bin/bash
# E36 v2: simpler alpha forms (stepO, rampO, linO) plus const and logO for reference. Outputs $P/data/alpha_by_length/grid_v2/
P=/broad/hptmp/bhaas/TRANSDECODER_PAPER; T=$P/TransDecoderPaper; S=$T/pipeline/alpha_by_length; G=$P/data/alpha_by_length/grid_v2
PY=$P/env/startclf/bin/python; mkdir -p $G; : > $G/grid.log
run() { nice -n 10 $PY $S/05_tune_grid_v2.py $T/data/alpha_by_length/manifest.tsv $P/data/alpha_by_length/labels $G/$1.tsv.gz "${@:2}" >> $G/grid.log 2>&1 || echo "FAILED $1" >> $G/grid.log; }
run lane1 const rampO linO &
run lane2 logO &
run lane3 stepO &
wait
echo ALLDONE >> $G/grid.log
