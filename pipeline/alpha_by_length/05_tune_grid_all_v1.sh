#!/bin/bash
# E36: 05_tune_grid_v1.py in 3 lanes (nice 10). Outputs $P/data/alpha_by_length/grid/<lane>.tsv.gz, log grid.log
P=/broad/hptmp/bhaas/TRANSDECODER_PAPER; T=$P/TransDecoderPaper; S=$T/pipeline/alpha_by_length; G=$P/data/alpha_by_length/grid
PY=$P/env/startclf/bin/python; mkdir -p $G; : > $G/grid.log
run() { nice -n 10 $PY $S/05_tune_grid_v1.py $T/data/alpha_by_length/manifest.tsv $P/data/alpha_by_length/labels $G/$1.tsv.gz "${@:2}" >> $G/grid.log 2>&1 || echo "FAILED $1" >> $G/grid.log; }
run const const &
run lane1 logL logO &
run lane2 OL uL &
run lane3 stepOL stepL &
wait
echo ALLDONE >> $G/grid.log
