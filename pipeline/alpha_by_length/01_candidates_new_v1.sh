#!/bin/bash
# E36: E35's 01_candidates_v1.py (unchanged) on the manifest rows whose candidate table is not E35's
# (human d2/d4/d16 ab initio, GM12878 Trinity ab initio and assist). 3 lanes, nice 10.
# Log: $P/data/alpha_by_length/candidates/extract.log
set -uo pipefail
P=/broad/hptmp/bhaas/TRANSDECODER_PAPER; T=$P/TransDecoderPaper
O=$P/data/alpha_by_length/candidates; mkdir -p $O
log=$O/extract.log; : > $log
export P O
tail -n +2 $T/data/alpha_by_length/manifest.tsv | awk -F'\t' -v O=$O 'index($11, O)==1' | cut -f1,6,7 | \
  xargs -P 3 -L 1 bash -c 'nice -n 10 env PYTHONPATH=$P/work/trees/PY_hstart $P/env/startclf/bin/python '"$T"'/pipeline/kozak_floor/01_candidates_v1.py "$2" "$1" "$0" $O/$0.tsv.gz 2>> $O/extract.log || echo "FAILED $0" >> $O/extract.log'
echo ALLDONE >> $log
