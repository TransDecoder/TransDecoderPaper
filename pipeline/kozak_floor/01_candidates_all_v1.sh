#!/bin/bash
# E35: run 01_candidates_v1.py on every manifest row; 4 lanes, nice 10. Log: $P/data/kozak_floor/candidates/extract.log
set -uo pipefail
P=/broad/hptmp/bhaas/TRANSDECODER_PAPER; T=$P/TransDecoderPaper
O=$P/data/kozak_floor/candidates; mkdir -p $O
log=$O/extract.log; : > $log
export P O
tail -n +2 $T/data/kozak_floor/manifest.tsv | cut -f1,6,7 | \
  xargs -P 4 -L 1 bash -c 'nice -n 10 env PYTHONPATH=$P/work/trees/PY_hstart $P/env/startclf/bin/python '"$T"'/pipeline/kozak_floor/01_candidates_v1.py "$2" "$1" "$0" $O/$0.tsv.gz 2>> $O/extract.log || echo "FAILED $0" >> $O/extract.log'
echo ALLDONE >> $log
