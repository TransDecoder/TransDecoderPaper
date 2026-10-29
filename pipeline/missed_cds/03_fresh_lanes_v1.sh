#!/bin/bash
# E34: two serial lanes of 02_fresh_runs_v1.sh (2 processes at a time on the shared host)
P=/broad/hptmp/bhaas/TRANSDECODER_PAPER; S=$P/TransDecoderPaper/pipeline/missed_cds; L=$P/work/missed_cds/logs
lane() { for job in "$@"; do bash $S/02_fresh_runs_v1.sh ${job%:*} ${job#*:} >> $L/lanes_v1.log 2>&1 || echo "FAILED $job" >> $L/lanes_v1.log; done; }
lane minlen30:human.d1 minlen30:athal.d1.A sensitive:athal.d1.A &
lane minlen30:human.d8 minlen30:human.d32 minlen30:athal.d16.A sensitive:human.d1 &
wait; echo ALLDONE >> $L/lanes_v1.log
