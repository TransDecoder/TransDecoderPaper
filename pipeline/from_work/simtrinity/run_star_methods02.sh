#!/bin/bash
# Run on methods02 (enough free memory; methods01 had ~29 GB available). Same 07_star.sh.
P=/broad/hptmp/bhaas/TRANSDECODER_PAPER
cd $P
echo "[$(date)] host $(hostname); starting STAR" >> $P/work/simtrinity/logs/star_gate.log
if bash $P/TransDecoderPaper/pipeline/simtrinity/07_star.sh > $P/work/simtrinity/logs/07_star.out 2>&1; then
  echo "[$(date)] STAR finished" >> $P/work/simtrinity/logs/star_gate.log
else
  echo "[$(date)] STAR failed" >> $P/work/simtrinity/logs/star_gate.log
fi
