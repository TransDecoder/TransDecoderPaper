#!/bin/bash
# E31 lane: for each <level>:<fraction> argument, in order: subsample, Trinity (CPU 4), label + run + score.
P=/broad/hptmp/bhaas/TRANSDECODER_PAPER
S=$P/TransDecoderPaper/pipeline/simtrinity
L=$P/work/simtrinity/depth/logs
for lf in "$@"; do
  level=${lf%%:*}; frac=${lf##*:}
  echo "[$(date)] $level ($frac) start"
  bash $S/11_depth_subsample.sh $level $frac > $L/$level.subsample.log 2>&1 || { echo "$level subsample failed"; exit 1; }
  CPU=4 bash $S/08b_trinity_gg_level.sh $level > $L/$level.trinity.log 2>&1 || { echo "$level trinity failed"; exit 1; }
  bash $S/12_depth_label_and_run.sh $level > $L/$level.run.log 2>&1 || { echo "$level label/run failed"; exit 1; }
  echo "[$(date)] $level LEVEL_DONE"
done
