#!/bin/bash
# E32 lane: for each <level>:<fraction> (d1:1 = all reads, no subsampling): subsample, Trinity, label + run + score.
P=/broad/hptmp/bhaas/TRANSDECODER_PAPER
A=$P/TransDecoderPaper/pipeline/simtrinity_athal
L=$P/work/simtrinity_athal/depth/logs
mkdir -p $L
for lf in "$@"; do
  level=${lf%%:*}; frac=${lf##*:}
  echo "[$(date)] $level ($frac) start"
  if [ $level != d1 ]; then bash $A/11_depth_subsample_athal.sh $level $frac > $L/$level.subsample.log 2>&1 || { echo "$level subsample failed"; exit 1; }; fi
  bash $A/08_trinity_gg_athal.sh $level > $L/$level.trinity.log 2>&1 || { echo "$level trinity failed"; exit 1; }
  bash $A/12_label_and_run_athal.sh $level > $L/$level.run.log 2>&1 || { echo "$level label/run failed"; exit 1; }
  echo "[$(date)] $level LEVEL_DONE"
done
