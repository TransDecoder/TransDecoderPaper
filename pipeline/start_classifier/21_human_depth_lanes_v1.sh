#!/bin/bash
# E33 / E31 human depth series: 3 lanes (<= 3 concurrent level jobs), each running its levels in order, niced.
# usage: 21_human_depth_lanes_v1.sh   (logs in $P/work/start_classifier/human_depth/lane_<i>.log)
P=/broad/hptmp/bhaas/TRANSDECODER_PAPER; S=$P/TransDecoderPaper/pipeline/start_classifier
O=$P/work/start_classifier/human_depth; mkdir -p $O
lane() { for level in "$@"; do nice -n 10 bash $S/20_human_depth_level_v1.sh $level || echo "[$(date)] $level FAILED"; done; echo LANE_DONE; }
lane d1 d32 > $O/lane_1.log 2>&1 &
lane d2 d16 > $O/lane_2.log 2>&1 &
lane d4 d8  > $O/lane_3.log 2>&1 &
wait
echo ALL_LANES_DONE
