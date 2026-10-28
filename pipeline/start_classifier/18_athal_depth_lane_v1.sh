#!/bin/bash
# E33 / E32 lane: levels in order; d32 only after E32's own d32 label/run/score has finished
# (its run log says LEVEL_DONE and its scores_overall.tsv exists).
P=/broad/hptmp/bhaas/TRANSDECODER_PAPER; S=$P/TransDecoderPaper/pipeline/start_classifier
for level in "$@"; do
  if [ $level = d32 ]; then
    until grep -q "d32: LEVEL_DONE" $P/work/simtrinity_athal/depth/logs/d32.run.log 2>/dev/null && \
          [ -s $P/work/simtrinity_athal/depth/d32/scores/scores_overall.tsv ]; do sleep 120; done
    echo "[$(date)] E32 d32 scored; starting"
  fi
  bash $S/17_athal_depth_level_v1.sh $level || echo "[$(date)] $level FAILED"
done
echo LANE_DONE
