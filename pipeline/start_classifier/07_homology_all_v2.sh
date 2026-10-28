#!/bin/bash
# E33 step 3 driver: one dataset at a time (3 threads).
set -euo pipefail
P=/broad/hptmp/bhaas/TRANSDECODER_PAPER; T=$P/TransDecoderPaper; O=$P/work/start_classifier/homology
S=$T/pipeline/start_classifier/07_homology_candidates_v2.sh
mkdir -p $O
bash $S bench.human.fragments $P/work/bench_v3/runs/human/fragments/pytd_v7 9606 $O/bench.human.fragments
bash $S sim.human.d1 $P/work/simtrinity/runs/pytd_v7 9606 $O/sim.human.d1
bash $S sim.human.d8 $P/work/simtrinity/depth/d8/runs/pytd_v7 9606 $O/sim.human.d8
bash $S sim.human.d32 $P/work/simtrinity/depth/d32/runs/pytd_v7 9606 $O/sim.human.d32
bash $S sim.athal.d1 $P/work/simtrinity_athal/depth/d1/runs/pytd_v7 3702 $O/sim.athal.d1
echo ALLDONE
