#!/bin/bash
# PyTD development branch: frame-contrast start refinement (first ATG / PWM-chosen ATG)
cd /broad/hptmp/bhaas/TRANSDECODER_PAPER
JOBS=12 TransDecoderPaper/benchmark_v2/scripts/run_all.sh data/bench_v3/inputs work/bench_v3/runs pytd_dev pytd_dev_pwm
echo ALL_DONE5
