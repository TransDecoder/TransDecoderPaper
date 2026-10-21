#!/bin/bash
# Minimum ORF length sweep with GENCODE lncRNA negatives (human, mouse)
cd /broad/hptmp/bhaas/TRANSDECODER_PAPER
VARIANTS="complete.wLnc fragments.wLnc" JOBS=4 TransDecoderPaper/benchmark_v2/scripts/run_all.sh \
  data/bench_v3/inputs work/bench_v3/runs pytd_dev_kozak_add pytd_dev_m90 pytd_dev_m80 pytd_dev_m70 pytd_dev_m60 pytd_dev_m50
echo ALL_DONE11
