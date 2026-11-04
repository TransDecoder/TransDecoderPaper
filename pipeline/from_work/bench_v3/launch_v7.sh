#!/bin/bash
# final benchmark at the release commits (PyTD 3a793ff, Perl 41fc5ee); fresh run directories
cd /broad/hptmp/bhaas/TRANSDECODER_PAPER
RA=/broad/hptmp/bhaas/TRANSDECODER_PAPER/TransDecoderPaper/benchmark_v2/scripts/run_all.sh
JOBS=6 $RA data/bench_v3/inputs work/bench_v3/runs pytd_v7 perl_td_v7
VARIANTS="complete.wRand fragments.wRand complete.wLnc fragments.wLnc" JOBS=6 $RA data/bench_v3/inputs work/bench_v3/runs pytd_v7_precise pytd_v7_sensitive
THREADS=4 JOBS=3 $RA data/bench_v3/inputs work/bench_v3/runs pytd_v7_homology
echo ALL_DONE
