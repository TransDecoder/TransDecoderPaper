#!/bin/bash
# Rerun the 7 Perl/PyTD mismatches after unifying the GC estimate for length-based auto-accept
cd /broad/hptmp/bhaas/TRANSDECODER_PAPER
S=TransDecoderPaper/benchmark_v2/scripts
for ov in athal/complete athal/complete_fullutr athal/fragments creinhardtii/complete creinhardtii/complete_fullutr creinhardtii/fragments mouse/complete_fullutr; do
  org=${ov%%/*}; v=${ov#*/}
  for t in pytd_dev_v2 perl_td_dev_v2; do
    od=work/bench_v3/runs/$org/$v/$t
    mkdir -p work/bench_v3/rerun_old/$org/$v; mv $od work/bench_v3/rerun_old/$org/$v/ 
    ( $S/run_tool.sh $t data/bench_v3/inputs/$org/$v.fa $od $org > work/bench_v3/runs/logs/$org.$v.$t.rerun.log 2>&1 && echo "done $org $v $t" || echo "FAILED $org $v $t" ) &
  done
done
wait
echo ALL_DONE15
