#!/bin/bash
# Balanced mixed inputs fragments_mixed_t{1,2,3} (25% complete / 5' partial / 3' partial / internal; E28).
# Same run_tool.sh commands as for `fragments` (the 'fragments' prefix selects CodAn *_partial models
# and RNAsamba partial weights); fresh run directories. Two lanes, 4 jobs at most: CodAn (2 at a
# time, longest inputs first) and all other tools (2 at a time). GeneMarkS-T, CodAn, CPC2, RNAsamba
# and Borf are not run on T. thermophila (no ciliate code).
cd /broad/hptmp/bhaas/TRANSDECODER_PAPER
RT=/broad/hptmp/bhaas/TRANSDECODER_PAPER/TransDecoderPaper/benchmark_v2/scripts/run_tool.sh
IN=data/bench_v3/inputs; RUNS=work/bench_v3/runs; LOGS=$RUNS/logs
VARS="fragments_mixed_t1 fragments_mixed_t2 fragments_mixed_t3"
NA="gmst codan cpc2 rnasamba borf"
run() {  # tool variant org
  local t=$1 v=$2 o=$3 log=$LOGS/$3.$2.$1.log
  if $RT $t $IN/$o/$v.fa $RUNS/$o/$v/$t $o > $log 2>&1; then echo "done $o $v $t"; else echo "FAILED $o $v $t (see $log)"; fi
}
export -f run; export RT IN RUNS LOGS
{
  for o in athal human mouse creinhardtii dmel pfalciparum spombe; do for v in $VARS; do echo "codan $v $o"; done; done |
    xargs -P 2 -L 1 bash -c 'run "$@"' _
  echo CODAN_LANE_DONE
} &
{
  for t in pytd_v7 perl_td_v7 perl_td_v600 prodigal gmst cpat cpc2 rnasamba borf; do
    for v in $VARS; do for o in athal creinhardtii dmel human mouse pfalciparum spombe tthermophila; do
      [ $o = tthermophila ] && [[ " $NA " == *" $t "* ]] && continue
      echo "$t $v $o"
    done; done
  done | xargs -P 2 -L 1 bash -c 'run "$@"' _
  echo OTHER_LANE_DONE
} &
wait
echo ALL_DONE_MIXED
