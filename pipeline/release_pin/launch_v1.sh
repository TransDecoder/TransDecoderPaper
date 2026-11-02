#!/bin/bash
# Release-pin reruns (E40), 4 processes at a time, normal priority, stage by stage.
# usage: launch_v1.sh    (run in the background; progress in $P/work/release_pin/launch_v1.log)
P=/broad/hptmp/bhaas/TRANSDECODER_PAPER; T=$P/TransDecoderPaper; R=$T/pipeline/release_pin; L=$P/work/release_pin/logs; mkdir -p $L
export PATH=$P/env/td/bin:$PATH
ORGS="athal creinhardtii dmel human mouse pfalciparum spombe tthermophila"
NEG20=$(for o in $ORGS; do for v in complete.wRand fragments.wRand complete.wLnc fragments.wLnc; do [ -s $P/data/bench_v3/inputs/$o/$v.fa ] && echo "$o $v"; done; done)
pool() { # run the lines of stdin as commands, 4 at a time; one log per stage
  xargs -d '\n' -P ${N:-4} -I{} bash -c 'if bash -c "$1"; then echo "ok $1"; else echo "FAILED $1"; fi' _ {} ; }
stage() { echo "=== stage $1 start $(date)"; }

stage A
{ echo "PYTHON=$P/env/td/bin/python bash $T/regression/run_regression.sh $P/work/trees/TD_rel_41fc5ee $P/work/trees/PY_rel_d3892b5 $P/work/release_pin/regression > $L/regression.log 2>&1"
  for f in $(cd $P/data/bench_v3/inputs; ls */*.fa); do o=${f%%/*}; v=$(basename $f .fa); echo "bash $R/run_bench_v1.sh $o $v rel > $L/$o.$v.rel.log 2>&1"; done
  echo "$NEG20" | while read o v; do for l in precise sensitive; do echo "bash $R/run_bench_v1.sh $o $v $l > $L/$o.$v.$l.log 2>&1"; done; done
} | pool

stage B
{ for a in lraa stringtie trinity; do echo "bash $R/run_gm_v1.sh $a rel > $L/gm.$a.rel.log 2>&1"; done; } | pool
{ for c in a0.4_w17 norefine a0.4_w0 precise sensitive m100_noq m50_noq a0.2_w0 a0.2_w17 a0.2_w35 a0.4_w35 a0.6_w0 a0.6_w17 a0.6_w35 a1_w0 a1_w17 a1_w35; do
    for a in lraa stringtie trinity; do echo "bash $R/run_gm_v1.sh $a $c > $L/gm.$a.$c.log 2>&1"; done; done
  echo "bash $R/run_gm_v1.sh trinity retonly > $L/gm.trinity.retonly.log 2>&1"
} | pool

stage C
{ echo "$NEG20" | while read o v; do
    for l in m100_noq m90_noq m80_noq m70_noq m60_noq m50_noq q0.02 q0.1 q0.3 q0.5 nolong_q0.05 nolong_q0.1; do echo "bash $R/run_bench_v1.sh $o $v $l > $L/$o.$v.$l.log 2>&1"; done; done
  for f in $(cd $P/work/bench_v3/runs; ls -d */*/pytd_v7h_homology); do o=${f%%/*}; x=${f#*/}; v=${x%%/*}; echo "bash $R/run_bench_v1.sh $o $v retonly > $L/$o.$v.retonly.log 2>&1"; done
} | pool

stage D
{ for o in human mouse athal dmel; do for k in A Amix; do echo "bash $R/run_preset_v1.sh ref $o $k > $L/presetref.$o.$k.log 2>&1"; done; done; } | pool
{ cat $P/work/bench_v3/presets/jobs.txt $P/work/bench_v3/presets/jobs2.txt | grep -v '^cliff_c200 ' | while read s fa out o; do
    n=self; [[ $s == preset_* ]] && { n=${s%_q0.05}; }; echo "bash $R/run_preset_v1.sh test $n $fa $o > /dev/null 2>&1"; done
} | pool

stage E
$P/env/td/bin/python $T/analyses/start_param_sensitivity/run_grid_v2.py --jobs 4
stage END
