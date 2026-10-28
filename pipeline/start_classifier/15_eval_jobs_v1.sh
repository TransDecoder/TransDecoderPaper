#!/bin/bash
# E33 adoption: all evaluation runs with the committed PY_hstart, at most 3 at a time.
#   GM12878 Trinity (13_*), 84 benchmark homology inputs (pytd_v7h_homology, hits from pytd_v7_homology),
#   4 simulated assemblies x {assist on, off} with the step-3 candidate searches.
set -euo pipefail
P=/broad/hptmp/bhaas/TRANSDECODER_PAPER; T=$P/TransDecoderPaper; S=$T/pipeline/start_classifier
J=$P/work/start_classifier/eval_jobs_v1.txt; : > $J
echo "bash $S/13_gm12878_trinity_v7h_v1.sh > $P/work/gm12878_td/logs/final_v7h_trinity_homology.log 2>&1" >> $J
for ds in sim.human.d1 sim.human.d8 sim.human.d32 sim.athal.d1; do
  run=$(awk -F'\t' -v d=$ds '$1==d{print $2}' $T/data/start_classifier/step2_manifest.tsv)
  for arm in on off; do
    ex=""; [ $arm = off ] && ex="--no-homology-start-assist"
    echo "bash $S/11_run_pytd_v7h_homology_v1.sh $P/work/start_classifier/hstart_runs/$ds/$arm $run/input.fa Universal 1 $P/work/start_classifier/homology/$ds $ex" >> $J
  done
done
for h in $P/work/bench_v3/runs/*/*/pytd_v7_homology; do
  v=$(basename $(dirname $h)); org=$(basename $(dirname $(dirname $h)))
  gcode=$(awk -F'\t' -v o=$org '$1==o{print $4}' $T/benchmark_v2/organisms.tsv)
  st=1; [[ $v == *unstranded* ]] && st=0
  echo "bash $S/11_run_pytd_v7h_homology_v1.sh $P/work/bench_v3/runs/$org/$v/pytd_v7h_homology $P/data/bench_v3/inputs/$org/$v.fa $gcode $st $h" >> $J
done
wc -l < $J
xargs -P 3 -I{} -d '\n' bash -c '{} || echo "FAILED: {}"' < $J
echo ALLDONE
