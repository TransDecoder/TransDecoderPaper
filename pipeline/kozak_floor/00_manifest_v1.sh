#!/bin/bash
# E35: manifest of the finished runs evaluated (all read only). Writes $T/data/kozak_floor/manifest.tsv.
# arm abinitio: release runs (benchmark pytd_v7 3a793ff; human sims pytd_v7 3a793ff; athal depth A a5b24c3)
# arm homology: start assist on (benchmark pytd_v7h_homology, human sims hstart_runs/<ds>/on, athal depth B; all a5b24c3)
set -euo pipefail
P=/broad/hptmp/bhaas/TRANSDECODER_PAPER; T=$P/TransDecoderPaper
out=$T/data/kozak_floor/manifest.tsv
printf "dataset\tarm\tgroup\torg\tinput\tgcode\trun_dir\ttruth\texcluded\tnegatives\n" > $out
while IFS=$'\t' read -r org _ _ gcode _; do
  [ "$org" = org ] && continue
  for v in complete complete_fullutr fragments fragments_mixed_t1 fragments_mixed_t2 fragments_mixed_t3; do
    g=bench_${v%_t[123]}
    printf "bench.$org.$v\tabinitio\t$g\t$org\t$v\t$gcode\t$P/work/bench_v3/runs/$org/$v/pytd_v7\t$P/data/bench_v3/inputs/$org/$v.truth.tsv\t-\t-\n" >> $out
  done
  for v in complete complete_fullutr fragments; do
    printf "bench.$org.$v.H\thomology\tbench_$v\t$org\t$v\t$gcode\t$P/work/bench_v3/runs/$org/$v/pytd_v7h_homology\t$P/data/bench_v3/inputs/$org/$v.truth.tsv\t-\t-\n" >> $out
  done
done < $T/benchmark_v2/organisms.tsv
for d in d1 d8 d32; do
  if [ $d = d1 ]; then run=$P/work/simtrinity/runs/pytd_v7; lab=$P/data/simtrinity/labels/sim_trinity
  else run=$P/work/simtrinity/depth/$d/runs/pytd_v7; lab=$P/data/simtrinity/depth/$d/labels/sim_trinity; fi
  printf "sim.human.$d\tabinitio\tsim_human\thuman\t$d\tUniversal\t$run\t$lab.truth.tsv\t$lab.excluded.txt\t$lab.negatives.txt\n" >> $out
  printf "sim.human.$d.H\thomology\tsim_human\thuman\t$d\tUniversal\t$P/work/start_classifier/hstart_runs/sim.human.$d/on\t$lab.truth.tsv\t$lab.excluded.txt\t$lab.negatives.txt\n" >> $out
done
for d in d1 d2 d4 d8 d16 d32; do
  lab=$P/data/simtrinity_athal/depth/$d/labels/sim_trinity
  printf "sim.athal.$d\tabinitio\tsim_athal\tathal\t$d\tUniversal\t$P/work/start_classifier/athal_depth/A/depth/$d/runs/pytd_v7\t$lab.truth.tsv\t$lab.excluded.txt\t$lab.negatives.txt\n" >> $out
  printf "sim.athal.$d.H\thomology\tsim_athal\tathal\t$d\tUniversal\t$P/work/start_classifier/athal_depth/B/depth/$d/runs/pytd_v7\t$lab.truth.tsv\t$lab.excluded.txt\t$lab.negatives.txt\n" >> $out
done
# every path must exist
tail -n +2 $out | while IFS=$'\t' read -r ds _ _ _ _ _ run truth ex neg; do
  for f in $run/preds.tsv $run/input.fa.transdecoder_dir/best_candidates.gff3.revised_starts.gff3 $truth; do [ -e $f ] || echo "MISSING $ds $f"; done
  for f in $ex $neg; do [ "$f" = - ] || [ -e $f ] || echo "MISSING $ds $f"; done
done
wc -l $out
