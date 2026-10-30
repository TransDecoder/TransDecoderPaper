#!/bin/bash
# E36: manifest of the finished runs evaluated (all read only). Writes $T/data/alpha_by_length/manifest.tsv.
# = the 90 E35 rows (data/kozak_floor/manifest.tsv; candidate tables reused from $P/data/kozak_floor/candidates)
# + human d2, d4, d16 ab initio (pytd_v7, 3a793ff; no homology runs exist at these depths)
# + GM12878 Trinity GG ab initio (final_v7_abinitio, 3a793ff; E34 link dir) and with the start assist
#   (final_v7h_homology, a5b24c3; link dir made here), truth = E29 transcript-level MANE truth (E34 copy with n_codons).
# Extra column cand = per-candidate table (E35's, or $P/data/alpha_by_length/candidates/<dataset>.tsv.gz made by 01).
set -euo pipefail
P=/broad/hptmp/bhaas/TRANSDECODER_PAPER; T=$P/TransDecoderPaper
out=$T/data/alpha_by_length/manifest.tsv
NEWC=$P/data/alpha_by_length/candidates
printf "dataset\tarm\tgroup\torg\tinput\tgcode\trun_dir\ttruth\texcluded\tnegatives\tcand\n" > $out
tail -n +2 $T/data/kozak_floor/manifest.tsv | awk -v C=$P/data/kozak_floor/candidates 'BEGIN{OFS="\t"}{print $0, C"/"$1".tsv.gz"}' >> $out
for d in d2 d4 d16; do
  lab=$P/data/simtrinity/depth/$d/labels/sim_trinity
  printf "sim.human.$d\tabinitio\tsim_human\thuman\t$d\tUniversal\t$P/work/simtrinity/depth/$d/runs/pytd_v7\t$lab.truth.tsv\t$lab.excluded.txt\t$lab.negatives.txt\t$NEWC/sim.human.$d.tsv.gz\n" >> $out
done
# GM12878 Trinity: link dirs so that the run looks like <dir>/input.fa + input.fa.transdecoder_dir + preds.tsv
G=$P/data/alpha_by_length/gm12878_inputs/v7h
if [ ! -e $G ]; then
  mkdir -p $G; src=$P/work/gm12878_td/trinity/final_v7h_homology
  ln -s $src/transcripts.fa $G/input.fa
  ln -s $src/transcripts.fa.transdecoder_dir $G/input.fa.transdecoder_dir
  ln -s $src/preds.tsv $G/preds.tsv
fi
gl=$P/data/missed_cds/gm12878_inputs/labels/mane
printf "gm12878.trinity\tabinitio\tgm12878\thuman\ttrinity\tUniversal\t$P/data/missed_cds/gm12878_inputs/run\t$gl.truth.tsv\t-\t-\t$NEWC/gm12878.trinity.tsv.gz\n" >> $out
printf "gm12878.trinity.H\thomology\tgm12878\thuman\ttrinity\tUniversal\t$G\t$gl.truth.tsv\t-\t-\t$NEWC/gm12878.trinity.H.tsv.gz\n" >> $out
tail -n +2 $out | while IFS=$'\t' read -r ds _ _ _ _ _ run truth ex neg cand; do
  for f in $run/preds.tsv $run/input.fa $run/input.fa.transdecoder_dir/best_candidates.gff3.revised_starts.gff3 $truth; do [ -e $f ] || echo "MISSING $ds $f"; done
  for f in $ex $neg; do [ "$f" = - ] || [ -e $f ] || echo "MISSING $ds $f"; done
  case $cand in $NEWC/*) ;; *) [ -e $cand ] || echo "MISSING $ds $cand";; esac
done
wc -l $out
