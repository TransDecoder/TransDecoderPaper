#!/bin/bash
# E33 / E32: PyTD a5b24c3 on one simulated Arabidopsis depth level, three fresh Predict conditions on one LongOrfs:
#   A ab initio; B homology (DIAMOND sprot_minus_3702 -k 1 1e-5 + hmmsearch --cut_ga Pfam-A on ALL LongOrfs
#   peptides, commands as run_tool.sh pytd_v7_homology), start assist on (default); C = B hits, --no-homology-start-assist.
# Reads E32's inputs/labels (read only); writes only under $P/work/start_classifier/athal_depth/.
# Each condition is scored with E32's scorer (score_simtrinity.py) in <cond>/depth/<level>/{runs/pytd_v7,scores}.
# usage: 17_athal_depth_level_v1.sh <level>
set -euo pipefail
level=$1
P=/broad/hptmp/bhaas/TRANSDECODER_PAPER; T=$P/TransDecoderPaper; H=$P/data/homology
PY=$P/env/td/bin/python; HS=$P/work/trees/PY_hstart; export PATH=$P/env/td/bin:$PATH
D=$P/data/simtrinity_athal/depth/$level; O=$P/work/start_classifier/athal_depth
fa=$D/inputs/fragments_simtrinity.fa
[ "$(git -C $HS rev-parse --short HEAD)" = a5b24c3 ] || { echo "PY_hstart not at a5b24c3" >&2; exit 1; }
L=$O/longorfs/$level
[ -e $L ] && { echo "exists: $L" >&2; exit 1; }
mkdir -p $L; cd $L; ln -s $fa input.fa
PYTHONPATH=$HS $PY -m pytransdecoder.cli longorfs -t input.fa -G Universal -S > longorfs.log 2>&1
pep=$L/input.fa.transdecoder_dir/longest_orfs.pep
echo "[$(date)] $level LongOrfs: $(grep -c '>' $pep) peptides"
predict() { # cond extra...
  local c=$1; shift
  local r=$O/$c/depth/$level/runs/pytd_v7
  mkdir -p $r; cd $r; ln -s $fa input.fa
  cp -r $L/input.fa.transdecoder_dir .   # the same LongOrfs for all conditions
  echo a5b24c3 > pytd_version.txt
  local t0=$(date +%s.%N)
  PYTHONPATH=$HS $PY -m pytransdecoder.cli predict -t input.fa -G Universal "$@" > predict.log 2>&1
  $PY $T/benchmark_v2/scripts/convert_predictions.py transdecoder_gff3 input.fa.transdecoder.gff3 > preds.tsv
  echo "$(echo "$(date +%s.%N) - $t0" | bc)" > runtime_sec.txt
  [ $c = A ] && $PY $T/benchmark_v2/scripts/convert_predictions.py transdecoder_gff3 --longest-per-transcript \
      input.fa.transdecoder_dir/longest_orfs.gff3 > longest_orf.preds.tsv
  echo "[$(date)] $level $c predict done: $(tail -n +2 preds.tsv | wc -l) ORFs"
}
predict A &
pa=$!
S=$O/searches/$level; mkdir -p $S; cd $S
diamond blastp --query $pep --db "$H/sprot_minus_3702.dmnd" --max-target-seqs 1 \
  --evalue 1e-5 --outfmt 6 --threads 3 --out blastp.outfmt6 > diamond.log 2>&1
mkdir -p pfam_chunks
awk -v d=pfam_chunks -v n=3 '/^>/{f=sprintf("%s/c%03d.pep",d,(i++)%n)} {print > f}' $pep
ls pfam_chunks/*.pep | xargs -P 3 -I{} hmmsearch --cpu 1 --cut_ga --noali -o /dev/null \
  --domtblout {}.domtblout "$H/Pfam-A.hmm" {} > hmmsearch.log 2>&1
cat pfam_chunks/*.domtblout > pfam.domtblout
echo "[$(date)] $level searches done: $(wc -l < blastp.outfmt6) blast, $(grep -vc '^#' pfam.domtblout) pfam lines"
wait $pa
predict B --retain-blastp-hits $S/blastp.outfmt6 --retain-pfam-hits $S/pfam.domtblout &
pb=$!
predict C --retain-blastp-hits $S/blastp.outfmt6 --retain-pfam-hits $S/pfam.domtblout --no-homology-start-assist &
pc=$!
wait $pb; wait $pc
for c in A B C; do
  $PY $T/pipeline/simtrinity/score_simtrinity.py $fa $D/labels/sim_trinity $O/$c/depth/$level/runs $O/$c/depth/$level/scores \
    > $O/$c/depth/$level/score.log 2>&1
done
echo "[$(date)] $level LEVEL_DONE"
