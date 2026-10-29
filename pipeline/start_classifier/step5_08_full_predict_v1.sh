#!/bin/bash
# E33 step 5, full homology mode: fresh PyTD a5b24c3 Predict with UniRef50 (or Swiss-Prot) DIAMOND hits on ALL
# LongOrfs peptides + the same Pfam (--cut_ga) hits as the Swiss-Prot runs; start assist on (default).
# The LongOrfs output of the existing a5b24c3 run is linked, not copied (Predict only reads longest_orfs.*,
# base_freqs.dat; checked in predict.py), so every arm starts from the same LongOrfs.
# usage: step5_08_full_predict_v1.sh <dataset> <arm: uniref50|sprot>
#   run dir: $P/work/start_classifier/step5_full/<arm>/<dataset>/
set -euo pipefail
ds=$1; arm=$2
P=/broad/hptmp/bhaas/TRANSDECODER_PAPER; T=$P/TransDecoderPaper; B=$P/data/start_classifier/step5
PY=$P/env/td/bin/python; HS=$P/work/trees/PY_hstart
[ "$(git -C $HS rev-parse --short HEAD)" = a5b24c3 ] || { echo "PY_hstart not at a5b24c3" >&2; exit 1; }
[ -z "$(git -C $HS status --porcelain --untracked-files=no)" ] || { echo "PY_hstart modified" >&2; exit 1; }
gcode=Universal
case $ds in
  bench.*)
    org=$(cut -d. -f2 <<< $ds); v=$(cut -d. -f3- <<< $ds); R=$P/work/bench_v3/runs/$org/$v
    fa=$(readlink -f $R/pytd_v7h_homology/input.fa); L=$R/pytd_v7h_homology/input.fa.transdecoder_dir
    pfam=$R/pytd_v7_homology/pfam.domtblout; sprot=$R/pytd_v7_homology/blastp.outfmt6
    gcode=$(awk -F'\t' -v o=$org '$1==o{print $4}' $T/benchmark_v2/organisms.tsv) ;;
  sim.athal.*)
    lv=${ds##*.}; fa=$P/data/simtrinity_athal/depth/$lv/inputs/fragments_simtrinity.fa
    L=$P/work/start_classifier/athal_depth/longorfs/$lv/input.fa.transdecoder_dir
    pfam=$P/work/start_classifier/athal_depth/searches/$lv/pfam.domtblout; sprot=$P/work/start_classifier/athal_depth/searches/$lv/blastp.outfmt6 ;;
  sim.human.d8)
    fa=$P/data/simtrinity/depth/d8/inputs/fragments_simtrinity.fa
    L=$P/work/start_classifier/hstart_runs/sim.human.d8/on/input.fa.transdecoder_dir
    pfam=$B/searches/sprot/sim.human.d8.full/pfam.domtblout; sprot=$B/searches/sprot/sim.human.d8.full/blastp.outfmt6 ;;
  *) echo "unknown dataset $ds" >&2; exit 1 ;;
esac
cmp -s <(readlink -f $L/../input.fa) <(echo $fa) || cmp -s $L/../input.fa $fa || { echo "input mismatch $ds" >&2; exit 1; }
blast=$sprot; [ $arm = uniref50 ] && blast=$B/hits/uniref50/$ds/blastp.outfmt6
[ -s $blast ] && [ -s $pfam ] || { echo "missing hits $ds" >&2; exit 1; }
run=$P/work/start_classifier/step5_full/$arm/$ds
[ -e $run ] && { echo "exists: $run" >&2; exit 1; }
mkdir -p $run/input.fa.transdecoder_dir; cd $run
ln -s $fa input.fa
for f in longest_orfs.pep longest_orfs.cds longest_orfs.gff3 base_freqs.dat longorfs.params.json; do ln -s $L/$f input.fa.transdecoder_dir/$f; done
echo "a5b24c3" > pytd_version.txt
printf "blast\t%s\npfam\t%s\nlongorfs\t%s\n" $blast $pfam $L > hits_used.tsv
t0=$(date +%s.%N)
PYTHONPATH=$HS nice -n 10 $PY -m pytransdecoder.cli predict -t input.fa -G "$gcode" \
  --retain-blastp-hits $blast --retain-pfam-hits $pfam > predict.log 2>&1
echo "$(echo "$(date +%s.%N) - $t0" | bc)" > runtime_sec.txt
$PY $T/benchmark_v2/scripts/convert_predictions.py transdecoder_gff3 input.fa.transdecoder.gff3 > preds.tsv
# keep the footprint small: compress the bulky outputs (preds.tsv, logs and the homology table stay readable)
gzip -f input.fa.transdecoder.cds input.fa.transdecoder.pep input.fa.transdecoder.bed input.fa.transdecoder.gff3 \
  input.fa.transdecoder_dir/best_candidates.gff3 input.fa.transdecoder_dir/best_candidates.gff3.revised_starts.gff3 \
  input.fa.transdecoder_dir/longest_orfs.cds.scores input.fa.transdecoder_dir/decoy_qvalues.tsv \
  input.fa.transdecoder_dir/start_refinement.homology.tsv 2>/dev/null || true
[ -e input.fa.transdecoder_dir/__checkpoints_TDpredict/top_training_orfs.cds ] && gzip -f input.fa.transdecoder_dir/__checkpoints_TDpredict/top_training_orfs.cds
echo "$ds $arm done $(tail -n +2 preds.tsv | wc -l) ORFs $(cat runtime_sec.txt)s"
