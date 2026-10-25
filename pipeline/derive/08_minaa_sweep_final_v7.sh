#!/bin/bash
# Notebook 08 (GM12878 parameters), Fig. 7b: post-hoc minimum-length filter on the RELEASE
# (final_v7, E25) ab initio predictions. Same procedure as E04's minaa_sweep.sh, which used the earlier
# final_abinitio predictions: eval_genome_cds.py --min-aa L keeps predictions encoding >= L aa
# (stop excluded) before collapsing to distinct CDSs.
#
# Step 1 (skipped when the summary exists): eval runs -> $P/work/gm12878_td/eval/minaa_v7/<asm>_min<L>.*
# Step 2: compact table -> data/gm12878/final_v7/08_minaa_sweep_v7.tsv (with header).
# usage: bash pipeline/derive/08_minaa_sweep_final_v7.sh [jobs=2]
set -euo pipefail
JOBS=${1:-2}
P=/broad/hptmp/bhaas/TRANSDECODER_PAPER; D=$P/data/gm12878; W=$P/work/gm12878_td
T=$(cd "$(dirname "$0")/../.." && pwd)
OUT=$W/eval/minaa_v7; mkdir -p $OUT
export P D W T OUT
export PYTHONPATH=$P/work/trees/PY_v7   # eval_genome_cds.py imports pytransdecoder.genome
export ST_lraa=$D/longread_lraa/lraa_denovo/LRAA.gtf ST_stringtie=$D/star_stringtie/stringtie.gtf \
       ST_trinity=$D/trinity/Trinity.minimap2.bam
run_one() {
  a=$1; L=$2; st=ST_$a
  [ -s $OUT/${a}_min$L.summary.tsv ] && { echo "have $a $L"; return 0; }
  $P/env/td/bin/python $T/genome_annotation/scripts/eval_genome_cds.py \
    --pred $W/$a/final_v7_abinitio/transcripts.fa.transdecoder.genome.gff3 --structures ${!st} \
    --mane $P/data/annot/human/MANE.GRCh38.v1.5.ensembl_genomic.gtf.gz \
    --gencode $D/ref/gencode.v50.primary_assembly.annotation.gtf \
    --label ${a}_v7_min$L --min-aa $L -o $OUT/${a}_min$L > $OUT/${a}_min$L.log 2>&1 \
    && echo "ok $a $L" || echo "FAILED $a $L"
}
export -f run_one
for a in lraa stringtie trinity; do for L in 50 60 70 80 90 100 150 200 300; do echo "$a $L"; done; done \
  | xargs -P $JOBS -L 1 bash -c 'run_one $0 $1'

keys="distinct_predicted_cds MANE_exact MANE_sens_exact MANE_sens_exact_given_assembled GENCODE_complete_sens_compatible pred_all_GENCODE_exact_frac pred_all_GENCODE_compatible_frac pred_all_no_GENCODE_CDS_overlap_frac pred_complete_GENCODE_exact_frac"
{ printf "assembly\tmin_aa"; for k in $keys; do printf "\t%s" $k; done; printf "\n"
  for a in lraa stringtie trinity; do for L in 50 60 70 80 90 100 150 200 300; do
    f=$OUT/${a}_min$L.summary.tsv; printf "%s\t%s" $a $L
    for k in $keys; do printf "\t%s" "$(awk -F'\t' -v k=$k '$1==k{print $2}' $f)"; done; printf "\n"
  done; done; } > $T/data/gm12878/final_v7/08_minaa_sweep_v7.tsv
echo "wrote $T/data/gm12878/final_v7/08_minaa_sweep_v7.tsv"
