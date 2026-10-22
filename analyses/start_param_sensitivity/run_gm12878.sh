#!/bin/bash
# GM12878 check of the start-refinement parameters on a sub-grid (+ --no-refine-starts):
# alpha in $GM_ALPHAS (default 0.2 0.4 0.6 1.0) x w in $GM_WEIGHTS (default 0 25 50), for the StringTie
# and LRAA assemblies (alpha 1.0 added because benchmark optima for two organisms lie at alpha >= 0.8).
# Reuses the LongOrfs output of the final-default runs ($P/work/gm12878_td/<asm>/final_abinitio;
# longest_orfs.{cds,gff3,pep} + base_freqs.dat symlinked, no checkpoints, no .cds.scores),
# then predict -> genome-map -> eval_genome_cds.py (MANE v1.5 / GENCODE v50), as in
# genome_annotation/scripts/gm12878_final_defaults.sh (Ribo-seq matching omitted).
# usage: run_gm12878.sh [jobs]      (writes $P/work/start_param_grid/gm12878/<asm>/<config>/)
set -euo pipefail
jobs=${1:-6}
P=/broad/hptmp/bhaas/TRANSDECODER_PAPER; D=$P/data/gm12878
PY=$P/env/td/bin/python; export PATH=$P/env/td/bin:$PATH; export PYTHONPATH=$P/work/trees/PY_dev
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
OUT=$P/work/start_param_grid/gm12878

one() {
  asm=$1 cfg=$2
  case $asm in
    stringtie) structs=$D/star_stringtie/stringtie.gtf;;
    lraa) structs=$D/longread_lraa/lraa_denovo/LRAA.gtf;;
  esac
  src=$P/work/gm12878_td/$asm/final_abinitio
  d=$OUT/$asm/$cfg; w=$d/transcripts.fa.transdecoder_dir
  [ -s $d/eval.summary.tsv ] && { echo "skip $asm $cfg"; return; }
  mkdir -p $w; cd $d
  ln -sfn $(readlink -f $src/transcripts.fa) transcripts.fa
  for f in longest_orfs.cds longest_orfs.gff3 longest_orfs.pep base_freqs.dat; do ln -sfn $src/transcripts.fa.transdecoder_dir/$f $w/$f; done
  [ -d $w/__checkpoints_TDpredict ] && mv $w/__checkpoints_TDpredict $w/__checkpoints_TDpredict.stale.$(date +%s)
  if [ $cfg = norefine ]; then opt=(--no-refine-starts)
  else a=${cfg#a}; a=${a%_w*}; wt=${cfg#*_w}; opt=(--start-contrast-alpha $a --start-kozak-weight $wt); fi
  t0=$(date +%s)
  $PY -m pytransdecoder.cli predict -t transcripts.fa "${opt[@]}" > predict.log 2>&1
  $PY -m pytransdecoder.cli genome-map --orfs transcripts.fa.transdecoder.gff3 --structures $structs -t transcripts.fa \
     -o transcripts.fa.transdecoder.genome.gff3 > genome_map.log 2>&1
  $PY $P/TransDecoderPaper/genome_annotation/scripts/eval_genome_cds.py --pred transcripts.fa.transdecoder.genome.gff3 \
     --structures $structs --mane $P/data/annot/human/MANE.GRCh38.v1.5.ensembl_genomic.gtf.gz \
     --gencode $D/ref/gencode.v50.primary_assembly.annotation.gtf --label ${asm}_$cfg -o eval.tmp > eval.log 2>&1
  for x in ref_status pred_status summary; do mv eval.tmp.$x.tsv eval.$x.tsv; done
  echo "$(( $(date +%s) - t0 ))" > runtime_sec.txt
  echo "done $asm $cfg $(cat runtime_sec.txt)s"
}
export -f one; export P D PY OUT

for asm in stringtie lraa; do
  echo "$asm a0.4_w25"; echo "$asm norefine"
  for a in ${GM_ALPHAS:-0.2 0.4 0.6 1.0}; do for w in ${GM_WEIGHTS:-0 25 50}; do [ $a/$w = 0.4/25 ] || echo "$asm a${a}_w${w}"; done; done
done | xargs -P $jobs -L 1 bash -c 'one $0 $1'
