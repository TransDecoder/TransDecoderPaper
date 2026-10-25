#!/bin/bash
# Fresh PyTransDecoder runs of the experimental UTR-length prior (evidence E29).
# PyTD tree: $P/work/trees/PY_utr (branch exp/utr-prior, local only; 3a793ff + --start-utr-prior-weight).
#
# usage: 03_utr_runs.sh bench <org> <variant> <lambda>     -> $P/work/utr_prior/runs/<org>/<variant>/utr_l<lambda>/
#        03_utr_runs.sh gm <lraa|stringtie> <lambda>        -> $P/work/utr_prior/gm12878_runs/<asm>/utr_l<lambda>/
# Benchmark inputs: stranded (-S), genetic code from organisms.tsv, preds.tsv via convert_predictions.py.
# GM12878: as genome_annotation/scripts/gm12878_final_v7.sh (-S, genome-map, eval_genome_cds.py), eval
# written to $P/work/utr_prior/gm12878_eval/<asm>_utr_l<lambda>.*; no Ribo-seq step.
# An existing output directory is refused (fresh runs only). lambda 0 = release behaviour.
set -euo pipefail
P=/broad/hptmp/bhaas/TRANSDECODER_PAPER; R=$P/TransDecoderPaper; D=$P/data/gm12878
PY=$P/env/td/bin/python; export PATH=$P/env/td/bin:$PATH PYTHONPATH=$P/work/trees/PY_utr
td() { $PY -m pytransdecoder.cli "$@"; }
mode=$1; shift
if [ $mode = bench ]; then
  org=$1; v=$2; lam=$3
  gcode=$(awk -F'\t' -v o=$org '$1==o{print $4}' $R/benchmark_v2/organisms.tsv)
  out=$P/work/utr_prior/runs/$org/$v/utr_l$lam
  [ -e $out ] && { echo "exists: $out"; exit 1; }
  mkdir -p $out; cd $out
  (cd $PYTHONPATH && git rev-parse --short HEAD) > pytd_version.txt
  ln -s $P/data/bench_v3/inputs/$org/$v.fa input.fa
  t0=$(date +%s)
  td longorfs -t input.fa -G $gcode -S > longorfs.log 2>&1
  td predict -t input.fa -G $gcode --start-utr-prior-weight $lam > predict.log 2>&1
  echo $(( $(date +%s) - t0 )) > runtime_sec.txt
  $PY $R/benchmark_v2/scripts/convert_predictions.py transdecoder_gff3 input.fa.transdecoder.gff3 > preds.tsv
  echo "DONE $org $v $lam"
elif [ $mode = gm ]; then
  asm=$1; lam=$2
  case $asm in
    lraa) fa=$D/longread_lraa/lraa_denovo/LRAA.transcripts.fa; st=$D/longread_lraa/lraa_denovo/LRAA.gtf ;;
    stringtie) fa=$D/star_stringtie/stringtie.transcripts.fa; st=$D/star_stringtie/stringtie.gtf ;;
    *) echo "unknown assembly $asm"; exit 1 ;;
  esac
  out=$P/work/utr_prior/gm12878_runs/$asm/utr_l$lam
  [ -e $out ] && { echo "exists: $out"; exit 1; }
  mkdir -p $out $P/work/utr_prior/gm12878_eval; cd $out
  (cd $PYTHONPATH && git rev-parse --short HEAD) > pytd_version.txt
  ln -s $fa transcripts.fa
  td longorfs -t transcripts.fa -S > longorfs.log 2>&1
  td predict -t transcripts.fa --start-utr-prior-weight $lam > predict.log 2>&1
  td genome-map --orfs transcripts.fa.transdecoder.gff3 --structures $st -t transcripts.fa \
     -o transcripts.fa.transdecoder.genome.gff3 > genome_map.log 2>&1
  $PY $R/genome_annotation/scripts/eval_genome_cds.py --pred transcripts.fa.transdecoder.genome.gff3 \
     --structures $st --mane $P/data/annot/human/MANE.GRCh38.v1.5.ensembl_genomic.gtf.gz \
     --gencode $D/ref/gencode.v50.primary_assembly.annotation.gtf --label ${asm}_utr_l$lam \
     -o $P/work/utr_prior/gm12878_eval/${asm}_utr_l$lam > eval.log 2>&1
  echo "DONE gm $asm $lam"
fi
