#!/bin/bash
# Derive compact GM12878 tables for notebooks/06_gm12878_application.Rmd (evidence E03, E05, E25).
# Inputs (not in git): $P/work/gm12878_td/eval/*final_v7_abinitio.{ref,pred}_status.tsv (release runs),
# $P/work/gm12878_td/eval/status/*final_abinitio.* (earlier runs), each run's genome_map.log and genome GFF3,
# and the StringTie -e quantification of GENCODE v50 ($P/data/gm12878/expr/, record E05).
# Outputs: data/gm12878/final_v7/06_*.tsv
#   06_gm12878_inputs.tsv  read and transcript counts of the upstream data and assemblies (E03, E11), from the
#       STAR, minimap2 and read-count logs and the assembly GTF/FASTA files under $P/data/gm12878/
#   06_gm12878_v7_expressed_sensitivity.tsv  expressed_sensitivity.py (E05 method) on the release runs
#   06_gm12878_class_by_length.tsv, 06_gm12878_genome_map.tsv, 06_gm12878_mane_not_exact.tsv
#       (06_gm12878_tables.py; earlier and release runs)
# The release eval summaries are already in data/benchmark/final_v7/eval.*_final_v7_abinitio.summary.tsv
# (identical to $P/work/gm12878_td/eval/), and the Ribo-seq by expression table in
# data/gm12878/final_v7/gm12878_v7_riboseq_stringtie_by_expression.tsv (pipeline/derive/07_gm12878_evidence_v7.sh).
# usage: bash pipeline/derive/06_gm12878_application.sh
set -euo pipefail
P=/broad/hptmp/bhaas/TRANSDECODER_PAPER; T=$P/TransDecoderPaper
W=$P/work/gm12878_td; E=$W/eval
O=$T/data/gm12878/final_v7; mkdir -p $O
PY=$P/env/td/bin/python
D=$P/data/gm12878

{
  printf "metric\tvalue\tsource\n"
  printf "short_read_pairs\t%s\tSTAR Log.final.out\n" \
    "$(grep 'Number of input reads' $D/star_stringtie/GM12878_ENCSR000AED_rep1.Log.final.out | cut -f2)"
  printf "short_read_unique_pct\t%s\tSTAR Log.final.out\n" \
    "$(grep 'Uniquely mapped reads %' $D/star_stringtie/GM12878_ENCSR000AED_rep1.Log.final.out | cut -f2 | tr -d %)"
  printf "long_reads_flnc\t%s\tlogs/readcount_long.txt\n" "$(awk '{print $NF}' $D/logs/readcount_long.txt)"
  printf "stringtie_transcripts\t%s\tstringtie.gtf\n" "$(awk '$3=="transcript"' $D/star_stringtie/stringtie.gtf | wc -l)"
  printf "lraa_transcripts\t%s\tLRAA.gtf\n" "$(awk '$3=="transcript"' $D/longread_lraa/lraa_denovo/LRAA.gtf | wc -l)"
  printf "trinity_transcripts\t%s\tTrinity.fasta\n" "$(grep -c '^>' $D/trinity/Trinity.fasta)"
  printf "trinity_primary_mapped\t%s\tTrinity.minimap2.flagstat\n" \
    "$(awk '/primary mapped/{print $1}' $D/trinity/Trinity.minimap2.flagstat)"
} > $O/06_gm12878_inputs.tsv

$PY -I $T/genome_annotation/scripts/expressed_sensitivity.py $P/data/gm12878/expr/gencode_v50.ref_quant.gtf \
  lraa=$E/lraa_final_v7_abinitio.ref_status.tsv stringtie=$E/stringtie_final_v7_abinitio.ref_status.tsv \
  trinity=$E/trinity_final_v7_abinitio.ref_status.tsv > $O/06_gm12878_v7_expressed_sensitivity.tsv

PYTHONPATH=$P/work/trees/PY_v7 $PY $T/pipeline/derive/06_gm12878_tables.py $W \
  $P/data/annot/human/MANE.GRCh38.v1.5.ensembl_genomic.gtf.gz $O/06_
ls -la $O/06_*
