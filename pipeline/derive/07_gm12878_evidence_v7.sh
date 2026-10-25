#!/bin/bash
# Derive compact GM12878 Ribo-seq / MS evidence tables for the release (final_v7) runs (evidence E01, E02, E25).
# Reads the full scoring outputs under $P/work/gm12878_td/final_v7_scoring/{ribo,proteomics}/ and the Ribo-seq
# processing logs; writes small TSVs to data/gm12878/final_v7/ for notebooks/07_gm12878_evidence.Rmd.
# The per-class support tables themselves are already in data/benchmark/final_v7/ (gm12878_v7_*).
# usage: bash pipeline/derive/07_gm12878_evidence_v7.sh
set -euo pipefail
P=/broad/hptmp/bhaas/TRANSDECODER_PAPER; T=$P/TransDecoderPaper
S=$P/work/gm12878_td/final_v7_scoring; R=$P/data/gm12878/riboseq_cenik2015
O=$T/data/gm12878/final_v7; mkdir -p $O
PY=$P/env/td/bin/python

# 1. Ribo-seq read processing counts per run (cutadapt, bowtie2 rRNA/tRNA removal, STAR unique).
{
  printf "run\treads\ttrimmed_25_34\ttoo_short\ttoo_long\tnon_rrna\tunique_mapped\n"
  for s in SRR1803147 SRR1803148; do
    c=$R/logs/$s.cutadapt.log; st=$R/star/$s.Log.final.out
    g() { grep "$1" $c | head -1 | sed -E 's/.*: +([0-9,]+).*/\1/; s/,//g'; }
    printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\n" $s "$(g 'Total reads processed')" "$(g 'Reads written')" \
      "$(g 'too short')" "$(g 'too long')" \
      "$(grep 'Number of input reads' $st | awk -F'|' '{gsub(/[ \t]/,"",$2); print $2}')" \
      "$(grep 'Uniquely mapped reads number' $st | awk -F'|' '{gsub(/[ \t]/,"",$2); print $2}')"
  done
} > $O/gm12878_v7_riboseq_read_processing.tsv

# 2. Footprint coverage by RNA expression of the source transcript (StringTie TPM), release default run.
{
  printf "class\ttpm_bin\tn\tfrac_ge10_psites\tfrac_supported\n"
  awk -f $T/analyses/riboseq/expr_matched_support.awk $P/data/gm12878/star_stringtie/stringtie.gtf \
    $S/ribo/gm12878_v7.stringtie.orf_support.tsv | sort -k1,1 -k2,2
} > $O/gm12878_v7_riboseq_stringtie_by_expression.tsv

# 3. MS peptide counts (filtering steps) and the unannotated ORFs with peptides (one row per ORF).
cp $S/proteomics/ms_peptide_support.peptide_counts.tsv $O/gm12878_v7_ms_peptide_counts.tsv
cp $S/proteomics/ms_peptide_support.novel_with_peptides.tsv $O/gm12878_v7_ms_novel_with_peptides.tsv

# 4. Non-reference peptides on unannotated (novel-class) ORFs only.
awk -F'\t' 'NR==1 || $3=="novel"' $S/proteomics/ms_peptide_support.novel_specific_hits.tsv \
  > $O/gm12878_v7_ms_novel_specific_hits_unannotated.tsv

# 5. Location of unannotated ORFs with peptides relative to the parent gene, and pseudogene overlap.
preds=()
for a in lraa lraa_precise lraa_sensitive stringtie stringtie_precise stringtie_sensitive trinity; do
  asm=${a%%_*}; sfx=${a#$asm}
  preds+=("$a=$P/work/gm12878_td/$asm/final_v7_abinitio$sfx/transcripts.fa.transdecoder.genome.gff3")
done
$PY -I $T/analyses/proteomics/novel_peptide_loci.py $S/proteomics/ms_peptide_support.novel_with_peptides.tsv \
  $P/data/gm12878/ref/gencode.v50.primary_assembly.annotation.gtf "${preds[@]}" \
  -o $O/gm12878_v7_ms_novel_peptide_loci
ls -la $O
