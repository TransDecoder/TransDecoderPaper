#!/bin/bash
# rRNA/tRNA contaminant index for Ribo-seq: human 45S rDNA repeat (GenBank U13369.1) plus GENCODE v50
# genes of type rRNA, rRNA_pseudogene, Mt_rRNA, Mt_tRNA (genomic spans; bowtie2 aligns both strands).
# Recorded from the commands run on 2026-10-21; output: data/gm12878/riboseq_cenik2015/rrna/rrna_trna.*
set -euo pipefail
P=/broad/hptmp/bhaas/TRANSDECODER_PAPER; export PATH=$P/env/riboseq/bin:$PATH
G=$P/data/gm12878/ref; cd $P/data/gm12878/riboseq_cenik2015; mkdir -p rrna; cd rrna
curl -s "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=U13369.1&rettype=fasta&retmode=text" > U13369.fa
awk -F'\t' '$3=="gene" && ($9 ~ /gene_type "(rRNA|rRNA_pseudogene|Mt_rRNA|Mt_tRNA)"/){match($9,/gene_name "[^"]+"/);
  n=substr($9,RSTART+11,RLENGTH-12); print $1"\t"$4-1"\t"$5"\t"n"\t0\t"$7}' $G/gencode.v50.primary_assembly.annotation.gtf > gencode_rrna.bed
awk '{print $1":"$2+1"-"$3}' gencode_rrna.bed | xargs samtools faidx $G/GRCh38.primary_assembly.genome.fa > gencode_rrna.fa
cat U13369.fa gencode_rrna.fa > rrna_trna.fa   # 575 sequences
bowtie2-build --threads 8 -q rrna_trna.fa rrna_trna
