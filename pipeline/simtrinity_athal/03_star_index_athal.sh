#!/bin/bash
# E32: STAR 2.7.11b index of TAIR12 (Chr1-5, ChrM, ChrC) with the TAIR12 exon GTF as junction database,
# overhang 100 (as the GM12878 index); --genomeSAindexNbases 12 for the 120 Mb genome (STAR manual).
set -euo pipefail
P=/broad/hptmp/bhaas/TRANSDECODER_PAPER
export PATH=$P/env/rnaseq_asm/bin:$PATH
R=$P/data/simtrinity_athal/ref
mkdir -p $R/star_index
[ -s $R/TAIR12.genome.fa ] || zcat $P/data/annot/athal/Athaliana_TAIR12_genomic_Chr_softmasked.fa.gz | awk '/^>/{print $1; next}{print toupper($0)}' > $R/TAIR12.genome.fa
samtools faidx $R/TAIR12.genome.fa
STAR --runMode genomeGenerate --runThreadN ${THREADS:-4} --genomeDir $R/star_index \
  --genomeFastaFiles $R/TAIR12.genome.fa --sjdbGTFfile $R/tair12.exons.gtf --sjdbOverhang 100 \
  --genomeSAindexNbases 12 --outTmpDir $R/star_index_tmp --outFileNamePrefix $R/star_index/
echo STAR_INDEX_DONE
