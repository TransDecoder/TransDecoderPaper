#!/bin/bash
# STAR alignment of the simulated reads with the index and parameters of the real GM12878 run
# (genome_annotation/scripts/data_prep/star_align_stringtie.sh; E03/E11), including --limitBAMsortRAM 40 GB
# (the reads are grouped by source transcript, so sorting needs ~30 GB; an 8 GB limit failed).
# Only the thread count differs (6 vs 12), which does not change the alignments.
set -euo pipefail
P=${P:-/broad/hptmp/bhaas/TRANSDECODER_PAPER}
export PATH=$P/env/rnaseq_asm/bin:$PATH
THREADS=${THREADS:-6}
F=$P/data/simtrinity/sim/fastq
O=$P/data/simtrinity/star
mkdir -p $O; cd $O
r1=$(ls $F/chunk_*_R1.fastq.gz | sort -V | paste -sd,)
r2=$(ls $F/chunk_*_R2.fastq.gz | sort -V | paste -sd,)
STAR --version
STAR --runThreadN $THREADS --genomeDir $P/data/gm12878/star_stringtie/star_index \
  --readFilesIn $r1 $r2 --readFilesCommand zcat \
  --twopassMode Basic \
  --outSAMtype BAM SortedByCoordinate --limitBAMsortRAM 40000000000 --outBAMsortingThreadN 4 \
  --outSAMstrandField intronMotif --outSAMattributes NH HI AS nM XS \
  --outTmpDir $O/star_tmp \
  --outFileNamePrefix sim_GM12878.
mv sim_GM12878.Aligned.sortedByCoord.out.bam sim_GM12878.STAR.bam
samtools index -@ 4 sim_GM12878.STAR.bam
samtools flagstat -@ 4 sim_GM12878.STAR.bam > sim_GM12878.STAR.flagstat
echo STAR_DONE
