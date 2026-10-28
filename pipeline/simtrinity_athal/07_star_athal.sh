#!/bin/bash
# E32: STAR 2.7.11b two-pass alignment of the simulated Arabidopsis reads, with the E30 parameters
# (07_star.sh) and the TAIR12 index (03_star_index_athal.sh).
set -euo pipefail
P=/broad/hptmp/bhaas/TRANSDECODER_PAPER
export PATH=$P/env/rnaseq_asm/bin:$PATH
THREADS=${THREADS:-6}
F=$P/data/simtrinity_athal/sim/fastq
O=$P/data/simtrinity_athal/star
mkdir -p $O; cd $O
r1=$(ls $F/chunk_*_R1.fastq.gz | sort -V | paste -sd,)
r2=$(ls $F/chunk_*_R2.fastq.gz | sort -V | paste -sd,)
STAR --version
STAR --runThreadN $THREADS --genomeDir $P/data/simtrinity_athal/ref/star_index \
  --readFilesIn $r1 $r2 --readFilesCommand zcat \
  --twopassMode Basic \
  --outSAMtype BAM SortedByCoordinate --limitBAMsortRAM 40000000000 --outBAMsortingThreadN 4 \
  --outSAMstrandField intronMotif --outSAMattributes NH HI AS nM XS \
  --outTmpDir $O/star_tmp \
  --outFileNamePrefix sim_athal.
mv sim_athal.Aligned.sortedByCoord.out.bam sim_athal.STAR.bam
samtools index -@ 4 sim_athal.STAR.bam
samtools flagstat -@ 4 sim_athal.STAR.bam > sim_athal.STAR.flagstat
echo STAR_DONE
