#!/bin/bash
# Route A: STAR (2-pass basic) alignment of ENCSR000AED rep1 -> StringTie de novo (no -G), stranded --rf
set -euo pipefail
export PATH=/broad/hptmp/bhaas/TRANSDECODER_PAPER/env/rnaseq_asm/bin:$PATH
THREADS=${THREADS:-12}
cd /broad/hptmp/bhaas/TRANSDECODER_PAPER/data/gm12878/star_stringtie
until grep -q STAR_INDEX_DONE /broad/hptmp/bhaas/TRANSDECODER_PAPER/data/gm12878/logs/star_index.log; do sleep 60; done
rm -rf star_align_tmp
STAR --runThreadN $THREADS --genomeDir star_index \
  --readFilesIn /broad/hptmp/bhaas/TRANSDECODER_PAPER/data/gm12878/reads/short/ENCFF001REK.fastq.gz /broad/hptmp/bhaas/TRANSDECODER_PAPER/data/gm12878/reads/short/ENCFF001REJ.fastq.gz --readFilesCommand zcat \
  --twopassMode Basic \
  --outSAMtype BAM SortedByCoordinate --limitBAMsortRAM 40000000000 --outBAMsortingThreadN 6 \
  --outSAMstrandField intronMotif --outSAMattributes NH HI AS nM XS \
  --outTmpDir /broad/hptmp/bhaas/TRANSDECODER_PAPER/data/gm12878/star_stringtie/star_align_tmp \
  --outFileNamePrefix GM12878_ENCSR000AED_rep1.
mv GM12878_ENCSR000AED_rep1.Aligned.sortedByCoord.out.bam GM12878_ENCSR000AED_rep1.STAR.bam
samtools index -@ 4 GM12878_ENCSR000AED_rep1.STAR.bam
echo STAR_ALIGN_DONE
touch STAR_ALIGN_DONE
stringtie GM12878_ENCSR000AED_rep1.STAR.bam --rf -p 4 -l STRG -o stringtie.gtf
gffread -w stringtie.transcripts.fa -g /broad/hptmp/bhaas/TRANSDECODER_PAPER/data/gm12878/ref/GRCh38.primary_assembly.genome.fa stringtie.gtf
samtools flagstat -@ 4 GM12878_ENCSR000AED_rep1.STAR.bam > GM12878_ENCSR000AED_rep1.STAR.flagstat
echo STRINGTIE_DONE
