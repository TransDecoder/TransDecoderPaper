#!/bin/bash
set -euo pipefail
export PATH=/broad/hptmp/bhaas/TRANSDECODER_PAPER/env/rnaseq_asm/bin:$PATH
cd /broad/hptmp/bhaas/TRANSDECODER_PAPER/data/gm12878/longread_lraa
minimap2 -t 4 -ax splice:hq -uf --secondary=no \
  /broad/hptmp/bhaas/TRANSDECODER_PAPER/data/gm12878/ref/GRCh38.primary_assembly.genome.fa /broad/hptmp/bhaas/TRANSDECODER_PAPER/data/gm12878/reads/long/ENCFF694DIE.fastq.gz 2> minimap2.log \
 | samtools sort -@ 2 -m 2G -T /broad/hptmp/bhaas/TRANSDECODER_PAPER/data/gm12878/longread_lraa/tmp/sort -o ENCFF694DIE.minimap2.bam -
samtools index ENCFF694DIE.minimap2.bam
samtools flagstat ENCFF694DIE.minimap2.bam > ENCFF694DIE.minimap2.flagstat
echo MINIMAP2_LONG_DONE
