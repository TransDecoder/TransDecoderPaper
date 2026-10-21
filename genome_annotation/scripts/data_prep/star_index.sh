#!/bin/bash
set -euo pipefail
export PATH=/broad/hptmp/bhaas/TRANSDECODER_PAPER/env/rnaseq_asm/bin:$PATH
cd /broad/hptmp/bhaas/TRANSDECODER_PAPER/data/gm12878/star_stringtie
STAR --runMode genomeGenerate --runThreadN 12 --genomeDir star_index \
  --genomeFastaFiles /broad/hptmp/bhaas/TRANSDECODER_PAPER/data/gm12878/ref/GRCh38.primary_assembly.genome.fa \
  --sjdbGTFfile /broad/hptmp/bhaas/TRANSDECODER_PAPER/data/gm12878/ref/gencode.v50.primary_assembly.annotation.gtf --sjdbOverhang 100 \
  --outTmpDir /broad/hptmp/bhaas/TRANSDECODER_PAPER/data/gm12878/star_stringtie/star_index_tmp
echo STAR_INDEX_DONE
