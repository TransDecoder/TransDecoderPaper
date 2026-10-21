#!/bin/bash
# Route C: LRAA v0.44.3 de novo (no --gtf), HiFi mode, on minimap2 splice:hq alignments of ENCFF694DIE
set -euo pipefail
P=/broad/hptmp/bhaas/TRANSDECODER_PAPER/data/gm12878/longread_lraa
until grep -q MINIMAP2_LONG_DONE /broad/hptmp/bhaas/TRANSDECODER_PAPER/data/gm12878/logs/longread_minimap2.log; do sleep 60; done
export PATH=/broad/hptmp/bhaas/tools/apptainer-1.5.3/bin:/broad/hptmp/bhaas/TRANSDECODER_PAPER/env/rnaseq_asm/bin:$PATH
export APPTAINER_CACHEDIR=$P/apptainer_cache APPTAINER_TMPDIR=$P/tmp
mkdir -p $P/lraa_denovo $P/tmp/lraa
cd $P/lraa_denovo
apptainer exec --containall --no-mount hostfs -B /broad/hptmp/bhaas/TRANSDECODER_PAPER/data/gm12878:/broad/hptmp/bhaas/TRANSDECODER_PAPER/data/gm12878 --pwd $P/lraa_denovo \
  --env TMPDIR=$P/tmp/lraa $P/image/lraa-core_0.44.3.sif \
  LRAA --genome /broad/hptmp/bhaas/TRANSDECODER_PAPER/data/gm12878/ref/GRCh38.primary_assembly.genome.fa --bam $P/ENCFF694DIE.minimap2.bam \
       --HiFi --cpu_budget ${CPU:-4} -o GM12878_ENCSR838WFC_rep2.LRAA
gffread -w GM12878_ENCSR838WFC_rep2.LRAA.transcripts.fa -g /broad/hptmp/bhaas/TRANSDECODER_PAPER/data/gm12878/ref/GRCh38.primary_assembly.genome.fa GM12878_ENCSR838WFC_rep2.LRAA.gtf
echo LRAA_DONE
