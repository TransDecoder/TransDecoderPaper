#!/bin/bash
# Route B: Trinity de novo (stranded RF, default in-silico normalization) on ENCSR000AED rep1, then minimap2 to genome.
# Gated on StringTie finishing (LRAA uses 4 cores alongside) so concurrent CPU stays within ~16 cores.
set -euo pipefail
unset PERL5LIB PERL_LOCAL_LIB_ROOT PERL_MB_OPT PERL_MM_OPT
export PATH=/broad/hptmp/bhaas/TRANSDECODER_PAPER/env/trinity215/bin:/usr/bin:/bin
cd /broad/hptmp/bhaas/TRANSDECODER_PAPER/data/gm12878/trinity
until grep -q STRINGTIE_DONE /broad/hptmp/bhaas/TRANSDECODER_PAPER/data/gm12878/logs/star_align_stringtie.log; do sleep 120; done
echo "[$(date)] starting Trinity"
Trinity --seqType fq --max_memory 100G --CPU ${CPU:-12} --SS_lib_type RF \
  --left /broad/hptmp/bhaas/TRANSDECODER_PAPER/data/gm12878/reads/short/ENCFF001REK.fastq.gz --right /broad/hptmp/bhaas/TRANSDECODER_PAPER/data/gm12878/reads/short/ENCFF001REJ.fastq.gz \
  --output /broad/hptmp/bhaas/TRANSDECODER_PAPER/data/gm12878/trinity/trinity_out_dir
F=/broad/hptmp/bhaas/TRANSDECODER_PAPER/data/gm12878/trinity/trinity_out_dir.Trinity.fasta; [ -s $F ] || F=/broad/hptmp/bhaas/TRANSDECODER_PAPER/data/gm12878/trinity/trinity_out_dir/Trinity.fasta
ln -sf $F /broad/hptmp/bhaas/TRANSDECODER_PAPER/data/gm12878/trinity/Trinity.fasta
echo "[$(date)] Trinity done; minimap2"
/broad/hptmp/bhaas/TRANSDECODER_PAPER/env/rnaseq_asm/bin/minimap2 -t ${CPU:-12} -ax splice:hq -uf --secondary=no /broad/hptmp/bhaas/TRANSDECODER_PAPER/data/gm12878/ref/GRCh38.primary_assembly.genome.fa /broad/hptmp/bhaas/TRANSDECODER_PAPER/data/gm12878/trinity/Trinity.fasta 2> /broad/hptmp/bhaas/TRANSDECODER_PAPER/data/gm12878/trinity/minimap2.log \
  | /broad/hptmp/bhaas/TRANSDECODER_PAPER/env/rnaseq_asm/bin/samtools sort -@ 4 -m 2G -T /broad/hptmp/bhaas/TRANSDECODER_PAPER/data/gm12878/trinity/sorttmp -o /broad/hptmp/bhaas/TRANSDECODER_PAPER/data/gm12878/trinity/Trinity.minimap2.bam -
/broad/hptmp/bhaas/TRANSDECODER_PAPER/env/rnaseq_asm/bin/samtools index /broad/hptmp/bhaas/TRANSDECODER_PAPER/data/gm12878/trinity/Trinity.minimap2.bam
/broad/hptmp/bhaas/TRANSDECODER_PAPER/env/rnaseq_asm/bin/samtools flagstat /broad/hptmp/bhaas/TRANSDECODER_PAPER/data/gm12878/trinity/Trinity.minimap2.bam > /broad/hptmp/bhaas/TRANSDECODER_PAPER/data/gm12878/trinity/Trinity.minimap2.flagstat
echo TRINITY_ROUTE_DONE
