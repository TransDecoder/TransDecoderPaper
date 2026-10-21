#!/bin/bash
# Route B (revised): Trinity genome-guided on the STAR alignments of ENCSR000AED rep1 (stranded RF),
# then minimap2 of the assembled transcripts to the genome. Replaces the de novo run, which was
# projected to take weeks on this 8-core machine.
set -euo pipefail
unset PERL5LIB PERL_LOCAL_LIB_ROOT PERL_MB_OPT PERL_MM_OPT
export PATH=/broad/hptmp/bhaas/TRANSDECODER_PAPER/env/trinity215/bin:/usr/bin:/bin
D=/broad/hptmp/bhaas/TRANSDECODER_PAPER/data/gm12878
E=/broad/hptmp/bhaas/TRANSDECODER_PAPER/env/rnaseq_asm/bin
cd $D/trinity
echo "[$(date)] starting Trinity genome-guided"
Trinity --genome_guided_bam $D/star_stringtie/GM12878_ENCSR000AED_rep1.STAR.bam \
  --genome_guided_max_intron 100000 --SS_lib_type RF \
  --max_memory 50G --CPU ${CPU:-8} --output $D/trinity/trinity_gg_out
F=$D/trinity/trinity_gg_out/Trinity-GG.fasta
ln -sf $F $D/trinity/Trinity.fasta
echo "[$(date)] Trinity done; minimap2"
$E/minimap2 -t ${CPU:-8} -ax splice:hq -uf --secondary=no $D/ref/GRCh38.primary_assembly.genome.fa $D/trinity/Trinity.fasta 2> $D/trinity/minimap2.log \
  | $E/samtools sort -@ 4 -m 2G -T $D/trinity/sorttmp -o $D/trinity/Trinity.minimap2.bam -
$E/samtools index $D/trinity/Trinity.minimap2.bam
$E/samtools flagstat $D/trinity/Trinity.minimap2.bam > $D/trinity/Trinity.minimap2.flagstat
echo TRINITY_ROUTE_DONE
