#!/bin/bash
# Genome-guided Trinity 2.15.2 on the STAR alignments of the simulated reads, then minimap2 of the
# assembled transcripts to GRCh38 (as genome_annotation/scripts/data_prep/trinity_gg.sh, E11; this run
# adds --min_kmer_cov 2 as specified for E30, CPU 8). Trinity's ParaFly is the xargs drop-in (E11).
set -euo pipefail
unset PERL5LIB PERL_LOCAL_LIB_ROOT PERL_MB_OPT PERL_MM_OPT
P=/broad/hptmp/bhaas/TRANSDECODER_PAPER
export PATH=$P/env/trinity215/bin:/usr/bin:/bin
E=$P/env/rnaseq_asm/bin
D=$P/data/simtrinity
G=$P/data/gm12878/ref/GRCh38.primary_assembly.genome.fa
mkdir -p $D/trinity; cd $D/trinity
echo "[$(date)] starting Trinity genome-guided"
Trinity --version || true
Trinity --genome_guided_bam $D/star/sim_GM12878.STAR.bam \
  --genome_guided_max_intron 100000 --SS_lib_type RF --min_kmer_cov 2 \
  --max_memory 50G --CPU 8 --output $D/trinity/trinity_gg_out
ln -sf $D/trinity/trinity_gg_out/Trinity-GG.fasta $D/trinity/Trinity.fasta
echo "[$(date)] Trinity done; minimap2"
$E/minimap2 --version
$E/minimap2 -t 4 -ax splice:hq -uf --secondary=no $G $D/trinity/Trinity.fasta 2> $D/trinity/minimap2.log \
  | $E/samtools sort -@ 2 -m 2G -T $D/trinity/sorttmp -o $D/trinity/Trinity.minimap2.bam -
$E/samtools index $D/trinity/Trinity.minimap2.bam
$E/samtools flagstat $D/trinity/Trinity.minimap2.bam > $D/trinity/Trinity.minimap2.flagstat
echo "[$(date)] TRINITY_ROUTE_DONE"
