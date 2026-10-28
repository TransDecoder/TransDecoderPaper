#!/bin/bash
# E32: genome-guided Trinity 2.15.2 + minimap2 for one Arabidopsis level (d1 = all simulated reads; dN = the
# subsample from 11_depth_subsample_athal.sh), with the E30/E31 parameters: --genome_guided_max_intron 100000
# --SS_lib_type RF --min_kmer_cov 2 --max_memory 50G; CPU from $CPU (default 4).
# usage: 08_trinity_gg_athal.sh <level>
set -euo pipefail
unset PERL5LIB PERL_LOCAL_LIB_ROOT PERL_MB_OPT PERL_MM_OPT
level=$1
P=/broad/hptmp/bhaas/TRANSDECODER_PAPER
export PATH=$P/env/trinity215/bin:/usr/bin:/bin
E=$P/env/rnaseq_asm/bin
D=$P/data/simtrinity_athal/depth/$level
G=$P/data/simtrinity_athal/ref/TAIR12.genome.fa
CPU=${CPU:-4}
mkdir -p $D; cd $D
[ $level = d1 ] && ln -sf $P/data/simtrinity_athal/star/sim_athal.STAR.bam $D/star.sub.bam && ln -sf $P/data/simtrinity_athal/star/sim_athal.STAR.bam.bai $D/star.sub.bam.bai
echo "[$(date)] $level: starting Trinity genome-guided (CPU $CPU)"
Trinity --genome_guided_bam $D/star.sub.bam \
  --genome_guided_max_intron 100000 --SS_lib_type RF --min_kmer_cov 2 \
  --max_memory 50G --CPU $CPU --output $D/trinity_gg_out
ln -sf $D/trinity_gg_out/Trinity-GG.fasta $D/Trinity.fasta
echo "[$(date)] $level: Trinity done; minimap2"
$E/minimap2 -t $CPU -ax splice:hq -uf --secondary=no $G $D/Trinity.fasta 2> $D/minimap2.log \
  | $E/samtools sort -@ 2 -m 2G -T $D/sorttmp -o $D/Trinity.minimap2.bam -
$E/samtools index $D/Trinity.minimap2.bam
$E/samtools flagstat $D/Trinity.minimap2.bam > $D/Trinity.minimap2.flagstat
echo "[$(date)] $level: TRINITY_ROUTE_DONE"
