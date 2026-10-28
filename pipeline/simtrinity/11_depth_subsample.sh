#!/bin/bash
# E31 read-depth ablation: subsample the E30 STAR BAM by read pair.
# usage: 11_depth_subsample.sh <level> <fraction>      e.g. 11_depth_subsample.sh d4 0.25
# samtools view --subsample hashes the read name with the seed, so mates stay together and, with one
# fixed seed, lower fractions are nested in higher ones. STAR is not re-run (its 2-pass junctions come
# from the full read set). Writes $P/data/simtrinity/depth/<level>/star.sub.bam(.bai) and pair_counts.tsv.
set -euo pipefail
level=$1; frac=$2
P=${P:-/broad/hptmp/bhaas/TRANSDECODER_PAPER}
S=$P/env/rnaseq_asm/bin/samtools
SEED=20261026
IN=$P/data/simtrinity/star/sim_GM12878.STAR.bam
O=$P/data/simtrinity/depth/$level
mkdir -p $O
$S view -@ 2 -b --subsample $frac --subsample-seed $SEED -o $O/star.sub.bam $IN
$S index $O/star.sub.bam
# read pairs = distinct templates; count primary first-mate records (STAR writes mapped reads only)
full=$($S view -c -f 64 -F 0x904 $IN)
sub=$($S view -c -f 64 -F 0x904 $O/star.sub.bam)
printf "level\tfraction\tseed\tpairs_full\tpairs_sub\tobserved_fraction\n%s\t%s\t%s\t%s\t%s\t%s\n" \
  $level $frac $SEED $full $sub $(awk -v a=$sub -v b=$full 'BEGIN{printf "%.5f", a/b}') > $O/pair_counts.tsv
cat $O/pair_counts.tsv
