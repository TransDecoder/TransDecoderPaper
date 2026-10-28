#!/bin/bash
# E32: subsample the Arabidopsis STAR BAM by read pair (as E31's 11_depth_subsample.sh; same seed).
# usage: 11_depth_subsample_athal.sh <level> <fraction>
set -euo pipefail
level=$1; frac=$2
P=/broad/hptmp/bhaas/TRANSDECODER_PAPER
S=$P/env/rnaseq_asm/bin/samtools
SEED=20261026
IN=$P/data/simtrinity_athal/star/sim_athal.STAR.bam
O=$P/data/simtrinity_athal/depth/$level
mkdir -p $O
$S view -@ 2 -b --subsample $frac --subsample-seed $SEED -o $O/star.sub.bam $IN
$S index $O/star.sub.bam
full=$($S view -c -f 64 -F 0x904 $IN)
sub=$($S view -c -f 64 -F 0x904 $O/star.sub.bam)
printf "level\tfraction\tseed\tpairs_full\tpairs_sub\tobserved_fraction\n%s\t%s\t%s\t%s\t%s\t%s\n" \
  $level $frac $SEED $full $sub $(awk -v a=$sub -v b=$full 'BEGIN{printf "%.5f", a/b}') > $O/pair_counts.tsv
cat $O/pair_counts.tsv
