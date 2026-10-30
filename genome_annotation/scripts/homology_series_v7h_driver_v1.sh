#!/bin/bash
# E06 release rerun driver: prep (LongOrfs + DIAMOND) per assembly, then 12 configurations per assembly
# (ab initio, BLAST x 5 clades, Pfam only, BLAST + Pfam x 5 clades), 3 at a time, niced.
# Log: $P/work/gm12878_td/logs/hseries_v7h/driver_v1.log
P=/broad/hptmp/bhaas/TRANSDECODER_PAPER; W=$P/work/gm12878_td; T=$P/TransDecoderPaper/genome_annotation/scripts
L=$W/logs/hseries_v7h; mkdir -p $L
echo "start $(date)"
for n in trinity stringtie lraa; do
  nice -n 10 bash $T/homology_series_v7h_prep.sh $n 4 > $L/prep.$n.log 2>&1 && echo "ok prep $n $(date)" || echo "FAILED prep $n"
done
job() { nice -n 10 bash $T/homology_series_v7h.sh "$@" > $L/$1.$2.$3.log 2>&1 && echo "ok $* $(date)" || echo "FAILED $*"; }
export -f job; export T L
CLADES="9606 9443 40674 7742 33208"
for n in trinity stringtie lraa; do
  echo "$n none nopfam"; echo "$n none pfam"
  for c in $CLADES; do echo "$n $c nopfam"; echo "$n $c pfam"; done
done | xargs -P 3 -L 1 bash -c 'job "$@"' _
echo "HSERIES_V7H_DONE $(date)"
