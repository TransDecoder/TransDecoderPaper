#!/bin/bash
# E06 release rerun driver, v2: as driver_v1 but without nice (Brian, 2026-10-30: normal priority).
# Runs the configurations that driver_v1 had not finished when it was stopped (prep for all three
# assemblies and 9 Trinity configurations were done by v1; the 3 Trinity runs in progress were
# killed and moved to $P/_trash/hseries_v7h_niced_killed). 3 at a time.
# Log: $P/work/gm12878_td/logs/hseries_v7h/driver_v2.log
P=/broad/hptmp/bhaas/TRANSDECODER_PAPER; W=$P/work/gm12878_td; T=$P/TransDecoderPaper/genome_annotation/scripts
L=$W/logs/hseries_v7h; mkdir -p $L
echo "start $(date)"
job() { bash $T/homology_series_v7h.sh "$@" > $L/$1.$2.$3.log 2>&1 && echo "ok $* $(date)" || echo "FAILED $*"; }
export -f job; export T L
CLADES="9606 9443 40674 7742 33208"
{ echo "trinity 7742 pfam"; echo "trinity 33208 nopfam"; echo "trinity 33208 pfam"
  for n in stringtie lraa; do
    echo "$n none nopfam"; echo "$n none pfam"
    for c in $CLADES; do echo "$n $c nopfam"; echo "$n $c pfam"; done
  done; } | xargs -P 3 -L 1 bash -c 'job "$@"' _
echo "HSERIES_V7H_DONE $(date)"
