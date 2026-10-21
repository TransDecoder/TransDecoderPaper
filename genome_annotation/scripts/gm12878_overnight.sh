#!/bin/bash
# 2026-10-20 overnight: LRAA + StringTie with final defaults now; Trinity (genome-guided) once its route finishes.
P=/broad/hptmp/bhaas/TRANSDECODER_PAPER; D=$P/data/gm12878; W=$P/work/gm12878_td
run() { bash $W/final_defaults.sh "$@" > $W/logs/final_$1_$4.log 2>&1 && echo "ok $1 $4" || echo "FAILED $1 $4"; }
mkdir -p $W/logs
run lraa $D/longread_lraa/lraa_denovo/LRAA.transcripts.fa $D/longread_lraa/lraa_denovo/LRAA.gtf abinitio &
run stringtie $D/star_stringtie/stringtie.transcripts.fa $D/star_stringtie/stringtie.gtf abinitio &
wait
until grep -q TRINITY_ROUTE_DONE $D/logs/trinity_gg.log; do
  pgrep -f 'bin/trinity_gg.sh' >/dev/null || { echo "FAILED trinity route (no TRINITY_ROUTE_DONE)"; echo OVERNIGHT_DONE; exit 1; }
  sleep 300
done
echo "trinity route done $(date)"
run trinity $D/trinity/Trinity.fasta $D/trinity/Trinity.minimap2.bam abinitio
run trinity $D/trinity/Trinity.fasta $D/trinity/Trinity.minimap2.bam homology 16
echo OVERNIGHT_DONE
