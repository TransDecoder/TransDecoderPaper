#!/bin/bash
D=/broad/hptmp/bhaas/TRANSDECODER_PAPER/data/gm12878; G=/broad/hptmp/bhaas/TRANSDECODER_PAPER/TransDecoderPaper/genome_annotation/scripts/gm12878_final_v7.sh; W=/broad/hptmp/bhaas/TRANSDECODER_PAPER/work/gm12878_td
run() { bash $G "$@" > $W/logs/final_v7_$1_$4_${6:-default}.log 2>&1 && echo "ok $1 $4 ${6:-default}" || echo "FAILED $1 $4 ${6:-default}"; }
run lraa $D/longread_lraa/lraa_denovo/LRAA.transcripts.fa $D/longread_lraa/lraa_denovo/LRAA.gtf abinitio 4 &
run stringtie $D/star_stringtie/stringtie.transcripts.fa $D/star_stringtie/stringtie.gtf abinitio 4 &
run trinity $D/trinity/Trinity.fasta $D/trinity/Trinity.minimap2.bam abinitio 4 &
wait
for q in precise sensitive; do
  run lraa $D/longread_lraa/lraa_denovo/LRAA.transcripts.fa $D/longread_lraa/lraa_denovo/LRAA.gtf abinitio 4 $q &
  run stringtie $D/star_stringtie/stringtie.transcripts.fa $D/star_stringtie/stringtie.gtf abinitio 4 $q &
done
wait
run trinity $D/trinity/Trinity.fasta $D/trinity/Trinity.minimap2.bam homology 8
echo ALL_DONE
