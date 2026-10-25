#!/bin/bash
# E26 frame-contrast experiment on GM12878; starts after the benchmark q-curve runs (frames/launch_q.sh) finish.
D=/broad/hptmp/bhaas/TRANSDECODER_PAPER/data/gm12878; G=/broad/hptmp/bhaas/TRANSDECODER_PAPER/TransDecoderPaper/genome_annotation/scripts/gm12878_frames.sh; W=/broad/hptmp/bhaas/TRANSDECODER_PAPER/work/gm12878_td
until grep -q ALL_DONE /broad/hptmp/bhaas/TRANSDECODER_PAPER/work/bench_v3/frames/logs/launch_q.out; do sleep 300; done
run() { bash $G "$@" > $W/logs/frames_$1_$4_${6:-default}.log 2>&1 && echo "ok $1 $4 ${6:-default}" || echo "FAILED $1 $4 ${6:-default}"; }
run lraa $D/longread_lraa/lraa_denovo/LRAA.transcripts.fa $D/longread_lraa/lraa_denovo/LRAA.gtf all 4 &
run stringtie $D/star_stringtie/stringtie.transcripts.fa $D/star_stringtie/stringtie.gtf all 4 &
run trinity $D/trinity/Trinity.fasta $D/trinity/Trinity.minimap2.bam all 4 &
wait
run lraa $D/longread_lraa/lraa_denovo/LRAA.transcripts.fa $D/longread_lraa/lraa_denovo/LRAA.gtf all 4 sensitive &
run stringtie $D/star_stringtie/stringtie.transcripts.fa $D/star_stringtie/stringtie.gtf all 4 sensitive &
wait
echo ALL_DONE
