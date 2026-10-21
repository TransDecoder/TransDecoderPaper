#!/bin/bash
# LRAA names final outputs <prefix>.LRAA.ref-free.*; lraa_denovo.sh's gffread line used the wrong name.
# This waits for lraa_denovo.sh (pid $1) to exit, then extracts transcript FASTA and writes the LRAA_DONE marker.
set -euo pipefail
while kill -0 $1 2>/dev/null; do sleep 60; done
cd /broad/hptmp/bhaas/TRANSDECODER_PAPER/data/gm12878/longread_lraa/lraa_denovo
G=$(ls GM12878_ENCSR838WFC_rep2.LRAA*.gtf | grep -v -E 'init|tmp' | head -1)
echo "final LRAA gtf: $G"
/broad/hptmp/bhaas/TRANSDECODER_PAPER/env/rnaseq_asm/bin/gffread -w ${G%.gtf}.transcripts.fa -g /broad/hptmp/bhaas/TRANSDECODER_PAPER/data/gm12878/ref/GRCh38.primary_assembly.genome.fa $G
ln -sf $G LRAA.gtf; ln -sf ${G%.gtf}.transcripts.fa LRAA.transcripts.fa
echo LRAA_DONE | tee -a /broad/hptmp/bhaas/TRANSDECODER_PAPER/data/gm12878/logs/lraa_denovo.log
