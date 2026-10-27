#!/bin/bash
# Quantify the real GM12878 short reads (ENCSR000AED rep 1, stranded dUTP = ISR) against the full
# GENCODE v50 transcriptome with salmon 1.10.3, to obtain per-transcript read counts and the
# fragment-length distribution used for the simulation. --keepDuplicates keeps identical
# sequences (otherwise salmon drops all but one, which could be a transcript of an excluded biotype).
# usage: 02_salmon.sh   (P defaults to the project root)
set -euo pipefail
P=${P:-/broad/hptmp/bhaas/TRANSDECODER_PAPER}
S=$P/env/trinity215/bin/salmon
R=$P/data/simtrinity/ref
O=$P/data/simtrinity/salmon
THREADS=${THREADS:-6}
mkdir -p $O
$S --version
[ -s $O/index/versionInfo.json ] && [ -e $O/index.done ] || {
  $S index -t $R/gencode.v50.transcripts.fa -i $O/index -k 31 --keepDuplicates -p $THREADS > $O/index.log 2>&1
  touch $O/index.done
}
$S quant -i $O/index -l ISR \
  -1 $P/data/gm12878/reads/short/ENCFF001REK.fastq.gz \
  -2 $P/data/gm12878/reads/short/ENCFF001REJ.fastq.gz \
  --seqBias --gcBias -p $THREADS --validateMappings -o $O/quant > $O/quant.log 2>&1
echo SALMON_DONE
