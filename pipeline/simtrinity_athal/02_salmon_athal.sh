#!/bin/bash
# E32: salmon 1.10.3 on the Arabidopsis Col-0 leaf reads (ENA ERR4140056; NEBNext Ultra Directional = ISR),
# full TAIR12 transcriptome (all transcript types; --keepDuplicates), as E30's 02_salmon.sh.
set -euo pipefail
P=${P:-/broad/hptmp/bhaas/TRANSDECODER_PAPER}
S=$P/env/trinity215/bin/salmon
R=$P/data/simtrinity_athal/ref
O=$P/data/simtrinity_athal/salmon
RD=$P/data/simtrinity_athal/reads
THREADS=${THREADS:-6}
mkdir -p $O
$S --version
[ -e $O/index.done ] || {
  $S index -t $R/all_transcripts.fa -i $O/index -k 31 --keepDuplicates -p $THREADS > $O/index.log 2>&1
  touch $O/index.done
}
[ "${INDEX_ONLY:-0}" = 1 ] && { echo INDEX_DONE; exit 0; }
$S quant -i $O/index -l ISR -1 $RD/ERR4140056_1.fastq.gz -2 $RD/ERR4140056_2.fastq.gz \
  --seqBias --gcBias -p $THREADS --validateMappings -o $O/quant > $O/quant.log 2>&1
echo SALMON_DONE
