#!/bin/bash
# E32: Polyester simulation of the Arabidopsis library, as E30's 06_simulate.sh (same R and conversion
# scripts and patches); chunk k uses seed 20261027 + k. Output: $O/fastq/chunk_<k>_R{1,2}.fastq.gz (RF).
set -euo pipefail
P=/broad/hptmp/bhaas/TRANSDECODER_PAPER
here=$P/TransDecoderPaper/pipeline/simtrinity
O=$P/data/simtrinity_athal/sim
JOBS=${JOBS:-4}
mean=$(awk -F'\t' '$1=="fld_mean"{print $2}' $O/sim_params.tsv)
sd=$(awk -F'\t' '$1=="fld_sd"{print $2}' $O/sim_params.tsv)
mkdir -p $O/fastq $O/logs
n=$(ls $O/chunk_*.counts.tsv | wc -l)
echo "[$(date)] $n chunks, fragment length $mean sd $sd, $JOBS at a time"
run_chunk() {
  k=$1
  $P/env/polyester/bin/Rscript $here/04_polyester_chunk.R $O/chunk_$k.fa $O/chunk_$k.counts.tsv \
    $O/polyester_chunk_$k $((20261027 + k)) $mean $sd > $O/logs/chunk_$k.polyester.log 2>&1
  grep -q POLYESTER_CHUNK_DONE $O/logs/chunk_$k.polyester.log
  $P/env/td/bin/python $here/05_polyester_to_fastq.py $O/polyester_chunk_$k c$k $O/fastq/chunk_$k > $O/fastq/chunk_${k}.done
  echo "[$(date)] chunk $k: $(cat $O/fastq/chunk_${k}.done) pairs"
}
export -f run_chunk; export P O here mean sd
seq 0 $((n - 1)) | xargs -P $JOBS -I{} bash -c 'run_chunk {}'
echo SIMULATE_DONE
