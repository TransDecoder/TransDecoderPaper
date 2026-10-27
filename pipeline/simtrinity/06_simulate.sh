#!/bin/bash
# Simulate the GM12878-like stranded paired-end library with Polyester, in parallel chunks (E30).
# usage: 06_simulate.sh            (after 03_make_counts.py; JOBS chunks at a time, default 5)
# Chunk k uses seed 20261025 + k. Each chunk: 04_polyester_chunk.R, then 05_polyester_to_fastq.py
# (R1 = antisense mate, RF layout). Output: $O/fastq/chunk_<k>_R{1,2}.fastq.gz
set -euo pipefail
P=${P:-/broad/hptmp/bhaas/TRANSDECODER_PAPER}
here=$(cd "$(dirname "$0")" && pwd)
O=$P/data/simtrinity/sim
JOBS=${JOBS:-5}
mean=$(awk -F'\t' '$1=="fld_mean"{print $2}' $O/sim_params.tsv)
sd=$(awk -F'\t' '$1=="fld_sd"{print $2}' $O/sim_params.tsv)
mkdir -p $O/fastq $O/logs
n=$(ls $O/chunk_*.counts.tsv | wc -l)
echo "[$(date)] $n chunks, fragment length $mean sd $sd, $JOBS at a time"
run_chunk() {
  k=$1
  [ -s $O/fastq/chunk_${k}.done ] && { echo "chunk $k done"; return 0; }
  rm -rf $O/polyester_chunk_$k
  $P/env/polyester/bin/Rscript $here/04_polyester_chunk.R $O/chunk_$k.fa $O/chunk_$k.counts.tsv \
    $O/polyester_chunk_$k $((20261025 + k)) $mean $sd > $O/logs/chunk_$k.polyester.log 2>&1
  grep -q POLYESTER_CHUNK_DONE $O/logs/chunk_$k.polyester.log
  $P/env/td/bin/python $here/05_polyester_to_fastq.py $O/polyester_chunk_$k c$k $O/fastq/chunk_$k > $O/fastq/chunk_${k}.done
  echo "[$(date)] chunk $k: $(cat $O/fastq/chunk_${k}.done) pairs"
}
export -f run_chunk; export P O here mean sd
seq 0 $((n - 1)) | xargs -P $JOBS -I{} bash -c 'run_chunk {}'
echo SIMULATE_DONE
