#!/bin/bash
# Larger-sample threshold simulation for AT-rich GC levels, where long random ORFs are rare.
# Submits seeded chunks of 100k random 5 kb sequences to UGER; pool the raw lengths with
#   pool_chunk_quantiles.py <rawdir>
# Chunks per GC level target ~100k ORFs >= 300 nt.
set -euo pipefail
PYTHON=${PYTHON:-$(command -v python3)}
here=$(cd "$(dirname "$0")" && pwd)
out=${1:?usage: $0 <outdir>}
mkdir -p "$out/raw" "$out/log"
declare -A CHUNKS=( [25]=50 [30]=20 [35]=10 [40]=5 [45]=3 )
for gc in "${!CHUNKS[@]}"; do
  for ((c=1; c<=CHUNKS[$gc]; c++)); do
    seed=$((gc*1000 + c))
    qsub -cwd -V -N thr${gc}_${c} -l h_rt=2:00:00 -l h_vmem=2G -o "$out/log" -j y -b y \
      "$PYTHON" "$here/simulate_random_orf_lengths.py" --gc "$gc" --n-seqs 100000 --seq-len 5000 \
        --partials --seed "$seed" --raw-out "$out/raw/gc${gc}_seed${seed}.tsv" > /dev/null
  done
done
