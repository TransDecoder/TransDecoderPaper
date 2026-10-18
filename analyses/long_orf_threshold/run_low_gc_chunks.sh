#!/bin/bash
# Larger-sample threshold simulation for AT-rich GC levels, where long random ORFs are rare.
# Runs seeded chunks of 100k random 5 kb sequences locally, JOBS (default 6) at a time; pool
# the raw lengths with
#   pool_chunk_quantiles.py <rawdir>
# Chunks per GC level target ~100k ORFs >= 300 nt.
set -euo pipefail
PYTHON=${PYTHON:-$(command -v python3)}
here=$(cd "$(dirname "$0")" && pwd)
out=${1:?usage: $0 <outdir>}
mkdir -p "$out/raw" "$out/log"
declare -A CHUNKS=( [25]=50 [30]=20 [35]=10 [40]=5 [45]=3 )
for gc in "${!CHUNKS[@]}"; do
  for ((c=1; c<=CHUNKS[$gc]; c++)); do echo "$gc $((gc*1000 + c))"; done
done | xargs -P "${JOBS:-6}" -L 1 sh -c \
  '"$0" "$1/simulate_random_orf_lengths.py" --gc "$3" --n-seqs 100000 --seq-len 5000 --partials --seed "$4" \
     --raw-out "$2/raw/gc$3_seed$4.tsv" > "$2/log/gc$3_seed$4.log" 2>&1' "$PYTHON" "$here" "$out"
