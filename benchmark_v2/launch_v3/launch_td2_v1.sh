#!/bin/bash
# E39: TD2 1.1.0 (ab initio) on all benchmark inputs, mixed inputs and simulated Trinity assemblies.
# Jobs: td2_jobs_v1.tsv (input, outdir, org; largest input first). 4 jobs at a time, 2 threads each.
P=/broad/hptmp/bhaas/TRANSDECODER_PAPER; W=$P/work/bench_v3
mkdir -p $W/runs/logs/td2
export P
xargs -P ${JOBS:-4} -L 1 bash -c 'n=$(echo "$1" | sed "s#$P/##; s#/#.#g"); bash $P/TransDecoderPaper/benchmark_v2/scripts/run_td2.sh "$0" "$1" "$2" > $P/work/bench_v3/runs/logs/td2/$n.log 2>&1; echo "[$(date)] $1 exit $?"' < $W/td2_jobs_v1.tsv
echo "[$(date)] TD2_ALL_DONE"
