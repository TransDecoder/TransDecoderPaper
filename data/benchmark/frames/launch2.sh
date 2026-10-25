#!/bin/bash
# runner: TransDecoderPaper/benchmark_v2/scripts/frames_run.sh (do not edit while this runs)
P=/broad/hptmp/bhaas/TRANSDECODER_PAPER; W=$P/work/bench_v3/frames; S=$P/TransDecoderPaper/benchmark_v2/scripts/frames_run.sh
cat $W/jobs.txt | xargs -P 3 -n 2 bash -c 'if '$S' '$P'/data/bench_v3/inputs/$1/$2.fa '$W'/runs/$1/$2/all $1 all > '$W'/logs/$1.$2.log 2>&1; then echo "done $1 $2"; else echo "FAILED $1 $2"; fi' _
echo ALL_DONE
