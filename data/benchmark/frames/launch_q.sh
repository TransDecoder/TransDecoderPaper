#!/bin/bash
# runner: TransDecoderPaper/benchmark_v2/scripts/frames_run.sh (do not edit while this runs); all-frames at q 0.1 / 0.2
P=/broad/hptmp/bhaas/TRANSDECODER_PAPER; W=$P/work/bench_v3/frames; S=$P/TransDecoderPaper/benchmark_v2/scripts/frames_run.sh
cat $W/jobs_q.txt | xargs -P 3 -n 3 bash -c 'if '$S' '$P'/data/bench_v3/inputs/$1/$2.fa '$W'/runs/$1/$2/all_q$3 $1 all --max-qvalue $3 > '$W'/logs/$1.$2.q$3.log 2>&1; then echo "done $1 $2 $3"; else echo "FAILED $1 $2 $3"; fi' _
echo ALL_DONE
