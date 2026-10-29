#!/bin/bash
# E34: fresh PyTD a5b24c3 runs for two counterfactuals (new directories only; the release runs are not touched).
#   minlen30: fresh LongOrfs -m 30 (-S, Universal) + default Predict, on the input of each listed dataset
#   sensitive: the release run's LongOrfs files copied to a new directory (no Predict checkpoints) + Predict --sensitive;
#              a check that the offline q 0.2 selection of 01_trace_missed_v1.py equals the program's
# usage: 02_fresh_runs_v2.sh <mode> <dataset> <manifest.tsv>   (v2: the manifest is an argument)
set -euo pipefail
mode=$1; ds=$2; manifest=$3
P=/broad/hptmp/bhaas/TRANSDECODER_PAPER; T=$P/TransDecoderPaper; HS=$P/work/trees/PY_hstart
PY=$P/env/td/bin/python
[ "$(git -C $HS rev-parse --short HEAD)" = a5b24c3 ] || { echo "PY_hstart not at a5b24c3" >&2; exit 1; }
run=$(awk -F'\t' -v d=$ds '$1==d{print $5}' $manifest)
[ -n "$run" ] || { echo "unknown dataset $ds" >&2; exit 1; }
fa=$(readlink -f $run/input.fa)
O=$P/work/missed_cds/$mode/$ds/runs/pytd_v7
[ -e $O ] && { echo "exists: $O" >&2; exit 1; }
mkdir -p $O; cd $O; ln -s $fa input.fa
echo a5b24c3 > pytd_version.txt
t0=$(date +%s)
case $mode in
  minlen30)
    PYTHONPATH=$HS nice -n 10 $PY -m pytransdecoder.cli longorfs -t input.fa -G Universal -S -m 30 > longorfs.log 2>&1
    PYTHONPATH=$HS nice -n 10 $PY -m pytransdecoder.cli predict -t input.fa -G Universal > predict.log 2>&1 ;;
  sensitive)
    mkdir input.fa.transdecoder_dir
    for f in longest_orfs.gff3 longest_orfs.cds longest_orfs.pep base_freqs.dat longorfs.params.json; do
      cp $run/input.fa.transdecoder_dir/$f input.fa.transdecoder_dir/; done
    PYTHONPATH=$HS nice -n 10 $PY -m pytransdecoder.cli predict -t input.fa -G Universal --sensitive > predict.log 2>&1 ;;
  *) echo "unknown mode $mode" >&2; exit 1 ;;
esac
$PY $T/benchmark_v2/scripts/convert_predictions.py transdecoder_gff3 input.fa.transdecoder.gff3 > preds.tsv
echo $(( $(date +%s) - t0 )) > runtime_sec.txt
echo "[$(date)] $mode $ds done: $(tail -n +2 preds.tsv | wc -l) ORFs"
