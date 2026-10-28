#!/bin/bash
# E33 adoption: fresh PyTD runs with the homology start assist (tool name pytd_v7h_homology), reusing the
# DIAMOND / hmmsearch outputs of an existing homology run as user-supplied inputs (as a user would).
# Same settings as run_tool.sh pytd_v7_homology otherwise (stranded -> -S; organism genetic code).
# usage: 11_run_pytd_v7h_homology_v1.sh <run_dir> <input.fa> <gcode> <stranded 0|1> <hits_dir> [extra Predict args]
#   hits_dir holds blastp.outfmt6 and pfam.domtblout searched on the same LongOrfs peptides
set -euo pipefail
P=/broad/hptmp/bhaas/TRANSDECODER_PAPER; T=$P/TransDecoderPaper; PY=$P/env/td/bin/python
HS=$P/work/trees/PY_hstart; here=$T/benchmark_v2/scripts
run=$1; input=$2; gcode=$3; stranded=$4; hits=$(readlink -f $5); shift 5; extra=("$@")
[ -e "$run" ] && { echo "exists: $run" >&2; exit 1; }
mkdir -p "$run"; cd "$run"
ln -s "$input" input.fa
git -C $HS rev-parse --short HEAD > pytd_version.txt
ss=(); [ "$stranded" = 1 ] && ss=(-S)
t0=$(date +%s.%N)
PYTHONPATH=$HS $PY -m pytransdecoder.cli longorfs -t input.fa -G "$gcode" ${ss[@]+"${ss[@]}"} > longorfs.log 2>&1
PYTHONPATH=$HS $PY -m pytransdecoder.cli predict -t input.fa -G "$gcode" \
  --retain-blastp-hits $hits/blastp.outfmt6 --retain-pfam-hits $hits/pfam.domtblout ${extra[@]+"${extra[@]}"} > predict.log 2>&1
$PY $here/convert_predictions.py transdecoder_gff3 input.fa.transdecoder.gff3 > preds.tsv
echo "$(echo "$(date +%s.%N) - $t0" | bc)" > runtime_sec.txt
echo "$run done $(tail -n +2 preds.tsv | wc -l) ORFs"
