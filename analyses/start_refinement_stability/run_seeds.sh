#!/bin/bash
# usage: run_seeds.sh <pytd_dir> <inputs_dir> <out_dir> "<datasets>" "<seeds>"
set -euo pipefail
PYTD=$(readlink -f $1); IN=$(readlink -f $2); OUT=$(mkdir -p $3 && readlink -f $3); DS=$4; SEEDS=$5
PY=${PYTHON:-python3}; HERE=$(dirname $(readlink -f $0))
declare -A LOFLAGS=( [trinity]="" [pasa]="" [stringtie]="-S" [supertx]="" [cufflinks]="" )
for ds in $DS; do for seed in $SEEDS; do
  d=$OUT/$ds/seed$seed; mkdir -p $d; cd $d; ln -sf $IN/$ds.fa .
  [ -s $ds.fa.transdecoder.gff3 ] && continue
  ( PYTHONPATH=$PYTD $PY $HERE/run_pytd_seed.py --seed $seed -- longorfs -t $ds.fa ${LOFLAGS[$ds]} > lo.log 2>&1 &&
    PYTHONPATH=$PYTD $PY $HERE/run_pytd_seed.py --seed $seed -- predict -t $ds.fa > pr.log 2>&1 ) &
done; done; wait
