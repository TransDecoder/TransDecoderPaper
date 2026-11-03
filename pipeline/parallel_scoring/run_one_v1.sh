#!/bin/bash
# One fresh PyTD run (LongOrfs + Predict) with the parallel-scoring branch at a given --cpu.
# usage: run_one_v1.sh <org|gm> <input|trinity|stringtie|lraa> <cpu>
# Output: $P/work/parallel_scoring/<org>/<input>/cpu<N>. Never resumes (an unfinished dir goes to $P/_trash).
# Same commands as pipeline/release_pin/run_bench_v1.sh (label rel) and run_gm_v1.sh (cfg rel), plus --cpu.
set -euo pipefail
org=$1; inp=$2; cpu=$3
P=/broad/hptmp/bhaas/TRANSDECODER_PAPER; T=$P/TransDecoderPaper
PYTD=${PYTD_PAR:-$P/work/trees/PY_par}; PY=$P/env/td/bin/python
export PATH=$P/env/td/bin:$PATH OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONPATH=$PYTD
out=$P/work/parallel_scoring/$org/$inp/cpu$cpu
[ -s $out/DONE ] && { echo "done already: $out"; exit 0; }
[ -e $out ] && { mkdir -p $P/_trash/parallel_scoring_partial; mv $out $P/_trash/parallel_scoring_partial/$org.$inp.cpu$cpu.$(date +%s); }
mkdir -p $out; cd $out
if [ $org = gm ]; then
  D=$P/data/gm12878
  case $inp in
    stringtie) fa=$D/star_stringtie/stringtie.transcripts.fa;;
    lraa) fa=$D/longread_lraa/lraa_denovo/LRAA.transcripts.fa;;
    trinity) fa=$D/trinity/Trinity.fasta;;
  esac
  ln -s $(readlink -f $fa) transcripts.fa; t=transcripts.fa; lo=(-S); g=()
else
  ln -s $P/data/bench_v3/inputs/$org/$inp.fa input.fa; t=input.fa
  gcode=$(awk -F'\t' -v o=$org '$1==o{print $4}' $T/benchmark_v2/organisms.tsv)
  g=(-G $gcode); lo=(-S); [[ $inp == *unstranded* ]] && lo=()
fi
(cd $PYTD && git rev-parse --short HEAD; git -C $PYTD status --porcelain | wc -l) > pytd_version.txt
t0=$(date +%s.%N)
$PY -m pytransdecoder.cli longorfs -t $t ${g[@]+"${g[@]}"} ${lo[@]+"${lo[@]}"} > longorfs.log 2>&1
t1=$(date +%s.%N)
$PY -m pytransdecoder.cli predict -t $t ${g[@]+"${g[@]}"} --cpu $cpu > predict.log 2>&1
t2=$(date +%s.%N)
echo "$t2 - $t0" | bc > runtime_sec.txt
echo -e "longorfs\t$(echo "$t1 - $t0" | bc)\npredict\t$(echo "$t2 - $t1" | bc)" > step_sec.tsv
date > DONE
