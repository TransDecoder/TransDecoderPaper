#!/bin/bash
# LongOrfs --cpu parity/timing (PyTD devel 725244f, worktree work/trees/PY_lo).
# run_lo_v3.sh bench <org> <input> <cpu> | gm <asm> <cpu> | pipe <org> <input> <cpu>
#   bench/gm: LongOrfs only, as the release runs (-G from organisms.tsv, -S unless unstranded);
#   pipe: pyTransdecoder (LongOrfs + Predict) with --cpu, release defaults.
# Output: work/longorfs_cpu/<mode>/<org|gm>/<input|asm>/cpu<N>/ with runtime_sec.txt, DONE.
set -eo pipefail
P=/broad/hptmp/bhaas/TRANSDECODER_PAPER; T=$P/TransDecoderPaper; PY=$P/env/td/bin/python
export PYTHONPATH=$P/work/trees/PY_lo
gc() { awk -F'\t' -v o=$1 '$1==o{print $4}' $T/benchmark_v2/organisms.tsv; }
mode=$1; shift
case $mode in
  bench|pipe) org=$1; inp=$2; cpu=$3; out=$P/work/longorfs_cpu/$mode/$org/$inp/cpu$cpu
              src=$P/data/bench_v3/inputs/$org/$inp.fa; t=input.fa
              ss=(-S); [[ $inp == *unstranded* ]] && ss=(); g=(-G $(gc $org));;
  gm) asm=$1; cpu=$2; out=$P/work/longorfs_cpu/gm/gm/$asm/cpu$cpu; D=$P/data/gm12878
      case $asm in
        stringtie) src=$D/star_stringtie/stringtie.transcripts.fa;;
        lraa) src=$D/longread_lraa/lraa_denovo/LRAA.transcripts.fa;;
        trinity) src=$D/trinity/Trinity.fasta;;
      esac
      t=transcripts.fa; ss=(-S); g=();;
esac
[ -e $out/DONE ] && { echo "exists: $out"; exit 0; }
mkdir -p $out; cd $out
ln -sf $(readlink -f $src) $t
git -C $P/work/trees/PY_lo rev-parse --short HEAD > pytd_version.txt
t0=$(date +%s.%N)
if [ $mode = pipe ]; then
  $PY -c "import sys; from pytransdecoder.cli import pipeline_main; sys.argv=['pyTransdecoder']+sys.argv[1:]; pipeline_main()" \
      -t $t ${g[@]+"${g[@]}"} ${ss[@]+"${ss[@]}"} --cpu $cpu > pipeline.log 2>&1
else
  $PY -m pytransdecoder.cli longorfs -t $t ${g[@]+"${g[@]}"} ${ss[@]+"${ss[@]}"} --cpu $cpu > longorfs.log 2>&1
fi
echo "$(date +%s.%N) - $t0" | bc > runtime_sec.txt
touch DONE
