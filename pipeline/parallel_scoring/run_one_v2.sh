#!/bin/bash
# One fresh run of the faster PyTD branch (dev/faster-serial; worktree work/trees/PY_fast) for parity
# with the release code (PyTD d3892b5 runs in work/release_pin and the same-code a5b24c3 homology runs).
# Never resumes: an unfinished output dir goes to $P/_trash.
# usage:
#   run_one_v2.sh bench <org> <input> <rel|retonly|hassist> <cpu>
#   run_one_v2.sh gm <lraa|stringtie|trinity> <cpu>
#   run_one_v2.sh presetbuild <org> <A|Amix>
#   run_one_v2.sh presettest <org> <test_set> <preset_name> <cpu>
# Commands as pipeline/release_pin/run_bench_v1.sh, run_gm_v1.sh, run_preset_v1.sh and
# pipeline/start_classifier/11_run_pytd_v7h_homology_v1.sh, plus --cpu.
set -euo pipefail
P=/broad/hptmp/bhaas/TRANSDECODER_PAPER; T=$P/TransDecoderPaper; W=$P/work/fast_scoring
PYTD=${PYTD_FAST:-$P/work/trees/PY_fast}; PY=$P/env/td/bin/python
export PATH=$P/env/td/bin:$PATH OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONPATH=$PYTD
td() { $PY -m pytransdecoder.cli "$@"; }
gc() { awk -F'\t' -v o=$1 '$1==o{print $4}' $T/benchmark_v2/organisms.tsv; }
mode=$1; shift
case $mode in
  bench) org=$1; inp=$2; label=$3; cpu=$4; out=$W/bench/$org/$inp/$label.cpu$cpu;;
  gm) asm=$1; cpu=$2; out=$W/gm/$asm/rel.cpu$cpu;;
  presetbuild) org=$1; kind=$2; name=$org; [ $kind = Amix ] && name=${org}mix; out=$W/presets/refruns/$name;;
  presettest) org=$1; set_=$2; pname=$3; cpu=$4; out=$W/presets/runs/$org/$set_/preset_$pname.cpu$cpu;;
  *) echo "bad mode" >&2; exit 1;;
esac
[ -s $out/DONE ] && { echo "done already: $out"; exit 0; }
[ -e $out ] && { mkdir -p $P/_trash/fast_scoring_partial; mv $out $P/_trash/fast_scoring_partial/$(echo $out | tr / _).$(date +%s); }
mkdir -p $out; cd $out
(git -C $PYTD rev-parse --short HEAD; git -C $PYTD status --porcelain | wc -l) > pytd_version.txt
t0=$(date +%s.%N)
case $mode in
  bench)
    ln -s $P/data/bench_v3/inputs/$org/$inp.fa input.fa; t=input.fa
    ss=(-S); [[ $inp == *unstranded* ]] && ss=()
    td longorfs -t input.fa -G $(gc $org) ${ss[@]+"${ss[@]}"} > longorfs.log 2>&1
    pr=()
    h=$P/work/bench_v3/runs/$org/$inp/pytd_v7_homology
    case $label in
      rel) ;;
      retonly) pr=(--retain-blastp-hits $h/blastp.outfmt6 --retain-pfam-hits $h/pfam.domtblout --no-homology-start-assist);;
      hassist) pr=(--retain-blastp-hits $h/blastp.outfmt6 --retain-pfam-hits $h/pfam.domtblout);;
    esac
    t1=$(date +%s.%N)
    td predict -t input.fa -G $(gc $org) --cpu $cpu ${pr[@]+"${pr[@]}"} > predict.log 2>&1;;
  gm)
    D=$P/data/gm12878
    case $asm in
      stringtie) fa=$D/star_stringtie/stringtie.transcripts.fa;;
      lraa) fa=$D/longread_lraa/lraa_denovo/LRAA.transcripts.fa;;
      trinity) fa=$D/trinity/Trinity.fasta;;
    esac
    ln -s $(readlink -f $fa) transcripts.fa
    td longorfs -t transcripts.fa -S > longorfs.log 2>&1
    t1=$(date +%s.%N)
    td predict -t transcripts.fa --cpu $cpu > predict.log 2>&1;;
  presetbuild)
    ln -s $P/work/bench_v3/presets/ref/$org.$kind.fa ref.fa
    td longorfs -t ref.fa -G $(gc $org) -S > longorfs.log 2>&1
    t1=$(date +%s.%N)
    td predict -t ref.fa -G $(gc $org) --cpu 8 > predict.log 2>&1
    td build-preset -t ref.fa --name $name -o $out/preset --source "E21 held-out half A of $org ($kind)" > build.log 2>&1;;
  presettest)
    ln -s $(readlink -f $P/work/release_pin/presets/runs/$org/$set_/preset_$pname/input.fa) input.fa
    td longorfs -t input.fa -G $(gc $org) -S > longorfs.log 2>&1
    t1=$(date +%s.%N)
    td predict -t input.fa -G $(gc $org) --cpu $cpu --preset $P/work/release_pin/presets/presets/$pname > predict.log 2>&1;;
esac
t2=$(date +%s.%N)
echo "$t2 - $t0" | bc > runtime_sec.txt
echo -e "longorfs\t$(echo "$t1 - $t0" | bc)\npredict\t$(echo "$t2 - $t1" | bc)" > step_sec.tsv
date > DONE
