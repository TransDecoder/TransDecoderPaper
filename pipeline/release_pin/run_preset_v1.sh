#!/bin/bash
# Held-out preset test at the pinned release code (E40; design E21).
# usage: run_preset_v1.sh ref <org> <A|Amix>                build the preset <org> or <org>mix from a fresh run on reference half A
#        run_preset_v1.sh test <setting> <test.fa> <org>    setting: self | preset_<name>
set -euo pipefail
P=/broad/hptmp/bhaas/TRANSDECODER_PAPER; T=$P/TransDecoderPaper; W=$P/work/release_pin/presets; OLD=$P/work/bench_v3/presets
PYTD=${PYTD_REL:-$P/work/trees/PY_rel_d3892b5}; PY=$P/env/td/bin/python
export PATH=$P/env/td/bin:$PATH PYTHONPATH=$PYTD OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
td() { $PY -m pytransdecoder.cli "$@"; }
gc() { awk -F'\t' -v o=$1 '$1==o{print $4}' $T/benchmark_v2/organisms.tsv; }
if [ $1 = ref ]; then
  org=$2; kind=$3; name=$org; [ $kind = Amix ] && name=${org}mix
  out=$W/refruns/$name
  [ -s $W/presets/$name/preset.json ] && exit 0
  [ -e $out ] && mv $out $P/_trash/release_pin_partial.preset.$name.$(date +%s)
  mkdir -p $out $W/presets; cd $out; ln -s $OLD/ref/$org.$kind.fa ref.fa
  (cd $PYTD && git rev-parse --short HEAD) > pytd_version.txt
  td longorfs -t ref.fa -G $(gc $org) -S > longorfs.log 2>&1
  td predict -t ref.fa -G $(gc $org) > predict.log 2>&1
  td build-preset -t ref.fa --name $name -o $W/presets/$name --source "E21 held-out half A of $org ($kind)" > build.log 2>&1
  exit 0
fi
setting=$2; fa=$3; org=$4; set_=$(basename $fa .fa)
out=$W/runs/$org/$set_/$setting
[ -s $out/DONE ] && exit 0
[ -e $out ] && { mkdir -p $P/_trash/release_pin_partial; mv $out $P/_trash/release_pin_partial/preset.$org.$set_.$setting.$(date +%s); }
mkdir -p $out; cd $out; ln -s $fa input.fa
pr=(); [ $setting != self ] && pr=(--preset $W/presets/${setting#preset_})
ss=(-S); (cd $PYTD && git rev-parse --short HEAD) > pytd_version.txt
td longorfs -t input.fa -G $(gc $org) "${ss[@]}" > longorfs.log 2>&1
td predict -t input.fa -G $(gc $org) ${pr[@]+"${pr[@]}"} > predict.log 2>&1
$PY $T/benchmark_v2/scripts/convert_predictions.py transdecoder_gff3 input.fa.transdecoder.gff3 > preds.tsv
date > DONE
