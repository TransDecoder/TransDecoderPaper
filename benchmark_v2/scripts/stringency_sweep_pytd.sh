#!/bin/bash
# Rerun PyTD Predict ORF selection (and start refinement) at one stringency setting, reusing the
# training, hexamer model and ORF scores of the final-default run (pytd_dev_v2).
#
# usage: stringency_sweep_pytd.sh org variant setting
#   setting: s<L>_c<C>[_noauto]
#     L        ORFs shorter than L codons need frame contrast >= C bases x the training ORFs'
#              median per-base contrast (--short-orf-length / --short-orf-min-contrast);
#              L = all applies the test to every ORF
#     _noauto  no ORF is kept on length alone (--retain-long-orfs-mode strict, length 10^9)
#   defaults are s100_c200.
# Output: $ROOT/work/bench_v3/stringency/<org>/<variant>/<setting>/preds.tsv
set -euo pipefail
org=$1 variant=$2 setting=$3
here=$(cd "$(dirname "$0")" && pwd)
ROOT=${ROOT:-/broad/hptmp/bhaas/TRANSDECODER_PAPER}
PYTHON=$ROOT/env/td/bin/python
PYTD_DEV=${PYTD_DEV:-$ROOT/work/trees/PY_dev}
src=$ROOT/work/bench_v3/runs/$org/$variant/pytd_dev_v2
base=$ROOT/work/bench_v3/stringency/$org/$variant/_base
out=$ROOT/work/bench_v3/stringency/$org/$variant/$setting
gcode=$(awk -F'\t' -v o="$org" '$1==o {print $4}' "$here/../organisms.tsv")

[[ $setting =~ ^s([0-9]+|all)_c([0-9.]+)(_noauto)?$ ]] || { echo "bad setting $setting" >&2; exit 1; }
L=${BASH_REMATCH[1]}; C=${BASH_REMATCH[2]}; [ "$L" = all ] && L=1000000000
args=(--short-orf-length "$L" --short-orf-min-contrast "$C")
[ -n "${BASH_REMATCH[3]}" ] && args+=(--retain-long-orfs-mode strict --retain-long-orfs-length 1000000000)

# one private copy of the scored work directory per input, so the benchmark run is never touched
if [ ! -e "$base/.ok" ]; then
  mkdir -p "$base"
  cp "$src/input.fa" "$base/"
  d=$src/input.fa.transdecoder_dir
  mkdir -p "$base/td"
  cp "$d"/{base_freqs.dat,hexamer.scores,longest_orfs.cds,longest_orfs.cds.scores,longest_orfs.gff3,longest_orfs.pep} "$base/td/"
  mkdir -p "$base/td/__checkpoints_TDpredict"
  cp "$d"/__checkpoints_TDpredict/{training_orfs.ok,top_training_orfs.cds,hexamer_scores.ok,cds_scores.ok} \
     "$base/td/__checkpoints_TDpredict/"
  touch "$base/.ok"
fi

[ -s "$out/preds.tsv" ] && exit 0
[ -e "$out" ] && mv "$out" "${TRASH:-$ROOT/work/bench_v3/stringency/_trash}/$org.$variant.$setting.$$"
mkdir -p "$out"; cd "$out"
ln -s "$base/input.fa" input.fa
mkdir -p input.fa.transdecoder_dir/__checkpoints_TDpredict
for f in "$base"/td/*.*; do ln -s "$f" input.fa.transdecoder_dir/; done
cp "$base"/td/__checkpoints_TDpredict/* input.fa.transdecoder_dir/__checkpoints_TDpredict/
PYTHONPATH=$PYTD_DEV $PYTHON -m pytransdecoder.cli predict -t input.fa -G "$gcode" "${args[@]}" > predict.log 2>&1
$PYTHON "$here/convert_predictions.py" transdecoder_gff3 input.fa.transdecoder.gff3 > preds.tsv
