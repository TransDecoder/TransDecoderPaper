#!/bin/bash
# Full, fresh PyTD run (LongOrfs + Predict) with an experimental short-ORF rule, for the
# stringency study (E18/E19). Never resumes: an output directory without preds.tsv is moved
# aside and the run starts from an empty directory, so no checkpoint of an earlier run is reused.
#
# usage: short_orf_rule_run.sh <rule>_c<C> <input.fa> <outdir> <org>
#   rule: cliff | fixed | inverse | taper   (Predict --short-orf-rule; see PY_taper predict.py)
#   C:    Predict --short-orf-min-contrast
# Tree: $PYTD_TAPER (default work/trees/PY_taper, branch exp/short-orf-taper, from PyTD 44e417a).
set -euo pipefail
setting=$1; input=$(readlink -f "$2"); outdir=$3; org=$4
here=$(cd "$(dirname "$0")" && pwd)
ROOT=${BENCH_ROOT:-/broad/hptmp/bhaas/TRANSDECODER_PAPER}
PYTHON=${PYTHON:-$ROOT/env/td/bin/python}
PYTD=${PYTD_TAPER:-$ROOT/work/trees/PY_taper}
export PATH=$ROOT/env/td/bin:$PATH
[[ $setting =~ ^(cliff|fixed|inverse|taper)_c([0-9.]+)$ ]] || { echo "bad setting $setting" >&2; exit 1; }
rule=${BASH_REMATCH[1]}; C=${BASH_REMATCH[2]}
gcode=$(awk -F'\t' -v o="$org" '$1==o {print $4}' "$here/../organisms.tsv")
[ -n "$gcode" ] || { echo "unknown org $org" >&2; exit 1; }
ss=(-S); [[ $(basename "$input") == *unstranded* ]] && ss=()

[ -s "$outdir/preds.tsv" ] && { echo "$outdir done; skipping"; exit 0; }
if [ -e "$outdir" ]; then  # incomplete earlier attempt: set aside, never resume from it
  trash=$(dirname "$outdir")/_trash; mkdir -p "$trash"; mv "$outdir" "$trash/$(basename "$outdir").$(date +%s)"
fi
mkdir -p "$outdir"; cd "$outdir"
ln -s "$input" input.fa
t0=$(date +%s.%N)
PYTHONPATH=$PYTD $PYTHON -m pytransdecoder.cli longorfs -t input.fa -G "$gcode" "${ss[@]}" > longorfs.log 2>&1
PYTHONPATH=$PYTD $PYTHON -m pytransdecoder.cli predict -t input.fa -G "$gcode" \
  --short-orf-rule "$rule" --short-orf-min-contrast "$C" > predict.log 2>&1
echo "$(date +%s.%N) - $t0" | bc > runtime_sec.txt
(cd "$PYTD" && git rev-parse --short HEAD && git diff --quiet || echo "+uncommitted") > pytd_version.txt
$PYTHON "$here/convert_predictions.py" transdecoder_gff3 input.fa.transdecoder.gff3 > preds.tsv
