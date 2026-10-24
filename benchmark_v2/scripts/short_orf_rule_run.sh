#!/bin/bash
# Full, fresh PyTD run (LongOrfs + Predict) with an experimental short-ORF rule, for the
# stringency study (E18/E19). Never resumes: an output directory without preds.tsv is moved
# aside and the run starts from an empty directory, so no checkpoint of an earlier run is reused.
#
# usage: short_orf_rule_run.sh <rule>_c<C> <input.fa> <outdir> <org>
#   rule: cliff | fixed | inverse | taper | curve   (Predict --short-orf-rule; see PY_taper predict.py)
#   curve needs a 5th argument: the log prior-odds curve file (--short-orf-curve)
#   C:    Predict --short-orf-min-contrast
#   or decoy_q<Q>[_<strata>][_long<gc|decoy>]: Predict --short-orf-rule decoy --decoy-q Q
#      [--decoy-strata strata] [--decoy-long-orfs gc|decoy]
# Tree: $PYTD_TAPER (default work/trees/PY_taper, branch exp/short-orf-taper, from PyTD 44e417a).
set -euo pipefail
setting=$1; input=$(readlink -f "$2"); outdir=$3; org=$4; curve=${5:-}
here=$(cd "$(dirname "$0")" && pwd)
ROOT=${BENCH_ROOT:-/broad/hptmp/bhaas/TRANSDECODER_PAPER}
PYTHON=${PYTHON:-$ROOT/env/td/bin/python}
PYTD=${PYTD_TAPER:-$ROOT/work/trees/PY_taper}
export PATH=$ROOT/env/td/bin:$PATH
if [[ $setting =~ ^decoy_q([0-9.]+)(_(len|len_tx|len_tx_type))?(_long(gc|decoy))?$ ]]; then
  rule=decoy; pargs=(--short-orf-rule decoy --decoy-q "${BASH_REMATCH[1]}")
  [ -n "${BASH_REMATCH[3]}" ] && pargs+=(--decoy-strata "${BASH_REMATCH[3]}")
  [ -n "${BASH_REMATCH[5]}" ] && pargs+=(--decoy-long-orfs "${BASH_REMATCH[5]}")
elif [[ $setting =~ ^(cliff|fixed|inverse|taper|curve)_c([0-9.]+)$ ]]; then
  rule=${BASH_REMATCH[1]}; pargs=(--short-orf-rule "$rule" --short-orf-min-contrast "${BASH_REMATCH[2]}")
else
  echo "bad setting $setting" >&2; exit 1
fi
gcode=$(awk -F'\t' -v o="$org" '$1==o {print $4}' "$here/../organisms.tsv")
[ -n "$gcode" ] || { echo "unknown org $org" >&2; exit 1; }
ss=(-S); [[ $(basename "$input") == *unstranded* ]] && ss=()

[ "$rule" = curve ] && [ ! -s "$curve" ] && { echo "curve rule needs a curve file" >&2; exit 1; }
[ -s "$outdir/preds.tsv" ] && { echo "$outdir done; skipping"; exit 0; }
if [ -e "$outdir" ]; then  # incomplete earlier attempt: set aside, never resume from it
  trash=$(dirname "$outdir")/_trash; mkdir -p "$trash"; mv "$outdir" "$trash/$(basename "$outdir").$(date +%s)"
fi
mkdir -p "$outdir"; cd "$outdir"
ln -s "$input" input.fa
[ -n "$curve" ] && cp "$curve" short_orf_curve.tsv
t0=$(date +%s.%N)
PYTHONPATH=$PYTD $PYTHON -m pytransdecoder.cli longorfs -t input.fa -G "$gcode" "${ss[@]}" > longorfs.log 2>&1
PYTHONPATH=$PYTD $PYTHON -m pytransdecoder.cli predict -t input.fa -G "$gcode" \
  "${pargs[@]}" ${curve:+--short-orf-curve "$curve"} > predict.log 2>&1
echo "$(date +%s.%N) - $t0" | bc > runtime_sec.txt
(cd "$PYTD" && git rev-parse --short HEAD && git diff --quiet || echo "+uncommitted") > pytd_version.txt
$PYTHON "$here/convert_predictions.py" transdecoder_gff3 input.fa.transdecoder.gff3 > preds.tsv
