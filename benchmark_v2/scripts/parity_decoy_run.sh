#!/bin/bash
# Perl TransDecoder vs PyTransDecoder on one input with the decoy-FDR ORF rule: runs both fresh
# (an existing outdir is moved aside) and compares decoy/preset q-values and predicted coordinates.
# usage: parity_decoy_run.sh <input.fa> <outdir> <genetic code> [stranded: S|U] [extra Predict args, Perl spelling]
#   extra args are passed to both (Perl spelling: --preset X, --max_qvalue Q, --precise, --sensitive)
set -uo pipefail
input=$(readlink -f "$1"); out=$2; gcode=$3; strand=${4:-S}; shift 4 2>/dev/null; extra=("$@")
ROOT=${BENCH_ROOT:-/broad/hptmp/bhaas/TRANSDECODER_PAPER}
TD=${TD_TREE:-$ROOT/work/trees/TD_fix}; PY=${PY_TREE:-$ROOT/work/trees/PY_taper}
export PATH=$ROOT/env/td/bin:$PATH PERL_HASH_SEED=0
here=$(cd "$(dirname "$0")" && pwd)
if [ -e "$out" ]; then trash=$(dirname "$out")/_trash; mkdir -p "$trash"; mv "$out" "$trash/$(basename "$out").$(date +%s)"; fi
mkdir -p "$out/pl" "$out/py"; ln -s "$input" "$out/pl/input.fa"; ln -s "$input" "$out/py/input.fa"
ss=(); [ "$strand" = S ] && ss=(-S)
pyextra=(); for a in "${extra[@]}"; do if [[ $a == --* ]]; then pyextra+=("${a//_/-}"); else pyextra+=("$a"); fi; done
( cd "$out/pl" && "$TD/util/TransDecoder.LongOrfs" -t input.fa -G "$gcode" "${ss[@]}" > lo.log 2>&1 \
  && "$TD/util/TransDecoder.Predict" -t input.fa -G "$gcode" "${extra[@]}" > pr.log 2>&1; echo $? > exit.txt ) &
( cd "$out/py" && export PYTHONPATH=$PY && python -m pytransdecoder.cli longorfs -t input.fa -G "$gcode" "${ss[@]}" > lo.log 2>&1 \
  && python -m pytransdecoder.cli predict -t input.fa -G "$gcode" "${pyextra[@]}" > pr.log 2>&1; echo $? > exit.txt ) &
wait
for d in pl py; do python "$here/convert_predictions.py" transdecoder_gff3 "$out/$d/input.fa.transdecoder.gff3" > "$out/$d/preds.tsv" 2>/dev/null; done
{
  echo -e "exit\tperl=$(cat $out/pl/exit.txt)\tpy=$(cat $out/py/exit.txt)"
  for f in decoy_qvalues.tsv; do
    a=$out/pl/input.fa.transdecoder_dir/$f; b=$out/py/input.fa.transdecoder_dir/$f
    if [ -e "$a" ] || [ -e "$b" ]; then cmp -s "$a" "$b" && echo -e "$f\tidentical\t$(grep -vc '^#' "$a")" || echo -e "$f\tDIFFER\t$(diff "$a" "$b" | grep -c '^<')"; fi
  done
  cmp -s "$out/pl/preds.tsv" "$out/py/preds.tsv" && echo -e "preds\tidentical\t$(($(wc -l < "$out/pl/preds.tsv") - 1))" \
    || echo -e "preds\tDIFFER\t$(diff "$out/pl/preds.tsv" "$out/py/preds.tsv" | grep -c '^<') of $(wc -l < "$out/pl/preds.tsv")"
} > "$out/parity.txt"
cat "$out/parity.txt"
