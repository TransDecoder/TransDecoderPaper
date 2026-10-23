#!/bin/bash
# Run one ORF predictor at a given minimum protein length (LongOrfs -m), for the lncRNA
# false-positive comparison across minimum lengths (as in the TD2 preprint, Table 4).
#
# usage: minlen_run.sh tool org variant min_aa
#   tool: pytd (final defaults, PY_dev), perl_v571 (TransDecoder v5.7.1), td2 (TD2 1.1.0, -m and -M)
# Output: $ROOT/work/bench_v3/minlen/<org>/<variant>/<tool>_m<min_aa>/preds.tsv
set -euo pipefail
tool=$1 org=$2 variant=$3 m=$4
here=$(cd "$(dirname "$0")" && pwd)
ROOT=${ROOT:-/broad/hptmp/bhaas/TRANSDECODER_PAPER}
PYTHON=$ROOT/env/td/bin/python
export PATH=$ROOT/env/td/bin:$PATH PERL_HASH_SEED=0
out=$ROOT/work/bench_v3/minlen/$org/$variant/${tool}_m$m
[ -s "$out/preds.tsv" ] && exit 0
mkdir -p "$out"; cd "$out"
ln -sf "$ROOT/data/bench_v3/inputs/$org/$variant.fa" input.fa
gcode=$(awk -F'\t' -v o="$org" '$1==o {print $4}' "$here/../organisms.tsv")
ncbi=$(awk -F'\t' -v o="$org" '$1==o {print $5}' "$here/../organisms.tsv")
s=$(date +%s.%N)
case $tool in
  pytd)
    PY=$ROOT/work/trees/PY_dev
    PYTHONPATH=$PY $PYTHON -m pytransdecoder.cli longorfs -t input.fa -m "$m" -S -G "$gcode" > longorfs.log 2>&1
    PYTHONPATH=$PY $PYTHON -m pytransdecoder.cli predict -t input.fa -G "$gcode" > predict.log 2>&1
    gff=input.fa.transdecoder.gff3 ;;
  perl_v571)
    TD=$ROOT/work/trees/TD_v571
    "$TD/TransDecoder.LongOrfs" -t input.fa -m "$m" -S -G "$gcode" > longorfs.log 2>&1
    "$TD/TransDecoder.Predict" -t input.fa -G "$gcode" > predict.log 2>&1
    gff=input.fa.transdecoder.gff3 ;;
  td2)
    export PATH=$ROOT/env/td2/bin:$PATH
    TD2.LongOrfs -t input.fa -m "$m" -M "$m" -S -G "$ncbi" -@ 2 > longorfs.log 2>&1
    TD2.Predict -t input.fa -G "$ncbi" > predict.log 2>&1
    gff=input.fa.TD2.gff3 ;;
  *) echo "unknown tool $tool" >&2; exit 1 ;;
esac
echo "$(date +%s.%N) - $s" | bc > runtime_sec.txt
$PYTHON "$here/convert_predictions.py" transdecoder_gff3 "$gff" > preds.tsv
