#!/bin/bash
# Frame-contrast experiment (E26): PyTD exp/contrast-frames (work/trees/PY_frames) with --contrast-frames <mode>.
# Fresh run in <outdir> (an existing <outdir> is refused). Same settings as run_tool.sh pytd_v7 otherwise.
# usage: frames_run.sh <input.fa> <outdir> <org> <forward|all> [extra predict args...]
set -euo pipefail
input=$(readlink -f "$1"); outdir=$2; org=$3; mode=$4; shift 4
here=$(cd "$(dirname "$0")" && pwd)
ROOT=/broad/hptmp/bhaas/TRANSDECODER_PAPER
PYTHON=$ROOT/env/td/bin/python
PYTD=${PYTD_FRAMES:-$ROOT/work/trees/PY_frames}
orgrow=$(awk -F'\t' -v o="$org" '$1==o' "$here/../organisms.tsv")
[ -n "$orgrow" ] || { echo "unknown org $org" >&2; exit 1; }
gcode=$(cut -f4 <<< "$orgrow")
ss=(); [[ $(basename "$input") == *unstranded* ]] || ss=(-S)
[ -e "$outdir" ] && { echo "$outdir exists; refusing to reuse" >&2; exit 1; }
mkdir -p "$outdir"; cd "$outdir"
ln -s "$input" input.fa
(cd "$PYTD" && git rev-parse --short HEAD) > pytd_version.txt
t0=$(date +%s.%N)
PYTHONPATH=$PYTD $PYTHON -m pytransdecoder.cli longorfs -t input.fa -G "$gcode" ${ss[@]+"${ss[@]}"} > longorfs.log 2>&1
PYTHONPATH=$PYTD $PYTHON -m pytransdecoder.cli predict -t input.fa -G "$gcode" --contrast-frames "$mode" "$@" > predict.log 2>&1
$PYTHON "$here/convert_predictions.py" transdecoder_gff3 input.fa.transdecoder.gff3 > preds.tsv
echo "$(date +%s.%N) - $t0" | bc > runtime_sec.txt
