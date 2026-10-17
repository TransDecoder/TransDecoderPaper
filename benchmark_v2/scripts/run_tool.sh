#!/bin/bash
# Run one predictor on one benchmark input and write standardized predictions.
#
# usage: run_tool.sh <tool> <input.fa> <outdir> <org>
#
# tools: pytd pytd_norefine perl_td perl_td_v600 perl_td_v571 prodigal codan
# Inputs named *.unstranded* are predicted on both strands; other inputs are treated
# as strand-specific (TransDecoder -S, Prodigal + strand only, CodAn --strand plus).
# Fragment inputs use CodAn's *_partial models, complete inputs its *_full models.
#
# Writes <outdir>/preds.tsv (transcript_id, start, end, strand; 1-based, CDS incl. stop)
# and <outdir>/runtime_sec.txt. For pytd also writes <outdir>/longest_orf.preds.tsv and
# <outdir>/all_orfs.preds.tsv (baselines from LongOrfs).
#
# Environment (defaults for the Broad workspace): BENCH_ROOT, PYTHON, CODAN_ENV, CODAN_MODELS,
# PRODIGAL, and TD_FIX/TD_V600/TD_V571/PYTD worktree paths.
set -euo pipefail

tool=$1; input=$(readlink -f "$2"); outdir=$3; org=$4
here=$(cd "$(dirname "$0")" && pwd)
ROOT=${BENCH_ROOT:-/broad/hptmp/bhaas/TRANSDECODER_PAPER}
PYTHON=${PYTHON:-$ROOT/env/td/bin/python}
CODAN_ENV=${CODAN_ENV:-$ROOT/env/codan}
CODAN_MODELS=${CODAN_MODELS:-$ROOT/tools/CodAn/models}
PRODIGAL=${PRODIGAL:-$ROOT/env/td/bin/prodigal}
TD_FIX=${TD_FIX:-$ROOT/work/trees/TD_fix}
TD_V600=${TD_V600:-$ROOT/work/trees/TD_base}
TD_V571=${TD_V571:-$ROOT/work/trees/TD_v571}
PYTD=${PYTD:-$ROOT/work/trees/PY_fix}
export PATH=$ROOT/env/td/bin:$PATH PERL_HASH_SEED=0

orgrow=$(awk -F'\t' -v o="$org" '$1==o' "$here/../organisms.tsv")
[ -n "$orgrow" ] || { echo "unknown org $org" >&2; exit 1; }
gcode=$(cut -f4 <<< "$orgrow"); table=$(cut -f5 <<< "$orgrow"); clade=$(cut -f6 <<< "$orgrow")

stranded=1; [[ $(basename "$input") == *unstranded* ]] && stranded=0
partial=0; [[ $(basename "$input") == fragments* ]] && partial=1

mkdir -p "$outdir"; outdir=$(readlink -f "$outdir"); cd "$outdir"
ln -sf "$input" input.fa
t0=$(date +%s.%N)

case $tool in
  pytd|pytd_norefine)
    ss=(); [ $stranded = 1 ] && ss=(-S)
    nr=(); [ $tool = pytd_norefine ] && nr=(--no-refine-starts)
    PYTHONPATH=$PYTD $PYTHON -m pytransdecoder.cli longorfs -t input.fa -G "$gcode" ${ss[@]+"${ss[@]}"} > longorfs.log 2>&1
    PYTHONPATH=$PYTD $PYTHON -m pytransdecoder.cli predict -t input.fa -G "$gcode" ${nr[@]+"${nr[@]}"} > predict.log 2>&1
    $PYTHON "$here/convert_predictions.py" transdecoder_gff3 input.fa.transdecoder.gff3 > preds.tsv
    if [ $tool = pytd ]; then
      lo=input.fa.transdecoder_dir/longest_orfs.gff3
      $PYTHON "$here/convert_predictions.py" transdecoder_gff3 $lo > all_orfs.preds.tsv
      $PYTHON "$here/convert_predictions.py" transdecoder_gff3 --longest-per-transcript $lo > longest_orf.preds.tsv
    fi
    ;;
  perl_td|perl_td_v600|perl_td_v571)
    case $tool in perl_td) TD=$TD_FIX;; perl_td_v600) TD=$TD_V600;; perl_td_v571) TD=$TD_V571;; esac
    LO=$TD/TransDecoder.LongOrfs; PR=$TD/TransDecoder.Predict
    [ -x "$LO" ] || { LO=$TD/util/TransDecoder.LongOrfs; PR=$TD/util/TransDecoder.Predict; }
    ss=(); [ $stranded = 1 ] && ss=(-S)
    "$LO" -t input.fa -G "$gcode" ${ss[@]+"${ss[@]}"} > longorfs.log 2>&1
    "$PR" -t input.fa -G "$gcode" > predict.log 2>&1
    $PYTHON "$here/convert_predictions.py" transdecoder_gff3 input.fa.transdecoder.gff3 > preds.tsv
    ;;
  prodigal)
    # single-genome mode trains on the input; -n forces a full motif scan (no SD trainer), as in 2018
    "$PRODIGAL" -i input.fa -f gff -o prodigal.gff -g "$table" -n -q > prodigal.log 2>&1
    sf=(); [ $stranded = 1 ] && sf=(--plus-only)
    $PYTHON "$here/convert_predictions.py" prodigal_gff ${sf[@]+"${sf[@]}"} prodigal.gff > preds.tsv
    ;;
  codan)
    model=$CODAN_MODELS/${clade}_full; [ $partial = 1 ] && model=$CODAN_MODELS/${clade}_partial
    strand=plus; [ $stranded = 0 ] && strand=both
    rm -rf codan_out
    PATH=$CODAN_ENV/bin:$PATH codan.py -t input.fa -m "$model" -s $strand -c 1 -o codan_out > codan.log 2>&1
    $PYTHON "$here/convert_predictions.py" codan_gtf codan_out/annotation.gtf > preds.tsv
    ;;
  *) echo "unknown tool $tool" >&2; exit 1;;
esac

echo "$(date +%s.%N) - $t0" | bc > runtime_sec.txt
