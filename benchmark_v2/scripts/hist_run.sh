#!/bin/bash
# Benchmark run of a historical TransDecoder release with its own defaults (supplementary history panel).
# Fresh directory per run; never resumes (an incomplete directory is moved aside).
#
# usage: hist_run.sh <version> <input.fa> <outdir> <org>
#   version: r2012   SourceForge transdecoder-r2012-08-07 (transcripts_to_best_scoring_ORFs.pl)
#            v2.0 | v4.1.0 | v5.0.0   GitHub tags (worktrees work/trees/TD_hist_<tag>)
# Releases up to v5.2.0 need DB_File: work/history/compat (links into the package cache; see
# work/history/VERSIONS.md). Genetic code: before v5.1 '-G universal' crashes Predict, so -G is passed
# only for non-standard codes (Tetrahymena); v2.0 Predict has no -G. Inputs are strand-specific (-S).
set -euo pipefail
ver=$1; input=$(readlink -f "$2"); outdir=$3; org=$4
here=$(cd "$(dirname "$0")" && pwd)
ROOT=${BENCH_ROOT:-/broad/hptmp/bhaas/TRANSDECODER_PAPER}
PYTHON=$ROOT/env/td/bin/python
export PATH=$ROOT/env/td/bin:$PATH PERL_HASH_SEED=0
export PERL5LIB=$ROOT/work/history/compat/perl5${PERL5LIB:+:$PERL5LIB}
export LD_LIBRARY_PATH=$ROOT/work/history/compat/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}
gcode=$(awk -F'\t' -v o="$org" '$1==o {print $4}' "$here/../organisms.tsv")
g=(); [ "$gcode" != Universal ] && g=(-G "$gcode")
ss=(-S); [[ $(basename "$input") == *unstranded* ]] && ss=()

[ -s "$outdir/preds.tsv" ] && { echo "$outdir done; skipping"; exit 0; }
if [ -e "$outdir" ]; then
  trash=$(dirname "$outdir")/_trash; mkdir -p "$trash"; mv "$outdir" "$trash/$(basename "$outdir").$(date +%s)"
fi
mkdir -p "$outdir"; cd "$outdir"; ln -s "$input" input.fa
t0=$(date +%s.%N)
case $ver in
  r2012)
    "$ROOT/work/history/sf_runs/r2012-08-07/transcripts_to_best_scoring_ORFs.pl" -t input.fa "${ss[@]}" "${g[@]}" > run.log 2>&1
    gff=best_candidates.eclipsed_orfs_removed.gff3 ;;
  v2.0|v4.1.0|v5.0.0)
    TD=$ROOT/work/trees/TD_hist_$ver
    "$TD/TransDecoder.LongOrfs" -t input.fa "${ss[@]}" "${g[@]}" > longorfs.log 2>&1
    pg=("${g[@]}"); [ "$ver" = v2.0 ] && pg=()
    "$TD/TransDecoder.Predict" -t input.fa "${pg[@]}" > predict.log 2>&1
    gff=input.fa.transdecoder.gff3 ;;
  *) echo "unknown version $ver" >&2; exit 1 ;;
esac
echo "$(date +%s.%N) - $t0" | bc > runtime_sec.txt
$PYTHON "$here/convert_predictions.py" transdecoder_gff3 "$gff" > preds.tsv
