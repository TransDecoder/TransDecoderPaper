#!/bin/bash
# Homology vs ab initio table for the released homology mode (PyTD a5b24c3, homology-assisted start
# refinement; runs pytd_v7h_homology) against the ab initio release runs (pytd_v7, 3a793ff).
# Same scoring (data/benchmark/final_v7/score_v7.sh -> score_predictions.py) and same pairing
# (benchmark_v2/scripts/compare_homology.py) as data/benchmark/homology_vs_abinitio.tsv (E25).
# Outputs:
#   data/benchmark/final_v7/v7h_homology_scores.tsv  pytd_v7 rows (copied from all_scores.tsv) + pytd_v7h_homology rows
#   data/benchmark/homology_v7h_vs_abinitio.tsv
# usage: bash pipeline/derive/05_homology_v7h_table.sh      (4 parallel scorers)
set -euo pipefail
P=/broad/hptmp/bhaas/TRANSDECODER_PAPER
T=$P/TransDecoderPaper
cd $T
ALL=data/benchmark/final_v7/all_scores.tsv
OUT=data/benchmark/final_v7/v7h_homology_scores.tsv
tmp=$(mktemp -d -p $P/tmp 2>/dev/null || mktemp -d)
# every input that has a pytd_v7 row
awk -F'\t' 'NR>1 && $3=="pytd_v7" {print $1"\t"$2}' $ALL > $tmp/inputs.txt
for f in $(awk -F'\t' '{print $1"/"$2}' $tmp/inputs.txt); do
  v=$(grep -c . $P/work/bench_v3/runs/$f/pytd_v7h_homology/pytd_version.txt >/dev/null && cat $P/work/bench_v3/runs/$f/pytd_v7h_homology/pytd_version.txt)
  [ "$v" = a5b24c3 ] || { echo "bad version $f: $v" >&2; exit 1; }
done
tr '\t' ' ' < $tmp/inputs.txt | xargs -P 4 -L 1 sh -c 'bash data/benchmark/final_v7/score_v7.sh "$0" "$1" pytd_v7h_homology > '"$tmp"'/"$0"__"$1".row'
{ head -1 $ALL
  awk -F'\t' 'NR>1 && $3=="pytd_v7"' $ALL
  cat $tmp/*.row | sort -k1,1 -k2,2; } > $OUT
n=$(awk -F'\t' '$3=="pytd_v7h_homology"' $OUT | wc -l)
[ "$n" -eq "$(wc -l < $tmp/inputs.txt)" ] || { echo "only $n v7h rows" >&2; exit 1; }
$P/env/td/bin/python benchmark_v2/scripts/compare_homology.py $OUT $P/data/bench_v3/inputs $P/work/bench_v3/runs \
  pytd_v7 pytd_v7h_homology > data/benchmark/homology_v7h_vs_abinitio.tsv
mkdir -p $P/_trash && mv $tmp $P/_trash/
echo "wrote $OUT ($n v7h rows) and data/benchmark/homology_v7h_vs_abinitio.tsv"
