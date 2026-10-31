#!/usr/bin/env bash
# E39: score the TD2 runs (benchmark_v2/scripts/run_td2.sh; launcher work/bench_v3/launch_td2_v1.sh) with the
# scorers used for the other comparators, and add them to the committed tables.
#
#   1. main benchmark (84 inputs): score_v7.sh (score_predictions.py) -> rows appended to
#      data/benchmark/final_v7/all_scores.tsv (tool td2; fails if td2 rows are already there)
#   2. mixed inputs (24): score_mixed.py -> mixed_scores.tsv, five_prime_by_distance.tsv;
#      02_partial_types.py -> partial_types.tsv; 02_mixed_table.py -> mixed_medians.tsv
#   3. simulated Trinity assemblies: score_simtrinity.py per level (human d1 into
#      $P/work/simtrinity/scores_e39, then copied to data/simtrinity/ and work/simtrinity/filtered_preds/;
#      human dN and A. thaliana levels in place), then 13_depth_collect.py / 13b_depth_collect.py and
#      14_negatives_by_class.py
#   4. 04_release_short_orfs.sh (lncRNA calls by GENCODE era now include td2)
# All other rows must be unchanged: the previous tables are copied to $P/_trash/e39_pre_<time>/ and
# compared after the rerun (non-td2 rows identical, else the script stops).
#
# usage: bash pipeline/derive/10_td2_tables.sh      (about 15 min, 1 core)
set -euo pipefail
P=/broad/hptmp/bhaas/TRANSDECODER_PAPER
T=$P/TransDecoderPaper
PY=$P/env/td/bin/python
F=$T/data/benchmark/final_v7
ST=$T/pipeline/simtrinity
PRE=$P/_trash/e39_pre_$(date +%Y%m%d_%H%M%S)
mkdir -p $PRE
cd $T

# non-td2 rows of a table (tool in column $2 of the header-named column 'tool') must be unchanged
same_without_td2() {   # old new
  local c; c=$(head -1 "$1" | tr '\t' '\n' | grep -nx tool | cut -d: -f1)
  diff <(awk -F'\t' -v c=$c '$c != "td2"' "$1" | sort) <(awk -F'\t' -v c=$c '$c != "td2"' "$2" | sort) > /dev/null \
    || { echo "CHANGED non-td2 rows: $2" >&2; exit 1; }
  echo "ok (non-td2 rows unchanged): ${2#$T/}"
}

# 1. main benchmark
ORGS="athal creinhardtii dmel human mouse pfalciparum spombe tthermophila"
if awk -F'\t' '$3 == "td2"' $F/all_scores.tsv | grep -q .; then echo "td2 rows already in all_scores.tsv" >&2; exit 1; fi
cp $F/all_scores.tsv $PRE/
n=0
for org in $ORGS; do
  for v in $(awk -F'\t' -v o=$org '$1 == o && $3 == "pytd_v7" {print $2}' $PRE/all_scores.tsv | sort); do
    [ -s $P/work/bench_v3/runs/$org/$v/td2/preds.tsv ] || { echo "missing td2 run $org $v" >&2; exit 1; }
    bash $F/score_v7.sh $org $v td2 >> $F/all_scores.tsv; n=$((n + 1))
  done
done
echo "all_scores.tsv: $n td2 rows appended"
same_without_td2 $PRE/all_scores.tsv $F/all_scores.tsv

# 2. mixed inputs and partial types
for f in mixed_scores.tsv five_prime_by_distance.tsv partial_types.tsv mixed_medians.tsv; do cp $F/$f $PRE/; done
$PY pipeline/benchmark/score_mixed.py $P/data/bench_v3/inputs $P/work/bench_v3/runs
$PY pipeline/derive/02_partial_types.py $P/data/bench_v3/inputs $P/work/bench_v3/runs > $F/partial_types.tsv
$PY pipeline/derive/02_mixed_table.py > /dev/null
for f in mixed_scores.tsv five_prime_by_distance.tsv partial_types.tsv mixed_medians.tsv; do same_without_td2 $PRE/$f $F/$f; done

# 3. simulated assemblies
D1=$T/data/simtrinity
mkdir -p $PRE/simtrinity_d1 $PRE/simtrinity_depth $PRE/simtrinity_athal_depth
cp $D1/*.tsv $PRE/simtrinity_d1/; cp $D1/depth/*.tsv $PRE/simtrinity_depth/; cp $T/data/simtrinity_athal/depth/*.tsv $PRE/simtrinity_athal_depth/
$PY $ST/score_simtrinity.py $P/data/simtrinity/inputs/fragments_simtrinity.fa $P/data/simtrinity/labels/sim_trinity \
  $P/work/simtrinity/runs $P/work/simtrinity/scores_e39
for f in scores_overall.tsv scores_by_type.tsv scores_by_expr.tsv negatives.tsv excluded_preds.tsv pred_types.tsv; do
  cp $P/work/simtrinity/scores_e39/$f $D1/$f; same_without_td2 $PRE/simtrinity_d1/$f $D1/$f
done
cp $P/work/simtrinity/scores_e39/_filtered_preds/td2.* $P/work/simtrinity/filtered_preds/
for l in d2 d4 d8 d16 d32; do
  d=$P/data/simtrinity/depth/$l
  $PY $ST/score_simtrinity.py $d/inputs/fragments_simtrinity.fa $d/labels/sim_trinity $P/work/simtrinity/depth/$l/runs \
    $P/work/simtrinity/depth/$l/scores
done
$PY $ST/13_depth_collect.py $D1/depth d1 d2 d4 d8 d16 d32
for f in $PRE/simtrinity_depth/*.tsv; do
  b=$(basename $f)
  if head -1 $f | tr '\t' '\n' | grep -qx tool; then same_without_td2 $f $D1/depth/$b; else cmp $f $D1/depth/$b; echo "ok (identical): $b"; fi
done
for l in d1 d2 d4 d8 d16 d32; do
  d=$P/data/simtrinity_athal/depth/$l
  $PY $ST/score_simtrinity.py $d/inputs/fragments_simtrinity.fa $d/labels/sim_trinity \
    $P/work/simtrinity_athal/depth/$l/runs $P/work/simtrinity_athal/depth/$l/scores
done
$PY $ST/13b_depth_collect.py --data-root $P/data/simtrinity_athal --work-root $P/work/simtrinity_athal \
  $T/data/simtrinity_athal/depth d1 d2 d4 d8 d16 d32
$PY pipeline/simtrinity_athal/14_negatives_by_class.py $T/data/simtrinity_athal/depth/negatives_by_class.tsv d1 d2 d4 d8 d16 d32
for f in $PRE/simtrinity_athal_depth/*.tsv; do
  b=$(basename $f)
  if head -1 $f | tr '\t' '\n' | grep -qx tool; then same_without_td2 $f $T/data/simtrinity_athal/depth/$b; else cmp $f $T/data/simtrinity_athal/depth/$b; echo "ok (identical): $b"; fi
done

# 4. lncRNA calls by GENCODE era (and the other release short-ORF tables, which must not change)
for f in by_length.v7.tsv lnc_by_era.human.v7.tsv lnc_call_lengths.human.tsv sens_by_length.v7.tsv lnc_riboseq_orfs.v7.tsv; do cp $F/$f $PRE/; done
bash pipeline/derive/04_release_short_orfs.sh
for f in by_length.v7.tsv lnc_call_lengths.human.tsv sens_by_length.v7.tsv lnc_riboseq_orfs.v7.tsv; do cmp $PRE/$f $F/$f; echo "ok (identical): $f"; done
same_without_td2 $PRE/lnc_by_era.human.v7.tsv $F/lnc_by_era.human.v7.tsv
echo "E39_TABLES_DONE (previous tables in $PRE)"
