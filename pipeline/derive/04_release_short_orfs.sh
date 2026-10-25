#!/usr/bin/env bash
# Notebook 04 (short ORFs and q-values): tables from the final release benchmark runs (E25).
#
# Writes into data/benchmark/final_v7/:
#   by_length.v7.tsv            score_by_length.py on pytd_v7, pytd_v7_precise, pytd_v7_sensitive
#                               (all *.wRand / *.wLnc inputs): sensitivity by true protein length and
#                               calls by predicted length (as data/benchmark/by_length.tsv, E04)
#   lnc_by_era.human.v7.tsv     lnc_calls_by_era.py (analyses/lnc_holdout) on human complete.wLnc and
#                               fragments.wLnc: release runs and the comparators of E04/E07, tidy format
#                               (variant column added)
#   lnc_call_lengths.human.tsv  human lncRNA calls by predicted ORF length ((end-start+1)/3 codons,
#                               stop included), for the release runs and the development q sweep of
#                               E22 (work/bench_v3/rules/human/<variant>/decoy_q<Q>_longgc)
#   sens_by_length.v7.tsv       sens_by_length.py on the release runs and cliff_c200 (codon bins, E20)
#   lnc_riboseq_orfs.v7.tsv     score_lnc_riboseq_orfs.py (GENCODE phase-2 comprehensive catalogue)
#                               on the release runs, plus cliff_c200 and decoy_q0.05 as a check against
#                               data/benchmark/rules/lnc_riboseq_orfs.tsv
#
# usage: bash pipeline/derive/04_release_short_orfs.sh      (from anywhere; about 2 min, 1 core)
set -euo pipefail
P=/broad/hptmp/bhaas/TRANSDECODER_PAPER
T=$P/TransDecoderPaper
PY=$P/env/td/bin/python
S=$T/benchmark_v2/scripts
OUT=$T/data/benchmark/final_v7
INPUTS=$P/data/bench_v3/inputs
RUNS=$P/work/bench_v3/runs
RULES=$P/work/bench_v3/rules
REL="pytd_v7 pytd_v7_precise pytd_v7_sensitive"

# 1. length-stratified scores of the release runs
$PY $S/score_by_length.py $INPUTS $RUNS $REL > $OUT/by_length.v7.tsv

# 2. lncRNA call rate by GENCODE era (human)
{
  printf "variant\ttool\tera\tn_lnc\tn_called\tfrac_called\n"
  for v in complete.wLnc fragments.wLnc; do
    $PY $T/analyses/lnc_holdout/lnc_calls_by_era.py $P/work/lnc/human_lnc.pick.gtf \
      $P/data/gencode_v28/v28_gene_ids.txt $RUNS/human $v \
      $REL pytd_dev_v2 cpat cpc2 rnasamba borf gmst prodigal | awk -v v=$v 'NR>1{print v"\t"$0}'
  done
} > $OUT/lnc_by_era.human.v7.tsv

# 3. lncRNA calls by predicted ORF length (human)
{
  printf "variant\tlabel\tsource\tn_lnc_calls\tge100\tge200\tge300\n"
  for v in complete.wLnc fragments.wLnc; do
    for t in $REL; do echo -e "$t\trelease\t$RUNS/human/$v/$t/preds.tsv"; done
    for q in 0.01 0.02 0.05 0.1 0.2 0.3 0.5; do echo -e "decoy_q${q}_longgc\tdevelopment\t$RULES/human/$v/decoy_q${q}_longgc/preds.tsv"; done
  done | while IFS=$'\t' read -r lab src f; do
    v=$(basename $(dirname $(dirname $f)))
    awk -F'\t' -v v=$v -v l=$lab -v s=$src 'NR>1 && $1 ~ /^lnc_/ {n++; c=($3-$2+1)/3; if(c>=100)a++; if(c>=200)b++; if(c>=300)d++}
      END{printf "%s\t%s\t%s\t%d\t%d\t%d\t%d\n", v, l, s, n, a, b, d}' $f
  done
} > $OUT/lnc_call_lengths.human.tsv

# 4. lncRNA calls ending at a catalogued Ribo-seq ORF stop (human)
$PY $S/score_lnc_riboseq_orfs.py $P/data/riboseq_orfs/phase2/Ribo-seq_ORFs.comprehensive.fna $INPUTS human \
  complete.wLnc,fragments.wLnc \
  pytd_v7=$RUNS/{org}/{variant}/pytd_v7/preds.tsv \
  pytd_v7_precise=$RUNS/{org}/{variant}/pytd_v7_precise/preds.tsv \
  pytd_v7_sensitive=$RUNS/{org}/{variant}/pytd_v7_sensitive/preds.tsv \
  cliff_c200=$RULES/{org}/{variant}/cliff_c200/preds.tsv \
  decoy_q0.05=$RULES/{org}/{variant}/decoy_q0.05/preds.tsv \
  > $OUT/lnc_riboseq_orfs.v7.tsv

# 5. stop sensitivity by true CDS length (codons) and negatives by predicted length, release runs
#    (as data/benchmark/rules/sens_by_length.decoy.tsv, E20); cliff_c200 from the rules runs for reference
$PY $S/sens_by_length.py $INPUTS complete.wRand,fragments.wRand,complete.wLnc,fragments.wLnc \
  cliff_c200=$RULES/{org}/{variant}/cliff_c200/preds.tsv \
  pytd_v7=$RUNS/{org}/{variant}/pytd_v7/preds.tsv \
  pytd_v7_precise=$RUNS/{org}/{variant}/pytd_v7_precise/preds.tsv \
  pytd_v7_sensitive=$RUNS/{org}/{variant}/pytd_v7_sensitive/preds.tsv \
  > $OUT/sens_by_length.v7.tsv
