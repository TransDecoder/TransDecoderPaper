#!/usr/bin/env bash
# Score the benchmark sweeps at the pinned release code (E40). Writes data/benchmark/release_pin/:
#   negatives_scores.tsv, negatives_medians.tsv   score_negatives.py: q settings, q sweep, no q-value test, long-ORF rule off,
#                                                 longest-ORF baselines (rel LongOrfs) and, as the historical comparator, cliff_c200
#   by_length.tsv, short_orf_summary.tsv          score_by_length.py on m*_noq, rel, precise, sensitive, q*
#   sens_by_length.tsv                            sens_by_length.py (codon bins)
#   lnc_call_lengths.human.tsv                    human lncRNA calls by predicted ORF length
#   lnc_by_era.human.tsv                          as data/benchmark/final_v7/lnc_by_era.human.v7.tsv, q sweep runs not included
set -euo pipefail
P=/broad/hptmp/bhaas/TRANSDECODER_PAPER; T=$P/TransDecoderPaper; PY=$P/env/td/bin/python; S=$T/benchmark_v2/scripts
OUT=$T/data/benchmark/release_pin; mkdir -p $OUT
IN=$P/data/bench_v3/inputs; R=$P/work/release_pin/bench; RULES=$P/work/bench_v3/rules
V=complete.wRand,fragments.wRand,complete.wLnc,fragments.wLnc
LAB="rel precise sensitive q0.02 q0.1 q0.3 q0.5 nolong_q0.05 nolong_q0.1 m100_noq m90_noq m80_noq m70_noq m60_noq m50_noq"
args=(); for l in $LAB; do args+=("$l=$R/{org}/{variant}/$l/preds.tsv"); done
# longest ORF per transcript (LongOrfs of the release run): >= 50 codons, and >= 100 codons (>= 303 nt)
for o in $(ls $R); do for v in ${V//,/ }; do d=$R/$o/$v/rel; [ -s $d/longest_orf.preds.tsv ] || continue
  awk -F'\t' 'NR==1 || $3-$2+1 >= 303' $d/longest_orf.preds.tsv > $d/longest_orf_m100.preds.tsv; done; done
$PY $S/score_negatives.py $IN $V "${args[@]}" longest_orf_m50=$R/{org}/{variant}/rel/longest_orf.preds.tsv \
  longest_orf_m100=$R/{org}/{variant}/rel/longest_orf_m100.preds.tsv cliff_c200=$RULES/{org}/{variant}/cliff_c200/preds.tsv > $OUT/negatives_scores.tsv
$PY - $OUT <<'PY'
import sys, pandas as pd
o = sys.argv[1]; d = pd.read_csv(f"{o}/negatives_scores.tsv", sep="\t")
d["fdr"] = 1 - d.precision_stop
m = d.groupby(["variant", "label"]).agg(sens_stop=("sens_stop", "median"), sens_exact=("sens_exact", "median"),
    precision_stop=("precision_stop", "median"), neg_rate=("neg_rate", "median"), n_orgs=("org", "size")).reset_index()
m.to_csv(f"{o}/negatives_medians.tsv", sep="\t", index=False, float_format="%.4f")
PY
$PY $S/score_by_length.py $IN $R rel precise sensitive q0.02 q0.1 q0.3 q0.5 m100_noq m90_noq m80_noq m70_noq m60_noq m50_noq > $OUT/by_length.tsv
awk -F'\t' 'NR==1{print; next} $2 ~ /\.(wRand|wLnc)$/ && $2 !~ /unstranded/' $OUT/by_length.tsv > $OUT/by_length.tmp && mv $OUT/by_length.tmp $OUT/by_length.tsv
{ printf "org\tv\tcfg\tsens_all\tsens_short\tneg\tfp_real\tprec\n"; awk -F'\t' 'NR>1{
  k=$1"\t"$2"\t"$3; K[k]=1
  if ($4=="truth"){tn[k]+=$6; ts[k]+=$7; if($5~/^[5-9]0-/){sn[k]+=$6; ss[k]+=$7}} else {tp[k]+=$7; fr[k]+=$9; fx[k]+=$10}}
  END{for(k in K) printf "%s\t%.17g\t%.17g\t%d\t%d\t%.17g\n",k,ts[k]/tn[k],(sn[k]?ss[k]/sn[k]:0),fx[k],fr[k],tp[k]/(tp[k]+fr[k]+fx[k])}' $OUT/by_length.tsv | sort; } > $OUT/short_orf_summary.tsv
args2=(); for l in rel precise sensitive q0.1 nolong_q0.05 nolong_q0.1; do args2+=("$l=$R/{org}/{variant}/$l/preds.tsv"); done
$PY $S/sens_by_length.py $IN $V cliff_c200=$RULES/{org}/{variant}/cliff_c200/preds.tsv "${args2[@]}" > $OUT/sens_by_length.tsv
{ printf "variant\tlabel\tn_lnc_calls\tge100\tge200\tge300\n"
  for v in complete.wLnc fragments.wLnc; do for l in rel precise sensitive q0.02 q0.1 q0.3 q0.5; do
    awk -F'\t' -v v=$v -v l=$l 'NR>1 && $1 ~ /^lnc_/ {n++; c=($3-$2+1)/3; if(c>=100)a++; if(c>=200)b++; if(c>=300)d++}
      END{printf "%s\t%s\t%d\t%d\t%d\t%d\n", v, l, n, a, b, d}' $R/human/$v/$l/preds.tsv; done; done
} > $OUT/lnc_call_lengths.human.tsv
# lncRNA calls at a catalogued Ribo-seq ORF stop (human), all q settings
args3=(); for l in rel precise sensitive q0.02 q0.1 q0.3 q0.5; do args3+=("$l=$R/{org}/{variant}/$l/preds.tsv"); done
$PY $S/score_lnc_riboseq_orfs.py $P/data/riboseq_orfs/phase2/Ribo-seq_ORFs.comprehensive.fna $IN human complete.wLnc,fragments.wLnc "${args3[@]}" > $OUT/lnc_riboseq_orfs.tsv
