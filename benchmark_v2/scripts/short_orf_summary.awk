# Summarise by_length.tsv (score_by_length.py) per organism x variant x configuration.
# RECONSTRUCTED 2026-10-21: short_orf_summary.tsv was made by a one-off command; this script
# reproduces its rows and values from by_length.tsv (identical as doubles; the original was
# written by Python, so the decimal text of some values differs in the last digits).
# usage: awk -F'\t' -f short_orf_summary.awk by_length.tsv | (read h; echo "$h"; sort) > short_orf_summary.tsv
# cfg: pytd_dev_kozak_add -> m100 (LongOrfs default -m 100 at the time); pytd_dev_m<N> -> m<N>;
#      pytd_dev_m50_sc<C> -> m50+g<C> (Predict --short-orf-min-contrast C)
# sens_all   truth stop-level sensitivity, all lengths
# sens_short truth stop-level sensitivity, true proteins of 50-99 aa
# neg        calls on negatives (rand_* or lnc_* sequences)
# fp_real    calls on real transcripts that do not match the true stop
# prec       TP / (TP + fp_real + neg), TP = calls matching the true stop
NR > 1 {
  c = $3; sub(/^pytd_dev_/, "", c); if (c == "kozak_add") c = "m100"; sub(/_sc/, "+g", c)
  k = $1 "\t" $2 "\t" c; K[k] = 1
  if ($4 == "truth") { tn[k] += $6; ts[k] += $7; if ($5 ~ /^[5-9]0-/) { sn[k] += $6; ss[k] += $7 } }
  else { tp[k] += $7; fr[k] += $9; fx[k] += $10 }
}
END {
  print "org\tv\tcfg\tsens_all\tsens_short\tneg\tfp_real\tprec"
  for (k in K) printf "%s\t%.17g\t%.17g\t%d\t%d\t%.17g\n", k, ts[k]/tn[k], ss[k]/sn[k], fx[k], fr[k], tp[k]/(tp[k]+fr[k]+fx[k])
}
