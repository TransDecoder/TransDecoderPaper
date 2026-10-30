#!/bin/bash
# E06 release rerun (PyTD a5b24c3): collect the GM12878 homology series as E06 did.
#   data/gm12878/gm12878_homology_series_v7h.tsv            (same columns as gm12878_homology_series.tsv)
#   data/gm12878/eval_summaries/hseries_v7h/*.summary.tsv
#   data/gm12878/gm12878_homology_primate_retention_v7h.txt  (homology_series_retention.py, minus human vs
#                                                             minus primates, BLAST only)
# The development-code files (48c870c) are kept unchanged for comparison.
# usage (from the paper repo root): bash pipeline/derive/08_homology_series_v7h.sh
set -euo pipefail
P=/broad/hptmp/bhaas/TRANSDECODER_PAPER; W=$P/work/gm12878_td; E=$W/eval/hseries_v7h
PY=$P/env/td/bin/python
mkdir -p data/gm12878/eval_summaries/hseries_v7h
CFG="abinitio blast_9606_nopfam blast_9443_nopfam blast_40674_nopfam blast_7742_nopfam blast_33208_nopfam pfam
     blast_9606_pfam blast_9443_pfam blast_40674_pfam blast_7742_pfam blast_33208_pfam"
COLS="distinct_predicted_cds MANE_sens_exact MANE_sens_compatible MANE_sens_exact_given_assembled
      GENCODE_complete_sens_compatible pred_complete_GENCODE_exact_frac pred_complete_no_GENCODE_CDS_overlap_frac"
out=data/gm12878/gm12878_homology_series_v7h.tsv
echo -e "assembly\tconfig\t$(echo $COLS | tr ' ' '\t')" > $out
for a in lraa stringtie trinity; do
  for c in $CFG; do
    [ "$(cat $W/$a/hseries_v7h/$c/pytd_version.txt)" = a5b24c3 ] || { echo "not a5b24c3: $a $c" >&2; exit 1; }
    s=$E/${a}_$c.summary.tsv
    cp $s data/gm12878/eval_summaries/hseries_v7h/
    row="$a\t$c"
    for k in $COLS; do row="$row\t$(awk -v k=$k '$1==k{print $2}' $s)"; done
    echo -e "$row" >> $out
  done
done
ret=data/gm12878/gm12878_homology_primate_retention_v7h.txt
: > $ret
for a in trinity stringtie lraa; do
  echo "== $a" >> $ret
  $PY analyses/homology_series_retention.py $E/${a}_blast_9606_nopfam.pred_status.tsv \
      $E/${a}_blast_9443_nopfam.pred_status.tsv $W/$a/hseries_v7h/blastp.minus_9606.outfmt6 >> $ret
done
echo "wrote $out $ret"
