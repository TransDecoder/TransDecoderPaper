#!/bin/bash
# Notebook 08, Fig. 7a release points: exact-start share of GENCODE-compatible complete predictions
# (complete & gencode_exact / complete & gencode_compatible; compatible includes exact), as in E13's
# pred_complete_start_precision, for the release runs (final_v7, E25) and the earlier final_abinitio runs.
# Counts from eval_genome_cds.py pred_status (distinct predicted CDSs).
# usage: bash pipeline/derive/08_start_share_v7.sh -> data/gm12878/final_v7/08_start_share.tsv
set -euo pipefail
P=/broad/hptmp/bhaas/TRANSDECODER_PAPER; E=$P/work/gm12878_td/eval
T=$(cd "$(dirname "$0")/../.." && pwd)
{ printf "assembly\trun\tn_complete\tcomplete_gencode_exact\tcomplete_gencode_compatible\n"
  for a in lraa stringtie trinity; do
    for r in final_abinitio final_v7_abinitio final_v7_abinitio_precise final_v7_abinitio_sensitive; do
      f=$E/${a}_$r.pred_status.tsv; [ -s $f ] || f=$E/status/${a}_$r.pred_status.tsv; [ -s $f ] || continue
      awk -F'\t' -v a=$a -v r=$r 'NR==1{for(i=1;i<=NF;i++)c[$i]=i; next}
        $c["orf_type"]=="complete"{n++; e+=$c["gencode_exact"]; m+=$c["gencode_compatible"]}
        END{printf "%s\t%s\t%d\t%d\t%d\n",a,r,n,e,m}' $f
    done
  done; } > $T/data/gm12878/final_v7/08_start_share.tsv
echo "wrote $T/data/gm12878/final_v7/08_start_share.tsv"
