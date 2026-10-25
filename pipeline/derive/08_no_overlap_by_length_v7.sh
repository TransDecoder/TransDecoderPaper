#!/bin/bash
# Notebook 08: share of release (final_v7, E25) ab initio predictions with no same-strand GENCODE CDS
# overlap, by predicted protein length (aa column of eval_genome_cds.py pred_status; bins as E04:
# <100 (includes the few < 50 aa), 100-149, 150-199, 200-299, >=300). Also the earlier final_abinitio
# runs (E04 table) for comparison.
# usage: bash pipeline/derive/08_no_overlap_by_length_v7.sh  -> data/gm12878/final_v7/08_no_overlap_by_length.tsv
set -euo pipefail
P=/broad/hptmp/bhaas/TRANSDECODER_PAPER; E=$P/work/gm12878_td/eval
T=$(cd "$(dirname "$0")/../.." && pwd)
{ printf "assembly\trun\tlength_bin\tn\tno_gencode_cds_overlap\n"
  for a in lraa stringtie trinity; do
    for r in final_v7_abinitio final_abinitio; do
      f=$E/${a}_$r.pred_status.tsv; [ -s $f ] || f=$E/status/${a}_$r.pred_status.tsv
      awk -F'\t' -v a=$a -v r=$r 'NR==1{for(i=1;i<=NF;i++)c[$i]=i; next}
        { L=$c["aa"]; b=(L<100)?"<100":(L<150)?"100-149":(L<200)?"150-199":(L<300)?"200-299":">=300";
          n[b]++; o[b]+=$c["no_gencode_cds_overlap"] }
        END{ split("<100 100-149 150-199 200-299 >=300",B," "); for(i=1;i<=5;i++) printf "%s\t%s\t%s\t%d\t%d\n",a,r,B[i],n[B[i]],o[B[i]] }' $f
    done
  done; } > $T/data/gm12878/final_v7/08_no_overlap_by_length.tsv
echo "wrote $T/data/gm12878/final_v7/08_no_overlap_by_length.tsv"
