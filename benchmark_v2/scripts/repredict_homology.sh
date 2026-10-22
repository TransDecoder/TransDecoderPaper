#!/bin/bash
# Re-run only PyTD Predict (and the conversion) in existing pytd_dev_v2_homology run dirs with a given PyTD tree,
# reusing their LongOrfs, DIAMOND and hmmsearch outputs. Used once (2026-10-22): the homology rerun (launch17)
# picked up PyTD 44e417a (Kozak in bits) for its last 4 runs; these were re-predicted with 48c870c so the whole
# set uses the same code as the ab initio pytd_dev_v2 runs. Previous Predict outputs are moved to predict_44e417a/.
# usage: repredict_homology.sh <pytd_tree> <run_dir> [...]
set -euo pipefail
P=/broad/hptmp/bhaas/TRANSDECODER_PAPER; PYTHON=$P/env/td/bin/python; here=$(dirname "$(readlink -f "$0")")
tree=$1; shift
for d in "$@"; do
  org=$(basename "$(dirname "$(dirname "$d")")"); variant=$(basename "$(dirname "$d")")
  gcode=$(awk -F'\t' -v o="$org" '$1==o{print $4}' "$here/../organisms.tsv")
  cd "$d"; mkdir -p predict_44e417a
  for f in predict.log preds.tsv input.fa.transdecoder.{gff3,bed,cds,pep} input.fa.transdecoder_dir/__checkpoints_TDpredict \
           input.fa.transdecoder_dir/{longest_orfs.cds.scores,longest_orfs.cds.top_longest_5000,longest_orfs.cds.top_500_longest,hexamer.scores,best_candidates.gff3,best_candidates.gff3.revised_starts.gff3} \
           input.fa.transdecoder_dir/start_refinement.kozak.{tsv,counts.tsv,logo.svg}; do
    [ -e "$f" ] && mv "$f" predict_44e417a/ || true
  done
  PYTHONPATH=$tree $PYTHON -m pytransdecoder.cli predict -t input.fa -G "$gcode" \
    --retain-blastp-hits blastp.outfmt6 --retain-pfam-hits pfam.domtblout > predict.log 2>&1
  $PYTHON "$here/convert_predictions.py" transdecoder_gff3 input.fa.transdecoder.gff3 > preds.tsv
  echo "repredicted $org $variant ($(tail -n +2 preds.tsv | wc -l) ORFs)"
done
