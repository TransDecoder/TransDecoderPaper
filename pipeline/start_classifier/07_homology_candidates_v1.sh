#!/bin/bash
# E33 step 3: DIAMOND + hmmsearch (parameters exactly as run_tool.sh pytd_v7_homology) on the
# longest-ORF peptides of the 5'-open start-decision candidates of existing pytd_v7 runs.
# usage: 07_homology_candidates_v1.sh <dataset> <run_dir> <taxid> <outdir>
# Per-query results do not depend on the other queries (DIAMOND -k 1; hmmsearch --cut_ga thresholds are
# bit scores), so searching the candidate subset gives the same hits for those ORFs as the full run
# (checked on bench.human.fragments against its pytd_v7_homology run).
set -euo pipefail
P=/broad/hptmp/bhaas/TRANSDECODER_PAPER; T=$P/TransDecoderPaper; H=$P/data/homology
export PATH=$P/env/td/bin:$PATH
THREADS=3
ds=$1; run=$2; taxid=$3; out=$4
mkdir -p $out; cd $out
zcat $T/data/start_classifier/features/$ds.tsv.gz | awk -F'\t' 'NR>1{print $3}' | sort -u > cand_ids.txt
awk 'NR==FNR{k[$1]=1; next} /^>/{split(substr($0,2),a," "); p=(a[1] in k)} p' cand_ids.txt \
  $run/input.fa.transdecoder_dir/longest_orfs.pep > cand.pep
echo "$ds: $(wc -l < cand_ids.txt) candidates, $(grep -c '>' cand.pep) peptides"
diamond blastp --query cand.pep --db "$H/sprot_minus_${taxid}.dmnd" --max-target-seqs 1 \
  --evalue 1e-5 --outfmt 6 --threads "$THREADS" --out blastp.outfmt6 > diamond.log 2>&1
mkdir -p pfam_chunks
awk -v d=pfam_chunks -v n="$THREADS" '/^>/{f=sprintf("%s/c%03d.pep",d,(i++)%n)} {print > f}' cand.pep
ls pfam_chunks/*.pep | xargs -P "$THREADS" -I{} hmmsearch --cpu 1 --cut_ga --noali -o /dev/null \
  --domtblout {}.domtblout "$H/Pfam-A.hmm" {} > hmmsearch.log 2>&1
cat pfam_chunks/*.domtblout > pfam.domtblout
echo "$ds: done $(wc -l < blastp.outfmt6) blast lines, $(grep -vc '^#' pfam.domtblout) pfam lines"
