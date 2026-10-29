#!/bin/bash
# E33 step 5: Swiss-Prot-side searches that step 5 needs and earlier steps lacked (serial, niced).
#  1. sim.athal.d16 candidates: DIAMOND fast -k 25 --max-hsps 0 --evalue 10 vs sprot_minus_3702 (as step4_01 v3)
#  2. sim.human.d8 full mode: DIAMOND -k 1 --evalue 1e-5 vs sprot_minus_9606 on ALL LongOrfs peptides, and
#     hmmsearch --cut_ga --noali vs Pfam-A in 7 chunks (--cpu 0 each), as run_tool.sh pytd_v7_homology
# usage: step5_03_sprot_searches_v1.sh <threads>
set -euo pipefail
th=${1:-8}
P=/broad/hptmp/bhaas/TRANSDECODER_PAPER; H=$P/data/homology; B=$P/data/start_classifier/step5
export PATH=$P/env/td/bin:$PATH
COLS="qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen"
O=$B/searches/sprot/sim.athal.d16; mkdir -p $O; cd $O
if [ ! -e blastp.fast.done ]; then
  zcat $B/queries/cand_map.tsv.gz | awk -F'\t' '$1=="sim.athal.d16"{print $3}' > cand_ids.txt
  awk 'NR==FNR{k[$1]=1; next} /^>/{split(substr($0,2),a," "); p=(a[1] in k)} p' cand_ids.txt \
    $P/work/simtrinity_athal/depth/d16/runs/pytd_v7/input.fa.transdecoder_dir/longest_orfs.pep > cand.pep
  echo "$(date +%T) athal d16 cand: $(grep -c '>' cand.pep)"
  /usr/bin/time -v nice -n 19 diamond blastp --query cand.pep --db $H/sprot_minus_3702.dmnd --max-target-seqs 25 \
    --max-hsps 0 --evalue 10 --outfmt 6 $COLS --threads $th --out blastp.fast.k25.e10.tsv > diamond.fast.log 2>&1
  gzip -f blastp.fast.k25.e10.tsv; touch blastp.fast.done
fi
O=$B/searches/sprot/sim.human.d8.full; mkdir -p $O; cd $O
pep=$P/work/start_classifier/hstart_runs/sim.human.d8/on/input.fa.transdecoder_dir/longest_orfs.pep
if [ ! -e blastp.done ]; then
  /usr/bin/time -v nice -n 19 diamond blastp --query $pep --db $H/sprot_minus_9606.dmnd --max-target-seqs 1 \
    --evalue 1e-5 --outfmt 6 --threads $th --out blastp.outfmt6 > diamond.log 2>&1
  touch blastp.done; echo "$(date +%T) human d8 sprot: $(wc -l < blastp.outfmt6) lines"
fi
if [ ! -e pfam.done ]; then
  mkdir -p pfam_chunks
  awk -v d=pfam_chunks -v n=7 '/^>/{f=sprintf("%s/c%03d.pep",d,(i++)%n)} {print > f}' $pep
  ls pfam_chunks/*.pep | nice -n 19 xargs -P 7 -I{} nice -n 19 hmmsearch --cpu 0 --cut_ga --noali -o /dev/null \
    --domtblout {}.domtblout "$H/Pfam-A.hmm" {} > hmmsearch.log 2>&1
  cat pfam_chunks/*.domtblout > pfam.domtblout; touch pfam.done
  echo "$(date +%T) human d8 pfam: $(grep -vc '^#' pfam.domtblout) lines"
fi
echo "$(date +%T) SPROT_DONE"
