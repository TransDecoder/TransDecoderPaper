#!/bin/bash
# E33 step 4: deeper / looser homology searches of the 5'-open start-decision candidates.
# One search at a time, niced, <= 8 threads:
#   DIAMOND blastp default (fast) mode   -k 25 --max-hsps 0 (all HSPs) --evalue 10
#   DIAMOND blastp --more-sensitive      -k 25 --max-hsps 0 --evalue 10
#   hmmsearch --domE 1 (instead of --cut_ga) against Pfam-A, --cpu 7 (one process)
# Queries: the longest-ORF peptides of the candidates (01_features.py tables), as step 3
# (07_homology_candidates_v2.sh); sim datasets reuse the step-3 cand.pep.
# usage: step4_01_deep_search_v1.sh   (datasets listed below; skips finished searches)
set -uo pipefail
P=/broad/hptmp/bhaas/TRANSDECODER_PAPER; H=$P/data/homology; O=$P/data/start_classifier/step4/searches
export PATH=$P/env/td/bin:$PATH
COLS="qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen"
run_ds() {  # dataset run_dir taxid
  local ds=$1 run=$2 taxid=$3 d=$O/$1
  mkdir -p $d; cd $d || return 1
  if [ ! -s cand.pep ]; then
    zcat $P/data/start_classifier/features/$ds.tsv.gz | awk -F'\t' 'NR>1{print $3}' | sort -u > cand_ids.txt
    awk 'NR==FNR{k[$1]=1; next} /^>/{split(substr($0,2),a," "); p=(a[1] in k)} p' cand_ids.txt \
      $run/input.fa.transdecoder_dir/longest_orfs.pep > cand.pep.tmp && mv cand.pep.tmp cand.pep
  fi
  echo "$(date +%T) $ds: $(grep -c '>' cand.pep) peptides"
  for mode in fast msens; do
    [ -e blastp.$mode.done ] && continue
    sens=(); [ $mode = msens ] && sens=(--more-sensitive)
    nice -n 19 diamond blastp --query cand.pep --db "$H/sprot_minus_${taxid}.dmnd" ${sens[@]+"${sens[@]}"} \
      --max-target-seqs 25 --max-hsps 0 --evalue 10 --outfmt 6 $COLS --threads 8 \
      --out blastp.$mode.k25.e10.tsv > diamond.$mode.log 2>&1 && gzip -f blastp.$mode.k25.e10.tsv && touch blastp.$mode.done
    echo "$(date +%T) $ds: diamond $mode rc=$?"
  done
  if [ ! -e pfam.domE1.done ]; then
    nice -n 19 hmmsearch --cpu 7 --domE 1 --noali -o /dev/null --domtblout pfam.domE1.domtblout \
      "$H/Pfam-A.hmm" cand.pep > hmmsearch.log 2>&1 && gzip -f pfam.domE1.domtblout && touch pfam.domE1.done
    echo "$(date +%T) $ds: hmmsearch rc=$?"
  fi
}
SIM=$P/work/start_classifier/homology
for ds in sim.human.d1 sim.human.d8 sim.human.d32 sim.athal.d1; do
  mkdir -p $O/$ds; [ -s $O/$ds/cand.pep ] || cp $SIM/$ds/cand.pep $SIM/$ds/cand_ids.txt $O/$ds/
done
run_ds sim.human.d1 - 9606
run_ds sim.human.d8 - 9606
run_ds sim.human.d32 - 9606
run_ds sim.athal.d1 - 3702
R=$P/work/bench_v3/runs
for org in human:9606 athal:3702 dmel:7227 tthermophila:312017 mouse:10090 spombe:284812 pfalciparum:36329 creinhardtii:3055; do
  o=${org%:*}; tx=${org#*:}
  for v in complete fragments; do run_ds bench.$o.$v $R/$o/$v/pytd_v7 $tx; done
done
echo "$(date +%T) ALLDONE"
