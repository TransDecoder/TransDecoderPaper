#!/bin/bash
# E33 step 5: DIAMOND 2.2.8 blastp (default = fast mode) of the step-5 query sets against UniRef50 2026_03,
# excluding the organism's taxon and its NCBI descendants (--taxon-exclude, on the representative's taxid).
# One search at a time, nice 19, 8 threads.
#   k25: cand.<org>.pep (5'-open candidate peptides)  --max-target-seqs 25 --max-hsps 0 --evalue 1e-5
#   k1 : k1.<org>.pep (full-mode LongOrfs + candidates) --max-target-seqs 1 --evalue 1e-5 (as run_tool.sh)
# usage: step5_04_uniref50_searches_v1.sh [org ...]   (skips finished searches)
set -uo pipefail
P=/broad/hptmp/bhaas/TRANSDECODER_PAPER; B=$P/data/start_classifier/step5
DB=$B/db/uniref50_2026_03.dmnd; Q=$B/queries; O=$B/searches/uniref50; mkdir -p $O $B/tmp
export PATH=$P/env/td/bin:$PATH
COLS="qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen staxids"
declare -A TAX=([human]=9606 [mouse]=10090 [athal]=3702 [dmel]=7227 [spombe]=284812 [creinhardtii]=3055 [pfalciparum]=36329 [tthermophila]=312017)
orgs=("$@"); [ ${#orgs[@]} -eq 0 ] && orgs=(human athal dmel tthermophila mouse spombe pfalciparum creinhardtii)
one() { # kind org
  local kind=$1 org=$2 out=$O/$1.$2
  [ -e $out.done ] && return 0
  local k=(--max-target-seqs 1); [ $kind = k25 ] && k=(--max-target-seqs 25 --max-hsps 0)
  echo "$(date +%T) start $kind $org: $(grep -c '>' $Q/$kind.$org.pep) queries"
  /usr/bin/time -v nice -n 19 diamond blastp --query $Q/$kind.$org.pep --db $DB --taxon-exclude ${TAX[$org]} \
    "${k[@]}" --evalue 1e-5 --outfmt 6 $COLS --threads 8 --tmpdir $B/tmp --out $out.tsv > $out.log 2>&1
  local rc=$?
  if [ $rc = 0 ]; then gzip -f $out.tsv && touch $out.done; fi
  echo "$(date +%T) end $kind $org rc=$rc $(grep -E 'Elapsed \(wall|Maximum resident' $out.log | tr -s ' ' | tr '\n' ' ')"
}
for org in "${orgs[@]}"; do
  one k25 $org
  one k1 $org
done
echo "$(date +%T) ALLDONE"
