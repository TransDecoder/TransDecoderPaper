#!/bin/bash
# E33 step 5: DIAMOND run time for one benchmark input as a user would run it (all LongOrfs peptides of
# bench.human.complete, 63,551 peptides; -k 1 --evalue 1e-5, 8 threads, nice 19), Swiss-Prot minus human vs
# UniRef50 minus human, back to back on the same (shared, loaded) machine. Also checks that the UniRef50 hits
# equal those of the pooled per-organism search (query-set independence).
set -euo pipefail
P=/broad/hptmp/bhaas/TRANSDECODER_PAPER; B=$P/data/start_classifier/step5; O=$B/timing; mkdir -p $O; cd $O
export PATH=$P/env/td/bin:$PATH
pep=$P/work/bench_v3/runs/human/complete/pytd_v7h_homology/input.fa.transdecoder_dir/longest_orfs.pep
/usr/bin/time -v nice -n 19 diamond blastp --query $pep --db $P/data/homology/sprot_minus_9606.dmnd --max-target-seqs 1 \
  --evalue 1e-5 --outfmt 6 --threads 8 --out sprot.human.complete.outfmt6 > sprot.log 2>&1
/usr/bin/time -v nice -n 19 diamond blastp --query $pep --db $B/db/uniref50_2026_03.dmnd --taxon-exclude 9606 \
  --max-target-seqs 1 --evalue 1e-5 --outfmt 6 --threads 8 --tmpdir $B/tmp --out uniref50.human.complete.outfmt6 > uniref50.log 2>&1
for f in sprot uniref50; do echo "$f $(grep -E 'Elapsed|User time|System time|Maximum resident' $f.log | tr -s ' ' | tr '\n' ' ')"; done
if cmp -s <(cut -f1-12 uniref50.human.complete.outfmt6 | sort) <(sort $B/hits/uniref50/bench.human.complete/blastp.outfmt6); then echo "pooled == single: identical"; else
  echo "pooled vs single differ: $(comm -3 <(cut -f1-12 uniref50.human.complete.outfmt6 | sort) <(sort $B/hits/uniref50/bench.human.complete/blastp.outfmt6) | wc -l) lines"; fi
gzip -f sprot.human.complete.outfmt6 uniref50.human.complete.outfmt6
