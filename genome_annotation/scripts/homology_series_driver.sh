#!/bin/bash
# Runs the full GM12878 homology-restriction series (see homology_series.sh). Log: logs/hseries.out
P=/broad/hptmp/bhaas/TRANSDECODER_PAPER; D=$P/data/gm12878; W=$P/work/gm12878_td; H=$P/data/homology
export PATH=$P/env/td/bin:$PATH
declare -A S=([lraa]=$D/longread_lraa/lraa_denovo/LRAA.gtf [stringtie]=$D/star_stringtie/stringtie.gtf
              [trinity]=$D/trinity/Trinity.minimap2.bam)
CLADES="9606 9443 40674 7742 33208"
mkdir -p $W/logs/hseries
job() { bash $W/homology_series.sh "$@" > $W/logs/hseries/$1.$3.$4.log 2>&1 && echo "ok $*" || echo "FAILED $*"; }
export -f job; export W

# 1. blast-only, all assemblies (DIAMOND per clade is shared via a lock)
for n in lraa stringtie trinity; do for c in $CLADES; do echo "$n ${S[$n]} $c nopfam"; done; done |
  xargs -P 5 -L 1 bash -c 'job "$@"' _ &

# 2. chunked Pfam for LRAA and StringTie (single-threaded hmmsearch per chunk)
for n in lraa stringtie; do
  c=$W/$n/pfam_chunks; mkdir -p $c
  [ -s $W/$n/pfam.domtblout ] && continue
  awk -v d=$c 'BEGIN{n=0} /^>/{f=sprintf("%s/chunk%02d.pep",d,(n++)%12)} {print > f}' \
      $W/$n/final_abinitio/transcripts.fa.transdecoder_dir/longest_orfs.pep
  ls $c/chunk*.pep | xargs -P 12 -I{} hmmsearch --cpu 1 --cut_ga --noali -o /dev/null --domtblout {}.domtblout $H/Pfam-A.hmm {}
  cat $c/chunk*.pep.domtblout > $W/$n/pfam.domtblout
  echo "pfam done $n"
done
for n in lraa stringtie; do echo "$n ${S[$n]} none pfam"; for c in $CLADES; do echo "$n ${S[$n]} $c pfam"; done; done |
  xargs -P 4 -L 1 bash -c 'job "$@"' _
wait

# 3. Trinity Pfam comes from the running final_homology hmmsearch
T=$W/trinity/final_homology/pfam.domtblout
until [ -s $T ] && ! pgrep -f "domtblout pfam.domtblout $H/Pfam-A.hmm" >/dev/null; do sleep 300; done
tail -1 $T | grep -q '^#' || { echo "FAILED trinity pfam incomplete"; echo HSERIES_DONE; exit 1; }
ln -sf $T $W/trinity/pfam.domtblout
{ echo "trinity ${S[trinity]} none pfam"; for c in $CLADES; do echo "trinity ${S[trinity]} $c pfam"; done; } |
  xargs -P 3 -L 1 bash -c 'job "$@"' _
echo HSERIES_DONE
