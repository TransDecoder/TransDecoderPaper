#!/bin/bash
# E06 release rerun (PyTD a5b24c3, PY_hstart): per-assembly preparation for the GM12878 homology series.
#   1. fresh LongOrfs (-S) in $W/<asm>/hseries_v7h/_longorfs; longest_orfs.pep must equal the release
#      ab initio run's (final_v7_abinitio, 3a793ff);
#   2. fresh DIAMOND blastp of all LongOrfs peptides against Swiss-Prot minus each clade
#      (--max-target-seqs 1 --evalue 1e-5 --outfmt 6, as E06 and gm12878_final_v7.sh);
#   3. Pfam: the existing hmmsearch --cut_ga domtblout on the identical peptide file is copied (checked
#      by the md5 of longest_orfs.pep recorded here; see E06 "Release rerun" for the verification).
# usage: homology_series_v7h_prep.sh <lraa|stringtie|trinity> [threads]
set -euo pipefail
name=$1; threads=${2:-4}
P=/broad/hptmp/bhaas/TRANSDECODER_PAPER; D=$P/data/gm12878; W=$P/work/gm12878_td; H=$P/data/homology
PY=$P/env/td/bin/python; export PATH=$P/env/td/bin:$PATH; export PYTHONPATH=$P/work/trees/PY_hstart
out=$W/$name/hseries_v7h
[ -e $out/_longorfs ] && { echo "exists: $out/_longorfs" >&2; exit 1; }
mkdir -p $out/_longorfs; cd $out/_longorfs
(cd $PYTHONPATH && git rev-parse --short HEAD) > pytd_version.txt
[ "$(cat pytd_version.txt)" = a5b24c3 ] || { echo "PyTD is not a5b24c3" >&2; exit 1; }
ln -s $(readlink -f $W/$name/final_v7_abinitio/transcripts.fa) transcripts.fa
$PY -m pytransdecoder.cli longorfs -t transcripts.fa -S > longorfs.log 2>&1
pep=$out/_longorfs/transcripts.fa.transdecoder_dir/longest_orfs.pep
cmp $pep $W/$name/final_v7_abinitio/transcripts.fa.transdecoder_dir/longest_orfs.pep
md5sum $pep > $out/longest_orfs.pep.md5
for c in 9606 9443 40674 7742 33208; do
  db=$H/sprot_minus_clade_$c.dmnd; [ $c = 9606 ] && db=$H/sprot_minus_9606.dmnd
  diamond blastp --query $pep --db $db --max-target-seqs 1 --evalue 1e-5 --outfmt 6 --threads $threads \
     --out $out/blastp.minus_$c.outfmt6 > $out/blastp.minus_$c.log 2>&1
done
src=$(readlink -f $W/$name/pfam.domtblout)
[ $name = trinity ] && src=$W/trinity/final_v7_homology/pfam.domtblout
cp $src $out/pfam.domtblout; echo "copied from $src" > $out/pfam.domtblout.source
echo "PREP_DONE $name"
