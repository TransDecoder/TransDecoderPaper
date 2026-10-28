#!/bin/bash
# E33 adoption: ab initio output of the PyTD homology-start-assist tree must be byte-identical to 3a793ff.
# (1) the five sample datasets (as regression/run_regression.sh): PY_v7 (3a793ff) and PY_hstart, both fresh;
# (2) four benchmark inputs: PY_hstart fresh vs the existing pytd_v7 (3a793ff) release runs.
# usage: 09_pytd_hstart_abinitio_regression_v1.sh <out_dir>
set -euo pipefail
P=/broad/hptmp/bhaas/TRANSDECODER_PAPER; T=$P/TransDecoderPaper; PY=$P/env/td/bin/python
V7=$P/work/trees/PY_v7; HS=$P/work/trees/PY_hstart; TD=$P/work/trees/TD_v7; SD=$TD/sample_data
OUT=$(mkdir -p "$1" && readlink -f "$1"); export PERL_HASH_SEED=0
echo "PY_v7 $(git -C $V7 rev-parse --short HEAD); PY_hstart $(git -C $HS rev-parse --short HEAD) $(git -C $HS status --short | wc -l) uncommitted paths"
mkdir -p $OUT/inputs; cd $OUT/inputs
gunzip -c $SD/simple_transcriptome_target/Trinity.fasta.gz > trinity.fa
gunzip -c $SD/pasa_example/pasa_assemblies.fasta.gz > pasa.fa
cp $SD/stringtie_example/stringtie_merged.transcripts.fasta stringtie.fa
cp $SD/supertranscripts_example/supertranscripts.fasta supertx.fa
gunzip -c $SD/cufflinks_example/test.genome.fasta.gz > cuff.genome.fa
gunzip -c $SD/cufflinks_example/transcripts.gtf.gz > cuff.gtf
$TD/util/gtf_genome_to_cdna_fasta.pl cuff.gtf cuff.genome.fa > cufflinks.fa
declare -A LOFLAGS=( [trinity]="" [pasa]="" [stringtie]="-S" [supertx]="" [cufflinks]="" )
FILES="input.fa.transdecoder.gff3 input.fa.transdecoder.bed input.fa.transdecoder.cds input.fa.transdecoder.pep input.fa.transdecoder_dir/best_candidates.gff3 input.fa.transdecoder_dir/best_candidates.gff3.revised_starts.gff3"
cmpall() { local a=$1 b=$2 tag=$3 n=0 d=0; for f in $FILES; do n=$((n+1)); cmp -s $a/$f $b/$f || { d=$((d+1)); echo "  DIFF $tag $f"; }; done; echo -e "$tag\tfiles_compared=$n\tdiffering=$d"; }
for ds in trinity pasa stringtie supertx cufflinks; do
  for impl in v7 hstart; do
    tree=$V7; [ $impl = hstart ] && tree=$HS
    d=$OUT/samples/$ds/$impl; mkdir -p $d; cd $d; ln -sf $OUT/inputs/$ds.fa input.fa
    PYTHONPATH=$tree $PY -m pytransdecoder.cli longorfs -t input.fa ${LOFLAGS[$ds]} > longorfs.log 2>&1
    PYTHONPATH=$tree $PY -m pytransdecoder.cli predict -t input.fa > predict.log 2>&1
  done
  cmpall $OUT/samples/$ds/v7 $OUT/samples/$ds/hstart sample:$ds
done
for oi in spombe/complete spombe/fragments dmel/complete pfalciparum/fragments_mixed_t1; do
  org=${oi%/*}; v=${oi#*/}
  gcode=$(awk -F'\t' -v o=$org '$1==o{print $4}' $T/benchmark_v2/organisms.tsv)
  d=$OUT/bench/$org.$v; mkdir -p $d; cd $d; ln -sf $P/data/bench_v3/inputs/$org/$v.fa input.fa
  PYTHONPATH=$HS $PY -m pytransdecoder.cli longorfs -t input.fa -G "$gcode" -S > longorfs.log 2>&1
  PYTHONPATH=$HS $PY -m pytransdecoder.cli predict -t input.fa -G "$gcode" > predict.log 2>&1
  cmpall $P/work/bench_v3/runs/$org/$v/pytd_v7 $d bench:$org.$v
done
echo ALLDONE
