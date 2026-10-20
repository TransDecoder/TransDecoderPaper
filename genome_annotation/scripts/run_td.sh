#!/bin/bash
# usage: run_td.sh <name> <transcripts.fa> <structures (gtf|bam)> <abinitio|homology> [threads]
# PyTransDecoder (dev tree) on assembled GM12878 transcripts, then genome-map to genome CDS models.
# homology: DIAMOND vs Swiss-Prot without human entries + hmmsearch vs Pfam-A (--cut_ga).
set -euo pipefail
name=$1; fa=$(readlink -f $2); structs=$(readlink -f $3); mode=$4; threads=${5:-8}
P=/broad/hptmp/bhaas/TRANSDECODER_PAPER
PY=$P/env/td/bin/python; export PATH=$P/env/td/bin:$PATH
export PYTHONPATH=$P/work/trees/PY_dev
out=$P/work/gm12878_td/$name/$mode; mkdir -p $out; cd $out
ln -sf $fa transcripts.fa
td() { $PY -m pytransdecoder.cli "$@"; }
[ -s transcripts.fa.transdecoder_dir/longest_orfs.pep ] || td longorfs -t transcripts.fa -S > longorfs.log 2>&1
extra=()
if [ $mode = homology ]; then
  pep=transcripts.fa.transdecoder_dir/longest_orfs.pep
  [ -s blastp.outfmt6 ] || diamond blastp --query $pep --db $P/data/homology/sprot_minus_9606.dmnd \
      --max-target-seqs 1 --evalue 1e-5 --outfmt 6 --threads $threads --out blastp.outfmt6 > diamond.log 2>&1
  [ -s pfam.domtblout ] || hmmsearch --cpu $threads --cut_ga --noali -o /dev/null --domtblout pfam.domtblout \
      $P/data/homology/Pfam-A.hmm $pep > hmmsearch.log 2>&1
  extra=(--retain-blastp-hits blastp.outfmt6 --retain-pfam-hits pfam.domtblout)
fi
td predict -t transcripts.fa ${extra[@]+"${extra[@]}"} > predict.log 2>&1
td genome-map --orfs transcripts.fa.transdecoder.gff3 --structures $structs -t transcripts.fa \
   -o transcripts.fa.transdecoder.genome.gff3 > genome_map.log 2>&1
echo DONE $name $mode
