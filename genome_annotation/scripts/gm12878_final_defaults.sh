#!/bin/bash
# GM12878, final PyTD defaults (PY_dev 48c870c: LongOrfs -m 50 + short-ORF guard, Kozak start rule):
# PyTD -> genome-map -> MANE/GENCODE eval -> Ribo-seq match, for one assembly.
# usage: final_defaults.sh <name> <transcripts.fa> <structures gtf|bam> <abinitio|homology> [threads]
set -euo pipefail
name=$1; fa=$(readlink -f $2); structs=$(readlink -f $3); mode=$4; threads=${5:-8}
P=/broad/hptmp/bhaas/TRANSDECODER_PAPER; D=$P/data/gm12878
PY=$P/env/td/bin/python; export PATH=$P/env/td/bin:$PATH; export PYTHONPATH=$P/work/trees/PY_dev
out=$P/work/gm12878_td/$name/final_$mode; mkdir -p $out; cd $out
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
$PY $P/TransDecoderPaper/genome_annotation/scripts/eval_genome_cds.py --pred transcripts.fa.transdecoder.genome.gff3 \
   --structures $structs --mane $P/data/annot/human/MANE.GRCh38.v1.5.ensembl_genomic.gtf.gz \
   --gencode $D/ref/gencode.v50.primary_assembly.annotation.gtf --label ${name}_final_$mode \
   -o $P/work/gm12878_td/eval/${name}_final_$mode > eval.log 2>&1
PYTHONPATH=$P/TransDecoderPaper/genome_annotation/scripts:$PYTHONPATH $PY $P/TransDecoderPaper/analyses/riboseq/match_riboseq.py \
   --pred transcripts.fa.transdecoder.genome.gff3 --phase1 $P/data/riboseq_orfs/phase1/Ribo-seq_ORFs.all.bed \
   --phase2 $P/data/riboseq_orfs/phase2/Ribo-seq_ORFs.comprehensive.bed > $P/work/gm12878_td/ribo/${name}_final_$mode.ribo.tsv
echo "DONE $name $mode"
