#!/bin/bash
# GM12878 homology-restriction series: one PyTD predict config per call, reusing the final_abinitio LongOrfs output.
# usage: homology_series.sh <name> <structures gtf|bam> <blast clade taxid|none> <pfam|nopfam> [threads]
#   blast DB: Swiss-Prot minus every entry in the clade (9606 = self, 9443 primates, 40674 mammals,
#   7742 vertebrates, 33208 metazoa). Pfam hits: $W/<name>/pfam.domtblout (must exist for 'pfam').
set -euo pipefail
name=$1; structs=$(readlink -f $2); clade=$3; pf=$4; threads=${5:-4}
P=/broad/hptmp/bhaas/TRANSDECODER_PAPER; D=$P/data/gm12878; W=$P/work/gm12878_td; H=$P/data/homology
PY=$P/env/td/bin/python; export PATH=$P/env/td/bin:$PATH; export PYTHONPATH=$P/work/trees/PY_dev
tag=blast_${clade}_$pf; [ $clade = none ] && tag=$pf
src=$W/$name/final_abinitio; out=$W/$name/hseries/$tag; mkdir -p $out; cd $out
ln -sf $(readlink -f $src/transcripts.fa) transcripts.fa
[ -d transcripts.fa.transdecoder_dir ] || rsync -a --exclude '__checkpoints*' $src/transcripts.fa.transdecoder_dir/ transcripts.fa.transdecoder_dir/
extra=()
if [ $clade != none ]; then
  db=$H/sprot_minus_clade_$clade.dmnd; [ $clade = 9606 ] && db=$H/sprot_minus_9606.dmnd
  hits=$W/$name/hseries/blastp.minus_$clade.outfmt6
  ( flock 9; [ -s $hits ] || diamond blastp --query $src/transcripts.fa.transdecoder_dir/longest_orfs.pep --db $db \
      --max-target-seqs 1 --evalue 1e-5 --outfmt 6 --threads $threads --out $hits > $hits.log 2>&1 ) 9>$hits.lock
  extra+=(--retain-blastp-hits $hits)
fi
[ $pf = pfam ] && extra+=(--retain-pfam-hits $W/$name/pfam.domtblout)
td() { $PY -m pytransdecoder.cli "$@"; }
td predict -t transcripts.fa ${extra[@]+"${extra[@]}"} > predict.log 2>&1
td genome-map --orfs transcripts.fa.transdecoder.gff3 --structures $structs -t transcripts.fa \
   -o transcripts.fa.transdecoder.genome.gff3 > genome_map.log 2>&1
mkdir -p $W/eval/hseries
$PY $P/TransDecoderPaper/genome_annotation/scripts/eval_genome_cds.py --pred transcripts.fa.transdecoder.genome.gff3 \
   --structures $structs --mane $P/data/annot/human/MANE.GRCh38.v1.5.ensembl_genomic.gtf.gz \
   --gencode $D/ref/gencode.v50.primary_assembly.annotation.gtf --label ${name}_$tag \
   -o $W/eval/hseries/${name}_$tag > eval.log 2>&1
echo "DONE $name $tag"
