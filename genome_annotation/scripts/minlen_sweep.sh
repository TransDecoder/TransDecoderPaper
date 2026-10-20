#!/bin/bash
# GM12878 ab initio annotation at minimum ORF lengths 90..50 aa (PyTD dev), then genome-map + eval
set -euo pipefail
P=/broad/hptmp/bhaas/TRANSDECODER_PAPER; export PYTHONPATH=$P/work/trees/PY_dev; PY=$P/env/td/bin/python
D=$P/data/gm12878
declare -A FA=([lraa]=$D/longread_lraa/lraa_denovo/LRAA.transcripts.fa [stringtie]=$D/star_stringtie/stringtie.transcripts.fa)
declare -A ST=([lraa]=$D/longread_lraa/lraa_denovo/LRAA.gtf [stringtie]=$D/star_stringtie/stringtie.gtf)
one() {
  a=$1; m=$2; d=$P/work/gm12878_td/$a/m$m; mkdir -p $d; cd $d; ln -sf ${FA[$a]} transcripts.fa
  $PY -m pytransdecoder.cli longorfs -t transcripts.fa -S -m $m > longorfs.log 2>&1
  $PY -m pytransdecoder.cli predict -t transcripts.fa > predict.log 2>&1
  $PY -m pytransdecoder.cli genome-map --orfs transcripts.fa.transdecoder.gff3 --structures ${ST[$a]} -t transcripts.fa \
     -o transcripts.fa.transdecoder.genome.gff3 > genome_map.log 2>&1
  $PY $P/TransDecoderPaper/genome_annotation/scripts/eval_genome_cds.py --pred transcripts.fa.transdecoder.genome.gff3 \
     --structures ${ST[$a]} --mane $P/data/annot/human/MANE.GRCh38.v1.5.ensembl_genomic.gtf.gz \
     --gencode $D/ref/gencode.v50.primary_assembly.annotation.gtf --label ${a}_m$m -o $P/work/gm12878_td/eval/${a}_m$m > eval.log 2>&1
  echo "done $a m$m"
}
for m in 90 80 70 60 50; do one lraa $m & one stringtie $m & wait; done
echo ALL_DONE
