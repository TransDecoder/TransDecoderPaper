#!/bin/bash
# Re-run Predict (current PyTD) with the existing DIAMOND/Pfam hits, then genome-map, eval, Ribo-seq match
set -euo pipefail
P=/broad/hptmp/bhaas/TRANSDECODER_PAPER; export PYTHONPATH=$P/work/trees/PY_dev; PY=$P/env/td/bin/python
D=$P/data/gm12878
declare -A ST=([lraa]=$D/longread_lraa/lraa_denovo/LRAA.gtf [stringtie]=$D/star_stringtie/stringtie.gtf)
one() {
  a=$1; src=$P/work/gm12878_td/$a/homology; d=$P/work/gm12878_td/$a/homology_v2; mkdir -p $d/transcripts.fa.transdecoder_dir; cd $d
  ln -sf $(readlink $src/transcripts.fa) transcripts.fa
  cp $src/transcripts.fa.transdecoder_dir/longest_orfs.* $src/transcripts.fa.transdecoder_dir/base_freqs.dat transcripts.fa.transdecoder_dir/
  $PY -m pytransdecoder.cli predict -t transcripts.fa --retain-blastp-hits $src/blastp.outfmt6 --retain-pfam-hits $src/pfam.domtblout > predict.log 2>&1
  $PY -m pytransdecoder.cli genome-map --orfs transcripts.fa.transdecoder.gff3 --structures ${ST[$a]} -t transcripts.fa \
     -o transcripts.fa.transdecoder.genome.gff3 > genome_map.log 2>&1
  $PY $P/TransDecoderPaper/genome_annotation/scripts/eval_genome_cds.py --pred transcripts.fa.transdecoder.genome.gff3 \
     --structures ${ST[$a]} --mane $P/data/annot/human/MANE.GRCh38.v1.5.ensembl_genomic.gtf.gz \
     --gencode $D/ref/gencode.v50.primary_assembly.annotation.gtf --label ${a}_homology -o $P/work/gm12878_td/eval/${a}_homology > eval.log 2>&1
  echo "done $a"
}
one lraa & one stringtie & wait
echo ALL_DONE
