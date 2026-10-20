#!/bin/bash
# GM12878: LongOrfs -m 50 with the short-ORF frame-contrast guard (Predict --short-orf-min-contrast C)
set -euo pipefail
P=/broad/hptmp/bhaas/TRANSDECODER_PAPER; export PYTHONPATH=$P/work/trees/PY_dev; PY=$P/env/td/bin/python
D=$P/data/gm12878
declare -A FA=([lraa]=$D/longread_lraa/lraa_denovo/LRAA.transcripts.fa [stringtie]=$D/star_stringtie/stringtie.transcripts.fa)
declare -A ST=([lraa]=$D/longread_lraa/lraa_denovo/LRAA.gtf [stringtie]=$D/star_stringtie/stringtie.gtf)
one() {
  a=$1; c=$2; d=$P/work/gm12878_td/$a/m50_sc$c; mkdir -p $d; cd $d; ln -sf ${FA[$a]} transcripts.fa
  # reuse the m50 LongOrfs output
  [ -d transcripts.fa.transdecoder_dir ] || { mkdir transcripts.fa.transdecoder_dir; cp $P/work/gm12878_td/$a/m50/transcripts.fa.transdecoder_dir/longest_orfs.* $P/work/gm12878_td/$a/m50/transcripts.fa.transdecoder_dir/base_freqs.dat transcripts.fa.transdecoder_dir/; }
  $PY -m pytransdecoder.cli predict -t transcripts.fa --short-orf-min-contrast $c > predict.log 2>&1
  $PY -m pytransdecoder.cli genome-map --orfs transcripts.fa.transdecoder.gff3 --structures ${ST[$a]} -t transcripts.fa \
     -o transcripts.fa.transdecoder.genome.gff3 > genome_map.log 2>&1
  $PY $P/TransDecoderPaper/genome_annotation/scripts/eval_genome_cds.py --pred transcripts.fa.transdecoder.genome.gff3 \
     --structures ${ST[$a]} --mane $P/data/annot/human/MANE.GRCh38.v1.5.ensembl_genomic.gtf.gz \
     --gencode $D/ref/gencode.v50.primary_assembly.annotation.gtf --label ${a}_m50_sc$c -o $P/work/gm12878_td/eval/${a}_m50_sc$c > eval.log 2>&1
  PYTHONPATH=$P/TransDecoderPaper/genome_annotation/scripts:$PYTHONPATH $PY $P/TransDecoderPaper/analyses/riboseq/match_riboseq.py \
     --pred transcripts.fa.transdecoder.genome.gff3 --phase1 $P/data/riboseq_orfs/phase1/Ribo-seq_ORFs.all.bed \
     --phase2 $P/data/riboseq_orfs/phase2/Ribo-seq_ORFs.comprehensive.bed > $P/work/gm12878_td/ribo/${a}_m50_sc$c.ribo.tsv
  echo "done $a sc$c"
}
one lraa 200 & one stringtie 200 & wait
echo ALL_DONE
