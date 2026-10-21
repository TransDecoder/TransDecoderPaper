#!/bin/bash
# Post-hoc minimum-length filter on the final-default predictions: eval at --min-aa L.
P=/broad/hptmp/bhaas/TRANSDECODER_PAPER; D=$P/data/gm12878; W=$P/work/gm12878_td
export PYTHONPATH=$P/work/trees/PY_dev
export ST_lraa=$D/longread_lraa/lraa_denovo/LRAA.gtf ST_stringtie=$D/star_stringtie/stringtie.gtf ST_trinity=$D/trinity/Trinity.minimap2.bam
mkdir -p $W/eval/minaa
for a in lraa stringtie trinity; do for L in 50 60 70 80 90 100 150 200 300; do echo "$a $L"; done; done | \
xargs -P 6 -L 1 bash -c '$0/env/td/bin/python $0/TransDecoderPaper/genome_annotation/scripts/eval_genome_cds.py \
  --pred $1/$2/final_abinitio/transcripts.fa.transdecoder.genome.gff3 --structures $(eval echo \${ST_$2}) \
  --mane $0/data/annot/human/MANE.GRCh38.v1.5.ensembl_genomic.gtf.gz --gencode $0/data/gm12878/ref/gencode.v50.primary_assembly.annotation.gtf \
  --label $2_min$3 --min-aa $3 -o $1/eval/minaa/$2_min$3 > /dev/null 2>&1 && echo "ok $2 $3" || echo "FAILED $2 $3"' $P $W
echo SWEEP_DONE
