#!/bin/bash
# E33 adoption: GM12878 Trinity GG with the PyTD homology start assist (PY_hstart), reusing the DIAMOND / hmmsearch
# outputs of the release homology run (trinity/final_v7_homology) as user inputs; otherwise as
# genome_annotation/scripts/gm12878_final_v7.sh (fresh LongOrfs -S, Predict, genome-map, MANE/GENCODE eval).
set -euo pipefail
P=/broad/hptmp/bhaas/TRANSDECODER_PAPER; D=$P/data/gm12878; PY=$P/env/td/bin/python
export PATH=$P/env/td/bin:$PATH; export PYTHONPATH=$P/work/trees/PY_hstart
src=$P/work/gm12878_td/trinity/final_v7_homology
structs=$D/trinity/Trinity.minimap2.bam
out=$P/work/gm12878_td/trinity/final_v7h_homology
[ -e $out ] && { echo "exists: $out" >&2; exit 1; }
mkdir -p $out; cd $out
(cd $PYTHONPATH && git rev-parse --short HEAD) > pytd_version.txt
ln -s $(readlink -f $src/transcripts.fa) transcripts.fa
td() { $PY -m pytransdecoder.cli "$@"; }
td longorfs -t transcripts.fa -S > longorfs.log 2>&1
cmp transcripts.fa.transdecoder_dir/longest_orfs.pep $src/transcripts.fa.transdecoder_dir/longest_orfs.pep
td predict -t transcripts.fa --retain-blastp-hits $src/blastp.outfmt6 --retain-pfam-hits $src/pfam.domtblout > predict.log 2>&1
td genome-map --orfs transcripts.fa.transdecoder.gff3 --structures $structs -t transcripts.fa \
   -o transcripts.fa.transdecoder.genome.gff3 > genome_map.log 2>&1
$PY $P/TransDecoderPaper/genome_annotation/scripts/eval_genome_cds.py --pred transcripts.fa.transdecoder.genome.gff3 \
   --structures $structs --mane $P/data/annot/human/MANE.GRCh38.v1.5.ensembl_genomic.gtf.gz \
   --gencode $D/ref/gencode.v50.primary_assembly.annotation.gtf --label trinity_final_v7h_homology \
   -o $P/work/gm12878_td/eval/trinity_final_v7h_homology > eval.log 2>&1
$PY $P/TransDecoderPaper/benchmark_v2/scripts/convert_predictions.py transdecoder_gff3 transcripts.fa.transdecoder.gff3 > preds.tsv
echo "DONE trinity v7h_homology"
