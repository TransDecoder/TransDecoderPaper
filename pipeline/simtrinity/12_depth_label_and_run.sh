#!/bin/bash
# E31: label one depth level's assembly, run the nine methods and score them, as in E30.
# usage: 12_depth_label_and_run.sh <level>
# Labels: label_assembly.py (E30 final rules). Input: <level>/inputs/fragments_simtrinity.fa (the 'fragments'
# prefix selects partial-length models; stranded; human). Runs: run_tool.sh (unmodified), 4 at a time,
# in $P/work/simtrinity/depth/<level>/runs/<tool>; longest ORF from pytd_v7. Scores: score_simtrinity.py.
set -euo pipefail
level=$1
P=/broad/hptmp/bhaas/TRANSDECODER_PAPER
here=$(cd "$(dirname "$0")" && pwd)
D=$P/data/simtrinity/depth/$level
W=$P/work/simtrinity/depth/$level
R=$P/data/simtrinity/ref
PY=$P/env/td/bin/python
export PATH=$P/env/rnaseq_asm/bin:$PATH
mkdir -p $D/labels $D/inputs $W/runs/logs
$PY $here/label_assembly.py --bam $D/Trinity.minimap2.bam --fasta $D/Trinity.fasta \
  --gencode $R/gencode_v50.tx.tsv.gz --kept $R/kept.tsv --salmon $P/data/simtrinity/salmon/quant/quant.sf \
  --out $D/labels/sim_trinity > $D/labels/label.log 2>&1
$PY $P/TransDecoderPaper/pipeline/derive/03_utr_gm12878_truth.py $D/Trinity.minimap2.bam $D/Trinity.fasta \
  $P/data/annot/human/MANE.GRCh38.v1.5.ensembl_genomic.gtf.gz $D/labels/sim_trinity_e29 > $D/labels/e29.log 2>&1
cp $D/Trinity.fasta $D/inputs/fragments_simtrinity.fa
cp $D/labels/sim_trinity.truth.tsv $D/inputs/fragments_simtrinity.truth.tsv
export P W D
printf "%s\n" pytd_v7 perl_td_v600 gmst prodigal cpat cpc2 rnasamba borf | \
  xargs -P ${JOBS:-4} -I{} bash -c 'THREADS=1 bash $P/TransDecoderPaper/benchmark_v2/scripts/run_tool.sh {} $D/inputs/fragments_simtrinity.fa $W/runs/{} human > $W/runs/logs/{}.log 2>&1; echo "[$(date)] {} exit $?"'
$PY $here/score_simtrinity.py $D/inputs/fragments_simtrinity.fa $D/labels/sim_trinity $W/runs $W/scores
echo "[$(date)] $level: LEVEL_DONE"
