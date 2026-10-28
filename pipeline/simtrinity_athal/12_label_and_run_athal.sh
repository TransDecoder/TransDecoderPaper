#!/bin/bash
# E32: label one Arabidopsis level's assembly (label_assembly.py, E30 final rules, TAIR12 tables), run the nine
# methods with organism 'athal' (run_tool.sh, unmodified; 4 at a time) and score (score_simtrinity.py).
# No E29 statistic (it needs MANE Select). usage: 12_label_and_run_athal.sh <level>
set -euo pipefail
level=$1
P=/broad/hptmp/bhaas/TRANSDECODER_PAPER
here=$P/TransDecoderPaper/pipeline/simtrinity
D=$P/data/simtrinity_athal/depth/$level
W=$P/work/simtrinity_athal/depth/$level
R=$P/data/simtrinity_athal/ref
PY=$P/env/td/bin/python
export PATH=$P/env/rnaseq_asm/bin:$PATH
mkdir -p $D/labels $D/inputs $W/runs/logs
$PY $here/label_assembly.py --bam $D/Trinity.minimap2.bam --fasta $D/Trinity.fasta \
  --gencode $R/tair12.tx.tsv.gz --kept $R/kept.tsv --salmon $P/data/simtrinity_athal/salmon/quant/quant.sf \
  --out $D/labels/sim_trinity > $D/labels/label.log 2>&1
cp $D/Trinity.fasta $D/inputs/fragments_simtrinity.fa
cp $D/labels/sim_trinity.truth.tsv $D/inputs/fragments_simtrinity.truth.tsv
export P W D
printf "%s\n" pytd_v7 perl_td_v600 gmst prodigal cpat cpc2 rnasamba borf | \
  xargs -P ${JOBS:-4} -I{} bash -c 'THREADS=1 bash $P/TransDecoderPaper/benchmark_v2/scripts/run_tool.sh {} $D/inputs/fragments_simtrinity.fa $W/runs/{} athal > $W/runs/logs/{}.log 2>&1; echo "[$(date)] {} exit $?"'
$PY $here/score_simtrinity.py $D/inputs/fragments_simtrinity.fa $D/labels/sim_trinity $W/runs $W/scores
echo "[$(date)] $level: LEVEL_DONE"
