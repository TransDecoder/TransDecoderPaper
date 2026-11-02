#!/bin/bash
# One fresh PyTD run on one benchmark input with the pinned release code (E40). Never resumes: an
# unfinished output directory is moved to $P/_trash first.
# usage: run_bench_v1.sh <org> <input> <label>
# labels:
#   rel               release defaults (q 0.05)
#   precise|sensitive --precise (q 0.01), --sensitive (q 0.2)
#   q<Q>              --max-qvalue Q
#   m<N>_noq          LongOrfs -m N; q-value test off (--max-qvalue 1000000; FDR estimates can exceed 1)
#   nolong_q<Q>       long-ORF rule off (--retain-long-orfs-mode strict --retain-long-orfs-length 1000000), q Q
#   retonly           homology hits for retention only (--no-homology-start-assist); hits of pytd_v7_homology,
#                     whose longest_orfs.pep must equal this run's (cmp)
set -euo pipefail
org=$1; inp=$2; label=$3
P=/broad/hptmp/bhaas/TRANSDECODER_PAPER; T=$P/TransDecoderPaper
PYTD=${PYTD_REL:-$P/work/trees/PY_rel_d3892b5}; PY=$P/env/td/bin/python; S=$T/benchmark_v2/scripts
export PATH=$P/env/td/bin:$PATH OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
out=$P/work/release_pin/bench/$org/$inp/$label
[ -s $out/DONE ] && { echo "done already: $out"; exit 0; }
[ -e $out ] && { mkdir -p $P/_trash/release_pin_partial; mv $out $P/_trash/release_pin_partial/$org.$inp.$label.$(date +%s); }
mkdir -p $out; cd $out
ln -s $P/data/bench_v3/inputs/$org/$inp.fa input.fa
gcode=$(awk -F'\t' -v o=$org '$1==o{print $4}' $T/benchmark_v2/organisms.tsv)
ss=(-S); [[ $inp == *unstranded* ]] && ss=()
lo=(); pr=()
case $label in
  rel) ;;
  precise|sensitive) pr=(--$label);;
  q*) pr=(--max-qvalue ${label#q});;
  m*_noq) n=${label#m}; lo=(-m ${n%_noq}); pr=(--max-qvalue 1000000);;
  nolong_q*) pr=(--retain-long-orfs-mode strict --retain-long-orfs-length 1000000 --max-qvalue ${label#nolong_q});;
  retonly) ;;
  *) echo "bad label $label" >&2; exit 1;;
esac
(cd $PYTD && git rev-parse --short HEAD) > pytd_version.txt
t0=$(date +%s.%N)
PYTHONPATH=$PYTD $PY -m pytransdecoder.cli longorfs -t input.fa -G $gcode "${ss[@]}" ${lo[@]+"${lo[@]}"} > longorfs.log 2>&1
if [ $label = retonly ]; then
  h=$P/work/bench_v3/runs/$org/$inp/pytd_v7_homology
  cmp input.fa.transdecoder_dir/longest_orfs.pep $h/input.fa.transdecoder_dir/longest_orfs.pep
  pr=(--retain-blastp-hits $h/blastp.outfmt6 --retain-pfam-hits $h/pfam.domtblout --no-homology-start-assist)
fi
PYTHONPATH=$PYTD $PY -m pytransdecoder.cli predict -t input.fa -G $gcode ${pr[@]+"${pr[@]}"} > predict.log 2>&1
echo "$(date +%s.%N) - $t0" | bc > runtime_sec.txt
$PY $S/convert_predictions.py transdecoder_gff3 input.fa.transdecoder.gff3 > preds.tsv
if [ $label = rel ]; then
  $PY $S/convert_predictions.py transdecoder_gff3 --longest-per-transcript input.fa.transdecoder_dir/longest_orfs.gff3 > longest_orf.preds.tsv
fi
date > DONE
