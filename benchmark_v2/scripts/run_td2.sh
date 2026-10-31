#!/bin/bash
# Run TD2 (Mao et al. 2025; pip package td2 1.1.0, $P/env/td2) on one benchmark input, ab initio
# (no homology), and write standardized predictions, as run_tool.sh does for the other comparators
# (E39). Kept separate from run_tool.sh so that script is not edited while other jobs may use it.
#
# usage: run_td2.sh <input.fa> <outdir> <org>
#
# Settings: TD2 defaults (TD2.LongOrfs minimum 90 aa; TD2.Predict PSAURON cutoff 0.5, dynamic
# long-ORF retention at FDR 0.10), except:
#   - strand: inputs named *unstranded* are predicted on both strands; all others with -S (as
#     TransDecoder -S, CodAn -s plus, GeneMarkS-T --strand direct);
#   - genetic code: -G <NCBI table> from organisms.tsv (6 for T. thermophila), for both steps;
#   - threads: LongOrfs -@ $TD2_THREADS (default 2); torch/OpenMP threads of PSAURON in Predict
#     limited to the same number (OMP_NUM_THREADS, MKL_NUM_THREADS).
# TD2.Predict reports at most one ORF per transcript by default (no --all-good); this is checked
# after conversion (the run fails if a transcript has more than one prediction).
#
# Writes <outdir>/preds.tsv (transcript_id, start, end, strand; 1-based, CDS incl. stop),
# <outdir>/runtime_sec.txt and <outdir>/td2_version.txt.
set -euo pipefail

input=$(readlink -f "$1"); outdir=$2; org=$3
here=$(cd "$(dirname "$0")" && pwd)
ROOT=${BENCH_ROOT:-/broad/hptmp/bhaas/TRANSDECODER_PAPER}
PYTHON=${PYTHON:-$ROOT/env/td/bin/python}
TD2_ENV=${TD2_ENV:-$ROOT/env/td2}
NT=${TD2_THREADS:-2}
export OMP_NUM_THREADS=$NT MKL_NUM_THREADS=$NT PATH=$TD2_ENV/bin:$PATH   # TD2.Predict calls psauron from PATH

orgrow=$(awk -F'\t' -v o="$org" '$1==o' "$here/../organisms.tsv")
[ -n "$orgrow" ] || { echo "unknown org $org" >&2; exit 1; }
table=$(cut -f5 <<< "$orgrow")

stranded=1; [[ $(basename "$input") == *unstranded* ]] && stranded=0

mkdir -p "$outdir"; outdir=$(readlink -f "$outdir"); cd "$outdir"
exec 9> .run.lock
flock -n 9 || { echo "$outdir is being run by another runner; skipping"; exit 0; }
[ -s preds.tsv ] && { echo "$outdir already has preds.tsv; skipping"; exit 0; }
ln -sf "$input" input.fa
$TD2_ENV/bin/python -c 'import importlib.metadata as m; print("td2", m.version("td2"), "psauron", m.version("psauron"))' > td2_version.txt
t0=$(date +%s.%N)

ss=(); [ $stranded = 1 ] && ss=(-S)
$TD2_ENV/bin/TD2.LongOrfs -t input.fa ${ss[@]+"${ss[@]}"} -G "$table" -@ "$NT" > longorfs.log 2>&1
$TD2_ENV/bin/TD2.Predict -t input.fa -G "$table" > predict.log 2>&1
$PYTHON "$here/convert_predictions.py" transdecoder_gff3 input.fa.TD2.gff3 > preds.tsv.tmp

echo "$(date +%s.%N) - $t0" | bc > runtime_sec.txt
dups=$(tail -n +2 preds.tsv.tmp | cut -f1 | sort | uniq -d | wc -l)
[ "$dups" = 0 ] || { echo "TD2 reported more than one ORF for $dups transcripts" >&2; exit 1; }
mv preds.tsv.tmp preds.tsv
