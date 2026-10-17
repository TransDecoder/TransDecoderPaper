#!/bin/bash
# Run Perl TransDecoder and PyTransDecoder on the shared sample datasets and
# compare their final ORF predictions.
#
# usage: run_regression.sh <perl_TD_dir> <pytd_dir> <out_dir> [extra args for both Predict steps]
#   e.g. run_regression.sh ../../repos/TransDecoder ../../repos/PyTransDecoder runs/baseline_v6.0.0
#
# Requires: perl, and a python with biopython+tqdm (PYTHON env var, default: python3).
set -euo pipefail

TD=$(readlink -f "$1"); PYTD=$(readlink -f "$2"); OUT=$(mkdir -p "$3" && readlink -f "$3"); shift 3
EXTRA=("$@")
PYTHON=${PYTHON:-python3}
HERE=$(dirname "$(readlink -f "$0")")
SD=$TD/sample_data
export PERL_HASH_SEED=0

mkdir -p "$OUT/inputs"
cd "$OUT/inputs"
[ -s trinity.fa ]    || gunzip -c $SD/simple_transcriptome_target/Trinity.fasta.gz > trinity.fa
[ -s pasa.fa ]       || gunzip -c $SD/pasa_example/pasa_assemblies.fasta.gz > pasa.fa
[ -s stringtie.fa ]  || cp $SD/stringtie_example/stringtie_merged.transcripts.fasta stringtie.fa
[ -s supertx.fa ]    || cp $SD/supertranscripts_example/supertranscripts.fasta supertx.fa
if [ ! -s cufflinks.fa ]; then
    gunzip -c $SD/cufflinks_example/test.genome.fasta.gz > cuff.genome.fa
    gunzip -c $SD/cufflinks_example/transcripts.gtf.gz > cuff.gtf
    $TD/util/gtf_genome_to_cdna_fasta.pl cuff.gtf cuff.genome.fa > cufflinks.fa
fi

# dataset -> extra LongOrfs flags (stringtie example is run strand-specific upstream)
declare -A LOFLAGS=( [trinity]="" [pasa]="" [stringtie]="-S" [supertx]="" [cufflinks]="" )

summary=$OUT/summary.tsv
echo -e "dataset\tmetric\tvalue" > "$summary"

for ds in trinity pasa stringtie supertx cufflinks; do
    fa=$OUT/inputs/$ds.fa
    for impl in perl py; do
        d=$OUT/$ds/$impl; mkdir -p "$d"; cd "$d"
        ln -sf "$fa" $ds.fa
        if [ -s $ds.fa.transdecoder.gff3 ]; then continue; fi
        start=$(date +%s.%N)
        if [ $impl == perl ]; then
            $TD/util/TransDecoder.LongOrfs -t $ds.fa ${LOFLAGS[$ds]} > longorfs.log 2>&1
            $TD/util/TransDecoder.Predict  -t $ds.fa "${EXTRA[@]}" > predict.log 2>&1
        else
            PYTHONPATH=$PYTD $PYTHON -m pytransdecoder.cli longorfs -t $ds.fa ${LOFLAGS[$ds]} > longorfs.log 2>&1
            PYTHONPATH=$PYTD $PYTHON -m pytransdecoder.cli predict  -t $ds.fa "${EXTRA[@]}" > predict.log 2>&1
        fi
        end=$(date +%s.%N)
        echo "$(echo "$end - $start" | bc)" > runtime_sec.txt
    done
    cd "$OUT/$ds"
    $PYTHON $HERE/compare_orfs.py perl/$ds.fa.transdecoder.gff3 py/$ds.fa.transdecoder.gff3 \
        --label-a perl --label-b py --tsv discordant.tsv > compare.txt
    sed "s/^/$ds\t/" compare.txt >> "$summary"
    echo -e "$ds\tperl_runtime_sec\t$(cat perl/runtime_sec.txt)" >> "$summary"
    echo -e "$ds\tpy_runtime_sec\t$(cat py/runtime_sec.txt)" >> "$summary"
done
column -t "$summary"
