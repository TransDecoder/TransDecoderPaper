#!/bin/bash
# Run one predictor on one benchmark input and write standardized predictions.
#
# usage: run_tool.sh <tool> <input.fa> <outdir> <org>
#
# tools: pytd pytd_norefine pytd_homology perl_td perl_td_v600 perl_td_v571 prodigal codan gmst
#        cpat cpc2 borf rnasamba angel
# Inputs named *.unstranded* are predicted on both strands; other inputs are treated
# as strand-specific (TransDecoder -S, Prodigal + strand only, CodAn --strand plus).
# Fragment inputs use CodAn's *_partial models, complete inputs its *_full models.
#
# Writes <outdir>/preds.tsv (transcript_id, start, end, strand; 1-based, CDS incl. stop)
# and <outdir>/runtime_sec.txt. For pytd also writes <outdir>/longest_orf.preds.tsv and
# <outdir>/all_orfs.preds.tsv (baselines from LongOrfs).
#
# Environment (defaults for the Broad workspace): BENCH_ROOT, PYTHON, CODAN_ENV, CODAN_MODELS,
# PRODIGAL, GMST_DIR (GeneMarkS-T; key in ~/.gm_key), HOMOLOGY_DIR, THREADS, and TD_FIX/TD_V600/TD_V571/PYTD worktree paths.
set -euo pipefail

tool=$1; input=$(readlink -f "$2"); outdir=$3; org=$4
here=$(cd "$(dirname "$0")" && pwd)
ROOT=${BENCH_ROOT:-/broad/hptmp/bhaas/TRANSDECODER_PAPER}
PYTHON=${PYTHON:-$ROOT/env/td/bin/python}
CODAN_ENV=${CODAN_ENV:-$ROOT/env/codan}
CODAN_MODELS=${CODAN_MODELS:-$ROOT/tools/CodAn/models}
PRODIGAL=${PRODIGAL:-$ROOT/env/td/bin/prodigal}
GMST_DIR=${GMST_DIR:-$HOME/software/genemarks_t/dist}  # licensed, not distributed here
CPAT_ENV=${CPAT_ENV:-$ROOT/env/cpat}; CPAT_MODELS=${CPAT_MODELS:-$ROOT/tools/models/cpat}
CPC2_ENV=${CPC2_ENV:-$ROOT/env/cpc2}; BORF_ENV=${BORF_ENV:-$ROOT/env/borf}; ANGEL_ENV=${ANGEL_ENV:-$ROOT/env/angel}
RNASAMBA_ENV=${RNASAMBA_ENV:-$ROOT/env/rnasamba}; RNASAMBA_MODELS=${RNASAMBA_MODELS:-$ROOT/tools/models/rnasamba}
HOMOLOGY_DIR=${HOMOLOGY_DIR:-$ROOT/data/homology}
THREADS=${THREADS:-4}
TD_FIX=${TD_FIX:-$ROOT/work/trees/TD_fix}
TD_V600=${TD_V600:-$ROOT/work/trees/TD_base}
TD_V571=${TD_V571:-$ROOT/work/trees/TD_v571}
PYTD=${PYTD:-$ROOT/work/trees/PY_fix}
PYTD_DEV=${PYTD_DEV:-$ROOT/work/trees/PY_dev}
export PATH=$ROOT/env/td/bin:$PATH PERL_HASH_SEED=0

orgrow=$(awk -F'\t' -v o="$org" '$1==o' "$here/../organisms.tsv")
[ -n "$orgrow" ] || { echo "unknown org $org" >&2; exit 1; }
gcode=$(cut -f4 <<< "$orgrow"); table=$(cut -f5 <<< "$orgrow"); clade=$(cut -f6 <<< "$orgrow")

stranded=1; [[ $(basename "$input") == *unstranded* ]] && stranded=0
partial=0; [[ $(basename "$input") == fragments* ]] && partial=1

mkdir -p "$outdir"; outdir=$(readlink -f "$outdir"); cd "$outdir"
# Several runners may share a runs directory: skip a run another runner holds or has finished.
exec 9> .run.lock
flock -n 9 || { echo "$outdir is being run by another runner; skipping"; exit 0; }
[ -s preds.tsv ] && { echo "$outdir already has preds.tsv; skipping"; exit 0; }
ln -sf "$input" input.fa
t0=$(date +%s.%N)

case $tool in
  pytd|pytd_norefine|pytd_dev|pytd_dev_pwm|pytd_dev_single|pytd_dev_kozak|pytd_dev_kozak_add|pytd_dev_m[0-9]*|pytd_dev_v2)
    # pytd_dev*: development tree (PYTD_DEV); pytd_dev_pwm chooses trimmed starts by PWM;
    # pytd_dev_single keeps one ORF per transcript (--single-best-only); pytd_dev_kozak: dev defaults
    # (Kozak scaled by segment length); pytd_dev_kozak_add: dev defaults with the additive Kozak term;
    # pytd_dev_m<N>: dev defaults with minimum ORF length N aa (LongOrfs -m); pytd_dev_m<N>_sc<C>: also
    # Predict --short-orf-min-contrast C; pytd_dev_v2: dev defaults after -m 50 + short-ORF test (2026-10-20)
    ss=(); [ $stranded = 1 ] && ss=(-S)
    nr=(); [ $tool = pytd_norefine ] && nr=(--no-refine-starts)
    [ $tool = pytd_dev_pwm ] && nr=(--start-atg-choice pwm)
    [ $tool = pytd_dev_single ] && nr=(--single-best-only)
    [[ $tool == pytd_dev* ]] && PYTD=$PYTD_DEV
    mm=(); [[ $tool =~ ^pytd_dev_m([0-9]+)(_sc([0-9.]+))?$ ]] && mm=(-m "${BASH_REMATCH[1]}")
    [[ $tool =~ _sc([0-9.]+)$ ]] && nr=(--short-orf-min-contrast "${BASH_REMATCH[1]}")
    PYTHONPATH=$PYTD $PYTHON -m pytransdecoder.cli longorfs -t input.fa -G "$gcode" ${ss[@]+"${ss[@]}"} ${mm[@]+"${mm[@]}"} > longorfs.log 2>&1
    PYTHONPATH=$PYTD $PYTHON -m pytransdecoder.cli predict -t input.fa -G "$gcode" ${nr[@]+"${nr[@]}"} > predict.log 2>&1
    $PYTHON "$here/convert_predictions.py" transdecoder_gff3 input.fa.transdecoder.gff3 > preds.tsv
    if [ $tool = pytd ]; then
      lo=input.fa.transdecoder_dir/longest_orfs.gff3
      $PYTHON "$here/convert_predictions.py" transdecoder_gff3 $lo > all_orfs.preds.tsv
      $PYTHON "$here/convert_predictions.py" transdecoder_gff3 --longest-per-transcript $lo > longest_orf.preds.tsv
    fi
    ;;
  pytd_homology)
    # homology-assisted: DIAMOND vs Swiss-Prot without the organism's own entries (-k 1, 1e-5, as
    # in the TransDecoder wrapper) and hmmsearch vs Pfam-A with gathering thresholds
    taxid=$(cut -f7 <<< "$orgrow")
    ss=(); [ $stranded = 1 ] && ss=(-S)
    PYTHONPATH=$PYTD $PYTHON -m pytransdecoder.cli longorfs -t input.fa -G "$gcode" ${ss[@]+"${ss[@]}"} > longorfs.log 2>&1
    pep=input.fa.transdecoder_dir/longest_orfs.pep
    diamond blastp --query $pep --db "$HOMOLOGY_DIR/sprot_minus_${taxid}.dmnd" --max-target-seqs 1 \
      --evalue 1e-5 --outfmt 6 --threads "$THREADS" --out blastp.outfmt6 > diamond.log 2>&1
    hmmsearch --cpu "$THREADS" --cut_ga --noali -o /dev/null --domtblout pfam.domtblout \
      "$HOMOLOGY_DIR/Pfam-A.hmm" $pep > hmmsearch.log 2>&1
    PYTHONPATH=$PYTD $PYTHON -m pytransdecoder.cli predict -t input.fa -G "$gcode" \
      --retain-blastp-hits blastp.outfmt6 --retain-pfam-hits pfam.domtblout > predict.log 2>&1
    $PYTHON "$here/convert_predictions.py" transdecoder_gff3 input.fa.transdecoder.gff3 > preds.tsv
    ;;
  perl_td|perl_td_dev|perl_td_dev_v2|perl_td_v600|perl_td_v571)
    # perl_td_dev: the TD_FIX tree with its current (development) defaults
    case $tool in perl_td|perl_td_dev|perl_td_dev_v2) TD=$TD_FIX;; perl_td_v600) TD=$TD_V600;; perl_td_v571) TD=$TD_V571;; esac
    LO=$TD/TransDecoder.LongOrfs; PR=$TD/TransDecoder.Predict
    [ -x "$LO" ] || { LO=$TD/util/TransDecoder.LongOrfs; PR=$TD/util/TransDecoder.Predict; }
    ss=(); [ $stranded = 1 ] && ss=(-S)
    "$LO" -t input.fa -G "$gcode" ${ss[@]+"${ss[@]}"} > longorfs.log 2>&1
    "$PR" -t input.fa -G "$gcode" > predict.log 2>&1
    $PYTHON "$here/convert_predictions.py" transdecoder_gff3 input.fa.transdecoder.gff3 > preds.tsv
    ;;
  prodigal)
    # single-genome mode trains on the input; -n forces a full motif scan (no SD trainer), as in 2018
    "$PRODIGAL" -i input.fa -f gff -o prodigal.gff -g "$table" -n -q > prodigal.log 2>&1
    sf=(); [ $stranded = 1 ] && sf=(--plus-only)
    $PYTHON "$here/convert_predictions.py" prodigal_gff ${sf[@]+"${sf[@]}"} prodigal.gff > preds.tsv
    ;;
  codan)
    model=$CODAN_MODELS/${clade}_full; [ $partial = 1 ] && model=$CODAN_MODELS/${clade}_partial
    strand=plus; [ $stranded = 0 ] && strand=both
    rm -rf codan_out
    # With -s plus, CodAn 1.1 leaves its predictions in ORFs.gtf (no annotation.gtf) and then
    # exits non-zero removing a FASTA index Biopython never created; accept the run only if
    # predict produced output.
    gtf=codan_out/annotation.gtf; [ $strand = plus ] && gtf=codan_out/ORFs.gtf
    # CodAn's predict indexes the FASTA with Bio::DB::Fasta, which rejects lines >65,535 chars
    # (and then silently yields no plus-strand predictions), so give it a line-wrapped copy.
    awk '/^>/ {print; next} {for (i = 1; i <= length($0); i += 80) print substr($0, i, 80)}' input.fa > codan_input.fa
    PATH=$CODAN_ENV/bin:$PATH codan.py -t codan_input.fa -m "$model" -s $strand -c 1 -o codan_out > codan.log 2>&1 || true
    rm -f codan_input.fa codan_input.fa.index*
    [ -s $gtf ] || { echo "CodAn produced no predictions (see codan.log)" >&2; exit 1; }
    $PYTHON "$here/convert_predictions.py" codan_gtf $gtf > preds.tsv
    ;;
  gmst)
    # GeneMarkS-T self-trains on the input. It supports genetic codes 1, 4 and 11 only, so
    # organisms with another nuclear code (Tetrahymena, table 6) are run with the standard code.
    gc=1; [ "$table" = 4 ] && gc=4
    strand=direct; [ $stranded = 0 ] && strand=both
    perl "$GMST_DIR/gmst.pl" --format GFF --gcode $gc --strand $strand --output gmst.gff input.fa > gmst.log 2>&1
    $PYTHON "$here/convert_predictions.py" gmst_gff gmst.gff > preds.tsv
    ;;
  cpat)
    # CPAT 3 with its nearest prebuilt model (human, mouse, fly; human for the others) and that
    # model's published cutoff; the ciliate code is given through --stop.
    case $org in mouse) m=Mouse; hx=Mouse; cut=0.44;; dmel) m=Fly; hx=fly; cut=0.39;; *) m=Human; hx=Human; cut=0.364;; esac
    stops=TAG,TAA,TGA; [ "$table" = 6 ] && stops=TGA
    as=(); [ $stranded = 0 ] && as=(--antisense)
    PATH=$CPAT_ENV/bin:$PATH cpat -g input.fa -o cpat -d "$CPAT_MODELS/${m}_logitModel.RData" \
      -x "$CPAT_MODELS/${hx}_Hexamer.tsv" --stop $stops ${as[@]+"${as[@]}"} --log-file cpat.run.log > cpat.log 2>&1
    $PYTHON "$here/convert_predictions.py" cpat_tsv --cutoff $cut --fasta input.fa cpat.ORF_prob.best.tsv > preds.tsv
    ;;
  cpc2)
    # CPC2 1.0.1: longest ATG ORF per transcript, kept when labelled coding (standard code only)
    rc=(); [ $stranded = 0 ] && rc=(-r)
    PATH=$CPC2_ENV/bin:$PATH CPC2.py -i input.fa -o cpc2 --ORF ${rc[@]+"${rc[@]}"} > cpc2.log 2>&1
    $PYTHON "$here/convert_predictions.py" cpc2_txt --fasta input.fa cpc2.txt > preds.tsv
    ;;
  borf)
    # Borf 1.2 defaults (>=100 aa; standard code only)
    ss=(); [ $stranded = 0 ] && ss=(-s)
    $BORF_ENV/bin/borf -f ${ss[@]+"${ss[@]}"} -o borf input.fa > borf.log 2>&1
    $PYTHON "$here/convert_predictions.py" borf_txt borf.txt > preds.tsv
    ;;
  rnasamba)
    # RNAsamba 0.2.5 with the authors' human-trained weights (full-length or partial); it scores the
    # forward strand only, so it is run on stranded inputs only.
    [ $stranded = 1 ] || { echo "rnasamba is forward-strand only; not run on unstranded inputs" >&2; exit 1; }
    w=$RNASAMBA_MODELS/full_length_weights.hdf5; [ $partial = 1 ] && w=$RNASAMBA_MODELS/partial_length_weights.hdf5
    PATH=$RNASAMBA_ENV/bin:$PATH rnasamba classify -p rnasamba.pep rnasamba.tsv input.fa "$w" > rnasamba.log 2>&1
    $PYTHON "$here/convert_predictions.py" rnasamba --fasta input.fa rnasamba.pep > preds.tsv
    ;;
  angel)
    # ANGEL 3.0 self-training workflow from its README: longest ORFs -> non-redundant training set ->
    # classifier -> prediction (best ORF per transcript; standard code only)
    rv=(); [ $stranded = 0 ] && rv=(--use_rev_strand)
    export PATH=$ANGEL_ENV/bin:$PATH PYTHONWARNINGS=ignore
    { dumb_predict.py input.fa dumb --min_aa_length 100 --cpus "$THREADS" ${rv[@]+"${rv[@]}"}
      angel_make_training_set.py dumb.final train --random --cpus "$THREADS"
      angel_train.py train.cds train.utr angel.pickle --cpus "$THREADS"
      angel_predict.py input.fa angel.pickle angel --output_mode best --min_angel_aa_length 100 \
        --min_dumb_aa_length 100 --cpus "$THREADS" ${rv[@]+"${rv[@]}"}
    } > angel.log 2>&1
    $PYTHON "$here/convert_predictions.py" angel_cds --fasta input.fa angel.ANGEL.cds > preds.tsv
    ;;
  *) echo "unknown tool $tool" >&2; exit 1;;
esac

echo "$(date +%s.%N) - $t0" | bc > runtime_sec.txt
