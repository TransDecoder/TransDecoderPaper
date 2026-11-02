#!/usr/bin/env bash
# Protein matches, Ribo-seq and MS support of the pinned-commit GM12878 runs without the q-value test, and of the
# default (E40). Methods of E01, E02, E25. Output: work/release_pin/gm_scoring/, tables in data/gm12878/release_pin/.
set -euo pipefail
P=/broad/hptmp/bhaas/TRANSDECODER_PAPER; T=$P/TransDecoderPaper; G=$P/work/release_pin/gm; W=$P/work/release_pin/gm_scoring; mkdir -p $W/ribo $W/proteomics
PY=$P/env/td/bin/python; D=$P/data/gm12878/proteome_sinitcyn2023; O=$T/data/gm12878/release_pin
ASM="lraa stringtie trinity"
# protein matches
args=(); for a in $ASM; do for c in rel precise sensitive m50_noq m100_noq; do [ -s $G/$a/$c/DONE ] || continue
  $PY $T/benchmark_v2/scripts/convert_predictions.py transdecoder_gff3 $G/$a/$c/transcripts.fa.transdecoder.gff3 > $W/$a.$c.preds.tsv
  args+=("$a.$c=$G/$a/$c/transcripts.fa:$W/$a.$c.preds.tsv"); done; done
$PY $T/benchmark_v2/scripts/score_protein_matches.py $D/refdb/gencode.v50.pc_translations.fa.gz "${args[@]}" > $O/protein_matches.tsv
# Ribo-seq and MS support: no q-value test (m50_noq), and rel (must reproduce E25)
ra=(); ma=(); for a in $ASM; do for c in rel m50_noq; do
  ra+=(--pred ${a}_$c=$G/$a/$c/transcripts.fa.transdecoder.genome.gff3,$G/$a/$c/eval.pred_status.tsv)
  ma+=(--pred ${a}_$c=$G/$a/$c/transcripts.fa.transdecoder.pep,$G/$a/$c/eval.pred_status.tsv); done; done
cd $P
env/riboseq/bin/python $T/analyses/riboseq/riboseq_orf_support.py --bam data/gm12878/riboseq_cenik2015/footprints.genome.bam \
  --mane data/annot/human/MANE.GRCh38.v1.5.ensembl_genomic.gtf.gz "${ra[@]}" -o $W/ribo/gm > $W/ribo/log 2>&1
env/proteo/bin/python -I $T/analyses/proteomics/ms_peptide_support.py --peptides $D/search/GM12878/peptides.txt \
  --ref $D/refdb/gencode.v50.pc_translations.fa.gz --ref $D/refdb/uniprot_human_UP000005640_with_isoforms.fasta.gz \
  "${ma[@]}" -o $W/proteomics > $W/proteomics/log 2>&1
cp $W/ribo/gm.summary.tsv $O/riboseq_summary.tsv; cp $W/proteomics/ms_peptide_support.summary.tsv $O/ms_summary.tsv
$PY $T/pipeline/release_pin/14_rule_support_v1.py > $O/rule_support.tsv
