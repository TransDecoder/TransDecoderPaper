#!/bin/bash
# GM12878 expression for the expressed-gene sensitivity analysis (evidence/E05):
# StringTie -e (reference transcripts only, GENCODE v50) on the STAR alignments of ENCSR000AED rep1.
# RECONSTRUCTED 2026-10-21 from the outputs: this command was run as a one-off, not from a script.
# The stringtie line is copied verbatim from the header of expr/gencode_v50.ref_quant.gtf (run from
# data/gm12878); expr/stringtie_e.log holds only the EXPR_DONE marker.
set -euo pipefail
cd /broad/hptmp/bhaas/TRANSDECODER_PAPER/data/gm12878
mkdir -p expr
../../env/rnaseq_asm/bin/stringtie star_stringtie/GM12878_ENCSR000AED_rep1.STAR.bam --rf -e -G ref/gencode.v50.primary_assembly.annotation.gtf -p 8 -A expr/gencode_v50.gene_abund.tsv -o expr/gencode_v50.ref_quant.gtf
echo EXPR_DONE
