# E05. Expressed-gene sensitivity (GM12878)

## Claims supported
Manuscript: Results, application to GM12878.

1. **Restricted to genes expressed in the sample, sensitivity is much higher.** With short-read RNA-seq from one cell line plus TransDecoder:
   - StringTie gives the exact MANE CDS for 74% (gene TPM ≥ 1) to 81% (TPM ≥ 10) of expressed MANE genes;
   - StringTie and Trinity GG give a correct full or partial (compatible) CDS for about 88–91%.
2. **Exact given assembled stays at 0.91–0.94** at every expression threshold. The gain from restricting to expressed genes comes from assembly completeness.
3. LRAA is lower (exact 0.51–0.58). Its long-read depth is lower, and the TPMs come from the short reads.

## Data
- **Alignments:** STAR alignments of ENCSR000AED rep 1 (`$P/data/gm12878/star_stringtie/GM12878_ENCSR000AED_rep1.STAR.bam`; records E03 and E11).
- **Annotation:** GENCODE v50 primary-assembly GTF (`$P/data/gm12878/ref/gencode.v50.primary_assembly.annotation.gtf`).
- **Per-reference status from E03:** `$P/work/gm12878_td/eval/{lraa,stringtie,trinity}_final_abinitio.ref_status.tsv`. These cover MANE Select v1.5 (19,229 on the primary assembly) and GENCODE complete (253,709).

## Software
- StringTie 3.0.3 (`$P/env/rnaseq_asm`).
- Python 3.11.17 (`$P/env/td`).

## Commands
Expression (2026-10-21 07:15, run from `$P/data/gm12878`). This was a one-off command, not a script. `$T/genome_annotation/scripts/data_prep/stringtie_expr.sh` reconstructs it, with the command copied verbatim from the header of the output GTF:
```bash
../../env/rnaseq_asm/bin/stringtie star_stringtie/GM12878_ENCSR000AED_rep1.STAR.bam --rf -e \
  -G ref/gencode.v50.primary_assembly.annotation.gtf -p 8 \
  -A expr/gencode_v50.gene_abund.tsv -o expr/gencode_v50.ref_quant.gtf
```
- `-e` restricts quantification to the reference transcripts.
- `--rf` matches the dUTP library.

Sensitivity:
```bash
cd $P; E=work/gm12878_td/eval
env/td/bin/python TransDecoderPaper/genome_annotation/scripts/expressed_sensitivity.py \
  data/gm12878/expr/gencode_v50.ref_quant.gtf \
  lraa=$E/lraa_final_abinitio.ref_status.tsv stringtie=$E/stringtie_final_abinitio.ref_status.tsv \
  trinity=$E/trinity_final_abinitio.ref_status.tsv > $E/expressed_sensitivity.tsv
```
Method (`expressed_sensitivity.py`):
- **Gene TPM** is the sum of the StringTie TPMs of the gene's transcripts.
- **Expressed:** a reference transcript counts as expressed when its gene's TPM is ≥ the threshold (default 1, 5, 10; 0 = all).
- Transcript IDs are matched without version suffixes.
- Assembled, exact, compatible and exact given assembled are as defined in E03 (`eval_genome_cds.py`), recomputed over the expressed subset.

## Outputs
- `$T/data/gm12878/gm12878_expressed_sensitivity.tsv` (identical copy: `eval_summaries/expressed_sensitivity.tsv`)
- Not in git: `$P/data/gm12878/expr/gencode_v50.ref_quant.gtf` (950 MB), `gencode_v50.gene_abund.tsv`, `stringtie_e.log` (contains only `EXPR_DONE`).

## Results
Expressed MANE Select genes: 11,076 of 19,229 at TPM ≥ 1; 9,360 at TPM ≥ 5; 7,823 at TPM ≥ 10.

MANE Select, gene TPM ≥ 1 (≥ 10 in parentheses):

| Assembly | Assembled | Exact | Compatible | Exact given assembled |
|---|---|---|---|---|
| LRAA | 0.544 (0.621) | 0.510 (0.580) | 0.655 (0.731) | 0.937 (0.935) |
| StringTie | 0.796 (0.874) | 0.741 (0.814) | 0.879 (0.898) | 0.931 (0.931) |
| Trinity GG | 0.688 (0.777) | 0.625 (0.707) | 0.892 (0.913) | 0.909 (0.909) |

All MANE (no filter): exact 0.306 / 0.430 / 0.362; compatible 0.395 / 0.569 / 0.547.

GENCODE complete (n = 253,709 all; 187,135 at TPM ≥ 1; 145,058 at TPM ≥ 10):

| Assembly | Exact, all → TPM ≥ 1 | Compatible, all → TPM ≥ 1 | Exact given assembled, TPM ≥ 1 |
|---|---|---|---|
| LRAA | 0.192 → 0.250 | 0.269 → 0.350 | 0.890 |
| StringTie | 0.234 → 0.314 | 0.352 → 0.415 | 0.888 |
| Trinity GG | 0.209 → 0.282 | 0.446 → 0.574 | 0.8775 |

The TPM ≥ 5 rows are in the table.

## Checks
- Re-running `expressed_sensitivity.py` on 2026-10-21 reproduced the committed table byte for byte.
- All reference transcripts in the `ref_status.tsv` files were found in the quantification (no "not in the quantification" warnings).
- The TPM-0 rows equal the E03 summary values (e.g. StringTie MANE exact 0.4300, exact given assembled 0.9270).

## Caveats
- **Gene-level TPM, not isoform-level.**
  - For MANE (one transcript per gene) this is adequate.
  - For GENCODE (all isoforms), an expressed gene does not mean every isoform is expressed. So the GENCODE rows rise far less than MANE (exact 0.23 → 0.31 for StringTie, vs 0.43 → 0.74 for MANE).
  - A transcript-level TPM filter would fit GENCODE better, but isoform quantification is itself uncertain.
- **The expression and the StringTie/Trinity assemblies come from the same reads.** The threshold therefore also selects genes with enough coverage to assemble.
- **Single sample;** the TPMs come from 97.5M read pairs.
- **The expression command was not scripted at run time.** `stringtie_expr.sh` is a reconstruction (see Commands).
