# Figure 7: parameter dependence on real assemblies (GM12878)

`make_fig7.py` → `fig7_parameters.{pdf,png}`. It reads committed tables only.

    cd $T/figures/fig7_parameters && $P/env/td/bin/python make_fig7.py

The PDF is 153 × 124 mm, with fonts embedded (DejaVu Sans). Evidence records: E04, E06, E13.

## Panels and sources

**(a) Start refinement.** Source: `analyses/start_param_sensitivity/results/gm12878.tsv`.
- Data: configs a{0.2,0.4,0.6,1.0}_w{25,0}, for LRAA and StringTie only, since the grid was not run on Trinity.
- x is `MANE_sens_exact_given_assembled`.
- y is `pred_complete_start_precision`: the exact-start share of complete calls compatible with a GENCODE CDS (E13).
- The weight in the grid is per nat. **w = 25 per nat is labelled "17/bit"**: 25 × ln 2 = 17.33 per bit. The shipped default, w = 17 per bit, differs from this by rounding only (E13 addendum; ≤ 0.13% of benchmark predictions differ).
- The circled point is the default (α 0.4).
- No refinement is off scale and goes in the caption: LRAA 0.591 / 0.817; StringTie 0.610 / 0.733.

**(b) Minimum length.**
- Lines: the post-hoc `--min-aa` filter on the final-default predictions. Source: `genome_annotation/results/gm12878_minaa_sweep.tsv`, which is headerless.
  - Column 5 is MANE exact given assembled.
  - Column 9 is the all-prediction no-GENCODE-CDS-overlap fraction.
  - Column definitions are in E04.
- Open markers: runs with a pre-defined LongOrfs minimum, for LRAA and StringTie only. Source: `genome_annotation/results/eval_summaries/{asm}_{m50,m50_sc200,kozak_add}.summary.tsv`.
  - Square: −m 50, no rule.
  - Star: −m 50 + rule; its summary equals the final defaults.
  - Triangle: −m 100 (`kozak_add`).
- Dotted line: 50 aa, the default.

**(c, d) Homology.** Source: `genome_annotation/results/gm12878_homology_series.tsv`.
- The value plotted is config minus `abinitio`. Panel c shows `MANE_sens_exact_given_assembled` × 100 (percentage points); panel d shows `distinct_predicted_cds` / 1000.
- Clade taxids: 9606 human, 9443 primates, 40674 mammals, 7742 vertebrates, 33208 animals (Metazoa).
- Line styles: solid = `blast_<taxid>_nopfam`; dashed = `blast_<taxid>_pfam`; diamond = `pfam` alone.
- The annotation in panel d gives the E06 count: 4,465 Trinity CDSs that are found with −human but not with −primates have chimpanzee UTY (Q6B4Z3) as best hit.

## Checks against the manuscript text

**(a)** `results_start_refinement.tex:53-64`:
- MANE exact given assembled, no refinement → default: 0.610 → 0.927 (StringTie), 0.591 → 0.935 (LRAA). ✓
- Exact-start share: 73.3 → 78.4% and 81.7 → 84.3%. ✓
- At α 1.0: 77.0% and 81.5%. ✓
- Kozak term (w 25 vs 0, α 0.4):
  - MANE: 0.9270 − 0.8957 = 3.1 points; 0.9345 − 0.9063 = 2.8 points. ✓
  - Exact start: 0.8 and 1.1 points. ✓

**(b)** `results_homology_length.tex:113-118, 126-128`:
- No-overlap with −m 50 and no rule: 10.7% and 34.7%. With the rule: 3.7% and 12.0%. ✓
- Post-hoc 150 aa:
  - no-overlap 1.1 / 4.3 / 3.3%; ✓
  - MANE exact given assembled 0.85 / 0.86 / 0.83, from 0.93 / 0.93 / 0.91. ✓

**(c)** `results_homology_length.tex:39-46`:
- −human, MANE exact given assembled: +1.83 / +1.10 / +1.41. ✓
- −human, compatible: +0.85 / +2.31 / +2.10. ✓
- −animals: +0.37 / +0.29 / +0.40 → "0.4, 0.3 and 0.4". ✓
- Pfam alone: 95.33 / 93.92 / 92.00 vs 95.28 / 93.80 / 91.96. ✓

## Discrepancies and caveats for the text

1. **Pre-defined −m 100 vs post-hoc 100 aa.** MANE exact given assembled is 0.928 (LRAA) and 0.920 (StringTie) with a pre-defined −m 100 (`kozak_add`), but 0.916 and 0.912 with the post-hoc 100 aa filter.
   - The figure shows this, and it is consistent with E04's caveat: training and per-transcript selection see the full candidate set.
   - However, `kozak_add` was run with older PyTD code and has no saved script (E04 caveats). The caption should say "earlier development code" or drop the triangles.
2. **Panel a is not at the shipped default.** It shows the grid at w = 25 per nat, not the shipped w = 17 per bit. The E13 addendum covers the equivalence. The caption says so.
3. **No numerical discrepancies with the text** in the numbers checked above.

## Draft caption

    \caption{Parameter dependence of \PyTD\ on the GM12878 assemblies.
    (a) Start refinement. Exact recovery of fully assembled MANE Select CDSs against the fraction of
    GENCODE-compatible complete predictions with the exact start, for $\alpha=0.2$, 0.4, 0.6 and 1.0
    with the Kozak weight $w=17$ per bit (solid; 25 per nat in the parameter grid) or $w=0$ (dashed).
    Circles: defaults. Without refinement the values are 0.59 and 0.82 (LRAA) and 0.61 and 0.73
    (StringTie), off scale. The start-parameter grid was not run on Trinity.
    (b) Minimum ORF length. Lines: final-default predictions filtered after prediction to a minimum
    protein length. Markers: runs with the minimum set in LongOrfs (square, 50 codons without the
    short-ORF rule; star, 50 codons with the rule, the default; triangle, 100 codons, run with earlier
    development code). Right: fraction of distinct predicted CDSs with no same-strand overlap with a
    GENCODE CDS.
    (c) Homology evidence. Change in exact recovery of fully assembled MANE CDSs over ab initio
    prediction when Swiss-Prot, with the indicated clade removed, is used alone (solid) or with Pfam
    (dashed), and with Pfam alone (diamonds).
    (d) Change in the number of distinct predicted CDSs for the same runs. Most of the Trinity increase
    with human proteins removed but other primates retained comes from short unannotated ORFs whose best
    hit is the Alu-encoded C-terminus of chimpanzee UTY.}
