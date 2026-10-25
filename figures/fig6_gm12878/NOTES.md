# Figure 6: PyTransDecoder on three GM12878 assemblies, tested against Ribo-seq and MS

Files:
- `make_fig6.py` → `fig6_gm12878.{pdf,png}`. It reads committed tables only, plus `data/`.
- `derive_class_by_length.py` → `data/pred_class_by_length.tsv`. Its input,
  `$P/work/gm12878_td/eval/status/<asm>_final_abinitio.pred_status.tsv`, is not in git.

Build:

    cd $T/figures/fig6_gm12878
    $P/env/td/bin/python derive_class_by_length.py > data/pred_class_by_length.tsv
    $P/env/td/bin/python make_fig6.py

All runs are final defaults, ab initio (PyTD 48c870c, E03). Evidence records: E01, E02, E03, E04, E05.
The PDF is 151 × 117 mm, with fonts embedded. It uses DejaVu Sans because Arial is not installed on this host.

## Panels and sources

**(a) MANE Select recovery.** Source: `data/gm12878/gm12878_expressed_sensitivity.tsv`, MANE rows, tpm_min 0, 1 and 10.
- Hollow bar: `assembled_frac`.
- Filled bar: `sens_exact`.
- Number above each bar: `sens_exact_given_assembled`.

**(b) Class composition of distinct predicted CDSs by predicted protein length.** Source: `data/pred_class_by_length.tsv`.
- The classes are the same as in the Ribo-seq/MS tables: exact > compatible > other overlap > unannotated (no same-strand GENCODE CDS overlap).
- The `<100` bin includes the 78 / 273 / 1,251 predictions under 50 aa, whose start refinement moved the start downstream (E04).

**(c) Ribo-seq.** Source: `gm12878_riboseq_support.tsv`. Light bar: with_min_psites / n (≥ 10 P-sites). Dark bar: supported / n (frame-0 test).

**(d) MS.** Source: `gm12878_ms_peptide_support.tsv`, `frac_ge1`. Bars are predictions of ≥ 100 aa ("long"); diamonds are < 100 aa ("short").

**(e) Ribo-seq by expression.** Source: `gm12878_riboseq_stringtie_by_expression.tsv`, `frac_ge10_psites` by source-transcript TPM bin (StringTie).

## Checks against the manuscript text

- **Class totals.** The class totals of `data/pred_class_by_length.tsv` equal the `n` of every class in `gm12878_riboseq_support.tsv`, for example LRAA exact 12,021, compatible 5,340, other 5,676, unannotated 885. The `<100` rows equal the `:short` rows. The per-bin unannotated fractions equal the E04 table "No-GENCODE-overlap fraction by predicted length" (18.2%, 20.5%, 4.5%, 1.8%, 0.4% for LRAA, and so on).
- `results_gm12878_application.tex:13-14`: 23,922 / 34,404 / 72,323 distinct CDSs; 21% / 21% / 28% partial. ✓
- `results_gm12878_application.tex:53-61`:
  - TPM ≥ 1: StringTie exact 0.741; compatible 0.879 (StringTie) and 0.892 (Trinity). ✓
  - TPM ≥ 10: StringTie exact 0.814; compatible 0.898 and 0.913. ✓
  - Exact given assembled is 0.909–0.937 at TPM ≥ 1. ✓
  - LRAA exact is 0.51–0.58. ✓
- `results_gm12878_application.tex:72-75`: unannotated 3.7 / 12.0 / 15.5%; partial 18% / 24%; complete 10% / 12%. ✓
- `results_gm12878_evidence.tex:87-94`:
  - Exact ≥ 10 P-sites 0.892 / 0.834 / 0.841 → 83–89%. ✓
  - Frame test 0.588 / 0.536 / 0.535 → 53–59%. ✓
  - Unannotated 4.3 / 1.0 / 0.6%. ✓ Frame pass 0.2–1.5%. ✓
  - At 5–20 TPM: 0.893 vs 0.034. ✓

## Discrepancies between data and text

1. `results_gm12878_evidence.tex:97`: "92–94% of exact-match predictions of at least 100 aa" — the data give 0.9438 / 0.9342 / 0.9406, i.e. **93–94%**.
2. `results_gm12878_evidence.tex:99`: "Only 10–25% of unannotated predictions shared any peptide" — the data give 0.235 / 0.206 / 0.091, i.e. **9–24%**.
3. E04's length table, cited at `results_homology_length.tex:124`: its "50–99 aa" bin includes predictions under 50 aa. For example, LRAA n = 822 = all predictions under 100 aa, of which 78 are under 50 aa. Label it "<100 aa" or say so.

## Draft caption

    \caption{\PyTD\ on three GM12878 assemblies (final defaults, no homology evidence), tested against
    ribosome profiling and mass spectrometry from the same cell line.
    (a) MANE Select CDSs assembled in full (hollow bars) and predicted exactly (filled bars), for all
    MANE genes and for genes expressed at a gene TPM of at least 1 and 10. Numbers: fraction of fully
    assembled CDSs predicted exactly.
    (b) Distinct predicted CDSs by predicted protein length and relation to GENCODE: identical CDS
    (exact), a correct partial CDS (compatible), overlap with an annotated CDS that is neither
    (other overlap), and no same-strand overlap with any GENCODE CDS (unannotated). The $<$100 bin
    includes the few predictions shortened below 50~aa by start refinement.
    (c) Fraction of predicted CDSs with at least ten Ribo-seq P-sites (light) and passing the frame-0
    enrichment test (dark); n is the number of unannotated predictions.
    (d) Fraction of predicted CDSs with at least one MS peptide at 1\% FDR, for predictions of at
    least 100~aa (bars) and shorter ones (diamonds).
    (e) Fraction of StringTie-derived predictions with at least ten P-sites, by the expression of the
    source transcript.}
