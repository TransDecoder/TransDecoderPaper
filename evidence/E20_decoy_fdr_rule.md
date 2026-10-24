# E20. Target-decoy FDR for candidate ORFs, with shuffled transcripts as decoys (2026-10-24)

## Claims supported (planned; would replace the short-ORF and long-ORF length rules)
- Decoy ORFs come from dinucleotide-shuffled copies of the input transcripts. They are found with the same ORF finder, genetic code, minimum length and strand setting as LongOrfs, and scored with the run's own trained model. Within strata of ORF length × transcript length × ORF type (complete or partial), each candidate gets a target-decoy q-value for its frame contrast. A single q cutoff then replaces the 100-codon short-ORF rule (C·mu) and needs no tuned length curve. Unlike the taper (E19), the rule works for any genetic code: T. thermophila needs no fallback.
- At q ≤ 0.05, plus the GC-table long-ORF auto-accept, against the release rule (medians):
  - higher stop sensitivity on complete transcripts (0.965 vs 0.957), with fewer calls on random sequences (0.36% vs 0.48%);
  - fragments: +5.7 points stop sensitivity (0.872 vs 0.815) and +6.5 points exact (0.773 vs 0.708); random calls 0.22% vs 0.11%;
  - lncRNA calls cut about 3.5-fold (complete 18.9% → 5.0%; fragments 18.7% → 5.8%), at −1.2 points stop sensitivity on complete transcripts and −0.1 on fragments.
- Stratifying by transcript length was Brian's suggestion: fragments cap the possible ORF length. It is the key ingredient. In the prototype, ORF-length strata alone only matched the taper; adding transcript length gave +5–6 points on fragments.
- The lncRNA calls are not mostly hidden translated ORFs. Only 1–2% of predictions on GENCODE lncRNAs end at the stop of a GENCODE/TransCODE phase-2 Ribo-seq ORF, at any stringency and for every tool (`lnc_riboseq_orfs.tsv`). Brian's caveat stands, though: lncRNAs can encode peptides, so the text should say "predictions on annotated lncRNAs".

## Code
- PyTD branch `exp/short-orf-taper`, local only, worktree `work/trees/PY_taper`.
  - 1b9a291: decoy rule (`pytransdecoder/core/decoys.py`). LongOrfs records its settings in `longorfs.params.json`.
  - e3b9e35: `--decoy-long-orfs {none,gc,decoy}`.
- Predict options: `--short-orf-rule decoy --decoy-q Q [--decoy-strata len|len_tx|len_tx_type] [--decoy-shuffles N] [--decoy-seed 42] [--decoy-long-orfs none|gc|decoy]`.
  - The number of shuffles is automatic: enough copies for about 20 Mb of decoy sequence, at most 50 per transcript; 1 for all benchmark inputs and GM12878.
  - q-values are written to `decoy_qvalues.tsv`.
- Dinucleotide shuffle: Altschul–Erickson, with the random Eulerian path drawn by Wilson's algorithm. Tested on 2,000 random sequences: dinucleotide counts and the first and last bases are preserved.
- Runs: fresh, never resumed (`short_orf_rule_run.sh`, settings `decoy_q<Q>[_long<gc|decoy>]`).
  - 126 runs: q 0.01–0.3 on 20 inputs, plus GM12878 at q 0.05/0.1/0.2.
  - 84 long-ORF variants.
  - All finished; no runs were set aside.
- Earlier steps:
  - The oracle "curve" taper (annotation-derived log prior odds; `oracle_logodds_curve.py`, `results_v3/rules/curves/`; commit 5799c35; 64 runs) was no better than the linear taper (E19).
  - Prototype: `decoy_fdr_prototype.py`, post-filter on fresh `cliff_c0` runs; results in `results_v3/rules/decoy_proto/`.

## Results: medians (stop sensitivity / exact sensitivity / % negatives with a prediction)
| setting | complete.wRand | fragments.wRand | complete.wLnc | fragments.wLnc |
|---|---|---|---|---|
| cliff_c200 (release) | 0.957 / 0.897 / 0.48 | 0.815 / 0.708 / 0.11 | 0.976 / 0.905 / 18.9 | 0.858 / 0.753 / 18.7 |
| taper_c150 (E19) | 0.952 / 0.892 / 0.48 | 0.810 / 0.711 / 0.16 | 0.972 / 0.901 / 10.7 | 0.848 / 0.747 / 10.4 |
| decoy_q0.05 | 0.958 / 0.894 / 0.14 | 0.871 / 0.772 / 0.19 | 0.960 / 0.890 / 4.6 | 0.856 / 0.754 / 5.5 |
| decoy_q0.1 | 0.962 / 0.897 / 0.15 | 0.884 / 0.787 / 0.26 | 0.965 / 0.895 / 5.6 | 0.883 / 0.776 / 7.3 |
| decoy_q0.2 | 0.967 / 0.899 / 0.23 | 0.893 / 0.796 / 0.37 | 0.971 / 0.900 / 7.4 | 0.908 / 0.796 / 10.3 |
| **decoy_q0.05_longgc** | **0.965 / 0.899 / 0.36** | **0.872 / 0.773 / 0.22** | **0.964 / 0.894 / 5.0** | **0.857 / 0.755 / 5.8** |
| decoy_q0.1_longgc | 0.968 / 0.903 / 0.37 | 0.888 / 0.788 / 0.31 | 0.969 / 0.898 / 5.9 | 0.884 / 0.777 / 7.5 |
| decoy_q0.1_longdecoy | 0.969 / 0.903 / 0.42 | 0.889 / 0.789 / 0.32 | 0.968 / 0.898 / 5.8 | 0.884 / 0.777 / 7.5 |

### Long-ORF auto-accept
Without it, the decoy rule loses about 0.4% of long true CDSs. On complete.wRand, stop sensitivity at ≥ 400 codons is 0.9982 (release) vs 0.9950 (decoy). These are long ORFs whose frame contrast is ≤ 0, the case Brian raised.
- The GC table restores them (0.9984).
- The decoy-derived length (0.999 quantile of decoy ORF lengths) also restores them (0.9991), but admits many long random ORFs in T. thermophila: its TGA-only code gives 1,183 nt, and random-sequence calls rise from 43 to 333. The GC table gives no auto-accept for non-standard codes.

So: `gc`.

### Per organism (q0.1 vs release, `rules_scores.tsv`)
- Most organisms improve on both axes.
- Fragments gain 9–13 points of sensitivity, with more random calls in P. falciparum (0 → 0.19%), D. melanogaster (0.08 → 0.42%) and S. pombe (0.10 → 0.51%).
- T. thermophila complete: 0.907 → 0.891, with random calls 1.77% → 0.16%.

### Calibration (`decoy_calibration.per_input.tsv`)
Observed FDR is 1 − stop precision of the final predictions, median per variant.
- Random negatives: conservative. Nominal q 0.1 gives an observed FDR of 0.006 (complete) and 0.009 (fragments).
- lncRNA negatives: roughly nominal. q 0.1 gives 0.089 (complete) and 0.121 (fragments).

The q-values are computed over all candidate ORFs, before one is chosen per transcript.

### Sensitivity by CDS length (`sens_by_length.decoy*.tsv`)
- There is no step at 100 codons.
- Fragments improve most for short ORFs. At 50–74 codons, stop sensitivity is 0.27 (release) vs 0.62 (q0.1).
- Complete transcripts at 100–124 codons: 0.909 (release) vs 0.742 (q0.1).

## Results: GM12878 (LRAA / StringTie; TP = complete + 30-mer match; `gm12878_protein_matches.tsv`)
| setting | TP | nonmatching |
|---|---|---|
| cliff_c200 | 33,587 / 36,547 | 1,005 / 4,119 |
| decoy_q0.05 | 32,726 / 35,325 | 275 / 1,108 |
| decoy_q0.1 | 33,099 / 36,114 | 386 / 1,672 |
| decoy_q0.2 | 33,472 / 36,850 | 577 / 2,783 |

ORFs dropped relative to cliff_c200 have a Ribo-seq frame-support rate of 0.32–0.40 (all ORFs: 0.61–0.65) and an MS rate of 0.16–0.35 (`gm12878_rule_support.tsv`).

## Caveats
- The GM12878 long-ORF variants (q0.1_longgc / longdecoy) are run but not yet scored.
- The Perl port is not yet done.
- Presets for small inputs are planned. In the history smoke test, a 300-transcript input gave devel 31 of 100 random sequences with calls; to be checked.
