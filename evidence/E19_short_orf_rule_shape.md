# E19. Shape of the short-ORF coding-contrast threshold: removing the 100-codon cliff (2026-10-23)

## Claims supported (planned; would change the default if adopted)
- The release rule has a cliff at 100 codons. Below 100 codons an ORF needs frame contrast ≥ C·mu (C = 200); from 100 codons up it needs only contrast > 0.
  - Real CDSs show a step in stop sensitivity at the cliff: 0.60 at 75–99 codons vs 0.91 at 100–124 (complete.wLnc, pooled).
  - lncRNA predictions pile up just above it: 4,971 at 100–124 codons vs 877 at 75–99.
- Replacing the cliff with a smooth threshold that falls with ORF length removes the step. The "taper" rule is C·mu at 100 codons, falling linearly to 0 at the GC-based long-ORF auto-accept length. At C = 150 it:
  - about halves lncRNA false positives (18.9% → 10.7%);
  - costs 0.4 points of stop sensitivity on complete transcripts and 1 point on fragments;
  - leaves predictions on random sequences unchanged;
  - on GM12878 removes 40–46% of predictions that match no GENCODE protein, for a 0.1–0.5% loss of true positives.
- The shape matters less than C. Taper, inverse (C·100/L) and fixed (C at all lengths) give nearly the same trade-off curve:
  - taper is ahead by 0.3–0.8 points on complete transcripts at a matched lncRNA rate;
  - fixed is ahead on fragments at the strict end (C ≥ 200 equivalents).
- Brian asked for this on 2026-10-23: no `--precise` flag; a smooth length-dependent C instead of the cliff. Fable (consulted the same day) gave the rationale for the linear taper. The contrast is approximately a log-likelihood ratio, so the per-ORF false-positive rate is length-independent. The prior log-odds of "coding" rise about linearly with ORF length. So the required contrast should fall linearly and reach 0 where length alone suffices, which is the auto-accept length.

## Code
- PyTD branch `exp/short-orf-taper`, commit 652d5e4, on top of the final 44e417a; local only, not pushed. Worktree `work/trees/PY_taper`.
- New Predict option `--short-orf-rule {cliff,fixed,inverse,taper}`; `cliff` is the release default and identical to it.
- With no auto-accept length (non-standard stop codons, e.g. T. thermophila, or GC > 80%), taper ends at 300 codons (`TAPER_FALLBACK_CODONS`).
- `benchmark_v2/scripts/short_orf_rule_run.sh`: full, fresh runs (LongOrfs + Predict) in new directories. It never resumes; an incomplete directory is moved aside first. Brian asked for no reuse of earlier checkpoints. All 172 runs were fresh: no `_trash` directories exist, and every run records 652d5e4 in `pytd_version.txt`.
- Runs: `work/bench_v3/rules/<org>/<variant>/<rule>_c<C>/` and `work/bench_v3/rules/gm12878/<asm>/<rule>_c<C>/`; job list `data/benchmark/rules/jobs.txt`.
- Settings: cliff_c200 (= release defaults at the final commit), taper_c150/200/250, fixed_c150/200, inverse_c200/250. Inputs: 8 organisms × complete/fragments wRand, human and mouse × complete/fragments wLnc, and GM12878 LRAA/StringTie (6 settings).

## Scoring (results in `data/benchmark/rules/`)
- `score_negatives.py` → `rules_scores.tsv`; `plot_rules.py` → `rules.{pdf,png,medians.tsv}`. Comparator points (v5.7.1, CPAT, RNAsamba, GMST) come from `work/bench_v3/runs`. TD2 points are in a private copy only (`work/td2_private/rules/`).
- `sens_by_length.py` → `sens_by_length.tsv`: stop sensitivity by true CDS length, and negatives with a prediction by predicted ORF length.
- `score_protein_matches.py` → `gm12878_protein_matches.tsv`: TD2-preprint-style complete / 30-mer match / nonmatching against GENCODE v50 proteins.
- `gm_rule_support.py` → `gm12878_rule_support.tsv`:
  - ORFs are matched by transcript, strand and stop to the final_abinitio ORFs, which have the same transcripts (checked with `cmp`), and inherit their Ribo-seq (E01; ≥ 10 P-sites) and MS (E02) support.
  - It compares ORFs kept, dropped and added relative to cliff_c200.
  - Added ORFs have no support data.

## Results: medians (stop sensitivity / exact sensitivity / % negatives with a prediction)
| setting | complete.wLnc | fragments.wLnc | complete.wRand | fragments.wRand |
|---|---|---|---|---|
| cliff_c200 (release) | 0.976 / 0.905 / 18.9 | 0.858 / 0.753 / 18.7 | 0.957 / 0.897 / 0.48 | 0.815 / 0.708 / 0.11 |
| **taper_c150** | **0.972 / 0.901 / 10.7** | **0.848 / 0.747 / 10.4** | **0.952 / 0.892 / 0.48** | **0.810 / 0.711 / 0.16** |
| taper_c200 | 0.965 / 0.895 / 7.1 | 0.807 / 0.712 / 7.1 | 0.939 / 0.882 / 0.20 | 0.778 / 0.678 / 0.06 |
| taper_c250 | 0.959 / 0.889 / 5.7 | 0.768 / 0.678 / 5.5 | 0.928 / 0.870 / 0.12 | 0.744 / 0.654 / 0.04 |
| fixed_c150 | 0.967 / 0.897 / 11.2 | 0.846 / 0.747 / 10.7 | 0.949 / 0.888 / 0.86 | 0.820 / 0.730 / 0.34 |
| fixed_c200 | 0.957 / 0.887 / 6.0 | 0.795 / 0.704 / 5.8 | 0.935 / 0.876 / 0.33 | 0.767 / 0.684 / 0.10 |
| inverse_c200 | 0.964 / 0.894 / 6.8 | 0.796 / 0.702 / 6.7 | 0.937 / 0.880 / 0.19 | 0.771 / 0.669 / 0.06 |
| inverse_c250 | 0.957 / 0.887 / 5.3 | 0.760 / 0.671 / 5.2 | 0.926 / 0.868 / 0.11 | 0.734 / 0.649 / 0.04 |
| perl v5.7.1 | 0.969 / 0.764 / 17.9 | 0.792 / 0.653 / 17.6 | 0.946 / 0.758 / 0.30 | 0.771 / 0.623 / 0.09 |

cliff_c200 at the final commit reproduces the stored pytd_dev_v2 medians to the third decimal. On GM12878 its prediction and match counts are identical to pytd_dev_v2.

Per organism, taper_c150 vs cliff_c200:
- complete.wRand stop sensitivity drops by 0.1 (C. reinhardtii) to 1.1 points (T. thermophila, which uses the 300-codon fallback);
- random-sequence calls change by at most +0.15 points.

## Results: stop sensitivity by true CDS length (complete.wLnc, human + mouse pooled)
| codons | 50–74 | 75–99 | 100–124 | 125–149 | 150–199 | 200–267 | ≥ 268 |
|---|---|---|---|---|---|---|---|
| cliff_c200 | 0.265 | 0.603 | 0.909 | 0.938 | 0.972 | 0.989 | 0.999 |
| taper_c150 | 0.272 | 0.606 | 0.789 | 0.887 | 0.962 | 0.989 | 0.999 |
| taper_c200 | 0.176 | 0.476 | 0.713 | 0.855 | 0.956 | 0.989 | 0.999 |

lncRNAs with a prediction, by predicted ORF length:

| codons | 75–99 | 100–124 | 125–149 | 150–199 |
|---|---|---|---|---|
| cliff_c200 | 877 | 4,971 | 2,227 | 1,605 |
| taper_c150 | 1,144 | 1,427 | 1,074 | 1,149 |

## Results: GM12878 (LRAA / StringTie)
| setting | TP (complete + 30-mer) | complete | nonmatching |
|---|---|---|---|
| cliff_c200 | 33,587 / 36,547 | 20,692 / 17,940 | 1,005 / 4,119 |
| taper_c150 | 33,431 / 36,495 | 20,601 / 17,884 | 602 / 2,236 |
| taper_c200 | 32,886 / 35,388 | 20,335 / 17,659 | 348 / 1,404 |
| taper_c250 | 32,370 / 34,395 | 20,056 / 17,436 | 233 / 1,053 |
| fixed_c200 | 32,507 / 35,221 | 20,095 / 17,500 | 304 / 1,143 |
| inverse_c200 | 32,770 / 35,104 | 20,279 / 17,615 | 326 / 1,294 |

Ribo-seq frame-support rate among ORFs with ≥ 10 P-sites, and MS rate, for ORFs dropped relative to cliff_c200:

| setting | dropped ORFs | Ribo-seq rate | MS rate |
|---|---|---|---|
| taper_c150 | 902 / 2,850 | 0.33 / 0.31 | 0.27 / 0.10 |
| taper_c200 | 1,444 / 4,122 | 0.39 / 0.36 | 0.34 / 0.16 |
| fixed_c200 | 1,896 / 4,611 | 0.35 / 0.35 | 0.38 / 0.19 |

All ORFs score 0.65 / 0.61 (Ribo-seq) and 0.63 / 0.64 (MS).

## Caveats
- The lncRNA negatives are not all noncoding.
- Dropped ORFs have about half the Ribo-seq support rate of kept ORFs, so some are real.
- The Perl implementation does not yet have the rule.
- Whether to change the default is Brian's decision.

## Status update (2026-10-28)
Superseded, not adopted. The target-decoy q-value rule (E20) replaced both the 100-codon cliff and these
smooth shapes. The experimental `--short-orf-rule` options were removed when the decoy rule became the
release default (PyTD be23241; E24).
