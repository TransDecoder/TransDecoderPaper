# E17. Published critiques of TransDecoder tested old vs new; bootstrap intervals (2026-10-23)

## Claims supported (planned Supplementary section "Earlier critiques"; not yet in the manuscript)
For each failure mode that earlier studies reported for TransDecoder (E16), these tests use our benchmark
inputs to compare the version those studies ran (v5.7.1; byte-identical to v6.0.0 on all 80 inputs, E07)
with the final defaults (v7.0.0 / PyTD 0.2.0, label `pytd_dev_v2`). They show which critiques the new
version addresses and which it does not. Brian's request (2026-10-23): address the earlier critiques in
the paper and show in the supplement that the new version fixes what the old one had.

## Sources
- `benchmark_v2/scripts/critiques_old_vs_new.py`: per input and tool, counts of multi-ORF transcripts, non-ATG (5'-open) predictions on coding transcripts, start upstream of the annotation (and running to the 5' end), wrong strand, coding transcripts without a prediction, and negatives with a prediction. The exact and stop TP counts agree with `score_predictions.py`.
- Results:
  - `benchmark_v2/results_v3/critiques/critiques_old_vs_new.interim.tsv`: 8 organisms × complete, complete.unstranded, fragments, fragments.unstranded, complete.wRand; tools v571, new, longest. Also has wLnc rows for `new` only.
  - `critiques_lnc.tsv`: human and mouse complete.wLnc and fragments.wLnc, tools v571, new, longest and longest_atg. The v571 wLnc runs were made 2026-10-23 with `run_tool.sh perl_td_v571` (tree `work/trees/TD_v571`, 5791fde). The longest baselines on wLnc come from v571's LongOrfs (identical to the pytd LongOrfs longest ORF on human complete; checked with `cmp`).
- `benchmark_v2/scripts/bootstrap_ci.py`: 1,000 paired bootstrap resamples of the sequences of each input (coding and negative), applied to all tools at once; percentile 95% CIs for each tool and for tool minus final defaults; medians over organisms per replicate. Point estimates equal `summary.tsv`. Results: `benchmark_v2/results_v3/bootstrap/interim_20261023.{per_input,median}.tsv`. These are interim: CodAn complete.wRand has 7 organisms; rerun when CodAn finishes.

## Results: critique → old (v5.7.1) vs new; median over 8 organisms unless noted
| Critique (source, E16) | Measure | v5.7.1 | new | longest ORF | Addressed? |
|---|---|---|---|---|---|
| 5'-extension, non-AUG starts (TranSuite; Borf; GMST "longer") | complete: predictions with no ATG start, % of predictions on coding transcripts | 15.7 (human 19.6) | 5.6 (human 5.2) | 26.3 | largely |
| same | complete: stop-TPs with start upstream of annotation, % | 17.4 (human 20.7) | 7.4 (human 6.7) | 29.6 | largely |
| same, fragments | fragments: start upstream of annotation, % of stop-TPs | 11.8 | 3.7 | 16.0 | yes |
| exact CDS (ORFanage 64%; GMST TIS 52%; CodAn) | complete exact, % of truth | 75.8 (human 74.5) | 89.7 (human 91.4) | 67.7 | yes, but the first-ATG rule is higher still (E07: F1 0.942 vs 0.910) |
| fragments / NoStart (CodAn: 0 correct) | fragments exact, % of truth | 62.4 | 70.8 | 62.1 | improved; CodAn's 0 looks like a scoring artifact (E16) |
| no ORF called (ORFanage: 16,484 of 72,909) | complete: coding transcripts without prediction, % | 4.6 | 4.1 | 3.5 | slightly; on fragments 22.5 → 17.7 |
| antisense picks (Borf: up to 15%) | unstranded: wrong-strand predictions, % of truth | complete 0.5, fragments 1.0 | 0.5, 1.0 | 2.2, 5.5 | already fine in v5.7.1; not reproduced |
| multiple ORFs per transcript (GMST Table S2; inGAP-CDG) | transcripts with > 1 prediction | 0 on coding inputs (1 in human complete.unstranded); 398 mouse / 916 human on lncRNA inputs | 0 | – | yes (one ORF per transcript by default) |
| worse than plain longest ORF (GeMoRNA) | complete exact F1 vs longest ORF | 0.779 vs 0.689 | 0.910 vs 0.689 | – | not reproduced on our inputs; but see first-ATG rule |
| lncRNA false positives (TD2) | lncRNAs with a prediction, % (human / mouse, complete.wLnc) | 19.6 / 16.1 | 20.6 / 17.2 | 28.7 / 24.1 | **no: slightly worse** (minimum ORF length 100 → 50 codons; E04) |
| random-sequence false positives | complete.wRand: random sequences with a prediction, % | 0.3 | 0.5 | 23.0 | already low; unchanged |
| speed (GMST, ORFanage, CodAn) | wall clock | Perl | PyTD 2.3× faster (E08) | – | partly |

Not testable on our inputs:
- Long downstream ORFs past PTCs in alternative isoforms (TranSuite): we use one transcript per gene.
- Genome-annotation completeness (Freedman & Sackton).
- Assembly redundancy (inGAP-CDG).

## Results: bootstrap 95% CIs (median over organisms)
- Intervals are narrow, about ±0.002–0.006, because each input has 5,000–27,000 transcripts. All differences between the final defaults and the other tools quoted in results_benchmark.tex exclude 0, except GeneMarkS-T fragments.wRand precision (diff 0.0001, CI −0.0008 to 0.0009).
- Examples:
  - complete exact F1: final 0.910 (0.906–0.914); longest ORF with first ATG 0.942 (0.940–0.945), diff +0.032 (0.028–0.037);
  - fragments exact F1: final 0.789 (0.783–0.795); Prodigal 0.782, diff −0.008 (−0.015 to −0.001);
  - fragments stop F1: CodAn 0.911, diff +0.016 (0.011–0.020).
- These intervals capture sampling of transcripts only, not the choice of organisms or annotation errors. With n = 8 organisms, a per-organism sign count is the more honest summary for claims across species.

## Caveats
- `critiques_old_vs_new.py` takes the first codon of each prediction from the input sequence; "non-ATG" therefore means a 5'-open ORF (the ciliate and standard codes both start at ATG in TransDecoder).
- The lncRNA rows cover human and mouse only (the only organisms with wLnc inputs).
