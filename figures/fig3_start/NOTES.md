# Figure 3 (start-codon refinement) and Fig. S2 (alpha x w grid)

`make_fig3.py` draws both figures from committed tables only. Run it with `$P/env/td/bin/python figures/fig3_start/make_fig3.py`.

It writes:
- `fig3_start.{pdf,png}`;
- `figS2_start_grid.{pdf,png}`;
- `fig3_values.tsv`, every plotted value;
- the checks listed below, printed to stdout.

Nothing under `$P/work` is read, so no `data/` directory is needed.

Evidence records: E13 (grid, open ends, GM12878, Kozak distributions), E07 (v6.0.0 rows of `summary.tsv`) and E09.

**Kozak weight units.** The grid was run with the weight in nats (0, 5, 10, 25, 50, 100). Axes and labels give it per bit, w_bits = w_nats * ln 2, which is 0, 3.5, 6.9, 17, 35 and 69. The default of 25 per nat is 17.3 per bit; the released code uses 17. The E13 addendum shows the two give identical exact F1 to within 0.0003 on every input, and at most 0.13% of predictions differ.

## Panels

### Fig. 3a: open 5' ends without refinement
- **Source:** `analyses/start_param_sensitivity/results/open_end_mix.tsv`.
- **Plotted:** the fraction truncated / (truncated + utr_extended), per organism, for complete and fragments. The `other` column (543 complete, 563 fragments) is left out.
- **Check:** complete inputs have 0 truncated in every organism (41,127 UTR-extended). Fragments are 0.642 (*T. thermophila*) to 0.955 (*P. falciparum*); ALL is 0.749. This matches results_start_refinement.tex:8-14 ("75%", "64% ... 96%").

### Fig. 3b: exact F1 per organism, complete and fragments
- **Sources:**
  - no refinement: `optima.tsv` `norefine_f1`;
  - α 0.4 with w 0, and the final defaults (a0.4_w25): `grid_scores.tsv`;
  - v6.0.0: `benchmark_v2/results_v3/summary.tsv`, tool `perl_td_v600`. Its f1_exact is identical in `$P/work/bench_v3/scores/summary.interim_20261022.tsv` (fixed scorer); the doubled `fp` column does not enter f1_exact.
- **Checks:** the script prints these, and they match results_start_refinement.tex:16-29.

  | Quantity | Text | Script |
  |---|---|---|
  | Refinement gain, complete | 3.1 to 38.7 points, median 20.1 | same |
  | Effect on fragments | −4.3 to +11.4, positive in 6 of 8 | same |
  | w=0 below no refinement on fragments | 6/8 | same |
  | Kozak term, fragments | +1.3 to +4.3 | same |
  | Kozak term, complete | −0.5 to +2.9, mean +0.6 | −0.5 to +2.85 (rounded to 2.8), mean +0.6 |
  | v6.0.0 vs final, complete | 1.7 to 27.3 | 1.65 to 27.25 from summary.tsv; the text rounds half up |
  | v6.0.0 vs final, fragments | 1.8 to 14.1 | 1.77 to 14.05 |
  | Median exact F1, complete | 0.779 → 0.910 | same |

- The final defaults in the grid (a0.4_w25) are byte-identical to pytd_dev_v2 (E13 checks).

### Fig. 3c: mean exact F1 against α
- **Source:** `plateau.tsv`, column mean_f1 for complete, fragments and both. "Both" is the per-organism mean of complete and fragments, averaged over the 8 organisms.
- **Plotted:** w=25 per nat (17 per bit) solid, w=0 dashed, 13 α values on a square-root x axis.
- **Check:** at the default, complete is 0.8838, fragments 0.7805 and both 0.8322, matching E13 Results.

### Fig. 3d: Kozak score distributions
- **Source:** `kozak_score_distribution.tsv`, sets `true_start`, `internal_ATG` and `frag_firstATG`, in bits.
- **Drawn:** boxes span p25–p75 with the median marked; whiskers span p10–p90. The file has quantiles only, so no other summary is possible.

### Fig. S2a: benchmark heatmaps
- **Source:** `plateau.tsv`, mean_f1 ×100, for complete, fragments and both.
- **Drawn:**
  - colour spans the best value down to 5 points below it;
  - the best cell is in bold;
  - the red box marks the default (α 0.4, w 17 per bit).
- **Check:** the values equal the first column of the old `heatmap_summary.png`, except that the old plot printed 0.7805 as 78.1 and this one prints 78.0 (formatting of a half-value).

### Fig. S2b: GM12878 heatmaps
- **Source:** `gm12878.tsv`. Columns `MANE_sens_exact_given_assembled` and `pred_complete_start_precision` (the exact-start share of GENCODE-compatible complete calls), ×100, for LRAA and StringTie. The no-refinement value is in each panel title.
- **Check:** the values equal the old `gm12878_heatmap.png`, and E13 Results: StringTie 90.69/92.70/93.79/94.92 and LRAA 91.80/93.45/94.47/95.47 at w 25 per nat. They also match results_start_refinement.tex:55-62: 61.0 → 92.7, 59.1 → 93.5, α 1.0 at 77.0 and 81.5.

## Draft captions

```latex
\begin{figure}[t]
\centering
\includegraphics[width=\textwidth]{figures/fig3_start/fig3_start.pdf}
\caption{Start-codon refinement on the eight-organism benchmark.
(a) Without refinement, predicted ORFs whose stop matches a reference CDS but whose 5$'$ end is
open were either UTR extensions of a complete CDS (light) or truly truncated CDSs (dark). Upper
bars, complete inputs; lower bars, fragments.
(b) Exact-CDS F1 per organism without refinement (open circle), with v6.0.0 (grey square), with
the final rule without the Kozak term ($\alpha=0.4$, $w=0$; triangle) and at the final defaults
($\alpha=0.4$, $w=17$ per bit; filled circle). On fragments, the difference between v6.0.0 and the
final defaults also includes the lower minimum ORF length.
(c) Mean exact-CDS F1 over the eight organisms as a function of $\alpha$, at $w=17$ per bit
(solid) and $w=0$ (dashed), for complete inputs, fragments and their mean. Dotted line, default.
(d) Kozak scores of the self-trained start-context model at annotated start codons, at internal
in-frame ATGs of the same CDSs, and at the first ATG of 5$'$-truncated fragments. Boxes,
25th--75th percentiles; whiskers, 10th--90th percentiles; line, median.}
\label{fig:start}
\end{figure}

\begin{figure}[t]
\centering
\includegraphics[width=\textwidth]{figures/fig3_start/figS2_start_grid.pdf}
\caption{Sensitivity of start refinement to $\alpha$ and the Kozak weight $w$.
(a) Mean exact-CDS F1 ($\times100$) over the eight benchmark organisms for complete inputs,
fragments and the per-organism mean of both. (b) GM12878 StringTie and LRAA assemblies: exact
recovery of fully assembled MANE CDSs, and the fraction of GENCODE-compatible complete calls with
the exact start ($\times100$); the value without refinement is given above each panel. The grid
was run with $w$ in nats (0, 5, 10, 25, 50, 100) and is labelled per bit ($\times\ln 2$); the
default, 25 per nat, equals 17.3 per bit (released code: 17). Red box, default; bold, best value;
colour spans the best value to 5 points below it.}
\label{fig:S2}
\end{figure}
```

## Discrepancies with the text
1. **methods_algorithm.tex:195-196:** "Eight values (α=0.3–0.5, w=7–35) were within 0.5 percentage points".
   - In `plateau.tsv` (both), the 8 cells within 0.5 points of 0.8322, including the default, are a0.2_w25, a0.3_w10, a0.3_w25, a0.3_w50, a0.4_w10, a0.4_w25, a0.4_w50 and a0.5_w25.
   - So α is 0.2–0.5, not 0.3–0.5. The w range of 7–35 per bit is right.
   - E13 Claims item 1 has the same error (0.3–0.5); FINDINGS.md lists 0.2.
2. **results_start_refinement.tex:27:** "$w=17$ bits" should read "$w=17$ per bit". K is in bits; w is UTR length per bit.
3. **results_start_refinement.tex:22:** the Kozak term's maximum on complete inputs (human) is +2.85 points.
   - The text has +2.9; rounding the 4-decimal F1s gives 2.8 or 2.9. This is not an error.
   - Similarly, the v6.0.0 ranges of 1.7 and 27.3 come from exact half-values (1.65, 27.25).
4. **No data problem found.** The v6.0.0 exact F1 values are identical in the committed and fixed-scorer summaries.
