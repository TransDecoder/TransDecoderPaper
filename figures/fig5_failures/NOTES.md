# Figure 5: true CDSs that fail the coding test, and the two fallbacks

Scripts (run with `$P/env/td/bin/python`; shared style `figures/style.py`):
- `make_fig5.py [homology_vs_abinitio.tsv] [out_dir]` writes `fig5_failures.{pdf,png}`.
- `derive_dmel_gc3_rescue.py > data/dmel_gc3_rescue.tsv` builds the panel c table from
  `data/long_orf_failures/per_cds/dmel.per_cds.tsv.gz`. That file is **not in git**; only
  the `per_cds/*.meta.tsv` files are committed. `data/dmel_gc3_rescue.tsv` should therefore be
  committed with the figure.

Evidence: E12 (panels a–c, f) and E07 (panels d, e).

**Panels d and e** use the release homology comparison (E25), built with

    python benchmark_v2/scripts/compare_homology.py data/benchmark/final_v7/all_scores.tsv $P/data/bench_v3/inputs \
        $P/work/bench_v3/runs pytd_v7 pytd_v7_homology > data/benchmark/homology_vs_abinitio.tsv

## Panels and data
| panel | content | source | computation |
|---|---|---|---|
| a | % of matched true CDSs failing the coding test, by annotated CDS length (6 bins, log y), one line per organism. Bold: all organisms. Dashed: all organisms, still lost after the length rule. Zeros are not drawn. | `data/long_orf_failures/summary_by_length.tsv` | `fail_coding_pct` (includes the short-ORF rule below 300 nt); `still_lost_pct` for org ALL |
| b | exclusive attribution of failures ≥300 nt, per organism and all; n at the right | `explanation_flags.tsv`, rows `attr:*` | `fail_pct` |
| c | *D. melanogaster* CDSs ≥300 nt failing the coding test (grey), and the part rescued by a GC3-matched Markov model (blue), by GC3 quintile and for chromosomes 4 and Y; counts as fail/n | `data/dmel_gc3_rescue.tsv` | fail_score as in `summarize.py`; rescued = `pass_gc3match500 == 1` |
| d | homology minus ab initio PyTD (release, 3a793ff; `pytd_v7_homology` − `pytd_v7`), stop-level sensitivity, stranded `complete` and `fragments` inputs, in percentage points | `data/benchmark/homology_vs_abinitio.tsv` | 100 × `d_sensitivity_stop` |
| e | the same for exact-CDS F1 | same | 100 × `d_f1_exact` |
| f | failing true CDSs retained by the length rule, split by whether they also have a homology hit; labels give length-only/total. *T. thermophila* is omitted because the rule does not apply there (0 retained). | `length_vs_homology_fallback.tsv` | `also_homology_hit`, `length_only` |

Check for panel c: n and fail counts per quintile (2,635/356, 2,634/8, 2,635/3, 2,634/1, 2,634/1)
equal `fail_rate_by_feature.tsv` (dmel, gc3). Chromosome 4 (80/40, 39 rescued) and Y (22/8, 8
rescued) equal `fail_by_chromosome.tsv`.

## Checks against the text (`manuscript/sections/results_coding_failures.tex`, `results_homology_length.tex`)
**Coding failures:**
- L12–13: 4.9% of 135,429 fail (`summary_by_org` ALL 4.89). ✓
- L15–19: by length, 67.0%, 9.35%, 2.32%, 0.60%, 0.30% and 0.41%. ✓
- L28–29: 447 retained by the length rule, 85 (19%) with no homology hit. ✓
- L34: 4.6% not predicted (`fail_lost_pct_of_matched` 4.56). ✓
- L40–41: "33%, 32%, 26%; 2.4% unexplained". These match all-organism atypical GC3 32.85, low complexity 32.05, no homology 26.24 and unexplained 2.38. See discrepancy 1.
- L46: *D. melanogaster* lowest GC3 quintile 13.5% vs ≤0.3% in the others (13.51; 0.30, 0.11, 0.04, 0.04). ✓
- L47: chromosome 4 40/80, Y 8/22. ✓
- L51–52: *T. thermophila* ribosomal proteins 21 of 60 (`rp_fail_score`; `rp_fail` = 22 including the short-ORF rule). Not plotted (see note).
- L53–54: GC3-matched model rescues 25% overall and 69% in *D. melanogaster* (`oracle_models.tsv` 25.37 / 68.72, of coding-score failures of all lengths); 39 of 40 on chromosome 4. ✓ Panel c uses ORFs ≥300 nt, where the lowest quintile is 326/356 rescued.

**Homology:**
- L13–19: complete inputs (48, excluding wLnc) stop sensitivity +0.08 to +2.46 points, i.e. 0.1–2.5; fragments (32) +1.25 to +9.85. ✓
- L20–22: stranded fragments, exact F1 human +3.77, mouse +4.65, *A. thaliana* +4.02, *T. thermophila* +0.56, *C. reinhardtii* +0.60. ✓

## Discrepancies and notes
1. **Attribution percentages (L40–41).** The text gives three causes (33/32/26%) plus 2.4% unexplained. That omits the fourth exclusive class, "stronger ORF in another frame" (6.5% overall, 11.8% in *T. thermophila*). The text folds it into cause 3 ("A quarter of them had a stronger ORF in another frame"), but the 26% figure is the no-homology class alone, so the four quoted numbers sum to 93.5%. Suggested wording: "(33%, 32% and 33%, the last comprising 26% without homology and 6.5% with a stronger ORF in another frame; 2.4% unexplained)", or name four classes.
2. **Panel c and the 69% in the text use different denominators.** The 69% is of all coding-score failures; the panel shows ORFs ≥300 nt. The caption says so.
3. The *T. thermophila* ribosomal-protein result (21/60) is not plotted; adding it as a second unit in panel c would mix axes. It stays text-only (`ribosomal_proteins.tsv`), and it could go into a supplementary table.
4. Panels d and e must be regenerated once the final `summary.tsv` exists (see above). The numbers will only change if the rescoring changes `pytd_dev_v2` or `pytd_dev_v2_homology` rows; the fixed-scorer change affected only the `fp` column.
5. The palette validator in the dataviz skill was not run, because `node` is not installed on this host.

## Draft caption
```latex
\begin{figure}[t]
\centering
\includegraphics[width=\textwidth]{figures/fig5_failures}
\caption{\textbf{True coding regions that fail the coding test, and the two fallbacks.}
\textbf{a}~Percentage of annotated CDSs (complete transcripts, matched to a LongOrfs candidate) that
fail the coding test, by CDS length, for each organism (thin lines) and for all organisms (bold).
Dashed: all organisms, failures not retained by the GC-calibrated length rule. Failures below 300\,nt
are mostly due to the short-ORF rule.
\textbf{b}~Exclusive attribution of failing CDSs of at least 300\,nt; $n$ at right.
\textbf{c}~\textit{D.~melanogaster} CDSs of at least 300\,nt that fail the coding test, by GC3
quintile and for chromosomes 4 and Y (failures/total above bars); blue, the part that passes with a
Markov model trained on true CDSs of matching GC3.
\textbf{d,e}~Change in stop-level sensitivity (\textbf{d}) and exact-CDS F1 (\textbf{e}) when
homology evidence (Swiss-Prot without the organism's own entries, and Pfam) is added, on stranded
complete transcripts and fragments.
\textbf{f}~Failing true CDSs retained by the length rule, with and without a homology hit (labels:
length rule only / total). The rule does not apply to \textit{T.~thermophila}.}
\label{fig:failures}
\end{figure}
```
