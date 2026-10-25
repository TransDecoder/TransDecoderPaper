# Figure 4: minimum ORF length and the short-ORF rule

Script: `make_fig4.py [results_dir] [out_dir]` (defaults `data/benchmark`, this dir), run with
`$P/env/td/bin/python`. Output: `fig4_short_orfs.{pdf,png}`. Shared style: `figures/style.py`.
Evidence: E04 (panels a–c), E07 (panel d).

## Panels and data
| panel | content | source | computation |
|---|---|---|---|
| a | human and mouse `complete.wLnc`: x = lncRNAs with a prediction (%), y = stop-level sensitivity for true proteins of 50–99 aa (%). Solid path: LongOrfs minimum 100, 90, …, 50 codons. Dashed path: 50 codons → 50 + rule ×150, ×200, ×250. Star: default (50 + rule ×200). | `data/benchmark/short_orf_summary.tsv` (cfg m100…m50, m50+g150/200/250) | x = `neg` / lncRNA count (33,668 human, 31,467 mouse); y = `sens_short` |
| b | the same for `fragments.wLnc` | same | same |
| c | predictions per negative sequence (log y) at 100 codons, 50 codons, 50 + rule ×200: GC-matched random sequences (`complete.wRand`, 8 organisms ordered by G+C) and lncRNAs (`complete.wLnc`, human and mouse) | same | `neg` / number of negatives. Random: one per truth transcript (Table 1 gene counts; checked: 19,167 `rand_` sequences in human, 5,281 in *P. falciparum*) |
| d | share of human lncRNAs with a prediction by method, lncRNAs present in GENCODE v28 vs added later (`complete.wLnc`) | `data/benchmark/lnc_by_era.human.tsv` (block `complete.wLnc`) | `n_called / n_lnc`, methods sorted by the overall rate |

Configurations: `m100` is `pytd_dev_kozak_add` (LongOrfs `-m 100`) and `m50+g200` is `pytd_dev_m50_sc200`
(see `benchmark_v2/scripts/short_orf_summary.awk`; E04 software table). Check: `m50+g200` reproduces the
final-default run `pytd_dev_v2`. lncRNAs with a prediction are 6,943 (human) and 5,423 (mouse), and random
calls are 467 (*T. thermophila* complete.wRand), the same as `neg_pred_abinitio` in
`homology_vs_abinitio.tsv`.

## Checks against the text (`manuscript/sections/results_homology_length.tex`)
All of these match the table:
- L74: 0.013–0.159 calls per random sequence at 50 codons, from *C. reinhardtii* 228/17,621 to *S. pombe* 813/5,100.
- L76–77: human short-protein sensitivity 0.31 → 0.86; lncRNAs with a prediction 19.6% (6,608) → 63.6% (21,416); precision 0.737 → 0.468.
- L87–90: default 0.564; 6,943 (20.6%); precision 0.729 vs 0.737.
- L92–93: mouse 0.202 → 0.474; 5,088 → 5,423.
- L93–94: fragments short sensitivity 0.04–0.16 at 100 codons, from *P. falciparum* 0.036 to *T. thermophila* 0.164 (`fragments.wRand`); 0.21–0.57 with the rule, from *P. falciparum* 0.205 to *C. reinhardtii* 0.569.
- L95–96: ×150 0.648 / 7,941; ×250 0.479 / 6,699.
- L97–98: *P. falciparum* complete 0.108 with the rule vs 0.030 at 100 codons.
- L103–107: TransDecoder 22.0% vs 19.7%; RNAsamba 6.2% vs 10.9%; overall TransDecoder 20.6%, CPAT 10.5%, CPC2 7.9%, RNAsamba 9.1%.

## Discrepancies and notes
- No numerical discrepancies were found.
- Terminology: E04 and the summary file call the rule the "guard" (`g200`, `--short-orf-min-contrast`). The manuscript and this figure call it the "short-ORF rule" with a "factor".
- `neg` counts predictions on negatives, and the x axis is labelled "lncRNAs with a prediction". The two are the same here: 6,943 predictions equals 6,943 distinct lncRNAs, as counted by `compare_homology.py`. They would diverge only if a run reported more than one ORF per lncRNA.
- Panel d: RNAsamba scores only the forward strand. Its weights may have seen pre-v28 lncRNAs, which the text discusses. Panel d could go to the supplement if space is tight.
- The palette validator in the dataviz skill was not run, because `node` is not installed on this host. The colours are the shared `style.py` ones.

## Draft caption
```latex
\begin{figure}[t]
\centering
\includegraphics[width=\textwidth]{figures/fig4_short_orfs}
\caption{\textbf{Minimum ORF length and the short-ORF rule.}
\textbf{a,b}~Sensitivity for true proteins of 50--99~aa against the share of GENCODE lncRNAs given a
prediction, for human (dark) and mouse (light) complete transcripts (\textbf{a}) and fragments
(\textbf{b}) with lncRNA negatives. Solid lines: LongOrfs minimum of 100 to 50 codons (labels).
Dashed lines: a 50-codon minimum with the short-ORF rule at factors of 150, 200 and 250 (labels);
star, the default (factor 200).
\textbf{c}~Predictions per negative sequence at a 100-codon minimum, a 50-codon minimum and the
default (50 codons with the rule), on GC-matched random sequences (complete transcripts, organisms
ordered by G+C) and on lncRNAs.
\textbf{d}~Share of human lncRNAs with a prediction for each method, for lncRNAs already present in
GENCODE v28 and those added later (complete transcripts with lncRNA negatives).}
\label{fig:short_orfs}
\end{figure}
```
