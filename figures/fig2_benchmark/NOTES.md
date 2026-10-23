# Figure 2: accuracy on eight eukaryotes

Script: `make_fig2.py <summary.tsv> [out_dir]`. Run with `$P/env/td/bin/python`.

**Current draft is drawn from the interim summary**
`$P/work/bench_v3/scores/summary.interim_20261022.tsv`. That file was produced by score_all.py with
the fixed scorer, where `fp` is no longer doubled. It also includes the patched tthermophila
fragments.wRand homology row. Once the CodAn runs finish, regenerate summary.tsv and rerun the
script on it:

    $P/env/td/bin/python make_fig2.py $T/benchmark_v2/results_v3/summary.tsv

Outputs:
- `fig2_benchmark.{pdf,png}`;
- `fig2_values.tsv` (every plotted cell);
- `fig2_medians.tsv` (median, min and max over organisms, and run times).

## Panels and sources
All panels use stranded inputs only.

| panel | variant | metric |
|---|---|---|
| a | complete | f1_stop |
| b | complete | f1_exact |
| c | complete.wRand | precision_stop |
| d | fragments | f1_stop |
| e | fragments | f1_exact |
| f | complete | runtime_sec (median and range over 8 organisms) |

**Rows (summary tool names).** pytd_dev_v2 (TransDecoder final defaults; identical to
perl_td_dev_v2 on all 84 inputs, E08), perl_td_v600, prodigal, gmst, codan, cpat, cpc2, rnasamba,
borf, longest_orf.

**Omitted:**
- From the heatmaps: pytd_dev_v2_homology, perl_td_v571 (identical to v6.0.0, E07), angel,
  all_orfs, and the development variants.
- From panel f: longest_orf and all_orfs have no runtime.

**Colours and markings:**
- Heatmap colour: single-hue sequential, clipped at 0.4.
- Cells in the T. thermophila column are hatched for tools run with the standard code: gmst,
  codan, cpc2, rnasamba, borf (methods_benchmark.tex l.122-127). Their values (≤0.06) are in
  fig2_values.tsv.
- Grey "pend." marks a missing row in the summary. In this draft that is only codan
  tthermophila complete.wRand.

**Panel f:**
- Perl = perl_td_dev_v2 and Python = pytd_dev_v2.
- ANGEL is shown from its one complete-input run (A. thaliana, 17.4 h). Its other two runs were on
  other variants.
- CodAn runs were 6-16 at a time on a shared host, so its run time in particular is indicative
  only.

## Checks against the text (results_benchmark.tex)
- Table tab:bench (l.23-31) matches fig2_medians.tsv to 3 decimals for every row and column,
  including the T. thermophila zeros in the medians. Examples:
  - TD final 0.978/0.910/0.895/0.789/0.993;
  - v6.0.0 0.972/0.779/0.869/0.715/0.995;
  - Prodigal 0.937/0.643/0.936/0.782/0.660;
  - CPC2 0.964/0.935/0.571/0.181/0.977.
- l.38-40 (TD stop F1 0.940-0.988; v6.0.0 0.948-0.985; longest ORF 0.954-0.984) match the min/max
  in fig2_medians.tsv.
- l.43-46 (exact F1 0.742-0.928, median 0.910; v6.0.0 0.567-0.912, median 0.779) match.
- l.62 (Prodigal precision 0.344-0.866, median 0.660), l.64 (CPAT 0.629 on C. reinhardtii),
  l.73-76 (T. thermophila 0.940/0.840; v6.0.0 0.948/0.567; Prodigal 0.856/0.395; CPAT 0.884/0.851;
  precision 0.954/0.519/0.344) all match.
- Methods l.118-119 (ANGEL 17-22 h) is consistent with the 17.4 h shown.

## Draft caption
\caption{\textbf{Accuracy on eight eukaryotes.} Stranded inputs; organisms ordered by
transcriptome G+C (\textit{P.~falciparum} 22\% to \textit{C.~reinhardtii} 68\%).
\textbf{a,b}, Stop-level and exact F1 on complete transcripts (UTRs $\le$250\,nt).
\textbf{c}, Stop-level precision when one GC-matched random sequence is added per transcript.
\textbf{d,e}, Stop-level and exact F1 on fragments, each transcript cut once within its CDS.
Hatched: tools that cannot use the ciliate genetic code, run with the standard code on
\textit{T.~thermophila} (F1 or precision $\le$0.06). \TD: final defaults, identical in the Perl and
Python implementations.
\textbf{f}, Wall-clock run time on complete inputs, median (dot) and range (line) over the eight
organisms; ANGEL was run on \textit{A.~thaliana} only. Runs shared a 28-core host, so times are
indicative. Values for every input are in Supplementary Table~\todo{S4}.}

## Discrepancies and notes
- **CodAn.** The interim summary has CodAn on all 8 organisms for complete and fragments, but
  misses tthermophila complete.wRand (still running), so panel c has one pending cell. CodAn's
  median complete run time here is 20,839 s (5.8 h). The Fable review quoted 18,322 s from the
  older committed summary.tsv. Use the regenerated file.
- **Perl run time.** The Fable review gave the median Perl run time as 408 s, but that is
  perl_td_v600. The final Perl defaults (perl_td_dev_v2) have a median of 421.5 s, and Python 185 s.
  The ratio of medians on complete inputs is 2.28×, consistent with the "2.3×" quoted (which is the
  median over all 84 inputs, E08).
- **CodAn on P. falciparum.** Low values (stop F1 0.71 complete, 0.59 fragments; exact 0.38/0.47)
  are real data from the invertebrate model. The text does not mention this yet; it is worth a
  sentence once CodAn is written up.
- **GeneMarkS-T and CPC2.** Their complete.wRand precision is ≥0.99 in most organisms, which
  could be read as "better than TD". Table tab:bench already shows this (0.996, 0.977).
