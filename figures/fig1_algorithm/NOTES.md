# Figure 1: the TransDecoder algorithm

Script: `make_fig1.py [out_dir]`. Run with `$P/env/td/bin/python`. It writes `fig1_algorithm.{pdf,png}`
at 170 mm width. Panel (a) is a draft schematic that Brian may redraw.

## Panels and sources
- **(a) Workflow schematic.** No data. Each statement was checked against
  `manuscript/sections/methods_algorithm.tex` as of 2026-10-22:
  - m = 50 codons (l.32);
  - 10T candidates with T = 500 and CDS ≤ 5,000 nt, then low-complexity and near-identical
    removal (l.50-57);
  - 3-periodic fifth-order Markov chain against the background (l.65);
  - six scores S1-S6 and C = S1 − max(S2,S3) (ORF selection);
  - rules 1-3 and C ≥ 200μ below 100 codons (l.88-93);
  - ranking by homology hits, S1, then length, with one ORF per transcript by default
    (l.126-134);
  - start refinement and outputs/genome projection (Start-codon refinement; Outputs).
  The footer states Perl/Python identity (E08).
- **(b) L\*(GC).** Source: `analyses/long_orf_threshold/final_thresholds.tsv`, column q0.999
  (E10). It is drawn as the lookup the code uses: G+C rounded up to the next 5% step, 25% value
  below 25%, no rule above 80% (methods_algorithm.tex l.114-117).
  - Organisms are placed at the G+C of their *complete* input (E07 Table 1, `style.ORG_GC`) with
    the looked-up L\*:
    - P. falciparum 22.1% → 519 nt;
    - S. pombe 37.9% → 615 nt;
    - A. thaliana 42.2% → 663 nt;
    - D. melanogaster, mouse, human 50.9-52.3% → 804 nt;
    - C. reinhardtii 68.3% → 1,347 nt.
  - T. thermophila is shown as "rule off" (ciliate code).
  - Check: 519 nt at 25% and 2,436 nt at 80% match l.115.
- **(c) Start rule.** The worked example of l.175-177: |u| = 90, α = 0.4 (baseline 36),
  coding(u) = 45, w = 17. K = +1 gives 53 (trim), K = 0 gives 36 (keep), K = −1 gives 19 (keep). All
  match the text.

## Draft caption
\caption{\textbf{The \TD\ algorithm.}
\textbf{a}, Workflow. LongOrfs enumerates ORFs of at least 50 codons, including ORFs that run
off either end of the transcript. Predict trains a three-periodic fifth-order Markov model and a
start-context (Kozak) model on the longest non-redundant ORFs of the input. It scores every
candidate in its own frame, the two shifted frames and the reverse complement, and keeps ORFs that
satisfy any of three rules. It reports the top-ranked ORF per transcript and refines open 5$'$
ends.
\textbf{b}, Length above which an ORF is kept without a coding test, $L^\ast$, as a function of
transcriptome G+C (0.999 quantile of ORF lengths in i.i.d.\ random sequence). Dots: benchmark
organisms at the G+C of their complete transcripts. The rule is off above 80\% G+C and under the
ciliate code.
\textbf{c}, Start-codon refinement for an ORF with an open 5$'$ end. The segment $u$ upstream of
the first in-frame ATG is trimmed if its UTR score, $\alpha|u|+wK$, exceeds its coding score,
$C(u)/\mu$. Example: $|u|=90$\,nt, coding score 45, $\alpha=0.4$, $w=17$ bases per bit; the ATG
is accepted only with a favourable context ($K=+1$).}

## Discrepancies and notes
- None between the figure and the text.
- The L\* actually used by a run depends on the G+C of the run's own input. With random negatives
  (.wRand) the G+C is the same by construction, and with lncRNA negatives (.wLnc) it is slightly
  different. Panel (b) uses the complete-input G+C only.
- Font: Arial is not installed on this host, so DejaVu Sans is used; the PDF embeds it (Type 42).
