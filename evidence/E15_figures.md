# E15. Manuscript figures (main Figs. 1-7, Supplementary Figs. S1-S2)

## Claims supported
The figures in `manuscript/sections/figures_main.tex` (PDFs copied to `manuscript/figures/`). Each figure
restates results that already have evidence records; this record says how each figure was drawn and checked.
Plan: Fable review of 2026-10-22 (main session scratchpad, `fable_review/fable_review.md`; not in git).

## Sources
| Figure | Script (`figures/`) | Data | Evidence | State |
|---|---|---|---|---|
| 1 algorithm | fig1_algorithm/make_fig1.py | methods_algorithm.tex; data/long_orf_threshold/final_thresholds.tsv | E10, E09, E13 | final |
| 2 benchmark | fig2_benchmark/make_fig2.py `<summary.tsv>` | `$P/work/bench_v3/scores/summary.interim_20261022.tsv` (fixed scorer) | E07 | **redraw** from the final summary.tsv when CodAn finishes (1 cell pending) |
| 3 start refinement | fig3_start/make_fig3.py | data/start_param_sensitivity/*; benchmark summary | E09, E13 | final |
| S2 start grid | fig3_start/make_fig3.py | same | E13 | final |
| 4 short ORFs | fig4_short_orfs/make_fig4.py | data/benchmark/{short_orf_summary,lnc_by_era.human}.tsv | E04 | final |
| 5 coding failures | fig5_failures/make_fig5.py (+ derive_dmel_gc3_rescue.py) | data/long_orf_failures/v7/*; homology_vs_abinitio.tsv (release) | E12, E25 | panels d,e from the release homology runs (2026-10-26) |
| 6 GM12878 | fig6_gm12878/make_fig6.py (+ derive_class_by_length.py) | data/gm12878/*; data/pred_class_by_length.tsv from `$P/work/gm12878_td/eval/status/` | E01-E05, E11 | final |
| 7 parameters | fig7_parameters/make_fig7.py | data/start_param_sensitivity/gm12878.tsv; data/gm12878/{gm12878_minaa_sweep,gm12878_homology_series}.tsv | E13, E04, E06 | final |
| S1 Kozak logos | kozak/make_kozak_figure.py (existing) | as in E09 | E09 | final (reused) |

Per-panel sources, computations and the checks against the manuscript text (with line numbers) are in
each directory's `NOTES.md`. Shared style: `figures/style.py`.

## Checks
- Fig. 2 medians reproduce Table tab:bench to 3 decimals and the ranges quoted in results_benchmark.tex.
- Each figure script prints or writes the plotted values (`fig2_values.tsv`, `fig3_values.tsv`, ...),
  and the figure agents compared them with every number the text quotes. Mismatches found were fixed in
  the text: MS 93-94% and 9-24%; plateau alpha 0.2-0.5; "w = 17 bases per bit"; GM12878 length bin
  "< 100 aa" (E04); failure attribution (33%, 32%, 33% = 26% no homology + 6.5% stronger other frame).
- Main session viewed every PNG.

## Caveats
- Fig. 7a and Fig. 3 grid values come from runs at w = 25 per nat (= 17.3 per bit); E13 addendum.
- Fig. 7b's -m 100 markers come from an earlier development run (`kozak_add`, no saved script; E04).
- Fig. 5c counts CDSs >= 300 nt; the 69% rescue in the text counts failures of all lengths.
- Fonts: DejaVu Sans (Arial is not installed on this host). Colour-blind check not run (no node).
- Cosmetic to-dos: Fig. 1a LongOrfs box text touches the border; Fig. 2f ANGEL label overlaps its marker;
  Fig. 5c legend crowds the y axis.
