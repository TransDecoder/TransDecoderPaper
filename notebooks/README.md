# Notebooks

One R Markdown notebook per part of the paper. Each notebook reads committed tables in `../data/`,
draws the figures of its section (written to `../figures/`), and computes every number the section
quotes. The knitted HTML is committed next to the Rmd.

| Notebook | Manuscript sections | Figures / tables | Evidence |
|---|---|---|---|
| `01_algorithm.Rmd` | Methods: algorithm; Results: algorithm | Fig. 1; presets table; long-ORF lengths; frame-contrast test | E08, E10, E20–E22, E24–E26 |
| `02_benchmark.Rmd` | Methods: benchmark; Results: benchmark | Fig. 2; Table 1; Perl/Python concordance; run times | E07, E16, E25 |
| `03_start_refinement.Rmd` | Results: start refinement | Fig. 3; Figs. S1, S2 | E09, E13 |
| `04_short_orfs_qvalues.Rmd` | Results: minimum ORF length and target-decoy acceptance | Fig. 4; q settings; presets on small inputs; lncRNAs | E04, E18–E22, E25 |
| `05_coding_failures.Rmd` | Results: coding-test failures | Fig. 5 | E12 |
| `06_gm12878_application.Rmd` | Methods: GM12878 assemblies; Results: GM12878 application | Fig. 6 | E03, E05, E11, E25 |
| `07_gm12878_evidence.Rmd` | Methods: experimental evidence; Results: Ribo-seq and MS support | Ribo-seq / MS tables | E01, E02, E25 |
| `08_gm12878_parameters.Rmd` | Results: GM12878 parameters | Fig. 7 | E04, E06, E13, E25 |
| `S1_release_history.Rmd` | Supplementary Note 1 | Fig. S3; Tables S7, S8 | E23 |

Knit with `../scripts/knit_notebooks.sh [notebook.Rmd ...]` (environment: `../environment/rnb.yml`).

## Conventions
- Start with `source(file.path("R", "common.R"))`. It provides `read_tab("benchmark", "summary.tsv")`
  (reads `data/benchmark/summary.tsv`), `save_fig(plot, "fig2_benchmark", width, height)`,
  the shared theme, colours and labels, number formatters (`f3`, `pct`, `num`) and `evidence("E25")`
  (a link to the evidence record).
- Notebooks only read `data/`. Tables are produced by the scripts in `pipeline/` from the full runs;
  a notebook never reads files outside the repository.
- Section headings follow the manuscript's paragraphs, so the LaTeX can point to them
  (`% notebook: notebooks/04_short_orfs_qvalues.html#one-setting`).
- Numbers quoted in the text are written with inline R (`` `r f3(x)` ``), and each one is also passed
  to `check_text(label, computed, quoted)`; `check_report()` at the end lists any that disagree
  with the manuscript.
- Each section links to the evidence record(s) that give provenance (software versions, run commands).
