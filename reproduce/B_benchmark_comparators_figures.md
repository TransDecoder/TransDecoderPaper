# Reproduction detail B: benchmark, comparators, Perl/Python concordance, version history, mixed inputs, TD2, figures and supplementary tables

> **Paths into `$P/work`.** This section was written on 2026-11-03, before `$P/work` was deleted. Files that were kept are now under `$P/data/work_keep/` at the same path relative to `work/`. For example, `$P/work/bench_v3/runs/human/complete/pytd_v7_homology/pfam.domtblout` is now `$P/data/work_keep/bench_v3/runs/human/complete/pytd_v7_homology/pfam.domtblout`. Everything else under `$P/work` has to be regenerated with the steps given here. The tables in Part 2 record what was kept and why; see `../REPRODUCE.md` for the full list.

`$P=/broad/hptmp/bhaas/TRANSDECODER_PAPER`, `$T=$P/TransDecoderPaper`. Records covered: E07, E08, E15, E16, E23,
E27, E28, E37, E39. E25 (release reruns) and E40 (pinned-commit reruns) are upstream of most final tables and are
cited here, not repeated. Read-only survey made 2026-11-03.

## 0. Common setup (applies to every section)

**Code versions.**
- PyTD: paper results = d3892b5 (release pin; identical `pytransdecoder/` tree to tag v0.2.0 = 275ff1b and to
  a5b24c3). The `pytd_v7` labels were run at 3a793ff, whose *ab initio* output is byte-identical to d3892b5
  (E40: 108/108 inputs). PyPI `pytransdecoder==0.4.0` (tag v0.4.0 = a285464) is byte-identical to v0.3.0/v0.2.x
  (E42). **Run times differ** (0.3.0+ is faster), so `runtime_sec`, Fig. 2f and the run-time ratios (2.7x) will
  not reproduce.
- Perl: release pin 41fc5ee (`devel`). GitHub tag `TransDecoder-v7.0.0` (47cedc9) differs from 41fc5ee only in
  the version string, README, Changelog and Dockerfile (checked with `git diff --stat`), so either works.
  v6.0.0 = 66d4712 (tag `TransDecoder-v6.0.0` = ec18d65 differs only in `Docker/Dockerfile`); v5.7.1 = tag
  `TransDecoder-v5.7.1` (5791fde).

**Using a pip-installed / cloned PyTD in the repo scripts.** `benchmark_v2/scripts/run_tool.sh` runs
`PYTHONPATH=$PYTD_V7 $PYTHON -m pytransdecoder.cli longorfs|predict` and writes
`(cd $PYTD_V7 && git rev-parse --short HEAD) > pytd_version.txt` under `set -e`, so `PYTD_V7` must be a git
checkout. Equivalent to the paper code:
```bash
git clone https://github.com/TransDecoder/PyTransDecoder $X/PyTD && git -C $X/PyTD checkout v0.4.0   # or d3892b5
git clone https://github.com/TransDecoder/TransDecoder  $X/TD   && git -C $X/TD   checkout 41fc5ee
export PYTD_V7=$X/PyTD TD_V7=$X/TD PYTHON=$P/env/td/bin/python   # env/td has numpy 2.4.6 (needed from 0.3.0)
```
(`pytransdecoder.cli` still has the `longorfs`/`predict` subcommands at a285464.) A plain `pip install
pytransdecoder==0.4.0` alone is not enough for the scripts without editing the `git rev-parse` line.
All defaults in `run_tool.sh` (`ROOT=$P`, `TD_FIX`, `TD_V600`, `TD_V571`, `PYTD`, `PYTD_DEV` = `$P/work/trees/...`)
must be overridden by env vars once `$P/work` is gone.

**Homology mode caveat.** `run_tool.sh pytd_v7_homology` at 3a793ff = homology for *retention only*. At
d3892b5/0.4.0 the default adds the homology start assist; the paper's homology rows are now `pytd_v7h_homology`
(E33/E40) = d3892b5 default. To reproduce the old retention-only rows add `--no-homology-start-assist`.

**Comparator software** (paths are `run_tool.sh` defaults; envs in `$P/env`, ymls in `$T/environment/conda_envs/`):
Prodigal 2.6.3 (`env/td`); CodAn 1.1 (`env/codan`, bioconda codan-1.1; models `$P/tools/CodAn/models`, git
43da738 of github.com/pedronachtigall/CodAn); GeneMarkS-T 5.1 (`~/software/genemarks_t/dist`, **licensed; key in
`~/.gm_key`; neither software nor key may enter any repo**); CPAT 3.0.5 (`env/cpat`, models
`$P/tools/models/cpat/{Human,Mouse,Fly}_*`); CPC2 1.0.1 (`env/cpc2`); Borf 1.2 (`env/borf`); RNAsamba 0.2.5
(`env/rnasamba`, weights `$P/tools/models/rnasamba/{full,partial}_length_weights.hdf5`); ANGEL 3.0
(`$P/tools/ANGEL`, git b805554, `env/angel`); TD2 1.1.0 + PSAURON 1.1.3 (`env/td2`, `pip install td2==1.1.0`,
freeze `benchmark_v2/td2_env.pip-freeze.txt`). DIAMOND 2.2.8, HMMER 3.4. `PERL_HASH_SEED=0` everywhere.
`$P/tools` is outside `$P/work` and must also be kept.

---

## PART 1. Experiments

### E07 Cross-species benchmark bench_v3 (status: in paper — Table 1, Fig. 2, Supp. Table S5; TD rows now from E25/E40)
- **Feeds**: Methods "Benchmark data/Compared methods/Accuracy measures"; Table 1; Fig. 2; Fig. S4 (ranks);
  Supp. Tables S5, S9; comparator rows of every benchmark table. The TD rows in the paper are the release labels
  (`pytd_v7`, `perl_td_v7`, `pytd_v7h_homology`; E25/E40); E07's own TD labels (`pytd`, `pytd_dev_v2`, ...) are
  development history.
- **Inputs**: 8 annotations, URLs/releases in `benchmark_v2/annotation_downloads.tsv` (= `$P/data/annot/sources.tsv`),
  downloaded 2026-10-18 to `$P/data/annot/<org>/`. Two URLs float (`MANE_human/current/...v1.5...`, TGD
  `downloads/tet/latest/...`), so `$P/data/annot` is the authoritative copy. Truth `$P/data/bench_v3/truth{,_fullutr}/`
  (md5s of `truth/<org>.truth.tsv` in E07); inputs `$P/data/bench_v3/inputs/<org>/<variant>.{fa,truth.tsv}` (10
  variants x 8 orgs, seed 1; + `complete.wLnc`, `fragments.wLnc` for human/mouse, from
  `analyses/lnc_holdout/pick_lnc.py` -> `$P/work/lnc/{human,mouse}_lnc.pick.fa`, see E04). Homology DBs
  `$P/data/homology/` (Swiss-Prot 2026_03 minus own taxid `sprot_minus_<taxid>.dmnd`; Pfam-A 38.2).
- **Software**: see §0. Comparators run once (2026-10-18..27); no TD version dependence. `longest_orf` baseline =
  longest LongOrfs ORF of the `pytd` run (27bfb5e, -m 100); E40 showed it equals the longest >=100-codon ORF of
  release LongOrfs (`m100_noq`), so 0.4.0 `longorfs -m 100` reproduces it. Dev labels need their commits
  (all on GitHub: PyTD 27bfb5e, 48c870c, 44e417a on origin/main; Perl 351634a, de3f57e on origin/devel).
- **Steps** (from `$P`; scripts `$T/benchmark_v2/scripts/`):
  ```bash
  bash TransDecoderPaper/benchmark_v2/scripts/download_annotations.sh data/annot      # needs data/annot/sources.tsv
  PYTHON=env/td/bin/python TransDecoderPaper/benchmark_v2/scripts/build_all_truth.sh data/annot data/bench_v3
  PYTHON=env/td/bin/python TransDecoderPaper/benchmark_v2/scripts/make_all_inputs.sh data/bench_v3
  # comparators (copies of the launchers are in benchmark_v2/launch_v3/; all hardcode cd $P)
  JOBS=16 TransDecoderPaper/benchmark_v2/scripts/run_all.sh data/bench_v3/inputs work/bench_v3/runs prodigal codan gmst
  JOBS=6  .../run_all.sh data/bench_v3/inputs work/bench_v3/runs cpat cpc2 borf                       # launch3.sh
  VARIANTS="complete complete.wRand complete_fullutr fragments fragments.wRand" JOBS=6 .../run_all.sh ... rnasamba
  VARIANTS="complete.wLnc fragments.wLnc" JOBS=4 .../run_all.sh ... cpat cpc2 borf rnasamba gmst prodigal  # launch14.sh
  # v6.0.0 comparator: TD_V600=<clone at 66d4712> .../run_all.sh ... perl_td_v600
  # release TD (E25 launcher, NOT in repo: $P/work/bench_v3/launch_v7.sh):
  JOBS=6 .../run_all.sh data/bench_v3/inputs work/bench_v3/runs pytd_v7 perl_td_v7
  VARIANTS="complete.wRand fragments.wRand complete.wLnc fragments.wLnc" JOBS=6 .../run_all.sh ... pytd_v7_precise pytd_v7_sensitive
  THREADS=4 JOBS=3 .../run_all.sh ... pytd_v7_homology
  env/td/bin/python TransDecoderPaper/benchmark_v2/scripts/score_all.py data/bench_v3/inputs work/bench_v3/runs > summary.tsv
  ```
  `run_all.sh` with no tool list runs the 2026-10-18 tool set incl. ANGEL (dropped, 17–22 h/run) — always name tools.
  `score_all.py` walks every run dir; filter out `fragments_mixed_t*` (scored by `score_mixed.py`, E28).
  `launch.sh`/`launch2.sh` (pytd, perl_td, v571 ...), `launch5/9–13/15–17` are dev runs (E04, E09, E10).
- **Outputs**: `data/benchmark/summary.tsv` (dev; `fp` doubled except CodAn rows), `summary.interim_20261025.tsv`
  (rescoring with fixed scorer, 79 CodAn rows), `summary.release.tsv` (= release + comparator rows of interim_20261025;
  read by notebook 02), `final_v7/all_scores.tsv` (via `data/benchmark/final_v7/score_v7.sh <org> <variant> <tool>`,
  which hardcodes `$P/work/bench_v3/runs`), `homology_vs_abinitio.dev_v2.tsv`, `homology_vs_abinitio.tsv`
  (`compare_homology.py`, E25), `final_v7/first_atg_baseline.tsv` (`pipeline/derive/02_first_atg_baseline.py`,
  reads `work/bench_v3/runs/*/*/pytd_v7/...longest_orfs.cds`, writes `work/bench_v3/first_atg_v7`). Consumers:
  notebook 02 (Table 1, Fig. 2, Fig. S4), `09_supplementary_tables.py` (S5, S9).
- **Cost** (sum of `runtime_sec.txt`, computed 2026-11-03): CodAn 1,629 h over 100 runs (max 112 h), Perl v7 27 h/108,
  PyTD v7 9.7 h/108, `pytd_v7_homology` 67 h/84, GMST 1.5 h, CPAT 1.5 h, Prodigal 1.3 h. ANGEL 17–22 h per run (3 runs).

### E08 Perl/Python concordance and runtime (status: in paper — Results "two implementations", Supp. Table S2)
- **Feeds**: Results "Two implementations with identical output"; Methods §Implementations; Supp. Table S2;
  Abstract/Results speed-up (2.7x; from E25 timings, not reproducible with 0.3.0+).
- **Inputs**: the 5 sample sets in `TransDecoder/sample_data` (extracted by `regression/run_regression.sh`);
  benchmark inputs (E07).
- **Software**: final: Perl 41fc5ee + PyTD d3892b5 (E40). Historical rows: v6.0.0 66d4712 vs PyTD 43f4293
  (`baseline_v6.0.0`), r4/r6 (8cf8d5f/694e204 vs dc7b514/f5b12e7), 48c870c/de3f57e, 44e417a/ab38ca6 — all on GitHub.
  With 0.4.0 + 41fc5ee the release rows reproduce (outputs byte-identical); timings will not.
- **Steps**:
  ```bash
  PYTHON=$P/env/td/bin/python bash $T/regression/run_regression.sh <TD clone @41fc5ee> <PyTD clone @v0.4.0> <outdir>
  $P/env/td/bin/python $T/regression/parity_check.py ...        # intermediates (args not recorded)
  cp <outdir>/summary.tsv $T/data/regression/release_d3892b5_41fc5ee.summary.tsv
  # benchmark parity (84/84): pipeline/derive/01_perl_pytd_identity.py (reads $P/work/bench_v3/runs, hardcoded)
  ```
- **Outputs**: `data/regression/*.summary.tsv` (7 files), `data/release_pin/parity_intermediates.txt`,
  `data/benchmark/final_v7/perl_pytd_identity.tsv` -> notebook 01, `09_supplementary_tables.py` (S2).
- **Cost**: minutes per set (Perl 3–285 s, Python 2–95 s).

### E16 How comparator papers benchmarked (status: background, literature only)
- **Feeds**: Background/Discussion; E16 table. No computation. Sources were read via paperclip/WebFetch into a
  session scratchpad (`papers/`, `codan_S2.tsv`; **not preserved**). DOIs/PMCs are in the record; re-read from
  those. Nothing in `$P/work`.

### E23 Version history of key releases (status: supplement — Supp. Note 1, Fig. S3, Tables S10/S11)
- **Feeds**: Supplementary Note 1; Fig. S3 (`figS3_history`); Tables S10 (releases), S11 (history); notebook S1.
- **Inputs**: 20 benchmark inputs (complete.wRand, fragments.wRand, complete.wLnc, fragments.wLnc x orgs; job list
  `data/benchmark/history/jobs.txt`, absolute `$P/...` paths, outdirs `$P/work/bench_v3/history/...`).
- **Software**: r2012-08-07 from SourceForge (archive `$P/work/history/sourceforge/archives/transdecoder-r2012-08-07.tgz`,
  Brian's download; script-only copy `$P/work/history/sf_runs/r2012-08-07/`); v2.0, v4.1.0, v5.0.0 = GitHub tags
  (c4ea24d, 8f0e552, 3b4f9b1; worktrees `work/trees/TD_hist_<tag>`). Releases <= v5.2 need Perl DB_File:
  `work/history/compat/{perl5,lib}` = symlinks into `$P/env/mamba/pkgs/.../perl-db_file-1.858...` and
  `.../libdb-6.2.32...` (recreate if lost; `data/benchmark/history/VERSIONS.md` §5). v2.0 needs cd-hit-est (env/td).
  v5.7.1 (tag), Perl devel, and the v7 rows come from the main benchmark / `rules/*/{cliff_c200,decoy_q0.05_longgc}` /
  release runs (E25). The `cliff_c200` row is an unreleased dev version (no release option; E40 exception 2).
- **Steps**:
  ```bash
  # hist_run.sh hardcodes $ROOT/work/history/sf_runs/r2012-08-07, $ROOT/work/history/compat, $ROOT/work/trees/TD_hist_$ver
  cd <history run dir>; xargs -P 6 -L 1 bash -c '$T/benchmark_v2/scripts/hist_run.sh "$0" "$1" "$2" "$3"' < jobs.txt
  #   (launcher $P/work/bench_v3/history/launch.sh — not in repo)
  # r2012 3'-partials trimmed: benchmark_v2/scripts/clip_preds_to_transcript.py -> preds.clipped.tsv
  # scores: score_negatives.py-style table -> history_scores.tsv (exact command NOT recorded)
  $P/env/td/bin/python benchmark_v2/scripts/plot_history.py data/benchmark/history/history_scores.tsv data/benchmark/history/history
  ```
- **Outputs**: `data/benchmark/history/{history_scores.tsv,history.medians.tsv,history.{pdf,png},VERSIONS.md,jobs.txt}`;
  `notebooks/S1_release_history.Rmd` -> `figures/figS3_history.pdf`, Tables S10/S11 (hand-entered in supplement.tex).
- **Cost**: 80 runs, 28 h single-thread total (v2.0 up to 6,005 s/run, cd-hit-est).

### E27 How the comparators work (status: in paper — Background, Supp. Table S1)
- **Feeds**: Background comparator paragraph; Supp. Table S1 (`tab:comparators_detail`); n/a treatment of
  *T. thermophila* for GMST/CodAn/CPC2/RNAsamba/Borf.
- **Inputs**: installed code (`$P/env/*`, `$P/tools/*`, `~/software/genemarks_t`), papers in `$P/literature_pdfs/`,
  per-transcript tallies from `work/bench_v3/runs/.../preds.tsv`.
- **Steps**: code reading (no runs). Tallies by scratch scripts `by_type.py`, `start_codons.py` (**not committed**).
  Table: `cd $T && $P/env/td/bin/python pipeline/derive/02_comparator_table.py > data/benchmark/comparator_table.tex`
  (reads `comparator_features.tsv`, `final_v7/table1.tsv`, `final_v7/partial_types.tsv`).
- **Outputs**: `data/benchmark/comparator_features.tsv` (hand-curated), `comparator_table.tex` -> supplement.tex S1.

### E28 Balanced mixed partial inputs (status: in paper — Table 1 "mixed", Supp. Tables S3/S4)
- **Feeds**: Table 1 mixed columns; Results "Balanced types of coding region"; Supp. Tables S3 (`tab:partials`),
  S4 (`tab:mixedtypes`); Discussion.
- **Inputs**: `$P/data/bench_v3/inputs/<org>/fragments_mixed_t{1,2,3}.{fa,truth.tsv}` (24; md5 prefixes in E28).
- **Software**: `run_tool.sh` labels pytd_v7 (3a793ff ≡ d3892b5 ≡ 0.4.0 output), perl_td_v7 (41fc5ee),
  perl_td_v600, prodigal, gmst, codan (`*_partial`), cpat, cpc2, rnasamba (partial weights), borf; TD2 (E39).
  E40 confirmed 24/24 mixed inputs byte-identical at d3892b5.
- **Steps**:
  ```bash
  cd $P
  env/td/bin/python TransDecoderPaper/pipeline/benchmark/make_mixed_inputs.py data/bench_v3/inputs   # seeds 1,2,3
  bash work/bench_v3/launch_mixed.sh      # NOT in repo; two lanes: codan xargs -P2, others xargs -P2 (see Part 3)
  env/td/bin/python TransDecoderPaper/pipeline/benchmark/score_mixed.py data/bench_v3/inputs work/bench_v3/runs
  env/td/bin/python TransDecoderPaper/pipeline/derive/02_partial_types.py data/bench_v3/inputs work/bench_v3/runs \
     > TransDecoderPaper/data/benchmark/final_v7/partial_types.tsv
  cd $T && $P/env/td/bin/python pipeline/derive/02_mixed_table.py
  ```
- **Outputs**: `data/benchmark/final_v7/{mixed_scores,partial_types,five_prime_by_distance,mixed_medians}.tsv` ->
  notebook 02; S3/S4 typed into supplement.tex from these tables (no generator).
- **Cost**: 225 runs; mean per run pytd_v7 3.3 min, perl_td_v7 8.7, v6.0.0 5.2, CodAn long (finished 2 days later).

### E37 Supplementary Tables S1–S9 (status: supplement)
- **Feeds**: `manuscript/supplement.tex` S1–S9. No runs.
- **Steps**: `cd $T && $P/env/td/bin/python pipeline/derive/09_supplementary_tables.py > data/supplementary_tables/latex_bodies.txt`
  then paste the rows (each table has a `% generated by` comment). S1 from `02_comparator_table.py`; S7 from
  `pipeline/derive/08_homology_series_v7h_table.py` (E06).
- **Outputs**: `data/supplementary_tables/{perl_python_sample_sets,benchmark_per_input,start_choice_legacy,
  long_orf_length_by_gc,benchmark_truth_sets}.tsv`, `latex_bodies.txt`. All inputs committed.

### E39 TD2 as a comparator (status: in paper — Table 1, Figs 2, 4d, S4–S7, Tables S1, S3–S5)
- **Inputs**: 84 main inputs, 24 mixed, 12 simulated assemblies (human and *A. thaliana* d1–d32 from E30–E32,
  `$P/work/simtrinity{,_athal}/...`). Job list `benchmark_v2/launch_v3/td2_jobs_v1.tsv` (absolute `$P/...` paths).
- **Software**: `td2==1.1.0`, `psauron` 1.1.3 (pip). No TD version dependence. Set `PYTHONHASHSEED` (not set for
  the published runs; unstranded ties can flip, 3/3,957 on S. pombe; TD2_REPRODUCE.md).
- **Steps**:
  ```bash
  bash $T/benchmark_v2/launch_v3/launch_td2_v1.sh        # hardcodes P=...; reads $P/work/bench_v3/td2_jobs_v1.tsv
  #   -> benchmark_v2/scripts/run_td2.sh <input.fa> <outdir> <org>   (BENCH_ROOT, PYTHON, TD2_ENV overridable)
  #   TD2.LongOrfs -t input.fa [-S] -G <table> -@ 2 ; TD2.Predict -t input.fa -G <table>; OMP/MKL threads 2
  bash $T/pipeline/derive/10_td2_tables.sh                # appends td2 rows; refuses if td2 rows already present
  $P/env/td/bin/python $T/pipeline/derive/10b_td2_mechanisms.py
  ```
  `10_td2_tables.sh` reads `$P/work/bench_v3/runs`, `$P/work/simtrinity/{runs,depth/*/runs,filtered_preds}`,
  `$P/work/simtrinity_athal/depth/*/runs`, and calls `04_release_short_orfs.sh` (reads `$P/work/lnc/human_lnc.pick.gtf`
  and `$P/work/bench_v3/rules`). A from-scratch rebuild needs the pre-td2 tables first.
- **Outputs**: td2 rows in `final_v7/{all_scores,mixed_scores,partial_types,five_prime_by_distance,mixed_medians,
  table1,lnc_by_era.human.v7,td2_fragments_by_cds_length,td2_random_calls_by_rule}.tsv`, `data/simtrinity*/`,
  `supplementary_tables/benchmark_per_input.tsv` -> notebooks 02, 04.
- **Cost**: 120 runs, 2026-10-31 12:01–15:05 at 4 jobs x 2 threads; 108 bench/mixed runs sum 11 h.

### E15 Manuscript figures (status: in paper) — see table below. The `figures/fig*/make_fig*.py` matplotlib
scripts are **superseded**: every manuscript figure is now drawn by a notebook (`save_fig` in
`notebooks/R/common.R`) into `$T/figures/`, then copied by hand to `manuscript/figures/` (S5–S7 renamed).

### Figures and supplementary tables
Regenerate all: `$T/scripts/knit_notebooks.sh [nb.Rmd ...]` (Rscript `$P/env/rnb`, env `environment/rnb.yml`);
then `cp figures/<name>.pdf manuscript/figures/<name>.pdf` and `cmp` (EXPERIMENT_LOG rule). Notebooks read only
`$T/data/` (README: "a notebook never reads files outside the repository"; confirmed by grep).

| Item | Producer | Main input tables (`data/...`) |
|---|---|---|
| Fig 1 algorithm | nb 01 `fig1_algorithm` | final_v7/all_scores, frames/*, release_pin/sens_by_length, long_orf_threshold/* |
| Fig 2 benchmark | nb 02 `fig2_benchmark` | benchmark/summary.release.tsv, final_v7/{all_scores,mixed_*,v7h_homology_scores,td2_*} |
| Fig 3 starts | nb 03 `fig3_start` | start_param_sensitivity/release/*, start_refinement/*, benchmark/summary.tsv, final_v7/all_scores |
| Fig 4 short ORFs | nb 04 `fig4_short_orfs` | final_v7/{by_length.v7,lnc_by_era.human.v7}, release_pin/*, rules/lnc_riboseq_orfs, short_orf_summary |
| Fig 5 failures | nb 05 `fig5_failures` | long_orf_failures/v7/*, homology_vs_abinitio, homology_v7h_vs_abinitio |
| Fig 6 GM12878 | nb 06 `fig6_gm12878` | gm12878/final_v7/06_*, final_v7/gm12878_v7_*; checks figures/fig6_gm12878/data/pred_class_by_length.tsv |
| Fig 7 parameters | nb 08 `fig7_parameters` | gm12878/release_pin/*, gm12878/final_v7/08_minaa_sweep_v7, gm12878_homology_series_v7h |
| Fig S1 Kozak logos | nb 03 `figS1_kozak_logos` | start_refinement/kozak_models.tsv |
| Fig S2 start grid | nb 03 `figS2_start_grid` | start_param_sensitivity/release/* |
| Fig S3 history | nb S1 `figS3_history` | benchmark/history/history_scores, final_v7/negatives_scores |
| Fig S4 ranks | nb 02 `figS4_rankings` | as Fig 2 |
| Fig S5 depth | nb 02 `figS_simtrinity_depth` -> `figS5_simtrinity_depth.pdf` | simtrinity/depth/*, start_classifier/adoption/human_depth_scores |
| Fig S6 depth sens | nb 02 `figS_simtrinity_depth_sens` -> `figS6_...` | as S5 |
| Fig S7 athal | nb 02 `figS_simtrinity_athal` -> `figS7_simtrinity_athal.pdf` | simtrinity_athal/*, start_classifier/adoption/athal_depth_scores |
| Table 1 | nb 02 (writes final_v7/table1.tsv); typed in results_benchmark.tex | summary.release.tsv, final_v7/* |
| S1 methods | `pipeline/derive/02_comparator_table.py > data/benchmark/comparator_table.tex` | comparator_features, table1, partial_types |
| S2 Perl vs Python | `09_supplementary_tables.py` | regression/release_d3892b5_41fc5ee, baseline_v6.0.0 |
| S3 partial types | hand-entered from nb 02 | final_v7/partial_types.tsv |
| S4 mixed types | hand-entered from nb 02 | final_v7/mixed_medians.tsv |
| S5 per input (TSV) | `09_supplementary_tables.py` | summary.release, final_v7/{all_scores,mixed_scores} |
| S6 start choice | `09_supplementary_tables.py` | summary.tsv, final_v7/all_scores, start_refinement/seed_start_changes, later_atg_kozak |
| S7 GM12878 homology | `pipeline/derive/08_homology_series_v7h_table.py` | gm12878/gm12878_homology_series_v7h.tsv (E06) |
| S8 L* by GC | `09_supplementary_tables.py` | long_orf_threshold/final_thresholds.tsv |
| S9 truth sets | `09_supplementary_tables.py` | benchmark_v2/annotation_downloads.tsv, truth_logs_v3/*.log |
| S10/S11 releases | hand-entered (nb S1) | history/VERSIONS.md, history_scores |

Figure-script reads of `$P/work`: only `figures/fig6_gm12878/derive_class_by_length.py` (default
`$P/work/gm12878_td/eval/status`, the pre-release "final" GM12878 runs). Its output
`figures/fig6_gm12878/data/pred_class_by_length.tsv` is committed and notebook 06 uses it only as a check of the
"earlier" rows; the release Fig. 6 does not need `$P/work`. Derive scripts (`pipeline/derive/*`) do read `$P/work`
(table below).

---

## PART 2. Things in `$P/work` that are inputs or expensive to regenerate

Sizes are `du -sh` (filesystem is compressed; apparent sizes are larger).

| path | size | what | needed by | copy elsewhere? | recommendation |
|---|---|---|---|---|---|
| `bench_v3/runs/*/*/pytd_v7_homology/{blastp.outfmt6,pfam.domtblout}` | 263 MB (168 files) | DIAMOND vs `sprot_minus_<taxid>` + Pfam hits on release LongOrfs peptides | E25, E33 (`start_classifier/11_run_pytd_v7h_homology_v1.sh`), E40 retonly (`release_pin/run_bench_v1.sh` hardcodes this path) | no; `pytd_dev_v2_homology` has identical peptides (84/84 `cmp`) and identical hits (human/complete checked) — also in work | **move** to `$P/data/bench_v3/homology_hits/<org>/<variant>/` (67 h to redo; the `sprot_minus_<taxid>` builder is not saved) |
| `bench_v3/runs/*/*/<comparator>/preds.tsv`, `runtime_sec.txt`, `*version.txt` (all labels) | ~1.35 GB apparent (5,902 files) | standardized predictions | rescoring any benchmark table (`score_all.py`, `score_v7.sh`, `score_mixed.py`, `02_partial_types.py`, `10_td2_tables.sh`, `01_perl_pytd_identity.py`, `02_first_atg_baseline.py`) | scores only (committed tables) | **keep at least codan, gmst, cpat, cpc2, borf, rnasamba, perl_td_v600, td2, prodigal, pytd(`longest_orf.preds.tsv`)**: CodAn = 1,629 h; GMST needs licence. Move a preds-only tree to `$P/data/bench_v3/preds/` |
| `bench_v3/runs/` full (all run dirs) | 161 GB | full outputs | as above | no | delete after the preds/hits extract (release TD runs regenerate in ~37 h total) |
| `bench_v3/launch_v7.sh`, `launch_mixed.sh`, `launch_internal.sh`, `history/launch.sh` | KB | launchers of E25, E28, superseded internal, E23 | reproduction | **not in repo** (others copied to `benchmark_v2/launch_v3/`) | **copy into `benchmark_v2/launch_v3/`** |
| `bench_v3/scores/` | 1.3 MB | `summary.interim_20261022.tsv` (source of `homology_vs_abinitio.dev_v2.tsv`, Fig.-era E07 checks), `longest_orf_atg.tsv`, `summary.partial.tsv`, err logs | E07 provenance | interim_20261025 and summary.tsv committed; 20261022 and longest_orf_atg not | copy the two uncommitted TSVs to `data/benchmark/` |
| `bench_v3/history/` | 7.2 GB | 80 historical-release runs | E23 | `history_scores.tsv` etc. committed | keep `*/*/*/preds{,.clipped}.tsv` + runtime (28 h to redo); delete rest |
| `history/sourceforge/` (archives 133 MB) + `history/sf_runs/` (1.7 MB) | 309 MB + 1.7 MB | SourceForge tarballs r2012-08-07 … r20140704 (`MD5SUMS.local`), r2012 script copy used by `hist_run.sh` | E23 (r2012 row) | no; SourceForge page now only says "moved to GitHub" | **move to `$P/data/history_releases/`** |
| `history/compat/` | 112 KB | DB_File + libdb symlinks into `$P/env/mamba/pkgs` | E23 (releases <= v5.2) | targets in env (kept) | trivially recreated; note in REPRODUCE |
| `lnc/{human_lnc.pick.gtf,human_lnc.pick.fa,mouse_lnc.pick.fa}` | 7 + 19 + 16 MB | lncRNA negative picks | building wLnc inputs; `lnc_calls_by_era.py` via `04_release_short_orfs.sh` (called by `10_td2_tables.sh`) | wLnc FASTAs in `$P/data/bench_v3/inputs` only | **move** to `$P/data/lnc/` (rest of `work/lnc`, 1.8 GB incl. mouse genome, is regenerable) |
| `gm12878_td/eval/status/` | 11 MB | pre-release GM12878 per-prediction status | `derive_class_by_length.py` default (Fig. 6 check only) | output committed | safe to delete (or keep 11 MB) |
| `bench_v3/first_atg_v7/` | 20 MB | first-ATG baseline preds | E07 addendum | `final_v7/first_atg_baseline.tsv` committed | safe (2 min to rebuild from release LongOrfs) |
| `bench_v3/rerun_old16/`, `rerun_old/` | 16 GB, 1.3 GB | superseded dev final-default runs | none | — | safe to delete |
| `bench_v3/rules/` | 47 GB | dev rule runs incl. `cliff_c200`, `decoy_q0.05_longgc` | E23 v7 rows, `04_release_short_orfs.sh`, Fig. 4 historical comparator (E19/E20/E22) | committed tables | keep preds of `cliff_c200` (unreleased code, cannot be rerun from a release) |
| `td2_private/` | 3.5 GB | private TD2 runs 2026-10-23, not reused | none | — | safe to delete |
| `bench_v2/` | 4.4 GB | superseded RefSeq-gbff benchmark | none (notes/04 early numbers, not quoted) | — | safe to delete |
| `regression/` | 1.3 GB | Perl/PyTD sample-set runs | E08 | summaries committed | safe (minutes to redo) |
| `release_pin/bench` (of 80 GB `release_pin/`) | 37 GB | E40 fresh d3892b5 runs | E40 tables | tables committed | out of scope here; keep the `retonly` dependency in mind (needs the hits above) |
| `simtrinity/`, `simtrinity_athal/` | (not measured) | E30–E32 assemblies + runs | `10_td2_tables.sh` reads them | — | see part A (E30–E32); Trinity assemblies are expensive |

**Comparator model files/containers:** none in `$P/work` — CodAn models, CPAT models, RNAsamba weights and ANGEL are
in `$P/tools` (keep it); GMST in `~/software/genemarks_t` + `~/.gm_key`.

**`$P/work/trees` (1.4 GB).** All clean (no modified/untracked files). `branch -r --contains HEAD` on the fetched refs:

| tree | HEAD | on remote? | role |
|---|---|---|---|
| PY_48c870c | 48c870c | origin/main | E07 homology re-predict |
| PY_base | 43f4293 | origin/main | E08 baseline |
| PY_dev | 44e417a | origin/main | dev (Kozak bits) |
| PY_docs, PY_rel_d3892b5 | d3892b5 | origin/main | release pin |
| PY_fix | 27bfb5e | origin/main | `pytd`, longest_orf baseline |
| PY_hstart | a5b24c3 | origin/main | E33 v7h |
| PY_taper, PY_v7 | 3a793ff | origin/main | E25 release runs |
| **PY_frames** | **83a4806** (local branch `exp/contrast-frames`) | **NO** (remote has only main, devel, tags) | E26 all-frame contrast (not adopted) |
| **PY_utr** | **c9d7875** (local branch `exp/utr-prior`) | **NO** | E29 UTR prior (not adopted) |
| TD_base | 66d4712 (v6.0.0) | origin/master, devel | perl_td_v600 |
| TD_fix, TD_v7, TD_rel_41fc5ee | 41fc5ee | origin/devel, master | release |
| TD_v571 | 5791fde (tag TransDecoder-v5.7.1) | yes | perl_td_v571 |
| TD_hist_v2.0 / v3.0.1 / v4.1.0 / v5.0.0 / TransDecoder-v5.5.0 | c4ea24d / 9b52e69 / 8f0e552 / 3b4f9b1 / 3c77d07 (tags) | yes | E23 |

Other commits cited by my records are on remotes: Perl de3f57e, 351634a, 694e204, 8cf8d5f, ab38ca6; PyTD efc1450,
f5b12e7, dc7b514. **Action: push `exp/contrast-frames` and `exp/utr-prior` to PyTransDecoder (or `git bundle` them
into `$P/data`) before deleting `work/trees`.** Remote-tag notes: PyTD v0.2.0 = 275ff1b (same `pytransdecoder/`
tree as d3892b5); Perl `TransDecoder-v7.0.0` = 47cedc9 (= 41fc5ee + version strings).

---

## PART 3. Gaps

1. **Launchers not in repo**: `work/bench_v3/launch_v7.sh` (E25 release benchmark), `launch_mixed.sh` (E28),
   `history/launch.sh` (E23). Copy them before deleting work.
2. **Hardcoded `$P/work` paths** in committed scripts: `run_tool.sh` defaults (trees), `hist_run.sh` (sf_runs, compat,
   TD_hist trees), `data/benchmark/final_v7/score_v7.sh`, `release_pin/run_bench_v1.sh` (PY_rel tree, homology hits),
   `start_classifier/11_run_pytd_v7h_homology_v1.sh` (PY_hstart), `10_td2_tables.sh`, `10b_td2_mechanisms.py`,
   `01_perl_pytd_identity.py`, `02_first_atg_baseline.py`, `04_release_short_orfs.sh` (lnc picks, rules),
   `launch_td2_v1.sh` and `td2_jobs_v1.tsv`/`history/jobs.txt` (absolute paths). REPRODUCE must say to set env vars or
   edit `P=`.
3. **No builder for `sprot_minus_<taxid>.dmnd`** (E07 gap); only the clade builder
   `genome_annotation/scripts/make_sprot_clade_excluded.py` exists. Hence keep the hit tables.
4. **History scoring command** producing `history_scores.tsv` (columns of `score_negatives.py`, with r2012 using
   `preds.clipped.tsv`) is not recorded.
5. **`summary.interim_20261025.tsv`** and `summary.release.tsv`: the exact `score_all.py` invocation/filter and
   the row selection are not recorded (E40 says release + comparator rows of interim, unchanged).
6. **E08**: `parity_check.py` arguments and the r1–r5 invocations not recorded; runtime ratios cannot be
   reproduced with PyTD >= 0.3.0 (faster).
7. **E27 tallies**: `by_type.py`, `start_codons.py` lived in a session scratchpad, not committed. E16 sources
   (paper copies) not kept.
8. **Model provenance**: download URLs of CPAT prebuilt models and RNAsamba weights are not recorded (files only in
   `$P/tools/models`); GMST download/version only as `gmst.pl --version`.
9. **Supplement tables S3, S4, S10, S11** and Table 1 are typed into LaTeX from notebook output (no generator).
10. **Figure copy step** `figures/ -> manuscript/figures/` (with S5–S7 renaming) is manual; `figures/fig*/make_fig*.py`
    are stale and should be labelled as such.
11. **Outputs only in `$P/work`, no committed table**: `scores/summary.interim_20261022.tsv` (source of
    `homology_vs_abinitio.dev_v2.tsv`), `scores/longest_orf_atg.tsv` (100-aa version, superseded by
    `first_atg_baseline.tsv`), all per-run `preds.tsv` (incl. CodAn 1,629 h and licensed GMST), homology hit tables,
    SourceForge archives, `rerun_old16/summary.pre16.tsv`.
12. **Unscripted moves**: the move to `rerun_old16/` before launch16 (E07); `score_all.py` will pick up any extra run
    dirs (e.g. mixed) — filter.
13. **Annotation URLs that float** (MANE `current/`, TGD `latest/`): reproduction depends on `$P/data/annot` (kept);
    truth md5s in E07/E28 allow verification.
14. **TD2**: `PYTHONHASHSEED` unset in published runs (tiny unstranded differences expected); `10_td2_tables.sh` is
    append-only and refuses to rerun on tables that already contain td2 rows.
15. **E07 ANGEL**: 3 completed + 3 partial athal runs; listed in S5 (3 rows) but not reproducible in reasonable time.
