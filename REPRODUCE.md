# Reproducing the TransDecoder paper from scratch

This guide explains how to regenerate every analysis in the paper, using the latest PyTransDecoder release
(**0.4.0**, PyPI `pytransdecoder==0.4.0`, tag `v0.4.0` = a285464) and the Perl TransDecoder release pin
(**41fc5ee** on `devel`; tag `TransDecoder-v7.0.0` differs only in version strings).

This page gives the setup, the order of runs, and an index of all 42 experiments. The step-by-step commands are in
four detail files:

| File | Experiments |
|---|---|
| [`reproduce/A_gm12878_and_simulated_assemblies.md`](reproduce/A_gm12878_and_simulated_assemblies.md) | E01–E03, E05, E06, E11 (GM12878); E30–E32 (simulated Trinity assemblies) |
| [`reproduce/B_benchmark_comparators_figures.md`](reproduce/B_benchmark_comparators_figures.md) | E07, E08, E15, E16, E23, E27, E28, E37, E39; **figure and supplementary-table map** |
| [`reproduce/C_orf_rules_start_refinement_presets.md`](reproduce/C_orf_rules_start_refinement_presets.md) | E04, E09, E10, E12, E13, E17–E22, E24–E26 |
| [`reproduce/D_release_pin_start_decision_speed.md`](reproduce/D_release_pin_start_decision_speed.md) | E40 (the release-pin backbone), E41, E29, E33–E36, E42, E14, E38 |

The evidence records (`evidence/E*.md`) remain the primary record of what was done and what was found. This guide
and the detail files are how to redo it.

Paths: `$P` = project root (`/broad/hptmp/bhaas/TRANSDECODER_PAPER` on the analysis machine), `$T` = this repo.
Every script hardcodes `P=`. To run elsewhere, change that line.

## 1. Why the latest PyTD reproduces the paper

- Every TransDecoder/PyTransDecoder number in the manuscript comes from one code version. That is the release
  pin PyTD **d3892b5** and Perl **41fc5ee** (E40), plus the E25 `pytd_v7` runs at 3a793ff, whose ab initio output
  is byte-identical to d3892b5.
- PyTD v0.2.0 has the same package tree as d3892b5; v0.2.1 differs only in the version string. v0.3.0 and v0.4.0
  were checked byte-identical on all output files, at any `--cpu` (E42: 167/167 and 226/226 comparisons). So
  **0.4.0 regenerates every prediction file**.
- **What 0.4.0 does not reproduce: run times.** 0.3.0+ is faster, so `runtime_sec.txt`, Fig. 2f, the
  run-time supplementary table and the "2.7× faster than Perl" statements will change. See
  `notes/09_handoff_2026-11-02.md` for the list of affected manuscript lines.
- **Homology-mode caveat.** At 3a793ff, `pytd_v7_homology` used hits for retention only. From a5b24c3
  (in d3892b5 and 0.4.0) the default also uses them for start refinement. The paper's homology rows are the
  latter (`pytd_v7h_homology`, E33/E40). Add `--no-homology-start-assist` to reproduce the retention-only rows
  (`retonly`).
- Development and not-adopted experiments need the commit they were run with. The index (section 6) says which.
  All of those commits are in the bundles of section 2. All are also on GitHub except the branches
  `exp/contrast-frames` (83a4806, E26), `exp/utr-prior` (c9d7875, E29) and `exp/short-orf-taper` (8dbe79b, E19;
  its commits are on `main`).

## 2. What exists and where (after `$P/work` was deleted)

| Location | Contents | Status |
|---|---|---|
| `$T` (this repo, GitHub TransDecoder/TransDecoderPaper) | all scripts, committed result tables (`data/`, `*/results/`), notebooks, evidence records, environment specs | in git |
| `$T/manuscript` | LaTeX (submodule TransDecoderPaper_LaTeX, Overleaf-synced) | in git |
| `$P/data` | inputs. Benchmark inputs and truth sets (`bench_v3/inputs`), annotations (`annot`; MANE/TGD URLs point to "latest", so this copy is authoritative), GM12878 assemblies, footprint BAM, MS search and StringTie quantification (`gm12878`), simulated transcript sets and assemblies (`simtrinity*`), homology databases (Swiss-Prot 2026_03, `sprot_minus_<taxid>.dmnd`, Pfam-A 38.2: `homology`), lncRNA sets, RefSeq, taxonomy | kept |
| `$P/data/work_keep/` | the parts of `$P/work` that are inputs, expensive to recompute, or the only record (section 5), at their original paths relative to `work/` | kept (copied 2026-11-03) |
| `$P/data/code_bundles/` | `git bundle --all` of `$P/repos/PyTransDecoder` and `$P/repos/TransDecoder`: every local branch and tag, including the unpushed `exp/*` branches | kept |
| `$P/env` | conda environments. Their exported specs are in `$T/environment/conda_envs/*.yml` (+ `td.pip-freeze.txt`, `td2.pip-freeze.txt`) | kept |
| `$P/tools` | CodAn (git 43da738) + models, CPAT models, RNAsamba weights, ANGEL (git b805554), ORFfinder. The CPAT/RNAsamba download URLs were not recorded, so keep this directory | kept |
| `$P/repos` | clones of PyTransDecoder, TransDecoder and others. The `work/trees/*` worktrees stored their git data here | kept |
| `~/software/genemarks_t`, `~/.gm_key` | GeneMarkS-T 5.1. **Licensed: neither the software nor the key may enter any repo** | outside project |
| `$P/work` | all run outputs | **deleted after this guide was committed** |

Raw reads were already deleted on 2026-10-29 (`notes/DISK_CLEANUP_2026-10-29.md`). These are the GM12878 short-read
FASTQs, STAR index and BAM, Ribo-seq FASTQs/BAMs, simulated FASTQs and STAR BAMs, and the ERR4140056 reads. Their
downstream products (assemblies, footprint BAM, labels) are in `$P/data`. Rebuilding the assemblies from re-downloaded
reads is described in detail A, but it has **not** been checked to give identical files. The paper's numbers depend
on the kept assemblies.

## 3. Setup

```bash
P=/broad/hptmp/bhaas/TRANSDECODER_PAPER; T=$P/TransDecoderPaper
# Code: git checkouts (the scripts run `git rev-parse` in the code tree under `set -e`, so a bare pip install is not enough)
mkdir -p $P/trees
git clone https://github.com/TransDecoder/PyTransDecoder $P/trees/PyTransDecoder     # or: git clone $P/data/code_bundles/PyTransDecoder.all-refs.2026-11-03.bundle
git -C $P/trees/PyTransDecoder worktree add --detach $P/trees/PY_v0.4.0 v0.4.0
git clone https://github.com/TransDecoder/TransDecoder $P/trees/TransDecoder
git -C $P/trees/TransDecoder worktree add --detach $P/trees/TD_41fc5ee 41fc5ee
export PYTD_REL=$P/trees/PY_v0.4.0 PYTD_V7=$P/trees/PY_v0.4.0 PYTD_TREE=$P/trees/PY_v0.4.0 \
       TD_V7=$P/trees/TD_41fc5ee PYTHON=$P/env/td/bin/python PERL_HASH_SEED=0
# Environments: recreate from the exported specs if $P/env is gone
micromamba create -p $P/env/td -f $T/environment/conda_envs/td.yml     # likewise codan, cpat, cpc2, borf, rnasamba, angel, trinity215, rnaseq_asm, ...
```

- **Env-var overrides the scripts honour:**
  - `PYTD_REL` (release_pin `run_bench_v1.sh`, `run_gm_v1.sh`, `run_preset_v1.sh`);
  - `PYTD_V7` / `TD_V7` (`benchmark_v2/scripts/run_tool.sh`, `genome_annotation/scripts/gm12878_final_v7.sh`);
  - `PYTD_TREE` (`pipeline/derive/05_long_orf_failures_v7.py`).
  - `run_tool.sh` also defaults `TD_FIX`, `TD_V600`, `TD_V571`, `PYTD`, `PYTD_DEV` to `$P/work/trees/...`. Set them,
    or create the trees at the same paths (detail B, section 0, lists each tree's commit).
- **Hardcoded tree paths** (edit, or symlink `$P/work/trees/PY_rel_d3892b5` and `$P/work/trees/PY_v7` to
  `$P/trees/PY_v0.4.0`):
  - `analyses/start_param_sensitivity/run_grid_v2.py`
  - `kozak_score_distribution_v2.py`
  - `later_atg_kozak_v1.py`
  - `pipeline/derive/08_minaa_sweep_final_v7.sh`
  - `pipeline/derive/03_start_refinement.py`
  - `pipeline/release_pin/launch_v1.sh` stage A (regression trees)
- **Using the pip package instead of a checkout:** replace `PYTHONPATH=$PYTD $PY -m pytransdecoder.cli` with the
  `pytransdecoder` entry point, and `git rev-parse` with `pytransdecoder --version`. Detail D, section 0, lists every
  place.
- **Offline re-derivation scripts** (E29, E33–E36) import private `TransDecoderPredict` methods that were
  refactored in 0.3.0. Run them with `pip install pytransdecoder==0.2.1` (same code as d3892b5).
- **`$P/work` paths inside scripts.** Most scripts write their outputs under `$P/work/...`, so a rerun recreates
  that tree. Scripts that *read* kept inputs from `$P/work` (homology hits, preset held-out sets, `cliff_c200`
  predictions, lncRNA picks, simulated-assembly labels) must be pointed at `$P/data/work_keep/...`. Or restore those
  files first: `rsync -a $P/data/work_keep/ $P/work/`. That restores them to the paths the scripts expect.
  The input symlinks `$P/data/{missed_cds,alpha_by_length}/gm12878_inputs/*` now point to the kept GM12878 Trinity
  runs `$P/data/work_keep/gm12878_td/trinity/{final_v7_abinitio,final_v7h_homology}`. These were repointed on
  2026-11-03.
- Run everything locally at normal priority (no UGER, no `nice`). A few E33–E36 scripts call `nice`; drop it.

## 4. Order of a full rebuild

| Step | What | How | Cost (recorded) |
|---|---|---|---|
| 0 | Inputs | `$P/data` as kept. From scratch: `benchmark_v2/{annotation_downloads,truth_sources,refseq_sources}.tsv`, `benchmark_v2/README.md`, E07 (benchmark truth sets, negatives), E11 (GM12878 reads, Trinity GG), E30/E32 (simulation) | days (assemblies) |
| 1 | Homology hits | reuse `work_keep` hits (DIAMOND 2.2.8 `blastp -k 1 --evalue 1e-5` vs `sprot_minus_<taxid>`; HMMER 3.4 `hmmsearch --cut_ga` vs Pfam-A 38.2), or redo via `run_tool.sh pytd_v7_homology` / `gm12878_final_v7.sh` homology mode. The `sprot_minus_<taxid>` builder was not saved; detail D reconstructs it | ~67 h CPU for the benchmark |
| 2 | Release-pin backbone (E40) | restore `work_keep`, then run `pipeline/release_pin/launch_v1.sh`. Its stages:<br>**A**: Perl/PyTD regression; `rel` on all benchmark inputs; `precise`/`sensitive` on the 20 negative-containing inputs.<br>**B**: GM12878 `rel` plus 17 configurations (α×w start grid, `norefine`, `precise`, `sensitive`, `m100_noq`, `m50_noq`) on 3 assemblies; Trinity `retonly`.<br>**C**: short-ORF (`m*_noq`), q and long-ORF-off sweeps on the 20 inputs; benchmark `retonly`. Its input list comes from the `work/bench_v3/runs/*/*/pytd_v7h_homology` directories of step 3.<br>**D**: presets. Reference builds, then tests from `work/bench_v3/presets/jobs*.txt` (kept).<br>**E**: benchmark start grid (`analyses/start_param_sensitivity/run_grid_v2.py`). Then score with `10_byte_identity_v1.py`, `11_score_sweeps_v1.sh`, `12_gm_tables_v1.py`, `13_gm_support_v1.sh`, `14_rule_support_v1.py`, `15_score_presets_v2.sh` | ~32 h wall at 4 processes with d3892b5; less with 0.4.0 `--cpu` |
| 3 | E25 release labels on the benchmark | `pipeline/from_work/bench_v3/launch_v7.sh` (→ `run_tool.sh pytd_v7*`, `perl_td_v7*`, `pytd_v7h_homology`). The derive scripts read these paths. They can be replaced by symlinks `work/bench_v3/runs/<o>/<i>/pytd_v7 → work/release_pin/bench/<o>/<i>/rel` (byte-identical) | hours |
| 4 | Comparators | `run_tool.sh <tool>` via `benchmark_v2/launch_v3/*` and `pipeline/from_work/bench_v3/launch_mixed.sh`. The kept `preds.tsv` of every comparator run (CodAn alone is 1,629 CPU-h; GeneMarkS-T is licensed) let you skip this and rescore | weeks of CPU for CodAn |
| 5 | Scoring the benchmark | `benchmark_v2/scripts/score_all.py`, `score_v7.sh`, `pipeline/benchmark/score_mixed.py`, `pipeline/derive/02_*` | minutes |
| 6 | GM12878 application | `genome_annotation/scripts/*` + `pipeline/from_work/gm12878_td/launch_v7.sh`; evidence scoring `pipeline/derive/06_*`, `07_*`, `08_*` (detail A; the Ribo-seq/MS scoring of the precise/sensitive rows is reconstructed there) | ~1 day |
| 7 | Simulated assemblies | `pipeline/simtrinity/` (human) and `pipeline/simtrinity_athal/` + `pipeline/from_work/simtrinity/*` | days (Trinity per depth level) |
| 8 | Derived analyses | `pipeline/derive/*`, `analyses/*`, `pipeline/{start_classifier,missed_cds,kozak_floor,alpha_by_length}` | hours |
| 9 | Figures | `scripts/knit_notebooks.sh` renders `notebooks/*.Rmd` from committed tables, then copy the figures to `manuscript/figures/`. S5–S7 are renamed (map in detail B). `figures/fig*/make_fig*.py` are superseded | minutes |
| 10 | Supplementary tables | E37 (S1–S9 from committed data). Table 1, S3, S4, S10, S11 were typed into the LaTeX by hand from committed tables | — |

## 5. What was kept from `work/` (`$P/data/work_keep/`)

Copied on 2026-11-03 by `pipeline/maintenance/keep_from_work_v1.sh` (28,477 files, 5.7 GB; log: `$P/data/work_keep.copy.log`), plus the two GM12878 Trinity runs `gm12878_td/trinity/{final_v7_abinitio,final_v7h_homology}` (0.7 GB) that the E34/E36 inputs link to; 6.5 GB in total. Copies only, original
relative paths.

| Kept | Why |
|---|---|
| `bench_v3/runs/*/*/pytd_v7_homology/{blastp.outfmt6,pfam.domtblout,…/longest_orfs.pep}` | benchmark homology hits (~67 CPU-h). E25, E33, E40 `retonly` read them, and `run_bench_v1.sh` checks the peptides with `cmp` |
| `gm12878_td/**/{blastp*.outfmt6,pfam.domtblout}` | GM12878 Trinity/LRAA/StringTie hits and the 15 clade-excluded series (E03, E06, E40 `retonly`) |
| `start_classifier/**/searches/*`, `start_classifier/homology/*` | per-depth homology searches of the simulated assemblies (E33–E35) |
| every `preds*.tsv`, `longest_orf.preds.tsv`, `runtime_sec.txt`, `*version*.txt` under `bench_v3/{runs,rules,history,rerun_old16,frames,stringency,minlen,presets,first_atg_v7}`, `simtrinity*` | standardized predictions of every benchmark run. These allow rescoring without rerunning: CodAn (1,629 CPU-h), GeneMarkS-T (licensed), unreleased-code rules such as `cliff_c200` (Fig. 4 comparator, needs 652d5e4), historical releases (E23) |
| `bench_v3/presets/{ref,tests}` | E21 held-out reference halves and test sets. Made by an uncommitted inline script, so they cannot be regenerated exactly |
| `history/{sourceforge,sf_runs}` | SourceForge release archives r2012…r20140704, which are no longer downloadable (E23) |
| `bench_v3/scores/*` | interim score summaries (two are now also in `data/benchmark/`) |
| `lnc/*.pick.*`, `lnc/*.py` | lncRNA negative picks (GTF needed by the era table) |
| `thresh_big/raw` | E10 low-GC simulation chunks |
| `utr_prior/annot`, `utr_prior/gm12878/*{truth,pairs}*` | E29 truth and UTR tables |
| `simtrinity/{realtrinity_label,filtered_preds}`, `simtrinity*/depth/*/scores` | real-Trinity labels and the full-depth scored predictions, which `09_summarize.py` and `13_depth_collect.py` read |
| `analysis/{pwm_audit,kozak,fork_findings.md}`, `start_grid` | E09 development notes and tables that existed only in work |
| `release_pin/**` small per-run files (preds, scores, eval, Ribo-seq, runtime, version, logs, presets) | re-derive E40 tables without rerunning |
| every `*.sh`, `*.py`, `*.R`, `*.awk`, `*.md`, `jobs*.txt` under `work/` | launchers and notes. The ones that were never committed are now also in `pipeline/from_work/` |

Deleted without a copy (all regenerable from the steps above, summarised in committed tables):
- full run directories of `bench_v3/runs` (161 GB), `rules`, `rerun_old16`, `frames`, `stringency`, `history`
  runs, `release_pin` (80 GB), `gm12878_td` runs, `start_param_grid`, `start_classifier` runs, `simtrinity*` runs;
- the superseded `bench_v2`, `td2_private` and `regression`;
- the code worktrees `work/trees` (their commits are in `$P/repos` and the bundles).

## 6. Experiment index

"0.4.0" = rerun with the latest PyTD (and Perl 41fc5ee where Perl is involved) reproduces the paper's values (except
run times). "commit X" = needs that historical code. "lit." = literature only.

| Rec. | Topic | In paper? | Code to reproduce | Detail |
|---|---|---|---|---|
| E01 | GM12878 Ribo-seq support | yes | 0.4.0 | A |
| E02 | GM12878 MS peptide support | yes | 0.4.0 | A |
| E03 | GM12878 three assemblies vs MANE/GENCODE | yes | 0.4.0 | A |
| E04 | Minimum ORF length sweep | yes | 0.4.0 (guard comparator: historical preds kept) | C |
| E05 | Expressed-gene sensitivity | yes | 0.4.0 (expression command reconstructed) | A |
| E06 | Homology restriction series | yes (Fig. 7c,d, Table S7) | 0.4.0 (release series `v7h`) | A |
| E07 | Cross-species benchmark | yes (Table 1, Fig. 2, S5) | 0.4.0 + comparators | B |
| E08 | Perl/Python concordance, run time | yes (Table S2) | 0.4.0 + Perl 41fc5ee; **run times change** | B |
| E09 | Start refinement and Kozak | yes (S1, S6) | 0.4.0 for release numbers; development history needs dev commits | C |
| E10 | GC-calibrated long-ORF rule | yes | code-independent simulation (part of the table merged by hand) | C |
| E11 | GM12878 Trinity genome-guided assembly | yes (Methods) | Trinity 2.15.2; reads deleted, assembly kept | A |
| E12 | True CDSs failing the coding test | yes (Fig. 5) | 0.4.0 | C |
| E13 | Start-parameter sensitivity | yes (Fig. 3, S2, Fig. 7) | 0.4.0 | C |
| E14 | References | lit. | — | D |
| E15 | Figures | yes | notebooks from committed tables | B |
| E16 | Comparator benchmarks in the literature | lit. | — | B |
| E17 | Critiques old vs new; bootstrap CIs | no | development code | C |
| E18 | Stringency sweep / `--precise` | not adopted | development code | C |
| E19 | Short-ORF threshold shape | not adopted | development code; `cliff_c200` comparator = commit 652d5e4 (preds kept) | C |
| E20 | Target-decoy FDR rule | adopted | 0.4.0; strata prototype = development code | C |
| E21 | Presets for small inputs | yes | 0.4.0 + kept held-out sets | C |
| E22 | q-value modes | yes | 0.4.0 | C |
| E23 | Version history | supplement (S3, S10, S11) | historical tags (in bundles) + kept SourceForge archives | B |
| E24 | Release decoy rule + presets, Perl = PyTD | yes (Methods) | 0.4.0 + Perl 41fc5ee | C |
| E25 | Final release reruns | yes | 0.4.0 | C |
| E26 | Frame contrast | yes (Methods) | commit 83a4806 (`exp/contrast-frames`, bundle; patch in `data/benchmark/frames/`) | C |
| E27 | How comparators work | yes (Background, S1) | lit. | B |
| E28 | Mixed partial inputs | yes (Table 1, S3, S4) | 0.4.0 + comparators | B |
| E29 | 5'-open starts; UTR prior | anatomy yes; prior not adopted | 0.4.0 (offline: 0.2.1); prior = commit c9d7875 (`exp/utr-prior`, bundle) | D |
| E30 | Simulated GM12878 Trinity benchmark | yes | 0.4.0 + comparators | A |
| E31 | Read-depth ablation | supplement (S5, S6) | 0.4.0 + comparators | A |
| E32 | Simulated Arabidopsis assembly + depth | supplement (S7) | 0.4.0 + comparators | A |
| E33 | Keep-open vs trim; homology start assist | assist yes; classifier no | 0.4.0 (offline parts: 0.2.1) | D |
| E34 | Missed-CDS anatomy | no (diagnostic) | 0.4.0 / 0.2.1 | D |
| E35 | Kozak floor | not adopted | 0.2.1 internals | D |
| E36 | Alpha by length | deferred | 0.2.1 internals | D |
| E37 | Supplementary Tables S1–S9 | supplement | from committed data | B |
| E38 | Background literature | lit. | — | D |
| E39 | TD2 comparator | yes | TD2 1.1.0 + PSAURON 1.1.3 (`benchmark_v2/TD2_REPRODUCE.md`) | B |
| E40 | Release-pin reruns (one code version) | yes (backbone) | 0.4.0 | D |
| E41 | Later in-frame ATGs by start context | yes (Table S6b) | 0.4.0 | D |
| E42 | Parallel and faster PyTD | engineering | 0.4.0 | D |

## 7. Known gaps

These steps have no committed command, or rely on results that were merged by hand. Details are in Part 3 of each
detail file.
- **GM12878:**
  - the release Ribo-seq/MS scoring (`final_v7_scoring`) has no script; a reconstruction in detail A matches the
    committed tables;
  - `13_gm_support_v1.sh` makes only the default and `m50_noq` rows;
  - the ERR4140056 download is not scripted;
  - the E05 expression command is a reconstruction.
- **Simulated assemblies:** the full-depth E30 labelling, scoring and summary steps are not scripted.
- **E10:** the 50–80% GC simulation runs have no command, and the final table was merged by hand.
- **E20 and E21:** the commands for the E20 prototype's shuffled-input runs were not recorded. The E21 held-out set
  builder was an inline script (its outputs are kept).
- **E26:** the scoring command for `negatives_scores.tsv` / `negatives_curve.tsv` was not recorded.
- **E40:**
  - no script makes `longest_orf_baseline_check.tsv`;
  - `10_byte_identity_v1.py` and the E42 parity summaries compare against runs that are now deleted, so the
    committed tables stand as the record.
- **E23:** the scoring command for `history_scores.tsv`, and the commands that built
  `summary.interim_20261025.tsv` and `summary.release.tsv`, were not recorded.
- **Homology databases:** the builder for `sprot_minus_<taxid>.dmnd` was not saved (reconstruction in detail D).
  The UniRef50 database (E33 step 5) was deleted on 2026-10-29.
- **Comparator models:** the CPAT model and RNAsamba weight download URLs were not recorded (the files are in
  `$P/tools`).
