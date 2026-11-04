# Reproduction detail C: ORF rules, start refinement, long-ORF rule, presets, release reruns (E04–E26)

> **Paths into `$P/work`.** This section was written on 2026-11-03, before `$P/work` was deleted. Files that were kept are now under `$P/data/work_keep/` at the same path relative to `work/`. For example, `$P/work/bench_v3/runs/human/complete/pytd_v7_homology/pfam.domtblout` is now `$P/data/work_keep/bench_v3/runs/human/complete/pytd_v7_homology/pfam.domtblout`. Everything else under `$P/work` has to be regenerated with the steps given here. The tables in Part 2 record what was kept and why; see `../REPRODUCE.md` for the full list.
(E04, E09, E10, E12, E13, E17, E18, E19, E20, E21, E22, E24, E25, E26; all read against E40)

`$P=/broad/hptmp/bhaas/TRANSDECODER_PAPER`, `$T=$P/TransDecoderPaper`. Python env `$P/env/td` (3.11.17).

## Common set-up (needed by every class-a recipe below)
The pinned release is PyTD d3892b5 (v0.2.0); v0.2.1 -> v0.3.0 -> v0.4.0 (a285464) are byte-identical in output
(E42). Every release-pin script runs `git rev-parse` in the PyTD tree, under `set -e`, and runs
`PYTHONPATH=$PYTD python -m pytransdecoder.cli`. So use a **git checkout**, not only `pip install`:
```bash
git clone https://github.com/TransDecoder/PyTransDecoder $P/trees/PyTransDecoder
git -C $P/trees/PyTransDecoder worktree add --detach $P/trees/PY_v0.4.0 v0.4.0     # a285464
git clone https://github.com/TransDecoder/TransDecoder $P/trees/TransDecoder
git -C $P/trees/TransDecoder worktree add --detach $P/trees/TD_41fc5ee 41fc5ee       # Perl v7.0.0 (origin/devel)
export PYTD_REL=$P/trees/PY_v0.4.0 PYTD_V7=$P/trees/PY_v0.4.0 PYTD_TREE=$P/trees/PY_v0.4.0 TD_V7=$P/trees/TD_41fc5ee
```
Env overrides honoured: `PYTD_REL` (run_bench_v1.sh, run_gm_v1.sh, run_preset_v1.sh), `PYTD_V7`/`TD_V7` (run_tool.sh,
gm12878_final_v7.sh), `PYTD_TREE` (05_long_orf_failures_v7.py). **Hardcoded, must be edited or symlinked**
(`ln -s $P/trees/PY_v0.4.0 $P/work/trees/PY_rel_d3892b5`; likewise `PY_v7`):
`analyses/start_param_sensitivity/run_grid_v2.py` (PYTD), `kozak_score_distribution_v2.py`, `later_atg_kozak_v1.py`
(sys.path), `pipeline/derive/08_minaa_sweep_final_v7.sh` (PYTHONPATH=work/trees/PY_v7),
`03_start_refinement.py` (documented PYTHONPATH=work/trees/PY_v7).
Development commits named below are all on `origin/main` of PyTransDecoder or `origin/devel` of TransDecoder
(checked with `git branch -a --contains`), **except 83a4806 (exp/contrast-frames, E26), which is local only**.

Two output roots are used by the derive scripts: `$P/work/release_pin/...` (E40 launcher) and the older
`$P/work/bench_v3/runs/<org>/<input>/pytd_v7*` and `$P/work/gm12878_td/<asm>/final_v7_*` (E25 labels; byte-identical
to release_pin, E40). Many derive scripts read the E25 paths, so a full rebuild needs both launchers (or symlinks
`work/bench_v3/runs/<o>/<i>/pytd_v7 -> work/release_pin/bench/<o>/<i>/rel`; note pytd_v7 has `longest_orf.preds.tsv`
and so does `rel`).

---------------------------------------------------------------------------------------------------------
# PART 1 — experiments

### E04 Minimum ORF length: benchmark sweep and GM12878 post-hoc filter  (status: in paper; class a, plus b for the guard comparator)
- **Feeds**: Results "short ORFs and the minimum length" (results_short_orfs.tex), Discussion; Fig. 4 (notebook 04),
  Fig. 7b and no-overlap-by-length (notebook 08); lncRNA era table (notebook 04).
- **Inputs**: `$P/data/bench_v3/inputs/<org>/{complete,fragments}.{wRand,wLnc}.{fa,truth.tsv}`; GM12878 assemblies in
  `$P/data/gm12878`; `$P/data/gencode_v28/v28_gene_ids.txt`; **`$P/work/lnc/human_lnc.pick.gtf`** (era table; work-only).
- **Software**: paper numbers = release (E40). `-m` sweep without q-test = LongOrfs `-m 100..50` +
  `predict --max-qvalue 1000000`. Development runs (PY_dev f5e6a39/7e9ffa2/48c870c, `-m`/`--short-orf-min-contrast`)
  are history (c). The "fixed short-ORF test" comparator `cliff_c200` (= -m 50 + guard 200) needs an unreleased tree (b, see E19).
- **Steps** (release):
  1. `bash $T/pipeline/release_pin/launch_v1.sh` stage A (`rel`) and C (`m100_noq..m50_noq` on 20 wRand/wLnc inputs) via
     `run_bench_v1.sh <org> <input> <label>` -> `$P/work/release_pin/bench/<org>/<input>/<label>/`.
  2. Stage B GM12878: `run_gm_v1.sh <asm> {rel,m100_noq,m50_noq}` (includes `match_riboseq.py` catalogue stop matches).
  3. `bash pipeline/release_pin/11_score_sweeps_v1.sh` -> `data/benchmark/release_pin/{by_length,short_orf_summary,negatives_*}.tsv`
     (also reads `$P/work/bench_v3/rules/{org}/{variant}/cliff_c200/preds.tsv`).
  4. `python pipeline/release_pin/12_gm_tables_v1.py`; `bash 13_gm_support_v1.sh` (-> protein_matches, riboseq_summary,
     ms_summary, rule_support in `data/gm12878/release_pin/`).
  5. Era table: `bash pipeline/derive/04_release_short_orfs.sh` (reads `work/bench_v3/runs/human/*/pytd_v7`,
     comparators in `work/bench_v3/runs`, `work/lnc/human_lnc.pick.gtf`, dev `work/bench_v3/rules/human/*/decoy_q*_longgc`)
     -> `data/benchmark/final_v7/lnc_by_era.human.v7.tsv`, `by_length.v7.tsv`, `sens_by_length.v7.tsv`, `lnc_riboseq_orfs.v7.tsv`.
  6. Post-hoc `--min-aa`: `bash pipeline/derive/08_minaa_sweep_final_v7.sh 2` (reads
     `work/gm12878_td/<asm>/final_v7_abinitio/...genome.gff3`; PYTHONPATH hardcoded to work/trees/PY_v7) ->
     `data/gm12878/final_v7/08_minaa_sweep_v7.tsv`; `bash pipeline/derive/08_no_overlap_by_length_v7.sh`.
- **Outputs**: `data/benchmark/release_pin/*`, `data/gm12878/release_pin/*`, `data/gm12878/final_v7/08_*` ->
  notebooks 04, 08 (Fig. 4, Fig. 7). Dev tables `data/benchmark/{by_length,short_orf_summary,lnc_by_era.human}.tsv`,
  `data/gm12878/gm12878_minaa_sweep.tsv` are history.
- **Cost**: part of E40 (stage C 6.3 h wall, 4 processes, all sweeps).

### E09 Start-codon refinement and Kozak context  (status: in paper + Supp. Fig. S1, Supp. Table S6; class a for release numbers, b for v6/PWM history)
- **Feeds**: Methods "Start-codon refinement"; Results start accuracy; Fig. 3, Supp. Fig. S1 (Kozak logos), Supp. Table S6
  (all drawn in notebook 03).
- **Inputs**: bench_v3 inputs; for PWM seed table, TransDecoder `sample_data` PASA/StringTie inputs in `$P/work/regression/*/inputs`.
- **Software**: release = defaults of v0.4.0 (alpha 0.4, w 17/bit); first-ATG column = `--start-contrast-alpha 0.4
  --start-kozak-weight 0` (grid cell `a0.4_w0`). v6.0.0 = Perl 66d4712 (`run_tool.sh perl_td_v600`, TD_V600 tree;
  `data/benchmark/summary.tsv`). Fixed/seeded PWM = PyTD 27bfb5e (`PY_fix`); dev labels efc1450/0886981/f5e6a39/eaea367
  (`run_tool.sh pytd_dev*`, `PYTD_DEV`) are history (c). PWM choice among later ATGs (`--start-atg-choice pwm`, efc1450) has no
  release option (E40 exception 4); S6 panel b now uses the release-model analysis `later_atg_kozak_v1.py` (E41).
- **Steps** (release):
  1. Release runs (E40 stage A) and the alpha x w grid (see E13).
  2. `PYTHONPATH=<PyTD tree> $P/env/td/bin/python pipeline/derive/03_start_refinement.py` -> `data/start_refinement/
     {kozak_models,training_mu,later_atg,seed_start_changes}.tsv` (reads `work/bench_v3/runs/*/complete/pytd_v7` and
     `$P/work/start_seed_runs`).
  3. `python analyses/start_param_sensitivity/later_atg_kozak_v1.py > data/start_refinement/later_atg_kozak.tsv`
     (reads `work/release_pin/bench`; tree path hardcoded).
  4. History (b): `bash analyses/start_refinement_stability/run_seeds.sh <PY_fix 27bfb5e tree> <regression inputs>
     $P/work/start_seed_runs "pasa stringtie" "1 2 3 4 5"`; `run_strategy_grid.sh <org> <variant> $P/work/start_grid/<org>.<variant>.tsv`
     (bench_v2 truth); `analyses/start_codon_dev/{pwm_audit,start_errors,kozak}/*.py`.
- **Outputs**: `data/start_refinement/*` -> notebook 03 (Fig. 3, S1, S6). Old `figures/kozak/*` (perl_td training sets) superseded.
- **Cost**: not recorded.

### E10 GC-calibrated long-ORF retention rule  (status: in paper, Supp. Table "Long-ORF length by G+C"; class a — code-independent simulation)
- **Feeds**: Methods "GC-calibrated long-ORF rule", "Background composition"; supplement table (notebook 01).
- **Inputs**: none (seeded i.i.d. simulation). Raw low-GC chunks: `$P/work/thresh_big/raw/gc<GC>_seed<seed>.tsv` (88 files).
- **Software**: `analyses/long_orf_threshold/simulate_random_orf_lengths.py` (stdlib only). The table L*(GC) is built into
  PyTD/Perl since f5b12e7/694e204; v0.4.0 = same table.
- **Steps**:
  1. `bash analyses/long_orf_threshold/run_low_gc_chunks.sh $P/work/thresh_big` (seeds GC*1000+c; 50/20/10/5/3 chunks of
     100,000 x 5 kb, `--partials`; was UGER, now local).
  2. `python analyses/long_orf_threshold/pool_chunk_quantiles.py $P/work/thresh_big/raw > lowgc_large_sample.tsv`
     (`--boot 200` for the CI).
  3. 50–80% rows: `sim_large/partial_L5000_gc*.tsv` — **command not recorded** (see Part 3); merge to
     `final_thresholds.tsv` by hand.
  4. Benchmark GC-estimate checks (launch15/16, Perl/PyTD parity) are history (c).
- **Outputs**: `data/long_orf_threshold/{final_thresholds,lowgc_large_sample,regenerated_thresholds,table_2017}.tsv`,
  `sim_runs/`, `sim_large/` -> notebook 01.
- **Cost**: not recorded (88 chunk jobs).

### E12 True CDSs that fail the coding criteria  (status: in paper, Results "True coding regions that fail...", Fig. 5; class a for v7 tables, c for the development version)
- **Feeds**: results_coding_failures.tex, Fig. 5a–c (notebook 05), Discussion.
- **Inputs**: `pytd_v7` runs of bench_v3 complete and complete_fullutr (+ wRand/wLnc for altrules), with `hexamer.scores`,
  `decoy_qvalues.tsv`, `longest_orfs.*`, `__checkpoints_TDpredict/top_training_orfs.cds`; homology flags from
  `work/bench_v3/runs/<org>/<input>/pytd_v7_homology/{blastp.outfmt6,pfam.domtblout}`; annotations `$P/data/annot`.
- **Software**: release runs (= v0.4.0 defaults). Analysis imports the PyTD tree (`PYTD_TREE`, default work/trees/PY_v7).
- **Steps**:
  1. Release runs at the E25 paths: `JOBS=6 $T/benchmark_v2/scripts/run_all.sh data/bench_v3/inputs work/bench_v3/runs pytd_v7`
     and `THREADS=4 JOBS=3 ... pytd_v7_homology` (run from `$P`; needs `$P/data/homology` DIAMOND/Pfam DBs).
  2. `bash pipeline/derive/05_long_orf_failures_v7.sh` (bench complete, complete_fullutr, altrules, summarize,
     `05_length_vs_homology_v7.py`, `05_fig5_tables_v7.py`; full outputs to `$P/work/long_orf_failures_v7/`).
- **Outputs**: `data/long_orf_failures/v7/` -> notebook 05 / Fig. 5. Dev tables in `data/long_orf_failures/*.tsv` and the
  GM12878 part (`gm12878_*.tsv`, final_abinitio, 48c870c) are history; the GM12878 part is not in the manuscript.
- **Cost**: 1–3 h, 2 processes (script header); plus the homology runs (DIAMOND + hmmsearch, 84 inputs).

### E13 Start-refinement parameter sensitivity (alpha, w)  (status: in paper, Fig. 3, Supp. Fig. S2, Fig. 7; class a)
- **Feeds**: Methods sec:start (choice of 0.4/17), Results start refinement, Fig. 3, S2 (notebook 03), Fig. 7 GM grid (notebook 08).
- **Inputs**: 16 stranded complete/fragments inputs; GM12878 LRAA/StringTie/Trinity. Reuses LongOrfs of the `rel` runs (symlinked).
- **Software**: release; flags `--start-contrast-alpha A --start-kozak-weight W` (bits), `--no-refine-starts`.
  Dev grid (PY_dev 48c870c, w per nat: `run_grid.py`, `run_gm12878.sh`, `summarise.py`, `compare_bits.py`) is history.
- **Steps**:
  1. E40 stage A (`rel`), then stage E: `python analyses/start_param_sensitivity/run_grid_v2.py --jobs 4`
     (13 alpha x 6 w + norefine = 1,264 Predict runs -> `$P/work/release_pin/start_grid/runs`; PYTD hardcoded).
  2. GM grid: `run_gm_v1.sh <asm> {a0.2,a0.4,a0.6,a1}_w{0,17,35} | norefine` (stage B).
  3. `python summarise_v2.py` (checks against `work/bench_v3/runs/*/pytd_v7/preds.tsv` too), `python open_end_mix_v2.py >
     data/start_param_sensitivity/release/open_end_mix.tsv`, `python kozak_score_distribution_v2.py >
     data/start_param_sensitivity/release/kozak_score_distribution.tsv` (usage line says results/; E40 says release/).
- **Outputs**: `data/start_param_sensitivity/release/{grid_scores,optima,plateau,mix,loo,gm12878,checks,open_end_mix,kozak_score_distribution}.tsv`
  + heatmaps; `data/gm12878/release_pin/gm_runs.tsv`. Dev `data/start_param_sensitivity/*.tsv` = history.
- **Cost**: stage E 17.6 h wall at 4 processes (2026-11-01 08:24 -> 11-02 01:58). Dev grid: 147,507 CPU-s.

### E17 Published critiques tested old vs new; bootstrap CIs  (status: not in manuscript (planned supplement); class c, interim)
- **Feeds**: nothing cited (no tex reference to E17). `bootstrap_ci.py` functions are reused by `02_first_atg_baseline.py`.
- **Inputs**: bench_v3; v5.7.1 runs (`run_tool.sh perl_td_v571`, tree TD_v571 = 5791fde, on origin); "new" = `pytd_dev_v2` (48c870c).
- **Steps**: `benchmark_v2/scripts/critiques_old_vs_new.py`, `benchmark_v2/scripts/bootstrap_ci.py` (exact command lines not recorded).
- **Outputs**: `data/benchmark/critiques/{critiques_old_vs_new.interim,critiques_lnc}.tsv`, `data/benchmark/bootstrap/interim_20261023.*`.
  If revived, rerun with `pytd_v7` as "new" (class a).

### E18 Stringency sweep for a contrast-based --precise mode  (status: not adopted; class c)
- **Software**: PY_dev 48c870c, checkpoint reruns of `pytd_dev_v2` (`--short-orf-length`, `--short-orf-min-contrast`;
  options removed in be23241). TD2-style minimum-length table: `minlen_run.sh {pytd,perl_v571,td2} <org> <variant> <m>`.
- **Steps**: `benchmark_v2/scripts/stringency_sweep_pytd.sh <org> <variant> s<L|all>_c<C>[_noauto]` (needs
  `work/bench_v3/runs/<org>/<variant>/pytd_dev_v2` checkpoints); `score_negatives.py`, `plot_stringency.py`, `prediction_features.py`.
- **Outputs**: `data/benchmark/stringency/*` (log only). Released `--precise` is q 0.01 (E22).

### E19 Shape of the short-ORF threshold (cliff, taper, fixed, inverse)  (status: not adopted; class c — but its `cliff_c200` runs are the paper's historical comparator, class b)
- **Feeds**: `cliff_c200` = "fixed short-ORF test, unreleased development version of v7": Fig. 4a–c,e, Supp. Table S11 /
  Fig. S3 ("v7 dev., fixed test"), Results 50–74-codon sensitivity 0.27 (notebooks 04, S1; via `11_score_sweeps_v1.sh`,
  `04_release_short_orfs.sh`).
- **Software**: PyTD 652d5e4 (branch exp/short-orf-taper; commit is on origin/main history). `short_orf_rule_run.sh`
  defaults to `work/trees/PY_taper`, **which is now checked out at 3a793ff** — set `PYTD_TAPER` to a worktree at 652d5e4.
- **Steps**: `PYTD_TAPER=<652d5e4 tree> bash benchmark_v2/scripts/short_orf_rule_run.sh cliff_c200 $P/data/bench_v3/inputs/<org>/<variant>.fa
  $P/work/bench_v3/rules/<org>/<variant>/cliff_c200 <org>` for the 20 wRand/wLnc inputs (job list `data/benchmark/rules/jobs.txt`).
  Scoring: `score_negatives.py`, `plot_rules.py`, `sens_by_length.py`, `score_protein_matches.py`, `gm_rule_support.py`.
- **Outputs**: `data/benchmark/rules/{rules_scores,rules.medians,sens_by_length}.tsv`, `gm12878_*` (log); the cliff_c200
  `preds.tsv` files themselves are only in `$P/work/bench_v3/rules` (see Part 2).

### E20 Target-decoy FDR rule  (status: adopted, in paper; class a for released rule, b for the strata prototype)
- **Feeds**: Methods (decoy rule), Results short ORFs, Discussion, Abstract; notebook 04. "ORF-length strata alone gained
  5–6 points less" = prototype `data/benchmark/rules/decoy_proto/` (E40 exception 3).
- **Software**: release = v0.4.0 default (`--max-qvalue 0.05`, GC-table long-ORF accept). Dev `decoy_q<Q>[_long..]` runs
  (PyTD 1b9a291/e3b9e35) are history. Long-ORF rule off in release: `--retain-long-orfs-mode strict
  --retain-long-orfs-length 1000000` (E40 `nolong_q0.05/0.1`).
- **Steps**: E40 stages A/C (`rel`, `q*`, `nolong_q*`), `11_score_sweeps_v1.sh` (calibration = 1 − precision_stop in
  `negatives_scores.tsv`; `sens_by_length.tsv`; `lnc_riboseq_orfs.tsv`). Prototype (b):
  `decoy_fdr_prototype.py <run_dir> <shuffled_run_dir> <truth.tsv> <gcode> <label>` on `cliff_c0` runs (PyTD 652d5e4).
- **Outputs**: `data/benchmark/release_pin/*`; dev `data/benchmark/rules/*` (decoy_proto and `lnc_riboseq_orfs.tsv` still read by notebook 04).

### E21 Presets for small inputs (held-out test)  (status: in paper, Results "Small inputs and presets"; class a)
- **Feeds**: results_short_orfs.tex §presetres, Methods, Discussion; notebook 04.
- **Inputs**: held-out split **in `$P/work/bench_v3/presets/ref/{org}.{A,Amix}.fa` and `tests/<org>/*.{fa,truth.tsv}`**;
  job lists `data/benchmark/presets/jobs{,2}.txt` (identical to the work copies; absolute paths into work/).
  `data/benchmark/presets/setup_presets.py` regenerates A refs and the 50/200/1000/5000 wRand tests only; the `Amix`
  refs, the 500-sequence and `codingonly` sets came from an uncommitted inline script.
- **Software**: release `build-preset` and `predict --preset DIR`; dev 6b024ca runs are history.
- **Steps**: E40 stage D: `run_preset_v1.sh ref <org> {A,Amix}` (LongOrfs -S, Predict, build-preset ->
  `work/release_pin/presets/presets/<name>`), then `run_preset_v1.sh test {self|preset_<name>} <test.fa> <org>` for every
  non-`cliff_c200` job line; `bash pipeline/release_pin/15_score_presets_v2.sh` (reads `work/bench_v3/presets/tests`).
- **Outputs**: `data/benchmark/presets_release/{preset_scores.tsv,preset_meta/humanmix/preset.json}` -> notebook 04.
- **Cost**: stage D 0.8 h wall (4 processes).

### E22 q-value modes (default 0.05, --precise 0.01, --sensitive 0.2)  (status: in paper; class a)
- **Feeds**: Abstract, Methods, Results short ORFs and GM12878 parameters, Discussion; notebooks 04, 08.
- **Software**: `--precise`, `--sensitive`, `--max-qvalue Q` (v0.4.0). Dev runs (8dbe79b) history.
- **Steps**: E40 stage A (`precise`, `sensitive`), stage C (`q0.02 q0.1 q0.3 q0.5`); longest-ORF baselines from
  `rel/longest_orf.preds.tsv` (>= 303 nt filter in 11_score_sweeps_v1.sh); GM: `run_gm_v1.sh <asm> {precise,sensitive}`;
  scoring `11_score_sweeps_v1.sh`, `13_gm_support_v1.sh` (protein matches via `score_protein_matches.py`).
- **Outputs**: `data/benchmark/release_pin/{negatives_scores,negatives_medians,lnc_call_lengths.human}.tsv`,
  `data/gm12878/release_pin/protein_matches.tsv`. Dev `data/benchmark/qmodes/*` history. GM12878 q 0.1/0.3 exist only as dev
  runs (not quoted in the tex as far as checked).

### E24 Release implementation: decoy rule + bundled presets, Perl/PyTD parity  (status: in paper (Methods, Implementations); class a)
- **Software**: Perl 41fc5ee, PyTD 3a793ff (= d3892b5 = v0.4.0 output). Bundled presets ship inside the packages.
- **Steps**: parity at the pins: `regression/run_regression.sh <TD tree> <PyTD tree> <outdir>` and
  `regression/parity_check.py` (E40: `data/regression/release_d3892b5_41fc5ee.summary.tsv`, `data/release_pin/parity_intermediates.txt`).
  Original 7-case parity: `benchmark_v2/scripts/parity_decoy_run.sh` with `data/benchmark/release/parity_jobs.txt`. Bundled preset
  builds: `make_bundled_preset_refs.py`, `data/benchmark/release/bundled_preset_run_ref.sh`; smoke `bundled_preset_smoke_jobs.txt`.
- **Outputs**: `data/benchmark/release/*`, `data/benchmark/presets/bundled_presets.tsv` (notebook 01).

### E25 Final release reruns: benchmark and GM12878  (status: in paper; class a — byte-identical to d3892b5, E40)
- **Feeds**: Table 1, Fig. 2/4/5/6, homology paragraph, GM12878 application; notebooks 01–08, S1; run-time ratio 2.71 (not re-timed).
- **Steps** (run from `$P`, with `PYTD_V7`/`TD_V7` set):
  1. Benchmark — `work/bench_v3/launch_v7.sh` (only in work/, content):
     `JOBS=6 $RA data/bench_v3/inputs work/bench_v3/runs pytd_v7 perl_td_v7`;
     `VARIANTS="complete.wRand fragments.wRand complete.wLnc fragments.wLnc" JOBS=6 $RA ... pytd_v7_precise pytd_v7_sensitive`;
     `THREADS=4 JOBS=3 $RA ... pytd_v7_homology` (`RA=$T/benchmark_v2/scripts/run_all.sh`).
  2. GM12878: `genome_annotation/scripts/gm12878_final_v7.sh <asm> <fa> <structures> abinitio 8 [precise|sensitive]`; Trinity homology.
  3. Scores: `data/benchmark/final_v7/score_v7.sh <org> <variant> <tool>`; `score_negatives.py ...`; `compare_homology.py`;
     `pipeline/derive/04_lnc_calls_by_autoaccept_length.py`, `05_homology_v7h_table.sh` (pytd_v7h_homology, E33), `07_gm12878_evidence_v7.sh`.
  4. Equivalent fresh path: E40 `launch_v1.sh` stages A–B + `10_byte_identity_v1.py` (only meaningful while old runs exist).
- **Outputs**: `data/benchmark/final_v7/*`, `data/benchmark/homology*_vs_abinitio.tsv`.
- **Cost**: run-time ratio data in `runtime_sec.txt`; PyTD v7 1.75x `pytd_dev_v2` (decoy step).

### E26 Frame contrast: forward frames only vs all five  (status: in paper (Methods, retest of v5 choice); class b)
- **Feeds**: methods_algorithm.tex lines ~105–117; notebook 01 (`data/benchmark/frames/*`).
- **Software**: PyTD **83a4806** (exp/contrast-frames = 3a793ff + `--contrast-frames {forward,all}`), **local only**;
  patch committed as `data/benchmark/frames/pytd_contrast_frames.patch` (`git am` onto 3a793ff). Comparison side = release.
- **Steps**: `PYTD_FRAMES=<83a4806 tree> bash data/benchmark/frames/launch2.sh` (84 inputs, q 0.05), `launch_q.sh` (36 inputs,
  q 0.1/0.2), `launch_gm12878_frames.sh` (`genome_annotation/scripts/gm12878_frames.sh`); per-run `frames_run.sh <fa> <outdir> <org> all [args]`;
  scoring `data/benchmark/frames/score_all.sh <org> <variant>` (+ release rows of E25); `score_negatives.py`, `score_protein_matches.py`,
  `match_riboseq.py` (exact negatives command not recorded).
- **Outputs**: `data/benchmark/frames/{all_scores,negatives_scores,negatives_curve,gm12878_frames_protein_matches}.tsv`.
- **Cost**: all-frames median 1.41x release run time; 156 benchmark runs + 5 GM runs.

### E40 (context) Release-pin reruns — the class-a engine for all of the above
`bash $T/pipeline/release_pin/launch_v1.sh` (stages A–E, 4 processes, ~32 h wall, ~80 GB in `$P/work/release_pin`), then
`11_score_sweeps_v1.sh`, `12_gm_tables_v1.py`, `13_gm_support_v1.sh`, `15_score_presets_v2.sh`, `summarise_v2.py`,
`open_end_mix_v2.py`, `kozak_score_distribution_v2.py`. Stage timings: A 3.4 h, B 4.2 h, C 6.3 h, D 0.8 h, E 17.6 h.
**Work-dir inputs read by the launcher**: `work/bench_v3/runs/*/*/pytd_v7h_homology` (list of retonly inputs) and
`pytd_v7_homology/{blastp.outfmt6,pfam.domtblout,…/longest_orfs.pep}`; `work/gm12878_td/trinity/final_v7_homology` hits;
`work/bench_v3/presets/{ref,tests,jobs*.txt}`; `11_score_sweeps_v1.sh` reads `work/bench_v3/rules/*/*/cliff_c200/preds.tsv`.

---------------------------------------------------------------------------------------------------------
# PART 2 — things in $P/work that are inputs or expensive to regenerate (for E04–E26 + E40)

| path (under $P/work) | size | what | needed by | copy elsewhere? | recommendation |
|---|---|---|---|---|---|
| bench_v3/presets/ref | 132M | E21 held-out reference halves A / Amix | run_preset_v1.sh ref (E40 D) | no | **keep** -> `$P/data/presets_heldout/ref` (Amix not regenerable: inline script) |
| bench_v3/presets/tests | 102M | E21 test sets (+truth) incl. 500 and codingonly | run_preset_v1 test, 15_score_presets_v2.sh | no | **keep** with ref; sed job lists' absolute paths |
| bench_v3/presets/jobs{,2}.txt | 24K+ | job lists | launch_v1 stage D | yes (data/benchmark/presets) | drop (committed copy identical) |
| bench_v3/runs/*/*/pytd_v7_homology/{blastp.outfmt6,pfam.domtblout} | 263M (hits only; dirs 9.6G) | DIAMOND/Pfam hits, 84 inputs | E40 retonly, E12 homology flags, E25 homology | no | **keep hits + longest_orfs.pep** (re-search = hours of DIAMOND+hmmsearch; DB versions fixed) |
| gm12878_td/trinity/final_v7_homology/{blastp.outfmt6,pfam.domtblout} | 18M | Trinity hits | run_gm_v1.sh retonly | no | **keep** |
| bench_v3/rules/*/*/cliff_c200/preds.tsv (+ gm12878/*/cliff_c200) | 6.6M (dirs 1.8G) | historical "fixed short-ORF test" predictions | 11_score_sweeps_v1.sh, 04_release_short_orfs.sh, notebooks 04/S1 | no | **keep preds.tsv** (regeneration needs 652d5e4 tree, 22 runs) |
| bench_v3/rules/human/*/decoy_q*_longgc/preds.tsv | small | dev q sweep | 04_release_short_orfs.sh lnc_call_lengths | no | keep preds.tsv or drop that block from the script |
| bench_v3/rules/decoy_proto, curves | 1.3M | E20 prototype/oracle | committed copies in data/benchmark/rules | yes | drop |
| lnc/human_lnc.pick.gtf (+mouse_lnc.pick.fa, pick_lnc.py) | 7M (+16M) | lncRNA negatives' GTF for the era table | lnc_calls_by_era.py (04_release_short_orfs.sh) | no (fa content is inside bench_v3 wLnc inputs) | **keep** -> `$P/data/lnc/` |
| thresh_big/raw | 4.5M | E10 low-GC simulation chunks | pool_chunk_quantiles.py | no | **keep** (cheap; regenerable by seed but slow) |
| start_seed_runs | 238M | v6-era PWM seed runs | 03_start_refinement.py seed table | table committed | keep only if seed table must be re-derived; else drop |
| start_grid, analysis/{pwm_audit,kozak,fork_findings.md} | 1.2M, 2M, 5.7M, 24K | E09 uncommitted findings/sim tables | E09 history | no | **copy into $T (e.g. data/start_codon_dev/)** — small, only record |
| analysis/start_errors | 99M | E09 start-error tables | history | no | optional archive |
| trees/PY_frames (branch exp/contrast-frames 83a4806) | 105M | E26 code | E26 rerun | patch committed | **push branch to GitHub** (or `git bundle`); patch is a fallback |
| trees/PY_taper (exp/short-orf-taper; now at 3a793ff) | — | E18–E22 code | cliff_c200 | commits on origin/main | push branch name exp/short-orf-taper (tip 8dbe79b) for convenience; else drop |
| trees/PY_rel_d3892b5, TD_rel_41fc5ee, PY_v7, TD_v7, TD_v571, TD_base, PY_fix, PY_dev, PY_48c870c | ~0.1G each | worktrees | defaults in scripts | commits on origin | drop; recreate with `git worktree add` (Part 1 set-up) |
| bench_v3/runs/*/*/pytd_v7{,_precise,_sensitive} | 8.5G (pytd_v7) | E25 release runs | derive scripts 02–05, 08, notebooks' tables | tables committed | drop; regenerable (launch_v7.sh ~hours) |
| release_pin | 80G | E40 reruns | committed tables | tables committed | drop; ~32 h at 4 processes to regenerate |
| long_orf_failures_v7 | 232M | E12 full outputs | compact tables committed | yes | drop (1–3 h) |
| start_param_grid | 4.4G | E13 dev grid | history | tables committed | drop |
| bench_v3/{stringency,minlen,frames,rules (rest)} | 9.2G, 3.1G, 14G, ~45G | E18/E19/E20/E22/E26 dev runs | history | scored tables committed | drop |
| bench_v3/presets/runs, refruns | 3.3G | E21 dev runs | history | tables committed | drop |
| gm12878_td/eval/minaa, gm12878_td/ribo | 95M, 16M | E04 dev GM | history | tables committed | drop |
| td2_private | 3.5G | private TD2 points (E19/E22 notes, e.g. qmodes/E22_td2_notes.md) | E39 owner | — | defer to E39 owner; keep notes .md |
| bench_v3/launch_v7.sh | <1K | E25 launcher | E25 rerun | no | **commit to benchmark_v2/launch_v3/** |

---------------------------------------------------------------------------------------------------------
# PART 3 — gaps

1. **E21 held-out sets not regenerable**: `Amix` refs, 500-sequence and `codingonly` test sets came from an uncommitted
   inline script; only `setup_presets.py` (A refs, 50/200/1000/5000 wRand tests) is committed. Preserve `work/bench_v3/presets/{ref,tests}`.
2. **E10**: commands for `sim_runs/` and `sim_large/` (the 50–80% GC rows of the adopted table) not recorded (n-seqs, seed,
   `--length-fasta`); `final_thresholds.tsv` merged by hand, no script; no committed bootstrap CI file.
3. **E09**: Trinity seed runs missing; coding-gate result (notes/04) has no table/command; additive-form/weight-25 choice and
   multiplicative 0.15 rest on commit messages; trim_all GM12878 runs (alpha 1000) have no saved command; `pwm_audit/FINDINGS.md`,
   `fork_findings.md`, `start_grid/`, `kozak/` sims only in work/ (fork_reg scripts lived in a session scratchpad);
   `start_grid/tthermophila.fragments.tsv` incomplete.
4. **E04**: GM12878 m100 baseline (`kozak_add`) had no saved script; dev runs were on uncommitted working trees (superseded by E40).
5. **E17**: exact command lines for `critiques_old_vs_new.py` / `bootstrap_ci.py` not recorded; tables interim (CodAn 7/8).
6. **E20**: prototype needs a `<shuffled_run_dir>`; how the shuffled inputs/runs (`cliff_c0`) were made is not recorded.
7. **E26**: `negatives_scores.tsv` / `negatives_curve.tsv` command not recorded (only "score_negatives.py"); code branch local-only.
8. **E25**: `launch_v7.sh` exists only in `work/bench_v3/`; GM12878 Trinity homology launcher not cited.
9. **E40**: `data/release_pin/longest_orf_baseline_check.tsv` has no generating script; `10_byte_identity_v1.py` compares against
   work/ runs and cannot be rerun after deletion; `kozak_score_distribution_v2.py` usage line names `results/` not
   `data/start_param_sensitivity/release/`.
10. **Stale figure scripts**: `figures/fig{3,4,5,7}_*/make_fig*.py` read dev tables (`short_orf_summary.tsv`, `data/start_param_sensitivity/*.tsv`);
    the manuscript figures are drawn by notebooks 03/04/05/08. REPRODUCE.md should point to the notebooks only.
11. **Hardcoded tree paths** (Common set-up) prevent a clean "latest PyTD" rerun without edits/symlinks.
12. **E22**: GM12878 q 0.1 / 0.3 rows exist only as dev runs (`work/bench_v3/rules/gm12878`) — fine if not quoted; confirm.
