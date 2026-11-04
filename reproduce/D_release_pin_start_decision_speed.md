# Reproduction detail D: release-pin reruns (E40), start-decision work (E29, E33–E36, E41), speed (E42)

> **Paths into `$P/work`.** This section was written on 2026-11-03, before `$P/work` was deleted. Files that were kept are now under `$P/data/work_keep/` at the same path relative to `work/`. For example, `$P/work/bench_v3/runs/human/complete/pytd_v7_homology/pfam.domtblout` is now `$P/data/work_keep/bench_v3/runs/human/complete/pytd_v7_homology/pfam.domtblout`. Everything else under `$P/work` has to be regenerated with the steps given here. The tables in Part 2 record what was kept and why; see `../REPRODUCE.md` for the full list.

`$P=/broad/hptmp/bhaas/TRANSDECODER_PAPER`, `$T=$P/TransDecoderPaper`. All scripts hardcode `P=` and most write to
`$P/work/...`. Rerunning elsewhere means changing `P=` at the top of each script; directories flagged **[work]**
below are outputs or inputs under `$P/work` that are deleted when `$P/work` goes.

## 0. Code versions and how to use pip-installed PyTransDecoder 0.4.0

Package-tree identity, checked with `git rev-parse <rev>:pytransdecoder` in `$P/repos/PyTransDecoder`:

| rev | package tree | note |
|---|---|---|
| a5b24c3 (homology start assist) | 4dde030 | `work/trees/PY_hstart` |
| d3892b5 (paper release pin) | 4dde030 | `work/trees/PY_rel_d3892b5` |
| tag v0.2.0 (275ff1b) | 4dde030 | same code as d3892b5 |
| v0.2.1 (PyPI) | cf4765f | differs only in `__init__.py` (version string) |
| v0.3.0 / v0.4.0 (PyPI) | 5d3a7bd / 2785a28 | refactored and faster. All output files byte-identical to d3892b5 (E42: 167/167 + 226/226 comparisons, any `--cpu`) |
| 3a793ff (`pytd_v7` runs, E25) | 9d9a84a | no homology assist. Ab initio output byte-identical to d3892b5 (E40) |

Perl: TransDecoder `devel` 41fc5ee (v7.0.0), on GitHub `origin/devel` and `origin/master`.

CLI entry points (`pyproject.toml`, unchanged d3892b5 -> v0.4.0): `pytransdecoder = pytransdecoder.cli:cli`
(subcommands `longorfs`, `predict`, `genome-map`, `build-preset`; `--version`) and `pyTransdecoder =
pytransdecoder.cli:pipeline_main`. All flags used by the E40 scripts exist in 0.4.0: `--precise`, `--sensitive`,
`--max-qvalue`, `--retain-long-orfs-mode/-length`, `--retain-blastp-hits`, `--retain-pfam-hits`,
`--no-homology-start-assist`, `--no-refine-starts`, `--start-contrast-alpha`, `--start-kozak-weight`, `--preset`.
0.4.0 adds `--cpu N` to longorfs and predict. Output does not depend on N.

**Option A, no script edits (recommended).** The workers take the tree from `PYTD_REL` and run `git rev-parse`
inside it. A git checkout of the tag therefore works as is:
```bash
git clone --branch v0.4.0 https://github.com/TransDecoder/PyTransDecoder.git $P/trees/PY_v0.4.0
export PYTD_REL=$P/trees/PY_v0.4.0      # read by run_bench_v1.sh, run_gm_v1.sh, run_preset_v1.sh
```
`$P/env/td` (Python 3.11.17, numpy 2.4.6, Biopython 1.88) satisfies the 0.4.0 dependencies (numpy >= 1.21).
- `analyses/start_param_sensitivity/run_grid_v2.py` hardcodes `PYTD = f"{P}/work/trees/PY_rel_d3892b5"` (line 28). Edit that line.
- `launch_v1.sh` stage A passes `$P/work/trees/PY_rel_d3892b5` to `run_regression.sh`. Edit that path too.
- The Perl tree for the regression: `git clone -b devel https://github.com/TransDecoder/TransDecoder; git checkout 41fc5ee`.

**Option B, the pip package.** `python -m venv $P/env/pytd040 && $P/env/pytd040/bin/pip install pytransdecoder==0.4.0`. Then, in each worker:
- replace `PYTHONPATH=$PYTD $PY -m pytransdecoder.cli` (run_bench_v1.sh) and `td() { $PY -m pytransdecoder.cli "$@"; }`
  (run_gm_v1.sh, run_preset_v1.sh, 13_gm12878_trinity_v7h_v1.sh, 11_run_pytd_v7h_homology_v1.sh) with
  `$P/env/pytd040/bin/pytransdecoder`;
- replace `(cd $PYTD && git rev-parse --short HEAD) > pytd_version.txt` with
  `pytransdecoder --version > pytd_version.txt`. Under `set -e`, `git rev-parse` on a non-git directory kills the job;
- in `run_grid_v2.py`, set `args = ["pytransdecoder", "predict", ...]`;
- in `later_atg_kozak_v1.py` and `kozak_score_distribution_v2.py`, drop the `sys.path.insert(... PY_rel_d3892b5)` line.
  They only use `pytransdecoder.core.kozak.read_model/score`, which 0.4.0 still has.

**Offline re-derivation scripts (E29, E33, E34, E35, E36).** These import `TransDecoderPredict` and private methods
(`_select_best_orfs`, `_revise_start_sites`), plus `core.homology_starts`, `core.decoys` and `core.prng`. `predict.py`
was refactored in 0.3.0/0.4.0 (73+/130- lines vs d3892b5). Use `pip install pytransdecoder==0.2.1` for them (package =
a5b24c3/d3892b5 apart from the version string), or `git checkout v0.2.0`. They were not tested against 0.4.0.

---

# PART 1. Experiments

### E40 One code version throughout: release-pin reruns (status: in paper, the backbone of every TD/PyTD number; reproducible with latest PyTD: yes, 0.4.0 output = d3892b5)
- **Feeds:** every TransDecoder/PyTD number, table and figure. Specifically: benchmark short-ORF/q-value/min-length
  paragraphs and Fig. 4; GM12878 start-refinement, min-length and q-value paragraphs (now including Trinity); held-out
  presets; the alpha x w grid (Supplementary heatmaps, LOO); Table S6 "first ATG" column (grid cell `a0.4_w0`); notebooks 01-08 and S1.
- **Inputs (all in `$P/data`, kept):**
  - benchmark: `data/bench_v3/inputs/<org>/<variant>.fa` and `.truth.tsv` (108 inputs) and `benchmark_v2/organisms.tsv` (genetic code, taxid);
  - GM12878: `data/gm12878/{star_stringtie,longread_lraa/lraa_denovo,trinity}/...` (FASTA plus GTF, or the Trinity minimap2 BAM),
    `data/annot/human/MANE.GRCh38.v1.5.ensembl_genomic.gtf.gz`, `data/gm12878/ref/gencode.v50...gtf`,
    `data/riboseq_orfs/phase{1,2}/`, the Ribo-seq BAM `data/gm12878/riboseq_cenik2015/footprints.genome.bam`, and
    `data/gm12878/proteome_sinitcyn2023/{search/GM12878/peptides.txt,refdb/*}`.
  - **[work] inputs:**
    - `work/bench_v3/presets/{ref/<org>.{A,Amix}.fa, tests/<org>/*.fa, jobs.txt, jobs2.txt}`: E21 held-out halves;
    - `work/bench_v3/runs/<org>/<inp>/pytd_v7_homology/{blastp.outfmt6,pfam.domtblout}`: 84 inputs, `retonly`;
    - `work/gm12878_td/trinity/final_v7_homology/{blastp.outfmt6,pfam.domtblout}`: GM `retonly`;
    - `work/bench_v3/rules/<org>/<v>/cliff_c200/preds.tsv`: historical comparator, unreleased dev code, cannot be rerun with any release.
- **Software:** PyTD d3892b5 (or 0.4.0), Perl TD 41fc5ee (regression only), `$P/env/td` (Python, DIAMOND 2.2.8,
  HMMER 3.4), `$P/env/riboseq`, `$P/env/proteo` (13_gm_support).
- **Steps:** `bash pipeline/release_pin/launch_v1.sh`. It runs 4 processes, local, normal priority, never resumes, and logs to `work/release_pin/launch_v1.log`. Stages:
  - **A:**
    - `regression/run_regression.sh <TD tree> <PyTD tree> work/release_pin/regression` (5 sample sets);
    - `run_bench_v1.sh <org> <inp> rel` on all 108 inputs;
    - `precise` and `sensitive` on the 20 negative-containing inputs ("NEG20": `{complete,fragments}.{wRand,wLnc}` that exist).
  - **B:** `run_gm_v1.sh <lraa|stringtie|trinity> rel`, then on each assembly the configs `a0.4_w17 norefine a0.4_w0 precise
    sensitive m100_noq m50_noq` and `a{0.2,0.4,0.6,1}_w{0,17,35}`, plus `run_gm_v1.sh trinity retonly`. Except for
    `rel`/`m*_noq`, the configs symlink rel's LongOrfs files. Each run then does predict, `genome-map`,
    `eval_genome_cds.py` and `analyses/riboseq/match_riboseq.py`.
  - **C:** on NEG20, labels `m{100,90,80,70,60,50}_noq` (LongOrfs `-m N`, `--max-qvalue 1000000`), `q{0.02,0.1,0.3,0.5}`
    and `nolong_q{0.05,0.1}` (`--retain-long-orfs-mode strict --retain-long-orfs-length 1000000`). Also `retonly` on
    the 84 inputs that have `pytd_v7_homology`: `--retain-blastp-hits/--retain-pfam-hits <hits> --no-homology-start-assist`,
    after `cmp` of `longest_orfs.pep` against the hits run.
  - **D:**
    - `run_preset_v1.sh ref <human|mouse|athal|dmel> <A|Amix>` (8 `build-preset` references);
    - `run_preset_v1.sh test <self|preset_X> <fa> <org>` on all jobs in `jobs{,2}.txt` except `cliff_c200` (728 tests).
  - **E:** `analyses/start_param_sensitivity/run_grid_v2.py --jobs 4`. 16 inputs (stranded complete and fragments x 8 orgs)
    x (13 alpha x w {0,3.5,7,17,35,69} + `norefine`) = 1,264 Predict runs on rel's symlinked LongOrfs, each scored by `score_predictions.py`.
  - **Outputs [work]:** `work/release_pin/{bench/<org>/<inp>/<label>, gm/<asm>/<cfg>, presets/{refruns,presets,runs}, start_grid/runs, regression, logs}`.
- **Scoring (commit outputs):**
  - `pipeline/release_pin/10_byte_identity_v1.py > data/release_pin/byte_identity.tsv`. Historical check against
    `work/bench_v3/runs/*/pytd_v7*` and `work/gm12878_td/*/final_v7_*`; it cannot be rerun once those are gone.
  - `11_score_sweeps_v1.sh` -> `data/benchmark/release_pin/{negatives_scores,negatives_medians,by_length,short_orf_summary,sens_by_length,lnc_call_lengths.human,lnc_riboseq_orfs}.tsv`.
    Needs `cliff_c200` preds [work].
  - `12_gm_tables_v1.py` -> `data/gm12878/release_pin/{gm_runs,no_overlap_by_length}.tsv`.
  - `13_gm_support_v1.sh`, which calls `14_rule_support_v1.py` -> `data/gm12878/release_pin/{protein_matches,riboseq_summary,ms_summary,rule_support}.tsv`.
    Intermediates go to `work/release_pin/gm_scoring`.
  - `15_score_presets_v2.sh` -> `data/benchmark/presets_release/preset_scores.tsv`, `preset_meta/humanmix/preset.json`.
    It reads the test truth from `work/bench_v3/presets/tests` [work].
  - `analyses/start_param_sensitivity/summarise_v2.py --grid work/release_pin/start_grid --gm work/release_pin/gm --out data/start_param_sensitivity/release`.
    Then run `open_end_mix_v2.py > .../release/open_end_mix.tsv` and `kozak_score_distribution_v2.py > .../release/kozak_score_distribution.tsv`.
  - `regression/parity_check.py` -> `data/release_pin/parity_intermediates.txt`, plus `data/regression/release_d3892b5_41fc5ee.summary.tsv`.
  - Consumers: notebooks 01-08 and S1 (`check_text` rows), `pipeline/derive/09_supplementary_tables.py`.
- **Cost:** stage A start 2026-10-31 17:40 to END 2026-11-02 01:58, about 32 h wall at 4 processes on methods01
  (2,528 jobs A-D plus 1,264 grid runs, 0 failures). 0.4.0 at `--cpu 1` is 1.75-3.45x faster per run (E42).
  Disk: 80 GB (bench 37, start_grid 34, gm 6.1, presets 3.1).

**How the homology hit tables are produced from scratch** (`retonly`, `pytd_v7h_homology`, GM Trinity homology):
- Databases (in `$P/data/homology`, kept, 3.8 GB):
  - Swiss-Prot 2026_03 (02-Sep-2026; `uniprot_sprot.fasta.gz`, `reldate.txt`). Archive:
    `ftp.uniprot.org/pub/databases/uniprot/previous_releases/release-2026_03/`;
  - `sprot_minus_<taxid>.dmnd` for 9606, 10090, 3702, 7227, 3055, 284812, 36329, 312017;
  - Pfam-A 38.2 (2026-01; `Pfam-A.hmm`, needs `hmmpress` for hmmscan, not for hmmsearch). Archive:
    `ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam38.2/`.
- `sprot_minus_<taxid>` has no saved builder (E07 gap). E33 step 5 verified that it is Swiss-Prot minus exactly the
  entries with `OX=<taxid>` (descendant taxa kept; human 575,748 -> 555,317). Reconstruction (not the original command):
  `zcat uniprot_sprot.fasta.gz | awk -v t=<taxid> '/^>/{keep = ($0 !~ (" OX=" t " "))} keep' > sprot_minus_<taxid>.fasta; diamond makedb --in sprot_minus_<taxid>.fasta -d sprot_minus_<taxid>`.
  The awk pattern needs care: an entry whose `OX=` is the last field has no trailing space. Check the counts against E06/E33.
- Benchmark hits: `benchmark_v2/scripts/run_tool.sh` tool `pytd_v7_homology` (THREADS=4). Run on each input's
  LongOrfs `longest_orfs.pep` (LongOrfs output is identical in every version above):
  `diamond blastp --query longest_orfs.pep --db sprot_minus_<taxid>.dmnd --max-target-seqs 1 --evalue 1e-5 --outfmt 6 --threads 4 --out blastp.outfmt6`.
  Then split the peptides round-robin into 4 chunks and run
  `hmmsearch --cpu 1 --cut_ga --noali -o /dev/null --domtblout cNNN.pep.domtblout Pfam-A.hmm cNNN.pep` on each.
  Concatenate the chunk outputs into `pfam.domtblout`. The search is hmmsearch, not hmmscan. GA thresholds are bit
  scores, so the domain set does not depend on chunking (E33 step 3 check); only the E-value columns do.
- GM12878 Trinity: `genome_annotation/scripts/gm12878_final_v7.sh trinity <fa> <bam> homology <threads>`
  (DIAMOND as above vs `sprot_minus_9606`; hmmsearch `--cpu $threads --cut_ga --noali` unchunked).
- A different Swiss-Prot/Pfam release or DIAMOND version changes the hits. The retention-only and assist numbers
  would then move slightly. To reproduce exactly, keep the hit tables (Part 2).

### E41 Later in-frame ATGs scored by start context (status: in paper, Methods "Start-codon refinement" plus Table S6b; reproducible with latest PyTD: yes)
- **Feeds:** Methods sentence after "Only the first in-frame ATG is considered"; Supplementary Table S6b; notebook 03 "Later ATGs scored by start context".
- **Inputs:** `data/bench_v3/inputs/<org>/complete.{fa,truth.tsv}`. **[work]** `work/release_pin/bench/<org>/complete/rel/input.fa.transdecoder_dir/start_refinement.kozak.tsv` (E40 stage A).
- **Software:** `pytransdecoder.core.kozak` from `work/trees/PY_rel_d3892b5` via `sys.path.insert` (0.4.0 has the same functions); `$P/env/td`.
- **Steps:** `python analyses/start_param_sensitivity/later_atg_kozak_v1.py > data/start_refinement/later_atg_kozak.tsv`.
- **Outputs:** `data/start_refinement/later_atg_kozak.tsv` -> notebook 03, S6b.
- **Cost:** about 2 min, one core (after E40 `rel` runs on the 8 complete inputs).

### E29 5'-open starts: error anatomy, real truncation, UTR-length prior (status: anatomy and truncation share in paper; the prior was not adopted, and the paper says so; reproducible with latest PyTD: anatomy yes; the prior needs commit c9d7875)
- **Feeds:** Results "Truncated ORFs and short UTRs are hard to tell apart" (about 5% left open, distance error rates,
  annotation UTR checks, pi 0.20-0.37 vs 0.76, prior not adopted); Discussion Limitations; notebook 03 (31 checks).
- **Inputs:**
  - release runs **[work]** `work/bench_v3/runs/<org>/<v>/pytd_v7` (stranded complete, complete_fullutr, fragments,
    fragments_mixed_t1-3) and `work/gm12878_td/<asm>/final_v7_abinitio`. They are byte-identical to E40 `rel`, so they can be
    replaced by `work/release_pin/{bench/<org>/<v>/rel, gm/<asm>/rel}`;
  - `data/annot/<org>/` (`benchmark_v2/truth_sources.tsv`), MANE v1.5 GTF, GM12878 assemblies/structures.
- **Software:** PyTD 3a793ff (`work/trees/PY_v7`) for `03_utr_features.py`, imported via PYTHONPATH. Use v0.2.1.
  The prior runs use branch `exp/utr-prior` **c9d7875**, worktree `work/trees/PY_utr`. **Local only: it exists only in
  `$P/repos/PyTransDecoder` and is not on GitHub.**
- **Steps (`pipeline/derive/`):**
  - `03_utr_features.py <run> <gcode> <out> [truth]` (`PYTHONPATH=<tree>`);
  - `03_utr_anatomy.py` -> `release.anatomy/decision_by_d` and `release_gm12878.*`;
  - `03_utr_annotation.py <org> <prefix>` and `03_utr_prior_table.py` -> `annotation_utr5_summary`, `trunc_d_summary`, `utr_prior_lr`
    (intermediates in **[work]** `work/utr_prior/annot/`);
  - `03_utr_gm12878_truth.py <structures> <fa> <mane> <prefix>` -> **[work]** `work/utr_prior/gm12878/<asm>.{truth,pairs}.tsv`.
    `03_utr_anatomy.py` and `03_utr_prior_sim.py` read them from that path, and E34 uses them too. Then run
    `03_utr_gm12878_summary.py` -> `gm12878_five_prime_truncation.tsv`, `missing_frac_hist.tsv`;
  - prior runs: `03_utr_runs.sh bench <org> <v> <lambda>` / `gm <asm> <lambda>` -> **[work]** `work/utr_prior/{runs,gm12878_runs,gm12878_eval}`;
    `03_utr_eval.py` -> `utr_prior_eval.{bench,bench_types,gm12878}.tsv`;
  - `03_utr_prior_sim.py` -> `offline_sim.release*.tsv`.
- **Outputs:** `data/start_refinement/utr_prior/*.tsv` (15 files) -> notebook 03, Results, Discussion.
- **Cost:** not recorded. The prior runs are 24 benchmark plus 24 mixed plus 4 GM Predict runs (4.1 GB + 0.5 GB under `work/utr_prior`).

### E33 Keep-open vs trim decision: classifier, homology assist, hit depth, UniRef50 (status: steps 1-2, 4 and 5 not adopted; step 3, the homology start assist, ADOPTED in PyTD (a5b24c3, in d3892b5/0.4.0) and in the paper; reproducible with latest PyTD: assist runs yes; offline steps need v0.2.1)
- **Feeds:**
  - in the paper: the homology rows (`pytd_v7h_homology`, 84 inputs; in `data/supplementary_tables/benchmark_per_input.tsv`),
    homology-mode text (5'-partial +11.4 points, retention-only comparison), GM12878 Trinity homology (7,063 MANE
    exact retention-only), and the simulated depth series A/B/C for Arabidopsis and human (`athal_depth_scores.tsv`, `human_depth_scores.tsv`);
  - not in the paper: steps 1-2, 4 and 5 (EXPERIMENT_LOG).
- **Inputs:**
  - release runs **[work]** `work/bench_v3/runs/<org>/<inp>/pytd_v7` (48);
  - simulated assemblies: `work/simtrinity/{runs,depth/<d>/runs}/pytd_v7` and `work/simtrinity_athal/depth/<d>/runs/pytd_v7`
    **[work]**; inputs and labels are in `$P/data/simtrinity{,_athal}`;
  - hits **[work]**: benchmark `pytd_v7_homology`; GM `trinity/final_v7_homology`;
    `work/start_classifier/homology/<ds>` (step-3 candidate searches); `work/start_classifier/{athal,human}_depth/searches/<d>` (full searches).
  - Kept in `$P/data`: step-4 deep searches `data/start_classifier/step4/searches`, step-5 hits, queries and features
    `data/start_classifier/step5/`, features and selftrain `data/start_classifier/{features,selftrain}`. The UniRef50
    2026_03 DIAMOND db (15.4 GB) was already deleted on 2026-10-29 (`notes/DISK_CLEANUP_2026-10-29.md`).
    Rebuild with `step5_01_uniref50_build_v3.sh`; 2026_03 may need the UniProt `previous_releases` area.
- **Software:**
  - PyTD 3a793ff (`PY_v7`, steps 1-2) and a5b24c3 (`PY_hstart`, assist runs, offline steps 3-5);
  - `$P/env/startclf` (Python 3.11.18, scikit-learn 1.9.1, pandas 3.1.0, numpy 2.4.6);
  - DIAMOND 2.2.8 and HMMER 3.4 (`$P/env/td`).
- **Steps (`pipeline/start_classifier/`):**
  - Step 1: `01_features_all_v1.sh`, `02_separation.py`, `03_ablation.py`.
  - Step 2: `04_selftrain_all_v1.sh`, `05_evaluate.py`, `06_summarize.py`.
  - Step 3: `07_homology_all_v2.sh` (DIAMOND `-k 1 1e-5` vs `sprot_minus_<taxid>`, hmmsearch `--cut_ga`, 5'-open candidates only), then `08_homology_rules.py`, `05_evaluate.py`.
  - Adoption runs: `15_eval_jobs_v1.sh` runs `11_run_pytd_v7h_homology_v1.sh <run> <fa> <gcode> <stranded> <hits_dir> [--no-homology-start-assist]`
    (84 benchmark inputs with `pytd_v7_homology` hits, writing **[work]** `work/bench_v3/runs/*/*/pytd_v7h_homology`;
    4 sims x on/off, writing `work/start_classifier/hstart_runs`) and `13_gm12878_trinity_v7h_v1.sh`
    (`work/gm12878_td/trinity/final_v7h_homology`). Score with `16_score_adoption.py`, `14_gm12878_e29_stat.py`,
    `12_check_inprogram_vs_offline.py` and `09_pytd_hstart_abinitio_regression_v1.sh`.
  - Depth series: `17_athal_depth_level_v1.sh` and `18_..._lane_v1.sh` -> `19_athal_depth_summary.py`. Human:
    `20_human_depth_level_v1.sh` and `21_human_depth_lanes_v1.sh` -> `22_human_depth_summary.py`. Each level gets one
    LongOrfs and conditions A (ab initio), B (hits plus assist) and C (`--no-homology-start-assist`), written to
    **[work]** `work/start_classifier/{athal,human}_depth`.
  - Step 4: `step4_01_deep_search_v3.sh` (DIAMOND `-k 25 --max-hsps 0 --evalue 10`, fast and `--more-sensitive`;
    hmmsearch `--domE 1 -Z 1 --domZ 1`), then `step4_02..06`.
  - Step 5: `step5_00..13`; `step5_08_full_predict_all_v{1,2}.sh` writes **[work]** `work/start_classifier/step5_full`.
  - With 0.4.0: set `HS=` (11_) or `PYTHONPATH=` (13_, 17_, 20_) to a v0.4.0 checkout, or use the Option B edits.
    The 84 `pytd_v7h_homology` runs and the GM `final_v7h_homology` run were not rerun under `release_pin`. They are
    a5b24c3, the same tree as d3892b5.
- **Outputs:** `data/start_classifier/{step1,step2,step3,step4,step5,adoption}/`. Paper consumers: `adoption/*` and
  `benchmark_per_input.tsv` (via `09_supplementary_tables.py`).
- **Cost:**
  - step 1 features about 25 min serial;
  - adoption jobs 3.3 h at 3 parallel (92 jobs);
  - Arabidopsis lane 2.6 h; human depth 2.4 h at 3 lanes (hmmsearch about 1 h 40 min at d1);
  - step 4 searches about 2.7 h;
  - step 5: UniRef50 build (streamed download, hours) plus searches about 2 h plus full predict 0.6 h.

### E34 Anatomy of missed CDSs and the 5' UTR-vs-coding decision (status: diagnostic, not in paper; reproducible with latest PyTD: fresh runs yes, offline trace needs v0.2.1)
- **Feeds:** internal (motivated E35); none in the manuscript.
- **Inputs:**
  - **[work]** runs: `work/simtrinity/runs/pytd_v7`, `work/simtrinity/depth/{d8,d32}/runs/pytd_v7`,
    `work/start_classifier/athal_depth/{A,B}/depth/{d1,d16}/runs/pytd_v7` (B hits in `.../athal_depth/searches/<d>`),
    and `work/gm12878_td/trinity/final_v7_abinitio`. The GM run is symlinked from `data/missed_cds/gm12878_inputs/run/`;
    those symlinks break when `work` is deleted.
  - GM truth: copied to `data/missed_cds/gm12878_inputs/labels/mane.truth.tsv` from `work/utr_prior/gm12878/trinity.truth.tsv`.
- **Software:** PyTD a5b24c3 (`PY_hstart`) imported; `$P/env/startclf`; `$P/env/td` for the fresh runs.
- **Steps (`pipeline/missed_cds/`):**
  - `01_trace_missed_v1.py <manifest> <ds> $P/data/missed_cds` for 7 sims plus GM;
  - `03_fresh_lanes_v1.sh` and `02_fresh_runs_v2.sh` (minlen30: LongOrfs `-m 30`; sensitive: Predict `--sensitive`)
    -> **[work]** `work/missed_cds/<mode>/<ds>/runs/pytd_v7`;
  - `01` again on minlen30, then `07_check_sensitive_v1.py`, `04_summarize_v1.py`, `05_source_annotation_v1.py`, `06_utr_coding_check_v1.py`.
  - The commands are in E34 "Commands". The scripts use `nice -n 10`; drop it per the project rule.
- **Outputs:** `data/missed_cds/*.tsv` (committed), per-CDS tables in `$P/data/missed_cds/<ds>/` (kept).
- **Cost:** not recorded (1.2 GB of fresh runs).

### E35 One-sided Kozak term ("Kozak floor"), offline (status: not adopted, not in paper; reproducible with latest PyTD: needs v0.2.1 / a5b24c3 internals)
- **Feeds:** internal decision record only.
- **Inputs:** 90 finished runs listed in `data/kozak_floor/manifest.tsv` (`00_manifest_v1.sh` checks every path), all **[work]**:
  - `work/bench_v3/runs/<org>/<inp>/{pytd_v7,pytd_v7h_homology}`;
  - `work/simtrinity/{runs,depth/<d>/runs}/pytd_v7`;
  - `work/start_classifier/hstart_runs/sim.human.<d>/on`;
  - `work/start_classifier/athal_depth/{A,B}/depth/<d>/runs/pytd_v7`.
  - Truth in `$P/data/{bench_v3,simtrinity,simtrinity_athal}`.
- **Software:** PyTD a5b24c3 (`PY_hstart`, `homology_start_call`, `_revise_start_sites` logic); `$P/env/startclf`.
- **Steps:** `00_manifest_v1.sh`, then `01_candidates_all_v1.sh` (4 lanes) -> `$P/data/kozak_floor/candidates`, then
  `02_evaluate_v1.py`, `03_check_official_v1.py`, `04_summarize_v1.py`, `05_changed_profile_v1.py` (E35 "Commands").
- **Outputs:** `data/kozak_floor/*.tsv` (12 tables); `$P/data/kozak_floor/{candidates,eval}` (46 MB, kept).
- **Cost:** about 3 min extraction plus minutes of evaluation. Only the input runs are expensive.

### E36 Per-transcript alpha from ORF/transcript length, offline (status: not adopted, deferred by BH 2026-10-30; not in paper; reproducible with latest PyTD: needs v0.2.1 internals)
- **Feeds:** internal; a possible future rule, alpha = 0.5 log10(O/200 nt).
- **Inputs:** E35's 90 runs and candidate tables, plus 5 new runs, all **[work]**:
  - `work/simtrinity/depth/{d2,d4,d16}/runs/pytd_v7`;
  - GM Trinity `final_v7_abinitio` and `final_v7h_homology`, symlinked from `data/{missed_cds,alpha_by_length}/gm12878_inputs`.
  - E29/E34 GM truth.
- **Software:** PyTD a5b24c3 via E35's `01_candidates_v1.py`; `$P/env/startclf`.
- **Steps (`pipeline/alpha_by_length/`):** `00_manifest_v1.sh`, `01_candidates_new_v1.sh`, `02_labels_v1.py`,
  `03_check_v1.py`, `04_descriptive_v1.py`, `05_tune_grid_all_v1.sh`, `06_select_v1.py`, `07_evaluate_v1.py`,
  `08_ceiling_v1.py`, `09_summarize_v1.py`, `10_figure_v1.py`. Addendum: `05_tune_grid_all_v2.sh`, `07_evaluate_v2.py`.
- **Outputs:** `data/alpha_by_length/*.tsv`, `simpler_v2/`, figure. Raw `$P/data/alpha_by_length/{candidates,labels,grid*,eval*}` (about 60 MB, kept).
- **Cost:** about 15 min in total (grid about 10 min at 3 lanes).

### E42 Parallel and faster PyTD (status: engineering, not a paper claim; it produced v0.3.0/v0.4.0; reproducible: yes, `pip install pytransdecoder==0.4.0`)
- **Feeds:** the claim that 0.4.0 equals d3892b5 byte for byte, which underpins this document. Not cited in the paper.
- **Inputs:** E40 release runs as references **[work]** `work/release_pin/{bench,gm,presets}`;
  `work/bench_v3/runs/*/*/pytd_v7h_homology`; `work/bench_v3/presets/ref`.
- **Software:** branches dev/parallel-scoring 98d77fb, dev/faster-serial 5859887, devel 725244f, all merged into
  v0.3.0 and v0.4.0. Their worktrees `PY_par`, `PY_fast` and `PY_lo` were deleted. Rerun with a tag checkout and `PYTD_FAST=` or `PYTHONPATH=`.
- **Steps (`pipeline/parallel_scoring/`):** `launch_v1.sh` and `run_one_v1.sh`, then `summarise_v1.py`; v2 the same with
  `launch_v2.sh`, `run_one_v2.sh`, `summarise_v2.py`; LongOrfs with `launch_lo_v3.sh`, `run_lo_v3.sh`, `summarise_lo_v3.py`.
- **Outputs:** `data/parallel_scoring/{v1,v2,lo_v3}.{parity,timing}.tsv`. The run outputs were already deleted
  (2026-11-03), so parity can only be recomputed by rerunning, and that needs the E40 release runs as references.
- **Cost:** about 156 + 226 runs of minutes each; timings are indicative (shared host).

### E14 / E38 (literature only)
- E14: 15 references added 2026-10-22, checked in paperclip and Crossref. No compute.
- E38: sources and passages for the Background, Abstract and Conclusions. No compute.

---

# PART 2. Things in `$P/work` that are inputs or expensive to regenerate (these experiments)

Sizes are `du -sh` on 2026-11-03. "Regen" = the cost to recreate with 0.4.0 if deleted.

| path | size | what | needed by | copy elsewhere? | recommendation |
|---|---|---|---|---|---|
| `work/bench_v3/runs/*/*/pytd_v7_homology/{blastp.outfmt6,pfam.domtblout}` (84 dirs) | 263 MB | DIAMOND vs sprot_minus and Pfam GA hits on the v7 LongOrfs | E40 `retonly`, E33 `pytd_v7h_homology`, E42, E33 step 5 | no | **Preserve.** Copy to `$P/data/homology_hits/bench/<org>/<inp>/`, plus `longest_orfs.pep` (or md5) for the `cmp` check. Regen is many CPU-hours of hmmsearch, and a newer DB release would not reproduce them |
| `work/gm12878_td/trinity/final_v7_homology/{blastp.outfmt6,pfam.domtblout}` | 18 MB | GM Trinity hits | E40 `gm trinity retonly`, E33 `final_v7h_homology`, E36 | no | **Preserve** to `$P/data/homology_hits/gm12878/trinity/` |
| `work/start_classifier/{athal,human}_depth/searches/` | 103 + 26 MB | full-mode hits per simulated depth level (Swiss-Prot/Pfam) | E33 depth series B/C (paper rows), E34, E35 | no | **Preserve** to `$P/data/homology_hits/sim_{athal,human}/`. Regen: about 1.7 h of hmmsearch per level at d1 |
| `work/start_classifier/homology/<ds>/` | 61 MB | step-3 candidate-only searches (4 sims) | E33 hstart_runs (adoption sim rows), E35 | no | Preserve (small) |
| `work/bench_v3/presets/{ref,tests,jobs.txt,jobs2.txt}` | 132 + 102 MB | E21 held-out halves and test subsets (with truth) | E40 stage D, `15_score_presets_v2.sh`, E42 | no (see E21 for builder) | **Preserve** to `$P/data/bench_v3/presets/`, and fix the `OLD=`/`TESTS=` paths |
| `work/bench_v3/rules/*/*/cliff_c200/preds.tsv` | 6.6 MB | fixed short-ORF test, unreleased dev code (historical comparator) | `11_score_sweeps_v1.sh` (negatives, sens_by_length); paper Fig. 4, S11 | no; cannot be regenerated with any release | **Preserve** (preds.tsv only) |
| `work/utr_prior/gm12878/*.{truth,pairs}.tsv*` | 2.9 MB | E29 transcript-level MANE truth for the 3 assemblies | E29 anatomy and offline sim, E34, E36 | trinity truth only (`data/missed_cds/gm12878_inputs/labels`) | **Preserve** (tiny; regen needs `03_utr_gm12878_truth.py`) |
| `work/utr_prior/annot/` | 5.6 MB | per-organism UTR and truncation tables | E29 `03_utr_prior_table.py` | no | Preserve (tiny) |
| `work/release_pin/` | 80 GB | all E40 runs (bench 37, start_grid 34, gm 6.1, presets 3.1, gm_scoring 0.1, regression 0.07) | E40 scoring 11-15, summarise_v2, E41, kozak_score_distribution_v2, E42 parity | committed tables only | Regen about 32 h (less with 0.4.0 `--cpu`). Optional: tar the small per-run files (`preds.tsv`, `score.tsv`, `eval.*`, `ribo.tsv`, `runtime_sec.txt`, logs, `start_refinement.kozak.tsv`, `pytd_version.txt`; likely < 2 GB) so tables can be re-derived without rerunning. Otherwise OK to delete |
| `work/release_pin/presets/presets/<name>/` | small | 8 built presets | stage D tests, E42 | only `humanmix/preset.json` in `data/benchmark/presets_release/preset_meta` | Optional preserve (deterministic regen in minutes) |
| `work/gm12878_td/trinity/{final_v7_abinitio,final_v7h_homology}` | 348 + 377 MB | GM runs symlinked from `data/{missed_cds,alpha_by_length}/gm12878_inputs` | E34, E36 | no | Symlinks break. Regen equals E40 `gm trinity rel` (abinitio) or `13_gm12878_trinity_v7h_v1.sh`. OK to delete if the hits above are kept |
| `work/bench_v3/runs/*/*/{pytd_v7,pytd_v7h_homology}` (rest of the 161 GB tree) | part of 161 GB | release runs read by E29/E33/E35/E36 | offline analyses | byte-identical to `release_pin/bench/*/rel` (pytd_v7); v7h needs a rerun | Delete; repoint manifests to `release_pin` or reruns |
| `work/simtrinity/`, `work/simtrinity_athal/` | 2.3 + 1.3 GB | pytd_v7 runs on the simulated assemblies (inputs are in `$P/data`) | E33, E34, E35, E36 (owned by E30-E32) | no | Delete OK (rerun minutes each). See the E30-E32 section |
| `work/start_classifier/{athal_depth,human_depth,hstart_runs,step5_full,hstart_regression}` | 1.3 + 2.1 GB + 688 + 1,200 + 161 MB | assist runs (A/B/C) | E33 tables (committed), E34-E36 | committed scores only | Delete after copying `searches/` (above). Regen is Predict-only, minutes per run |
| `work/missed_cds/`, `work/utr_prior/{runs,gm12878_runs}` | 1.2 GB, 4.1 + 0.5 GB | fresh minlen30/sensitive runs; UTR-prior runs (c9d7875) | E34; E29 (`utr_prior_eval`) | committed tables | Delete OK. UTR-prior regen needs c9d7875 (see Part 3) |
| `work/trees/PY_rel_d3892b5`, `TD_rel_41fc5ee`, `PY_hstart`, `PY_v7`, `PY_utr` | about 105 MB each (TD 37) | git worktrees | all scripts (hardcoded) | git objects in `$P/repos/{PyTransDecoder,TransDecoder}` | Delete OK **only if `$P/repos` is kept** (worktree metadata lives there; `git worktree prune` afterwards) |
| `work/parallel_scoring`, `fast_scoring`, `longorfs_cpu` | 5.5, 35, 38 MB | logs and timings only (outputs already deleted) | none (tables committed) | committed tsvs | Delete OK |
| DBs: Swiss-Prot, `sprot_minus_*.dmnd`, Pfam-A 38.2 | 3.8 GB | in `$P/data/homology` | all homology runs | n/a (not in work) | Safe (data kept). UniRef50 db already deleted |

---

# PART 3. Gaps

1. **`exp/utr-prior` c9d7875 (E29 prior) is local only**, in `$P/repos/PyTransDecoder`. The other unmerged
   experimental branches there, `exp/contrast-frames` (E26) and `exp/short-orf-taper`, are local too. `$P/repos` is not
   on the "kept" list (`data`, `env`). Push these branches to GitHub (or `git bundle` them into `$P/data`) before any cleanup
   of `$P/repos`. Deleting `$P/work/trees/*` alone is safe.
2. **`sprot_minus_<taxid>.dmnd` builder not saved** (E07 gap). The reconstruction in E40 above is verified only by counts
   (E33 step 5). The `.dmnd` files are kept in `$P/data/homology`, so this matters only for a from-scratch rebuild.
3. **Homology hit tables live only in `$P/work`.** Without them the paper's retention-only and `pytd_v7h_homology`
   numbers can be regenerated only approximately, if Swiss-Prot/Pfam/DIAMOND versions differ. Copy them first (Part 2).
4. **Run directories of 3a793ff/a5b24c3 are referenced by absolute `$P/work` paths** in the committed manifests
   (`data/kozak_floor/manifest.tsv`, `data/alpha_by_length/manifest.tsv`, `data/missed_cds/manifest*.tsv`,
   `data/start_classifier/step{2,3}_manifest.tsv`). Symlinks under `data/{missed_cds,alpha_by_length}/gm12878_inputs` point
   into `work/gm12878_td`. After deletion these analyses need reruns and manifests repointed, e.g. to `release_pin/*/rel`
   (byte-identical ab initio). The `pytd_v7h_homology` and depth A/B/C runs have no `release_pin` twin, so rerun them with
   `11_`/`17_`/`20_` and the preserved hits.
5. **Offline re-derivation scripts are untested on 0.3.0/0.4.0.** They call private `TransDecoderPredict` methods and
   `predict.py` changed. Use `pytransdecoder==0.2.1` (code identical to d3892b5) for E29 features and E33-E36.
6. **`10_byte_identity_v1.py` and the E42 parity summaries cannot be recomputed** after `work` is deleted. Their reference
   runs (`bench_v3/runs/*/pytd_v7*`, `gm12878_td/*/final_v7_*`, the deleted E42 outputs) are gone. The committed tsvs stand as the record.
7. **E40 "run-time ratio" (2.7x) comes from E25 timings**, not re-timed. 0.4.0 timings differ (faster), so rerunning
   the timing benchmark with 0.4.0 would not reproduce 2.71. State the version in REPRODUCE.md.
8. **Hardcoded paths:** `run_grid_v2.py` (PYTD); `later_atg_kozak_v1.py` and `kozak_score_distribution_v2.py` (`sys.path`);
   `launch_v1.sh` (regression trees, `work/bench_v3/runs` glob for `retonly`, `work/bench_v3/presets/jobs*.txt`);
   `run_preset_v1.sh` (`OLD=$P/work/bench_v3/presets`); `15_score_presets_v2.sh` (`TESTS=`);
   `run_bench_v1.sh`/`run_gm_v1.sh` (`retonly` hit dirs). All need repointing to the preserved copies.
9. **Cost not recorded** for E29 and E34 runs.
10. **`nice` in E33-E36 scripts** (`nice -n 10/19`) conflicts with the current "do not nice" rule. Drop it on rerun; outputs are unaffected.
11. **UniRef50 2026_03 DB already deleted** (E33 step 5, a negative result). A rebuild depends on UniProt still serving
    2026_03 (`previous_releases`) and takes hours plus 15 GB.
12. **GM12878 LRAA/StringTie have no v7 homology searches** (E33). Only Trinity has `retonly` and `v7h` runs. This is consistent with the paper but worth stating.
