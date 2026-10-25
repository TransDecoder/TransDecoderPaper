# E09. Start-codon refinement and Kozak context

## Claims supported
Manuscript: Methods, "Start-codon refinement" (`methods_algorithm.tex`, `sec:start`); Results on start accuracy (benchmark, record E07). The Methods text still describes the old PWM procedure; see Caveats.

1. **The v6.0.0 start refinement had a scale bug.** The threshold was learned on length-normalised PWM scores but applied to raw sums, on a window shifted 1 nt upstream. Both implementations were fixed (Perl 8cf8d5f, PyTD 0525f5a).
2. **The PWM refinement, even when fixed, was unstable and added little.**
   - The chosen PWM sub-window and the resulting starts depend on the training RNG seed.
   - The PWM learned mainly the coding/non-coding transition, not the start context.
   - Choosing among in-frame ATGs by the PWM did not beat "first in-frame ATG" in 7 of 8 organisms. *C. reinhardtii* was the exception.
3. **The adopted rule is the final default.** An open 5′ end of length L is trimmed to its first in-frame ATG when the segment's frame contrast is below (α·L + w·K) × the training ORFs' median per-base contrast. Here α = 0.4, w = 25, and K is the Kozak log-odds of that ATG.
   - The Kozak model covers positions −6..−1 and +4..+5. It compares training-ORF starts against in-frame ATGs inside the training ORFs, with +1 pseudocounts.
   - Commits: PyTD efc1450 → 0886981 → f5e6a39 → eaea367; Perl f7e75a1 → 7e30234 → 313be29.
4. **Exact-start accuracy rises substantially on complete transcripts** (8-organism benchmark).
5. **Self-trained Kozak models agree with annotation-trained ones** (Pearson r 0.85–0.98 over 8 positions × 4 bases). They recover the known −3 A/G preference in all eight organisms.

## Data
- **Benchmark v3:** inputs and truth in `$P/data/bench_v3/inputs/<org>/<variant>.{fa,truth.tsv}`; record E07.
- **Benchmark v2 (superseded RefSeq truth):** `$P/data/bench_v2/inputs`. Used only by the strategy grid and the early fix comparison.
- **Regression sample data** (PASA, StringTie) from the TransDecoder `sample_data`, for the seed runs. Inputs are under `$P/work/regression/*/inputs`.
- **GM12878 LRAA and StringTie assemblies, labelled by MANE Select v1.5:** record E03 and `$P/work/gm12878_td/`.

## Software
- **PyTransDecoder:**
  - `PY_fix` 27bfb5e: `pytd`, `pytd_norefine` (fixed PWM).
  - `PY_dev` (branch `dev/frame-contrast-starts`): the development labels below. Final: 48c870c.
- **Perl TransDecoder:**
  - `TD_base` 66d4712 (v6.0.0; `perl_td_v600`).
  - `TD_fix` (branch `fix/start-refinement`): final de3f57e.
- **Environment:** `$P/env/td` (Python 3.11.17, Biopython 1.88, matplotlib 3.10.9; perl 5.32.1).
- **Repo commits:** this repo's scripts are at 4090a7e (start_labels/start_features, start_kozak_form), 501bbae (figures/kozak), 7a7b10e (start_codon_dev) and cbfd0d0..7750fc7 (start_refinement_stability).

Start-refinement history, with the benchmark label and when each was run (from `preds.tsv` mtimes):

| Benchmark label | Code | Rule | Run |
|---|---|---|---|
| `perl_td_v600` | TD 66d4712 | v6.0.0 PWM (scale bug) | 2026-10-18 |
| `pytd_norefine` | PyTD 27bfb5e `--no-refine-starts` | none | 2026-10-18 |
| `pytd` | PyTD 27bfb5e | fixed PWM (0525f5a); multi-ORF output | 2026-10-18 |
| `pytd_dev` | PyTD efc1450 | trim open 5′ end to first ATG by frame contrast (α = 0.4) | 2026-10-19 15:09–15:51 |
| `pytd_dev_pwm` | PyTD efc1450 `--start-atg-choice pwm` | as above, but the ATG is chosen by the PWM | 2026-10-19 15:09–15:51 |
| `pytd_dev_single` | PyTD efc1450 `--single-best-only` | as `pytd_dev`, one ORF per transcript | 2026-10-19 16:24–16:51 |
| `pytd_dev_kozak` | PyTD 0886981 (+8412016 single-best default) | multiplicative Kozak: (α + 0.15·K)·L | 2026-10-19 17:41–18:03 |
| `pytd_dev_kozak_add` | PyTD f5e6a39 | additive Kozak: α·L + 25·K | 2026-10-20 07:43–08:44 |
| `pytd_dev_v2` / `perl_td_dev_v2` | PyTD 48c870c / TD de3f57e | final defaults: additive Kozak, plus `-m 50` with short-ORF guard and GC/tie fixes | 2026-10-20 12:41–15:32 |

## Commands
Benchmark runs: `$T/benchmark_v2/scripts/run_tool.sh`. The labels are mapped in its `pytd*` case; the launchers are `launch_v3/launch.sh`, `launch2.sh`, `launch5.sh`, `launch9.sh`, `launch13.sh` and `launch16.sh`. Scoring is `score_all.py` (record E07).

Seed sensitivity of the PWM (PyTD 27bfb5e-era code; seeds 1–5):
```bash
bash $T/analyses/start_refinement_stability/run_seeds.sh <PY tree> <regression inputs> $P/work/start_seed_runs "pasa stringtie" "1 2 3 4 5"
python $T/regression/compare_orfs.py $ds/seed1/$ds.fa.transdecoder.gff3 $ds/seed$s/$ds.fa.transdecoder.gff3
```

PWM strategy grid (legacy / kfold / fixed; bench_v2 `pytd_norefine` runs):
```bash
bash $T/analyses/start_refinement_stability/run_strategy_grid.sh <org> <variant> $P/work/start_grid/<org>.<variant>.tsv
```

PWM audit, start-error classification and offline rule simulations. These are read-only with respect to the code trees; outputs are in `$P/work/analysis/{pwm_audit,start_errors,kozak}`.
- `$T/analyses/start_codon_dev/pwm_audit/{audit_signal.py,simulate_rules.py,tabulate.py} <org> <variant>`
- `$T/analyses/start_codon_dev/start_errors/{classify_starts.py,upstream_segment_features.py,simulate_trim.py,simulate_advance.py,nc_model_eval.py}`
- `$T/analyses/start_codon_dev/kozak/sim_kozak.py` (ATG choice by Kozak PWM)
- `trim_kozak.py <org> > trim.<org>.tsv` (α, β grid for the multiplicative form)
- `trim_features.py <org>` (per-row r, K, L for comparing rule forms; same file as `analyses/start_kozak_form/trim_features.py`)

GM12878 start labels:
- Run PyTD with every open 5′ end trimmed (`trim_all`; the predict log shows α = 1000).
- `python $T/genome_annotation/scripts/start_labels.py $P/work/gm12878_td/<asm>/trim_all <labels.tsv>`
- `python $T/genome_annotation/scripts/start_features.py <run> <labels.tsv> $P/work/gm12878_td/eval/<asm>.start_features.tsv`

Kozak figure (`PYTHONPATH` = a PyTD tree with `core/kozak.py`, ≥ fdc6bbe). It uses the training sets of the `perl_td` runs on `complete`.
```bash
python $T/figures/kozak/make_kozak_figure.py $P/data/bench_v3/inputs $P/work/bench_v3/runs $T/figures/kozak
```

## Outputs
- `$T/data/benchmark/summary.tsv` (exact-level columns `tp_exact`, `f1_exact`).
- `$T/figures/kozak/kozak_logos{,_bits,_self_vs_truth}.{pdf,png}` and `kozak_stats.tsv`.
- Not in git:
  - `$P/work/start_seed_runs/{pasa,stringtie}/seed{1..5}`
  - `$P/work/start_grid/*.tsv`
  - `$P/work/analysis/pwm_audit/{FINDINGS.md,out/}`
  - `$P/work/analysis/fork_findings.md`
  - `$P/work/analysis/kozak/{sim_*.tsv,trim.*.tsv,trim_features.*.tsv}`
  - `$P/work/analysis/start_errors/`
  - `$P/work/gm12878_td/eval/{lraa,stringtie}.start_features.tsv`
- Per-run Kozak models: `<run>/input.fa.transdecoder_dir/start_refinement.kozak{,.counts}.tsv` and `.kozak.logo.svg`.

## Results
**Benchmark exact-start TP, summed over the 8 organisms** (`summary.tsv`). The stop-level TP is shown for reference.

| Label | complete exact | complete_fullutr exact | fragments exact | complete stop TP | complete exact F1, median [range] |
|---|---|---|---|---|---|
| pytd_norefine | 83,307 | 84,895 | 83,642 | 128,376 | 0.687 [0.455–0.897] |
| pytd (fixed PWM) | 90,163 | 91,501 | 84,032 | 128,376 | 0.728 [0.484–0.906] |
| perl_td_v600 (bugged PWM) | 96,957 | 97,933 | 83,691 | 128,562 | 0.779 [0.567–0.912] |
| pytd_dev (frame contrast, first ATG) | 114,400 | 114,165 | 85,103 | 128,017 | 0.893 [0.740–0.928] |
| pytd_dev_pwm | 113,742 | 113,357 | 84,357 | 128,017 | 0.887 [0.744–0.928] |
| pytd_dev_kozak (multiplicative) | 115,106 | 114,480 | 86,678 | 128,017 | 0.903 [0.737–0.929] |
| pytd_dev_kozak_add | 115,406 | 114,816 | 88,434 | 128,017 | 0.904 [0.735–0.927] |
| pytd_dev_v2 (final) | 116,384 | 115,781 | 95,520 | 129,058 | 0.910 [0.742–0.928] |

- **Per organism, complete, `pytd_norefine` → `pytd_dev_v2`:**

  | Organism | Exact TP |
  |---|---|
  | human | 10,590 → 17,514 |
  | mouse | 14,691 → 19,223 |
  | athal | 17,707 → 23,518 |
  | dmel | 10,055 → 12,481 |
  | creinhardtii | 10,324 → 12,702 |
  | spombe | 3,537 → 4,568 |
  | pfalciparum | 4,659 → 4,830 |
  | tthermophila | 11,744 → 21,548 |

- **PWM choice vs first ATG** (`pytd_dev_pwm` vs `pytd_dev`, complete): lower in 6 organisms, higher in creinhardtii (12,661 vs 12,592), equal within 1 in pfalciparum (4,824 vs 4,823).
- **Additive vs multiplicative Kozak** (`kozak_add` − `kozak`), exact F1 per organism:
  - fragments: +0.005 to +0.023;
  - complete: −0.006 to +0.015, mean +0.001.
- **Fragments, `pytd_dev_v2`:** the gain on fragments also includes `-m 50` (stop TP 102,662 → 110,333), so it is not attributable to the start rule alone.
- **Fixed PWM vs v6.0.0** (`fork_findings.md`): the fix lowered complete-transcript exact accuracy. v6.0.0's scale bug acted as a permissive threshold.
  - Revised starts, complete (fix vs v6.0.0): human 3,540 vs 4,207; athal 2,241 vs 3,371; tthermophila 615 vs 4,593.

**PWM seed sensitivity** (`work/start_seed_runs`). Starts differing from seed 1; stop sets are identical in all cases.

| Dataset | ORFs | Seeds 2–5 |
|---|---|---|
| PASA | 792 | 53, 30, 61, 29 |
| StringTie | 3,037 | 144, 47, 111, 54 |

Best-AUC (u,d) window, PASA, seeds 1–5: (15,2), (20,2), (3,7), (20,9), (7,2).

**Strategy grid** (`work/start_grid`, bench_v2 truth). Exact TP:

| Dataset | none | legacy, seeds 1–5 | kfold f1_balanced | kfold f1_raw |
|---|---|---|---|---|
| athal complete | 21,934 | 22,384–22,815 | 23,444–23,530 | 22,270–22,412 |
| creinhardtii complete | 10,324 | 10,521–11,674 | 12,160–12,184 | 10,982–11,246 |
| athal fragments | 17,055 | 16,490–17,044 | 15,110–15,311 | 16,814–16,950 |

kfold picks (19,10) in athal and (20,9) in creinhardtii. This line was superseded by the frame-contrast rule and not adopted.

**Where the true start lies among in-frame ATGs** (`pwm_audit/FINDINGS.md` Q4; annotated 5′ UTRs, bench_v3):
- The true start is a later in-frame ATG in 0.7–5.0% of transcripts in seven organisms, and in 19.6% in creinhardtii.
- Open-ORF exact starts, "first" vs oracle (complete): human 6,249 / 6,366; creinhardtii 2,267 / 3,097.

**Multiplicative Kozak weight grid** (`work/analysis/kozak/trim.<org>.tsv`, α = 0.4). Correct keep/trim decisions gained by β = 0.15 over β = 0:

| Organism | Gain |
|---|---|
| human | +665 |
| athal | +565 |
| mouse | +485 |
| dmel | +264 |
| tthermophila | +250 |
| spombe | +75 |
| creinhardtii | +59 |
| pfalciparum | +51 |

**Kozak logos** (`figures/kozak/kozak_stats.tsv`; self-trained = starts of the 500 training ORFs with a usable context, 253–476 per organism):

| org | self starts | truth starts | r(self, truth) | −3 A/G log2-odds, self / truth | +4 G, self / truth |
|---|---|---|---|---|---|
| human | 325 | 18,816 | 0.921 | 0.54 / 0.51 | 0.51 / 0.57 |
| mouse | 362 | 20,575 | 0.964 | 0.53 / 0.52 | 0.48 / 0.55 |
| dmel | 465 | 13,367 | 0.975 | 0.76 / 0.83 | 0.08 / 0.02 |
| athal | 404 | 25,725 | 0.908 | 0.19 / 0.35 | 0.88 / 0.71 |
| spombe | 402 | 4,501 | 0.883 | 0.41 / 0.70 | 0.17 / 0.12 |
| tthermophila | 253 | 26,177 | 0.845 | 0.10 / 0.21 | −0.31 / 0.12 |
| pfalciparum | 476 | 4,406 | 0.958 | 0.39 / 0.60 | 0.44 / 0.46 |
| creinhardtii | 395 | 17,621 | 0.961 | 0.58 / 0.57 | 0.24 / 0.10 |

**GM12878 start labels** (5′-open ORFs after trim_all; `eval/<asm>.start_features.tsv`):

| Label | LRAA | StringTie |
|---|---|---|
| utr_first_atg | 4,874 | 3,742 |
| truncated | 1,712 | 1,103 |
| upstream_atg | 214 | 216 |
| start_downstream_in_tx | 1,931 | 1,127 |
| other | 401 | 327 |

## Checks
- `pytd_dev_v2` and `perl_td_dev_v2` give identical predictions on all 84 benchmark inputs (record E08).
  - eaea367 removed the last start difference: a Kozak-score rounding near-tie, 1 of 4,099 ORFs on S. pombe fragments (per the commit message).
- The kozak.py model files are byte-identical to Perl `train_kozak_model.pl` (per fdc6bbe's message; not re-verified here).
- The `pytd_dev` vs `pytd_dev_pwm` numbers quoted in `FINDINGS.md` match `summary.tsv` exactly.
- **The fixed-PWM drop was reproduced.** Reverting 8cf8d5f restores v6.0.0-level exact F1: spombe complete 0.764 vs v600 0.761, athal 0.757 vs 0.755 (`fork_findings.md`).
- The seed counts (PASA, StringTie) and the PASA window swings were re-derived from `work/start_seed_runs` for this record.
- The strategy-grid ranges were re-derived from `work/start_grid` for this record.

## Gaps
- **Trinity seed runs are missing.** The Trinity seed range in notes/03 (14–26 differing starts) has no run directory in `work/start_seed_runs`, which holds PASA and StringTie only. It cannot be verified.
- **The coding-gate result is in notes only.** notes/04 says the gate blocks 76% of correct moves in athal and 16% in creinhardtii. No saved table or command reproduces this.
- **Two choices rest only on commit messages:**
  - the choice of the additive form and of weight 25 (f5e6a39: "on GM12878 … trims more true UTR extensions and keeps more truly truncated ORFs open");
  - the multiplicative weight 0.15 and "a later ATG by start context lowered accuracy … except Chlamydomonas" (0886981).

  The feature tables exist (`work/analysis/kozak/trim_features.*.tsv`, `work/gm12878_td/eval/*.start_features.tsv`), but the grid or summary comparing rule forms, and the command that produced it, are not saved.
- **The trim_all GM12878 runs have no saved command.** α = 1000 is known only from `predict.log`.
- **Analysis findings are not committed.** `pwm_audit/FINDINGS.md`, `fork_findings.md`, `start_grid/` and the `kozak/` simulation tables are only in `$P/work`.
  - `fork_findings.md` says its v6.0.0-revert comparison used scripts in a session scratchpad (`fork_reg/cmp_starts.py`), which are not in this repo.
  - `start_grid/tthermophila.fragments.tsv` is incomplete (`.tmp`).
- **Code state for two runs is not exactly identifiable.**
  - `pytd_dev_kozak_add` runs started at 07:43, before f5e6a39 was committed (08:19).
  - The final-default runs (12:41–15:32) predate 48c870c / de3f57e (15:37–15:38).

  The code that ran was uncommitted working-tree state, presumed equal to those commits. A full rerun at the final commits would confirm it.
  Partial check (2026-10-21, E07): rerunning `pytd_dev_v2` and `perl_td_dev_v2` at the committed trees (PY_dev 48c870c, TD_fix de3f57e) on spombe/complete and creinhardtii/complete gave byte-identical `preds.tsv`. The other 82 inputs per tool were not rerun.

## Caveats
- **The manuscript Methods (`methods_algorithm.tex` §Start-codon refinement) still describes the PWM procedure**: noise depletion, (u,d) AUC search, F1 threshold, choice among ATGs within max(30 nt, 15%). It needs rewriting for frame contrast + additive Kozak.
  - §Implementations still says Perl/Python "differ only in start codons chosen by the randomised PWM training". This is no longer true; they are identical (E08).
- **The headline gain is confounded.** Most of the gain from `pytd_norefine` to `pytd_dev` comes from trimming open 5′ ends into UTRs (frame contrast). The Kozak term adds about 1% exact TP on complete transcripts and more on fragments.
- **The `pytd_dev_kozak` label also includes the single-best-ORF default** (8412016). The `pytd_dev_v2` label also includes `-m 50`, the short-ORF guard and the GC/tie fixes.
- **The strategy grid and the v2 fix comparison use bench_v2 (RefSeq) truth**, not bench_v3.
- **The Kozak figure is trained on `perl_td` (TD_fix, 2026-10-18) training sets**, not on final-default runs. The training sets do not depend on the start rule.
- **The `start_downstream_in_tx` label is large** (1,931 LRAA, 1,127 StringTie), although `start_labels.py` documents it as "should not happen". It has not been investigated.
- **C. reinhardtii is the one organism where later in-frame starts are common (19.6%).** No tested rule recovers most of them.

## Addendum (2026-10-22): units and wording changed
The rule is now written as two scores, with the Kozak score in bits and the weight at 17 per bit (25 per nat
before). Numbers in this record are per nat (multiply K by 1.443 for bits). Behaviour and the equivalence
check are in E13 (Addendum); Kozak score distributions in bits are in
`data/start_param_sensitivity/kozak_score_distribution.tsv` (`kozak_score_distribution.py`).
`genome_annotation/scripts/start_labels.py` reads both the old (upstream_frame_contrast, nats) and new
(upstream_coding_score, bits) GFF3 attributes.
