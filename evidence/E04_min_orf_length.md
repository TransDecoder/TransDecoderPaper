# E04. Minimum ORF length: benchmark sweep and GM12878 post-hoc filter

## Claims supported
Manuscript: Results, short ORFs and the minimum length; Discussion, short ORFs and the dark proteome.

1. **Lowering LongOrfs `-m` from 100 to 50 recovers short proteins but adds many false calls on lncRNAs.**
   - Human: 50–99 aa sensitivity 0.31 → 0.86, while lncRNAs with a call rise 6,608 → 21,416 (of 33,668), and precision falls 0.74 → 0.47.
   - At `-m 50` there are 14,880 lncRNA calls of 50–99 aa vs 430 correct calls of that length (~35:1).
2. **Random sequence is a much easier negative than lncRNA.** Calls per GC-matched random sequence at `-m 50` are 0.013–0.159 across 8 organisms.
3. **The short-ORF frame-contrast guard (`-m 50` + guard 200, the final default) keeps most of the short-protein gain at near `-m 100` false-call levels.**
   - Human complete: short sensitivity 0.31 → 0.56; lncRNA calls 6,608 → 6,943; precision 0.737 → 0.729.
4. **The guard also removes most Ribo-seq-supported noncanonical ORFs.** On GM12878, catalogue matches drop 458 → 178 (LRAA) and 578 → 217 (StringTie).
5. **GM12878, post-hoc length filter: 50 vs 100 aa barely changes precision.** The large step is between 100 and 150 aa, where unannotated predictions fall to 1–4% but MANE exact-given-assembled falls to ~0.85.

## Data
- **Benchmark v3 inputs:** `$P/data/bench_v3/inputs/<org>/<variant>.fa` with truth `.truth.tsv`. Truth sets and inputs are described in record E07.
  - **Variants:**
    - `complete.wRand`, `fragments.wRand`: true transcripts (complete or fragmented) plus i.i.d. GC-matched random sequences (`rand_*`);
    - `complete.wLnc`, `fragments.wLnc`: human and mouse true transcripts plus GENCODE lncRNAs (`lnc_*`).
  - **Random negatives per organism:** athal 26,867; creinhardtii 17,621; dmel 13,884; human 19,167; mouse 21,409; pfalciparum 5,281; spombe 5,100; tthermophila 26,456.
- **lncRNA negatives** (`$T/analyses/lnc_holdout/pick_lnc.py`): one transcript per gene (the longest, ≥ 300 nt) whose exons overlap no same-strand GENCODE CDS, renamed `lnc_<id>`.
  - Human: 33,668 GENCODE v50 lncRNAs (`$P/work/lnc/human_lnc.pick.{gtf,fa}`). Mouse: `$P/work/lnc/mouse_lnc.pick.fa`.
- **Era split:** GENCODE v28 gene IDs, `$P/data/gencode_v28/v28_gene_ids.txt` (58,336 version-stripped `gene` IDs from `gencode.v28.annotation.gtf.gz`).
- **GM12878:** final-default predictions and references from record E03.

## Software
- PyTransDecoder dev tree `$P/work/trees/PY_dev`, run by `run_tool.sh` (`PYTD_DEV`). The runs used the working tree of 2026-10-20; per-run commits were not recorded. From timestamps:

| Runs | Time (2026-10-20) | Code |
|---|---|---|
| bench `pytd_dev_kozak_add` (m100), wRand; GM12878 `kozak_add` | 07:43–07:56 | working tree before commit f5e6a39 (08:19, additive Kozak) |
| bench and GM12878 `m90`..`m50`; bench `kozak_add` on wLnc | 08:26–08:58 | f5e6a39 |
| bench and GM12878 `m50_sc*` (guard) | 09:18–09:26 | working tree before commit 7e9ffa2 (09:43, short-ORF guard) |
| bench `pytd_dev_v2` (lncRNA era table) | 13:47–13:58 | eaea367 (Kozak scores to 3 decimals) |
| GM12878 final defaults and `--min-aa` sweep | 21:05 Oct 20 – 06:52 Oct 21 | 48c870c |

- In these runs `-m 100` was the LongOrfs default. It became 50 at 429a89a.
- Python 3.11.17 (`$P/env/td`).

## Commands
**(a) Benchmark sweep** (launchers in `$T/benchmark_v2/launch_v3/`, run from `$P`):
```bash
# launch9.sh  (m100 baseline, all variants)
JOBS=8 TransDecoderPaper/benchmark_v2/scripts/run_all.sh data/bench_v3/inputs work/bench_v3/runs pytd_dev_kozak_add
# launch10.sh (random negatives)
VARIANTS="complete.wRand fragments.wRand" JOBS=8 TransDecoderPaper/benchmark_v2/scripts/run_all.sh \
  data/bench_v3/inputs work/bench_v3/runs pytd_dev_m90 pytd_dev_m80 pytd_dev_m70 pytd_dev_m60 pytd_dev_m50
# launch11.sh (lncRNA negatives, human and mouse)
VARIANTS="complete.wLnc fragments.wLnc" JOBS=4 TransDecoderPaper/benchmark_v2/scripts/run_all.sh \
  data/bench_v3/inputs work/bench_v3/runs pytd_dev_kozak_add pytd_dev_m90 pytd_dev_m80 pytd_dev_m70 pytd_dev_m60 pytd_dev_m50
# launch12.sh (guard)
VARIANTS="complete.wRand fragments.wRand complete.wLnc fragments.wLnc" JOBS=8 \
  TransDecoderPaper/benchmark_v2/scripts/run_all.sh data/bench_v3/inputs work/bench_v3/runs \
  pytd_dev_m50_sc200 pytd_dev_m50_sc150 pytd_dev_m50_sc250
```
`run_tool.sh` maps `pytd_dev_m<N>` to `longorfs -m <N>` and `_sc<C>` to `predict --short-orf-min-contrast <C>` (plus `-S` for stranded variants, `-G <gcode>`).

Scoring:
```bash
cd $P
env/td/bin/python TransDecoderPaper/benchmark_v2/scripts/score_by_length.py data/bench_v3/inputs work/bench_v3/runs \
  pytd_dev_kozak_add pytd_dev_m90 pytd_dev_m80 pytd_dev_m70 pytd_dev_m60 pytd_dev_m50 \
  pytd_dev_m50_sc150 pytd_dev_m50_sc200 pytd_dev_m50_sc250 > work/bench_v3/scores/by_length.tsv
awk -F'\t' -f TransDecoderPaper/benchmark_v2/scripts/short_orf_summary.awk work/bench_v3/scores/by_length.tsv \
  | (read h; echo "$h"; sort) > work/bench_v3/scores/short_orf_summary.tsv
for v in complete.wLnc fragments.wLnc; do echo "== $v"
  env/td/bin/python TransDecoderPaper/analyses/lnc_holdout/lnc_calls_by_era.py work/lnc/human_lnc.pick.gtf \
    data/gencode_v28/v28_gene_ids.txt work/bench_v3/runs/human $v pytd_dev_v2 cpat cpc2 rnasamba borf gmst prodigal
done > work/bench_v3/scores/lnc_by_era.human.tsv
```
- `score_by_length.py`:
  - truth is binned by true protein length (CDS nt / 3, less the stop);
  - predictions are binned by predicted length ((end − start + 1) // 3, stop included);
  - a TP is a stop match within 2 nt.
- `short_orf_summary.awk`:
  - `sens_short` = stop-level sensitivity for true proteins of 50–99 aa;
  - `neg` = calls on `rand_*`/`lnc_*` sequences;
  - `prec` = TP / (TP + fp_real + neg).

**(b) GM12878.** The pre-defined-minimum runs used `$T/genome_annotation/scripts/minlen_sweep.sh` (`longorfs -S -m 90..50`, LRAA and StringTie) and `guard_sweep.sh` (`-m 50` + `--short-orf-min-contrast 200`). The post-hoc filter on the final-default predictions:
```bash
bash $T/genome_annotation/scripts/minaa_sweep.sh
#   eval_genome_cds.py ... --label <asm>_min<L> --min-aa <L> -o $P/work/gm12878_td/eval/minaa/<asm>_min<L>
#   for L in 50 60 70 80 90 100 150 200 300, asm in lraa stringtie trinity
```
`--min-aa L` keeps predictions whose CDS encodes ≥ L aa (stop excluded), applied before collapsing to distinct CDSs. The headerless table `gm12878_minaa_sweep.tsv` was assembled from the 27 summaries by a one-off command. This reconstruction reproduces it exactly:
```bash
cd $T/genome_annotation/results/eval_summaries/minaa
for a in lraa stringtie trinity; do for L in 50 60 70 80 90 100 150 200 300; do f=${a}_min$L.summary.tsv
  g(){ awk -v k=$1 '$1==k{print $2}' $f; }
  printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" $a $L $(g distinct_predicted_cds) $(g MANE_sens_exact) \
    $(g MANE_sens_exact_given_assembled) $(g GENCODE_complete_sens_compatible) $(g pred_all_GENCODE_exact_frac) \
    $(g pred_all_GENCODE_compatible_frac) $(g pred_all_no_GENCODE_CDS_overlap_frac) $(g pred_complete_GENCODE_exact_frac)
done; done > ../../gm12878_minaa_sweep.tsv
```
Columns: assembly, min aa, distinct CDSs, MANE sens exact, MANE exact given assembled, GENCODE complete sens compatible, all-prediction GENCODE exact frac, GENCODE compatible frac, no-GENCODE-overlap frac, complete-prediction GENCODE exact frac.

## Outputs
- `$T/benchmark_v2/results_v3/by_length.tsv`, `short_orf_summary.tsv`, `lnc_by_era.human.tsv` (identical to `$P/work/bench_v3/scores/`)
- `$T/genome_annotation/results/gm12878_minaa_sweep.tsv` (= `eval_summaries/minaa_sweep.summary.tsv`)
- `$T/genome_annotation/results/eval_summaries/minaa/<asm>_min<L>.summary.tsv`
- `$T/genome_annotation/results/eval_summaries/{lraa,stringtie}_{kozak_add,m90,m80,m70,m60,m50,m50_sc200}.summary.tsv`
- Not in git:
  - benchmark predictions `$P/work/bench_v3/runs/<org>/<variant>/<tool>/preds.tsv`;
  - GM12878 runs `$P/work/gm12878_td/<asm>/{kozak_add,m*,m50_sc200}/`;
  - `$P/work/gm12878_td/eval/minaa/*.{ref,pred}_status.tsv`;
  - `$P/work/gm12878_td/ribo/*.ribo.tsv`.

## Results
**(a) Benchmark, human, lncRNA negatives** (`short_orf_summary.tsv`; 33,668 lncRNAs):

| Config | Input | Sens (all) | Sens 50–99 aa | lncRNA calls | Precision |
|---|---|---|---|---|---|
| `-m 100` | complete | 0.972 | 0.311 | 6,608 (19.6%) | 0.737 |
| `-m 90` | complete | 0.978 | 0.502 | 8,400 | 0.689 |
| `-m 80` | complete | 0.982 | 0.652 | 10,812 | 0.634 |
| `-m 70` | complete | 0.985 | 0.771 | 13,840 (41.1%) | 0.576 |
| `-m 60` | complete | 0.987 | 0.841 | 17,416 | 0.520 |
| `-m 50` | complete | 0.988 | 0.859 | 21,416 (63.6%) | 0.468 |
| `-m 50` + guard 150 | complete | 0.982 | 0.648 | 7,941 | 0.702 |
| **`-m 50` + guard 200** | complete | 0.980 | 0.564 | 6,943 (20.6%) | 0.729 |
| `-m 50` + guard 250 | complete | 0.977 | 0.479 | 6,699 | 0.735 |
| `-m 100` | fragments | 0.801 | 0.127 | 6,508 | 0.700 |
| `-m 50` | fragments | 0.934 | 0.913 | 21,100 | 0.455 |
| **`-m 50` + guard 200** | fragments | 0.865 | 0.521 | 6,833 | 0.706 |

- Mouse is similar: complete `-m 100` → `-m 50` + guard 200 gives short sensitivity 0.202 → 0.474 and lncRNA calls 5,088 → 5,423.
- Guard 150 reaches the short sensitivity of `-m 80` (0.648 vs 0.652) with 2,871 fewer lncRNA calls.
- **By predicted length, human complete.wLnc** (`by_length.tsv`):
  - `-m 50`: 14,880 lncRNA calls of 50–99 aa vs 430 correct calls in those bins.
  - `-m 100`: 4,139 lncRNA calls of 100–149 aa vs 958 correct.

**lncRNA call rate by GENCODE era** (`lnc_by_era.human.tsv`, complete.wLnc; TransDecoder = `pytd_dev_v2`, i.e. `-m 50` + guard 200):

| Tool | Genes in v28 (n = 12,987) | Added after v28 (n = 20,681) | All |
|---|---|---|---|
| TransDecoder | 0.220 | 0.197 | 20.6% |
| CPAT | 0.109 | 0.102 | 10.5% |
| CPC2 | 0.079 | 0.079 | 7.9% |
| RNAsamba | 0.062 | 0.109 | 9.1% |

RNAsamba calls fewer pre-v28 lncRNAs than later ones, which is consistent with training leakage. CPAT and CPC2 show no era effect. The other comparators and the sensitivity of each tool belong to record E07.

**Random negatives** (calls on `rand_*` per random sequence; 50–99 aa sensitivity):

| Organism | Complete: calls/seq m100 → m50 → m50+g200 | Fragments: short sens m100 → m50 → m50+g200 |
|---|---|---|
| A. thaliana | 0.0030 → 0.121 → 0.0054 | 0.076 → 0.886 → 0.303 |
| C. reinhardtii | 0.0033 → 0.013 → 0.0033 | 0.098 → 0.823 → 0.569 |
| D. melanogaster | 0.0019 → 0.070 → 0.0023 | 0.065 → 0.908 → 0.445 |
| Human | 0.0058 → 0.065 → 0.0080 | 0.124 → 0.908 → 0.514 |
| Mouse | 0.0026 → 0.053 → 0.0042 | 0.098 → 0.885 → 0.475 |
| P. falciparum | 0.0002 → 0.045 → 0.0002 | 0.036 → 0.947 → 0.205 |
| S. pombe | 0.0045 → 0.159 → 0.0061 | 0.055 → 0.945 → 0.438 |
| T. thermophila | 0.0138 → 0.146 → 0.0177 | 0.164 → 0.742 → 0.393 |

- On fragments, 50–99 aa sensitivity rises from 0.04–0.16 (`-m 100`) to 0.74–0.95 (`-m 50`).
- The guard is strict in AT-rich P. falciparum: complete short sensitivity is 0.108 with the guard vs 0.030 at `-m 100`.

**(b) GM12878, pre-defined minimum** (`eval_summaries/<asm>_<cfg>.summary.tsv`; m100 = `kozak_add`):

| | LRAA m100 → m50 → m50+g200 | StringTie m100 → m50 → m50+g200 |
|---|---|---|
| Distinct CDSs | 23,546 → 26,544 → 23,922 | 32,672 → 49,601 → 34,404 |
| GENCODE complete exact (refs) | 48,086 → 49,899 (+1,813) → 48,765 (+679) | 58,228 → 60,014 (+1,786) → 59,229 (+1,001) |
| MANE exact (refs) | 5,838 → 5,946 (+108) → 5,882 (+44) | 8,201 → 8,313 (+112) → 8,267 (+66) |
| No GENCODE CDS overlap | 3.6% → 10.7% → 3.7% | 11.8% → 34.7% → 12.0% |

- The m50+g200 summaries equal the final-default runs (E03), except StringTie GENCODE exact 59,229 vs 59,231 and MANE exact 8,267 vs 8,268. The code changed between the two runs (7e9ffa2 → 48c870c).
- **Ribo-seq ORF catalogue stop matches** (prediction models; `$P/work/gm12878_td/ribo/`):

| | m100 | m50 | m50+g200 |
|---|---|---|---|
| LRAA: all matches | 184 | 458 | 178 |
| LRAA: matches with any phase-1 peptide tier | 39 | 148 | 37 |
| StringTie: all matches | 221 | 578 | 217 |
| StringTie: matches with any phase-1 peptide tier | 53 | 157 | 53 |

  - The peptide tiers are HLA_tier1B, HLA_tier2B, peptidein_HLA, peptidein_nonHLA and nonHLA_peptide.
  - HLA tiers only (HLA_tier1B + HLA_tier2B + peptidein_HLA) give LRAA 27 → 132 → 24.

**GM12878, post-hoc `--min-aa`** (`gm12878_minaa_sweep.tsv`):

| Asm | Min aa | Distinct CDSs | MANE exact given assembled | GENCODE exact frac (all) | No GENCODE overlap |
|---|---|---|---|---|---|
| LRAA | 50 | 23,844 | 0.934 | 0.504 | 3.6% |
| LRAA | 100 | 23,100 | 0.916 | 0.506 | 3.2% |
| LRAA | 150 | 20,652 | 0.851 | 0.519 | 1.1% |
| LRAA | 200 | 18,380 | 0.768 | 0.519 | 0.7% |
| StringTie | 50 | 34,131 | 0.927 | 0.404 | 11.6% |
| StringTie | 100 | 31,448 | 0.912 | 0.427 | 10.3% |
| StringTie | 150 | 25,903 | 0.859 | 0.478 | 4.3% |
| StringTie | 200 | 22,221 | 0.792 | 0.503 | 2.4% |
| Trinity | 50 | 71,072 | 0.905 | 0.194 | 14.5% |
| Trinity | 100 | 64,017 | 0.889 | 0.207 | 12.3% |
| Trinity | 150 | 48,628 | 0.829 | 0.247 | 3.3% |
| Trinity | 200 | 40,306 | 0.755 | 0.263 | 1.5% |

The 60–90 and 300 aa rows are in the table. No-GENCODE-overlap fraction by predicted length, final defaults (from `$P/work/gm12878_td/eval/status/<asm>_final_abinitio.pred_status.tsv`):

| Asm | 50–99 aa | 100–149 | 150–199 | 200–299 | ≥ 300 |
|---|---|---|---|---|---|
| LRAA | 18.2% (n=822) | 20.5% (2,448) | 4.5% (2,272) | 1.8% (4,320) | 0.4% (14,060) |
| StringTie | 29.6% (2,956) | 38.6% (5,545) | 15.6% (3,682) | 6.3% (5,472) | 1.2% (16,749) |
| Trinity | 40.9% (8,306) | 40.7% (15,389) | 11.6% (8,322) | 3.6% (11,443) | 0.7% (28,863) |

## Checks
- Reproduced exactly on 2026-10-21 from the run outputs:
  - `by_length.tsv`, by re-running `score_by_length.py` (byte-identical);
  - `lnc_by_era.human.tsv`, by re-running `lnc_calls_by_era.py` (byte-identical);
  - `short_orf_summary.tsv`, from `by_length.tsv` with `short_orf_summary.awk` (identical rows and values as doubles);
  - `gm12878_minaa_sweep.tsv`, from the minaa summaries (byte-identical).
- `v28_gene_ids.txt` equals the version-stripped `gene` IDs of `gencode.v28.annotation.gtf.gz` (set comparison).
- The current `pytd_dev_v2` and `pytd_dev_m50_sc200` outputs give the same lncRNA call counts as the table (6,943 complete, 6,833 fragments). `pytd_dev_v2` was rerun at 13:47, after the table was written at 12:25.

## Caveats
- **Predicted-length bins include the stop codon** (`(end − start + 1) // 3`). So `-m 100` runs have calls in the "<50" and 50–99 bins (e.g. 331 and 603 on human complete.wLnc), from partial ORFs and start refinement.
- **"Calls per random sequence"** counts calls on `rand_*` sequences divided by the number of random sequences. A sequence with two calls counts twice.
- **The lncRNA negatives are not perfectly clean.** Some contain translated sORFs (E01, Ribo-seq catalogues).
- **The `--min-aa` filter is post hoc.** It differs from running with a different `-m`, because Predict's training and per-transcript selection see the full candidate set.
  - Even `--min-aa 50` removes 78 LRAA distinct CDSs from the final defaults (23,922 → 23,844; 77 of them complete). These encode < 50 aa, presumably after start refinement moved the start downstream.
- **The GM12878 m100 baseline (`kozak_add`) was run without a saved script.** Its LRAA directory has no LongOrfs log, so it presumably reused the earlier ab initio LongOrfs output.
- **The benchmark and GM12878 sweeps used uncommitted or intermediate PyTD code** (Software table). Only the final defaults (E03) and the `--min-aa` sweep are at 48c870c.
