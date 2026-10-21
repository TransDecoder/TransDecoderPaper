# E10. GC-calibrated long-ORF retention rule

## Claims supported
Manuscript: Methods, "ORF selection" → "GC-calibrated long-ORF rule" (`methods_algorithm.tex`); Methods, "Background composition".

1. An ORF whose CDS length is at least $L^\ast(\mathrm{GC})$ is kept whatever its coding score. $L^\ast$ is the 0.999 quantile of the lengths of ORFs that LongOrfs extracts from i.i.d. random sequence of the same G+C content. The sequence is scanned on both strands, partial ORFs are allowed, and the minimum is 300 nt.
2. $L^\ast$ was estimated by simulation in 5% GC steps from 25% to 80%, with at least 85,000 random ORFs per level. It rises from 519 nt at 25% GC to 2,436 nt at 80% GC.
3. The rule is disabled above 80% GC. `--retain_long_orfs_mode strict` sets a fixed length instead.
4. The table is new. The 2017 table (TransDecoder `33247f5`) had no surviving simulation code. It is reproduced within about 3% at ≥50% GC, but at 25–45% GC it was 3–10% too low because its sample was small (notes/01 §6, notes/03, notes/04).
5. Implementation details not yet in the manuscript (see Caveats):
   - GC is now (G+C)/(A+C+G+T) from base counts.
   - The rule applies only to genetic codes whose stops are TAA/TAG/TGA.
   - Perl and PyTD apply it identically.

## Data
No external data are used. Random sequences are simulated with seeds:
- Low-GC chunks (25–45% GC):
  - Seeds `gc*1000 + c`; chunk counts are 50/20/10/5/3 for 25/30/35/40/45% GC.
  - Each chunk has 100,000 sequences of 5 kb.
  - Raw lengths (88 files, 2.4 MB): `$P/work/thresh_big/raw/gc<GC>_seed<seed>.tsv`.
  - UGER job logs: `$P/work/thresh_big/log/thr<GC>_<c>.o<jobid>`.
- Single-run simulations (seed 1): `$T/analyses/long_orf_threshold/sim_runs/` and `sim_large/`.
- Benchmark inputs used for the GC estimate checks are bench_v3 (record E07).

## Software
- Python 3.11.17 (`$P/env/td`), standard library only.
- This repo's simulation code was added in `a015f4a`. The tables were added in `0c83303` and `3c617bf`.
- Table in code:

  | Repo | Commit | Change |
  |---|---|---|
  | PyTD | `f5b12e7` (on origin/main) | regenerated table |
  | PyTD | `efc1450` | length auto-accept only for TAA/TAG/TGA codes |
  | PyTD | `48c870c` (branch dev/frame-contrast-starts, local only) | GC over A/C/G/T only; no length retention above 80% GC (it previously used 2,436 nt); deterministic ORF tie-break |
  | Perl TD | `694e204` | regenerated table |
  | Perl TD | `f7e75a1` | dynamic retention only for TAA/TAG/TGA codes |
  | Perl TD | `de3f57e` (branch fix/start-refinement, local only) | GC from G and C counts in `base_freqs.dat`, previously 2×freq(C); deterministic tie-break |

- The table is identical in `$P/work/trees/PY_dev/pytransdecoder/predict.py` (`gc_to_min_length`, line ~140) and `$P/work/trees/TD_fix/util/TransDecoder.Predict` (`@GC_TO_MIN_LONG_ORF_LENGTH`, line ~152). Both equal the q0.999 column of `final_thresholds.tsv`. The Perl header comment also carries the 0.95/0.99/0.9999 columns, rounded.
- Lookup: take the first table row with GC% ≤ the row's GC. So the threshold rounds up to the next 5% step, and GC < 25% uses 519 nt. Above 80%: no length retention (Perl returns 10⁶ nt, PyTD returns ∞). For a code with other stops (e.g. Tetrahymena, TGA only): no length retention.

## Commands
Simulator: `$T/analyses/long_orf_threshold/simulate_random_orf_lengths.py`.
- Bases are i.i.d. with P(G)=P(C)=gc/2.
- For each stop, the ORF runs from the most upstream ATG after the previous in-frame stop.
- Both strands are scanned. With `--partials`, ORFs may be open at the 5′ end and/or the 3′ end.
- Lengths include the stop. ORFs must be ≥ 300 nt.
- Quantiles use linear interpolation (R type 7).

Low-GC levels (25–45%). The exact `qsub` form is in git `3c617bf`. The local equivalent is:
```bash
bash $T/analyses/long_orf_threshold/run_low_gc_chunks.sh $P/work/thresh_big
# each chunk: simulate_random_orf_lengths.py --gc <GC> --n-seqs 100000 --seq-len 5000 --partials \
#             --seed <GC*1000+c> --raw-out raw/gc<GC>_seed<seed>.tsv
python $T/analyses/long_orf_threshold/pool_chunk_quantiles.py $P/work/thresh_big/raw > lowgc_large_sample.tsv
python $T/analyses/long_orf_threshold/pool_chunk_quantiles.py $P/work/thresh_big/raw --boot 200   # 95% CI of q0.999
```
`final_thresholds.tsv` = the 25–45% rows of `lowgc_large_sample.tsv` + the 50–80% rows of `sim_large/partial_L5000_gc*.tsv`.

Benchmark checks of the GC estimate: `$T/benchmark_v2/launch_v3/launch15.sh` reruns the 7 Perl/PyTD-mismatched inputs after the GC estimate was unified. The earlier runs were moved to `$P/work/bench_v3/rerun_old/`. `launch16.sh` reruns all inputs after the tie-break fix (earlier runs in `$P/work/bench_v3/rerun_old16/`).

## Outputs
- `$T/analyses/long_orf_threshold/final_thresholds.tsv`: the adopted table (n_orfs, q0.95/0.99/0.999/0.9999).
- `lowgc_large_sample.tsv`: pooled low-GC chunks.
- `regenerated_thresholds.tsv`: single 5 kb run, seed 1, with the 2017 q0.999 column.
- `table_2017.tsv`: the 2017 table, from the comment in TransDecoder `33247f5`.
- `sim_runs/{complete,partial}_{L2000,L5000,mouseLen}.tsv`: sensitivity to sequence length and partials.
- `sim_large/partial_L{5000,10000}_gc*.tsv`: per-GC large runs.

## Results
q0.999 ORF length (nt):

| GC% | 25 | 30 | 35 | 40 | 45 | 50 | 55 | 60 | 65 | 70 | 75 | 80 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| Adopted ($L^\ast$) | 519 | 552 | 573 | 615 | 663 | 726 | 804 | 927 | 1083 | 1347 | 1746 | 2436 |
| Random ORFs (n) | 85,915 | 85,056 | 101,458 | 115,742 | 150,612 | 103,782 | 200,546 | 357,795 | 584,858 | 860,669 | 1,121,478 | 1,264,679 |
| 95% CI (bootstrap, 200) | 510–525 | 546–561 | 564–582 | 606–621 | 651–672 | – | – | – | – | – | – | – |
| 2017 | 465 | 511 | 555 | 590 | 645 | 749 | 797 | 927 | 1086 | 1358 | 1743 | 2423 |
| 2017 vs adopted | −10.3% | −7.4% | −3.2% | −4.1% | −2.8% | +3.1% | −0.9% | 0.0% | +0.3% | +0.8% | −0.2% | −0.5% |

- The CIs were recomputed for this record with `pool_chunk_quantiles.py --boot 200`. The rerun reproduced `lowgc_large_sample.tsv` exactly. The CI half-widths are 1–2%.
- At ≥50% GC the 2017 values lie within 3.1% of the adopted table.
- At 25% GC, the 2017 q0.999 (465) is close to its own q0.9999 (468). That fits a small sample.

Sensitivity of q0.999 to the simulation set-up (single runs):
- Sequence length 10 kb vs 5 kb (`sim_large`): within 1% at ≤65% GC, +1.3% at 70%, +1.5% at 75%, and +3.1% at 80% (2,511 vs 2,436).
- Partials included (`sim_runs`, 5 kb, ~10× fewer ORFs): within 0.5–4.2% of `final` at ≥55% GC (largest gap at 70%: 1,290 vs 1,347). The low-GC rows have few ORFs (146–4,082), so they are noisy.
- Mouse transcript lengths (`*_mouseLen`) give lower values, e.g. 2,184 vs 2,436–2,455 at 80% and 904 vs 914–927 at 60%. So the threshold depends weakly on the transcript length distribution.

Effect of the GC estimate (bench_v3, stranded inputs, Perl `PCT_GC` vs PyTD log):

| Input | Perl before de3f57e (2×C) | After (G+C) | Threshold before → after |
|---|---|---|---|
| athal complete | 39.4% | 42.17% | 615 → 663 |
| creinhardtii complete | 64.4% | 68.29% | 1083 → 1347 |
| creinhardtii complete_fullutr | 61.4% | 66.16% | 1083 → 1347 |
| mouse complete_fullutr | 48.2% | 48.35% | 726 → 726 |

Perl/PyTD differences on the 7 rerun inputs are differing prediction lines (`comm -3` on transcript/start/end/strand):

| Input | Before launch15 | After launch15 (GC unified) | Final (launch16, + tie-break) |
|---|---|---|---|
| athal complete / complete_fullutr / fragments | 9 / 7 / 19 | 0 / 0 / 0 | 0 / 0 / 0 |
| creinhardtii complete / complete_fullutr / fragments | 7 / 11 / 1 | 0 / 0 / 0 | 0 / 0 / 0 |
| mouse complete_fullutr | 2 | 2 | 0 |

In the final runs:
- Thresholds applied: P. falciparum 519 nt (GC 20–22%); athal 663; human/mouse complete 804, complete_fullutr 726; C. reinhardtii 1,347 (GC 66–68%).
- T. thermophila: no length retention (TGA-only code).
- Perl and PyTD logs agree to 0.1%.

## Checks
- **Code tables:** the Perl and PyTD tables equal `final_thresholds.tsv` q0.999.
- **Reproduction:** the bootstrap rerun of `pool_chunk_quantiles.py` on `$P/work/thresh_big/raw` reproduced `lowgc_large_sample.tsv` exactly (n and all quantiles). The ≥50% rows of `final_thresholds.tsv` are identical to `sim_large/partial_L5000_gc*.tsv` and `regenerated_thresholds.tsv`.
- **2017 regeneration:** the simulator reproduces the 2017 table within about 3% at ≥50% GC, e.g. 60%: 927 vs 927; 80%: 2,436 vs 2,423.
- **Perl/PyTD parity:** with the unified GC estimate, Perl and PyTD give identical predictions on all 84 bench_v3 inputs (E08).

## Gaps
- **The 2017 table's origin is not reproducible.** Only the pasted quantiles exist (TransDecoder `33247f5`). Its simulation code, sample size and settings are unknown.
- **Commands not saved for `sim_runs/` and `sim_large/`.**
  - Not recorded: `--n-seqs`, the seeds (seed 1 is assumed from the default), and the `--length-fasta` file used for `*_mouseLen`.
  - So the 50–80% rows of the adopted table come from runs whose exact command line is not recorded. The sample sizes (n_orfs) are in the files.
- **The final table was merged by hand.** `final_thresholds.tsv` is 25–45% from `lowgc_large_sample.tsv` plus 50–80% from `sim_large`, with no script.
- **Bootstrap CI.** The "±1–2%" figure appears in notes/04 and in the commit messages of `694e204`/`f5b12e7`, but no CI output file was saved. The CI above was recomputed for this record; it is not a committed file. There is no CI for ≥50% GC.
- **The low-GC chunks ran on UGER** (`thr*.o<jobid>` logs). `run_low_gc_chunks.sh` now runs them locally, with the same arguments.
- **No analysis supports the stated motivation.** The claim is that the rule "recovers genes with atypical composition, e.g. horizontally transferred genes or mixed-organism transcriptomes". No analysis in the repo shows this, and none counts how many benchmark or GM12878 predictions are retained by length alone.
- **The final-default benchmark runs predate their commits.** `launch16`, 2026-10-20 12:38–15:32, ran before `de3f57e` (15:37:58) and `48c870c` (15:38:00) were committed. The trees were uncommitted at run time, so equivalence to the commits is inferred, not recorded. Partial check (2026-10-21, E07): rerunning `pytd_dev_v2` and `perl_td_dev_v2` at the committed trees (PY_dev 48c870c, TD_fix de3f57e) on spombe/complete and creinhardtii/complete gave byte-identical `preds.tsv`. The other 82 inputs per tool were not rerun.
- 48c870c is on PyTransDecoder origin/main and de3f57e on TransDecoder origin/devel (verified 2026-10-21 with git rev-list against the fetched remotes).

## Caveats
- **The manuscript is out of date** (`methods_algorithm.tex`):
  - "Background composition" says G+C = 2q(C). Since `de3f57e`/`48c870c` it is (G+C)/(A+C+G+T) from counts. 2q(C) is wrong for stranded data, where q(C) ≠ q(G) (athal 39.4% vs 42.2%).
  - It does not say that the rule is skipped for genetic codes whose stops are not TAA/TAG/TGA, nor how the 5% table is looked up (round up to the next step).
- **What the threshold models:** a single i.i.d. sequence model with universal stops, 5 kb sequences and a 300-nt LongOrfs minimum. It is a per-ORF 0.999 quantile, not a transcriptome-wide FDR; the old usage text said "1% FDR".
- **LongOrfs minimum:** the threshold was simulated with a 300-nt minimum. The final default LongOrfs `-m 50` changes the candidate set, but not the upper tail that the q0.999 describes.
- **Not covered here:** `$T/analyses/short_orf_threshold/collect_calls.py` and `results_v3/short_orf_summary.tsv` concern the short-ORF frame-contrast guard (`-m 50`, notes/05). That is a separate rule (E04), not the long-ORF retention.
- **Minor discrepancies between the notes and the files:**
  - notes/03 says "10 kb sequences give the same values". The files show up to +3.1% at 80% GC (2,511 vs 2,436).
  - notes/03 gives "2–10% lower" for the 2017 table at ≤45% GC, but against the adopted table the difference is 2.8–10.3%. Against the 5 kb single run (`regenerated_thresholds.tsv`) it is 1.2–10.4%, at 40% GC 590 vs 597.
  - The 2017 value at 80% GC is written as 2423 in notes/03 and as 2422 in notes/04; the file has 2422.656.
