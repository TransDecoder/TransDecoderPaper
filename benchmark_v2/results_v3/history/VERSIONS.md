# TransDecoder version history: release table, key versions, and smoke tests

Prepared 2026-10-24 for the supplementary "historical perspective" in the Genome Biology paper.
Sources: git tags and `git log`/`git diff` between tags in `repos/TransDecoder`; `Changelog.txt`;
the wiki checkout (which has no per-version history); the SourceForge archives that Brian downloaded to
`work/history/sourceforge/`; Trinity 2.1.1 `Release.Notes` (`env/trinity/opt/trinity-2.1.1`).
Statements marked *(code)* come from reading the source, not from the changelog.

## 1. Lineage before GitHub

| Release (source) | Date | Key features / changes |
|---|---|---|
| Trinity `Analysis/Coding/transcripts_to_best_scoring_ORFs.pl` | 2011 to 2012 | The ORF predictor began as a Trinity utility (Haas et al. 2013 Nat Protoc). Trinity Release.Notes: 2011-11-26 "bugfix: updated handling of partial genes on reverse strand"; 2012-06-08 "find orfs on both strands by default, use -S for top-strand only". Later Trinity bundled TransDecoder as a plugin: "updated to release 11-10-2013", then "updated TransDecoder to 20140704 version". This repo has no copy of the Trinity-era script (Trinity 2.1.1 no longer ships it). The 2012 SourceForge release keeps the same script name, so the Trinity script was very likely the same algorithm (inferred; not verified). |
| transdecoder-r2012-08-07 (SF) | 2012-08-07 | First standalone release. A single script, `transcripts_to_best_scoring_ORFs.pl`, writes fixed-name outputs (`best_candidates.eclipsed_orfs_removed.{gff3,pep,cds,bed}`) to the current directory. ORFs: longest ORF per stop codon, partials allowed, `-m 100` aa, `-S`, `-G` (universal, Euplotes, Tetrahymena, Candida, Acetabularia, mitochondrial codes). **Markov model** *(code)*: framed 5th-order (hexamer) log-likelihood against background base frequencies from the transcripts (both strands unless -S), with a pentamer initial-state term. **Training set:** the 500 longest ORFs with **no redundancy filter**. **Selection:** keep an ORF if it has a Pfam hit (built-in `hmmscan --search_pfam`), or length ≥ 900 nt (`--retain_long_orfs`), or frame-1 score > 0 and greater than all 5 other frames. Then remove "eclipsed" ORFs, meaning ORFs contained inside a longer selected ORF. *(code)* Here the gene features span only the CDS, so this filter really removes nested ORFs. Known bug: `-T` is declared without `=i`, so it cannot be set. |
| transdecoder_r2012-08-15 (SF) | 2012-08-15 | The same algorithm. Adds docs, `cufflinks_gtf_to_bed.pl`, and a BED `track name` header. |
| TransDecoder_r20131110 (SF) | 2013-11-10 | Script renamed `TransDecoder`. Adds a work directory (`transdecoder.tmp.<pid>`, `--workdir`), `--reuse`, `--train`, and a Pfam runner (`pfam_runner.pl`, MPI/ffindex, `--prepare_pfam`). Outputs are now named `<input>.transdecoder.{gff3,pep,cds,bed}`. **Training set:** if cd-hit-est is in PATH, the 4×T longest ORFs are clustered with cd-hit-est at its default 90% identity and the top 500 are kept. **The gene/exon features now span the whole transcript** (UTRs included), and model IDs become `<acc>|m.N`. *(code)* Side effect: `remove_eclipsed_ORFs.pl` compares gene spans, which are now identical for every ORF on a transcript, so from this release through v4.1.0 the eclipsed-ORF filter does nothing (confirmed empirically for v2.0, see §5). |
| TransDecoder_r20131117 (SF) | 2013-11-17 | cd-hit-est becomes **required** (`check_program`). Adds POD usage, ffindex/ParaFly builds via Makefile, and a `pfam/` directory. The prediction algorithm is unchanged. |
| transdecoder_rel13JAN2014 / rel16JAN2014 (SF) | 2014-01-13 / 2014-01-16 | Documentation and build changes, plus an MPI LD_LIBRARY_PATH check. The algorithm is unchanged. |
| TransDecoder_r20140704 (SF; last SF release) | 2014-07-04 | Release.Notes: "make simple" builds only parafly and cd-hit, and the `cds.` prefix is removed from pep/cds accessions. Writes an `.mRNA` output. The algorithm is the same as r20131110. |

r20140704 differs from v2.0 only in the training-set redundancy step (cd-hit-est 90% on the top 2,000 becomes 80%
on the top 5,000), the two-phase split, and Pfam/BLAST now being supplied as precomputed results. The Markov
model, scoring and selection rule (Pfam/BLAST, ≥ 900 nt, frame rule) are the same.

## 2. GitHub tags (repos/TransDecoder)

| Tag | Date | Substantive changes (algorithmic or user-facing) |
|---|---|---|
| v2.0 | 2015-01-26 | Moved to GitHub, with two phases: **TransDecoder.LongOrfs** (`-m 100`, `-S`, `-G`; writes `<input>.transdecoder_dir/longest_orfs.*`) and **TransDecoder.Predict**. Active Pfam/MPI search is removed; Predict instead takes precomputed `--retain_pfam_hits` (hmmscan domtbl) and, **new**, `--retain_blastp_hits` (outfmt6). **Training:** top 5,000 longest ORFs, then cd-hit-est `-c 0.80`, then top 500 (`-T`). The Markov model and frame rule are as in 2012 (6-frame rule, 900 nt long-ORF retention, eclipsed filter, ineffective as noted above). Outputs `.gff3/.pep/.cds/.bed/.mRNA`. Predict has **no `-G`**. Needs cd-hit-est (bundled 4.6.1 tgz) and the Perl modules DB_File and URI::Escape. |
| 2.0.1 | 2015-01-27 | Packaging only (keeps the util/bin directory). |
| v2.1.0 | 2016-03-11 | Predict `--cpu` (passed to cd-hit-est). LongOrfs `-p <pct>` optionally trims 5′ partials to an internal Met. Keeps Cufflinks gene IDs. Adds sample data and `get_longest_ORF_per_transcript.pl`. Defaults unchanged. |
| v3.0.0 | 2016-04-26 | LongOrfs `--gene_trans_map`. Predict **`--single_best_orf`** (opt-in; ranks by homology, then ORF length). Predict `-G`. Phase in the GFF3. The `.mRNA` output is dropped. Kept-ORF category counts are logged. Default predictions are unchanged. |
| v3.0.1 | 2016-10-02 | Better alternative genetic code support; `which` check fix. (Changelog "v3.0.2" corresponds to this tag.) Smoke predictions are **identical to v2.0**. |
| v4.0.0 | 2017-08-04 | Retracted ("all v4 releases pre-v4.1 were fairly quickly retracted due to bugs"). Introduced everything listed for v4.1.0, plus **GC-dynamic long-ORF retention** (first table: 25% GC → 449 nt … 80% GC → 1,728 nt). |
| v4.0.1 / v4.0.2 / v4.0.3 | 2017-08-13 / 08-16 / 08-16 | Fixes for non-GATC characters in the start PWM, the start-refinement workdir, and `--retain_long_orfs_mode strict`. Retracted. |
| **v4.1.0** | 2017-08-17 | **Markov model rewritten** *(code)*: a variable-order chain (orders 0 to 5) with +1/+4 pseudocounts, scoring from position 0, and the terminal stop codon excluded from scoring. **Training:** cd-hit is replaced by the in-house `exclude_similar_proteins.pl` (spaced 5-mer protein comparison at 80% proxy identity; low-complexity skip), and training ORFs are capped at 5,000 nt. **Single best ORF per transcript by default** (`--all_good_orfs` turns it off; ranking is homology, then frame-1 score, then length). **Start-codon refinement on by default** (`--no_refine_starts` turns it off): a +/− start-site PWM is trained on the training ORFs and used to move 5′-partial ORF starts to an internal ATG. **No long-ORF auto-retention by default** (`--retain_long_orfs_length`, default 1,000,000). The selection rule is still frame-1 > 0 and > the other 5 frames. Adds Pipeliner checkpoints (`<workdir>.__checkpoints`). Optional R seqLogo plots. |
| **v5.0.0** | 2017-08-26 | Selection moves to `select_best_ORFs_per_transcript.pl`. **Frame rule: frame-1 score > 0 and > frames 2 and 3 only** (the reverse frames are no longer competitors). **GC-dynamic long-ORF retention returns as the default** (`--retain_long_orfs_mode dynamic`, "1% FDR in random sequence of the same GC"; table 25% → 465 nt … 50% → 749 nt … 80% → 2,422 nt; `strict` mode uses `--retain_long_orfs_length`). Candidates are ranked by homology count, then frame-1 score, then length, then kept greedily if they overlap selected ORFs by ≤ 10%. `--all_good_orfs` is removed and **`--single_best_only`** is added (opt-in). Score and BLAST/Pfam info go into the GFF3 Name. *(code)* **Bug:** the overlap test used gene span (the whole transcript), so v5.0.0 through v5.3.0 report **effectively one ORF per transcript** (the top-ranked one) even without `--single_best_only`. The v5.4.0 changelog acknowledges this. |
| v5.0.1 | 2017-09-08 | Copes with non-GATC characters in start adjustment; GMAP-aligned transcripts supported. |
| TransDecoder-v5.0.2 / v5.0.3 | 2017-10-20 | Protein IDs become `<acc>.p1, .p2, …`. Gzipped input accepted. In-frame stop validation. Supertranscripts example. LongOrfs `-p` is gone. |
| TransDecoder-v5.1.0 | 2018-03-27 | **NCBI genetic code tables** (names such as `Universal`, `Tetrahymena`, `Ciliate`, …). No UTR records when the CDS starts at a partial codon. `--version`. |
| TransDecoder-v5.2.0 | 2018-04-21 | Predict can be rerun with BLAST/Pfam results (redoes only the selection). The pipeliner is reworked. Test data include mini Pfam/SwissProt. |
| TransDecoder-v5.3.0 | 2018-05-11 | **DB_File requirement removed** (GFF3_utils2 in-memory indexing). gene-GTF output. |
| TransDecoder-v5.4.0 | 2018-10-17 | **Bug fix: the overlap test uses the CDS span**, so multiple non-overlapping ORFs per transcript are reported again unless `--single_best_only` is given. Run-mode tests; less verbose. |
| **TransDecoder-v5.5.0** | 2018-10-22 | `--output_dir/-O` for both phases. The long-lived release (2018 to 2023). |
| (v5.6.0, untagged) | 2022-04-07 | Genome propagation fixes. `compute_base_probs` speed-up. With `--single_best_only`, selection happens before overlap removal. **All BLAST-supported ORFs are retained regardless of overlap** (homology-supported ORFs bypass the overlap filter). |
| TransDecoder-v5.7.0 | 2023-01-27 | Pfam domtbl from **hmmsearch or hmmscan**. LongOrfs `--complete_orfs_only`. GFF3 attributes no longer URI-encoded, and the **URI::Escape dependency is dropped**. Lengths in FASTA headers corrected. `srand(1234)`. Checkpoints reorganized. |
| TransDecoder-v5.7.1 | 2023-07-16 | Docker support. The workdir name is fixed within `-O`. |
| TransDecoder-v6.0.0 | 2026-03-20 | New `TransDecoder` full-pipeline wrapper (transcripts or genome+GTF; optional BLAST/DIAMOND `--blast_search_pep` and Pfam `--pfam-search-db`). Phase scripts move to `util/TransDecoder.{LongOrfs,Predict}`. BED output deduplicated. *(code)* **The phase scripts are algorithmically identical to v5.7.1** (diffs touch only paths and `--version`), and smoke predictions are identical. |

### v6.0.0 → current development (summary)

`git log TransDecoder-v6.0.0..` gives 15 commits on `origin/devel` (de3f57e, 2026-10-20). The TD_fix worktree is on branch
`fix/start-refinement` (ab38ca6), which is origin/devel plus one commit. Note that there is no local `devel` branch;
use `origin/devel`.

- **Start refinement redesigned:** scored by frame contrast plus Kozak context (f7e75a1, 8cf8d5f); Kozak training is its own step (`train_kozak_model.pl`, with a logo; 7e30234); Kozak evidence enters the keep-or-trim cutoff additively (313be29); refinement is skipped when there are too few training starts (351634a); and in TD_fix, coding and UTR are scored as two scores with Kozak in bits (ab38ca6).
- **Single best ORF per transcript by default** again (f7e75a1).
- **Minimum ORF length 50 aa** (was 100), with a frame-contrast test for ORFs under 100 codons (8e99201).
- **Long-ORF retention thresholds** regenerated from a new random-sequence simulation (694e204). Retention applies only to codes whose stops are TAA/TAG/TGA (f7e75a1). GC is computed from G+C counts (de3f57e).
- **Genetic code:** the terminal stop of the selected code is excluded from Markov training and scoring (b81ffd0), and the training redundancy filter uses the selected code (351634a).
- **Determinism:** deterministic tie-break in ORF selection (de3f57e), GFF3 indexing in coordinate order (4405a85), and Fasta_reader ignores text before the first header (108e0fb).
- **Genome mode:** BED output for genome-propagated ORFs (c2f1741, d437a6c), plus gene-ID and usage fixes (49e78f6).

## 3. Recommended key versions to benchmark

Each one marks a distinct change in method. Versions whose default predictions are identical to a neighbour were dropped.

| # | Version | Represents | Why it is a key point | Tree |
|---|---|---|---|---|
| 1 | **r2012-08-07** (SF) | 2012 standalone, which is effectively the Trinity-era predictor | Original method: 5th-order Markov, training on the 500 longest ORFs with no redundancy filter, 6-frame rule, ≥ 900 nt retention, and a working eclipsed-ORF filter. The baseline for the whole history. | `work/history/sf_runs/r2012-08-07` (script-only copy) |
| 2 | **v2.0** | r20131110 … v3.0.1 ("classic" TransDecoder, 2013 to 2017) | The same model and rule as 2012, but with a redundancy-filtered training set (cd-hit-est) and the eclipsed filter made a no-op. Two-phase interface, BLAST input. Default output is identical to v3.0.1 in the smoke test. r20140704 differed from v2.0 by 1 ORF out of 196. | `work/trees/TD_hist_v2.0` |
| 3 | **v4.1.0** | August 2017 model rewrite | Variable-order Markov with pseudocounts and stop excluded; in-house redundancy filter; single best ORF by default; PWM start refinement by default; no long-ORF retention. Short-lived (9 days before v5.0.0) but cleanly isolates these changes from the v5 selection changes. Drop it if only 5 versions are wanted. | `work/trees/TD_hist_v4.1.0` |
| 4 | **v5.0.0** | v5.0 to v5.3 | 3-frame (same-strand) rule; GC-dynamic long-ORF retention; greedy overlap-based selection (de facto one ORF per transcript because of the span bug). | `work/trees/TD_hist_v5.0.0` |
| 5 | **v6.0.0** (≡ v5.7.1) | v5.4 to v6.0 (2018 to 2026; the widely used releases) | Multiple non-overlapping ORFs per transcript (v5.4 fix); homology ORFs bypass the overlap filter (v5.6+). Already in the benchmark as `perl_td_v600` and `perl_td_v571`. v5.7.1 and v6.0.0 give identical predictions, so run only one of them for the history panel. | `work/trees/TD_base` (v6.0.0+Dockerfile commit), `work/trees/TD_v571` |
| 6 | **current devel / v7** | this paper | Kozak/frame-contrast start refinement, single best ORF by default, 50 aa minimum with the short-ORF test, new retention thresholds, code-aware stops. | `work/trees/TD_fix` |

Optional variants, worth adding only if a specific question needs them:
- **v5.5.0** (`work/trees/TD_hist_TransDecoder-v5.5.0`) adds nothing over v6.0.0 in ab initio mode (smoke-identical). It differs only in homology mode, where homology-supported ORFs were still subject to the overlap filter before v5.6/5.7. Use it if the history panel includes BLAST/Pfam runs.
- **v3.0.1 `--single_best_orf`** (`work/trees/TD_hist_v3.0.1`) is the first single-best option (homology, then length). Its default output equals v2.0.

## 4. Per-version working command lines and smoke results

Smoke input: `work/history/smoke/input/subset300.fa`, a deterministic 300-transcript subset of
`data/bench_v3/inputs/spombe/complete.wRand.fa` (200 S. pombe transcripts taken every ~25th, plus 100 `rand_*`). The input is stranded, so
`-S` is used throughout. Runners: `work/history/scripts/run_smoke.sh` (two-phase versions) and `run_smoke_sf.sh`
(SF releases). Both set `PATH=env/td/bin:$PATH` and `PERL_HASH_SEED=0`, and symlink the input as `input.fa` in the run directory.
Outputs are in `work/history/smoke/<label>/`. Each run executes in its own directory, because old versions write into the current directory.

`COMPAT` means: `export PERL5LIB=/broad/hptmp/bhaas/TRANSDECODER_PAPER/work/history/compat/perl5:$PERL5LIB LD_LIBRARY_PATH=/broad/hptmp/bhaas/TRANSDECODER_PAPER/work/history/compat/lib`
(see §5). `$T=/broad/hptmp/bhaas/TRANSDECODER_PAPER/work/trees`, `$R=/broad/hptmp/bhaas/TRANSDECODER_PAPER/work/history/sf_runs`.

| Version | Commands (run inside a fresh directory containing input.fa) | rc | LongOrfs s | Predict s | Final ORFs (rand tx) | Converter |
|---|---|---|---|---|---|---|
| r2012-08-07 | COMPAT; `$R/r2012-08-07/transcripts_to_best_scoring_ORFs.pl -t input.fa -S` → `best_candidates.eclipsed_orfs_removed.gff3` | 0 | 5.0 (single step) | – | 196 (3) | OK |
| r20140704 | COMPAT; `$R/r20140704/TransDecoder -t input.fa -S --CPU 2` → `input.fa.transdecoder.gff3` | 0 | 5.8 (single step) | – | 196 (3) | OK |
| v2.0 | COMPAT; `$T/TD_hist_v2.0/TransDecoder.LongOrfs -t input.fa -S`; `$T/TD_hist_v2.0/TransDecoder.Predict -t input.fa` | 0 | 0.7 | 4.5 | 195 (3) | OK |
| v3.0.1 | COMPAT; `$T/TD_hist_v3.0.1/TransDecoder.LongOrfs -t input.fa -S`; `…/TransDecoder.Predict -t input.fa` (optionally add `--single_best_orf`) | 0 | 1.5 | 5.7 | 195 (3), = v2.0 | OK |
| v4.1.0 | COMPAT; `$T/TD_hist_v4.1.0/TransDecoder.LongOrfs -t input.fa -S`; `…/TransDecoder.Predict -t input.fa` | 0 | 1.5 | 71.5 | 195 (3) | OK |
| v5.0.0 | COMPAT; `$T/TD_hist_v5.0.0/TransDecoder.LongOrfs -t input.fa -S`; `…/TransDecoder.Predict -t input.fa` (logs `PCT_GC: 37.8`, auto-accept ≥ 590 nt) | 0 | 1.3 | 74.1 | 195 (3) | OK |
| v5.5.0 | `$T/TD_hist_TransDecoder-v5.5.0/TransDecoder.LongOrfs -t input.fa -S -G Universal`; `…/TransDecoder.Predict -t input.fa -G Universal` | 0 | 1.1 | 83.3 | 195 (3) | OK |
| v5.7.1 | `$T/TD_v571/TransDecoder.LongOrfs -t input.fa -S -G Universal`; `…/TransDecoder.Predict -t input.fa -G Universal` | 0 | 1.3 | 81.8 | 195 (3) | OK |
| v6.0.0 | `$T/TD_base/util/TransDecoder.LongOrfs -t input.fa -S -G Universal`; `$T/TD_base/util/TransDecoder.Predict -t input.fa -G Universal` | 0 | 1.4 | 84.6 | 195 (3), = v5.7.1 | OK |
| devel (TD_fix) | `$T/TD_fix/util/TransDecoder.LongOrfs -t input.fa -S -G Universal`; `$T/TD_fix/util/TransDecoder.Predict -t input.fa -G Universal` | 0 | 1.4 | 7.3 | 231 (31) | OK |

Agreement on (transcript, start, end, strand):
- v2.0 = v3.0.1 (195/195).
- v4.1.0 = v5.0.0 = v5.5.0 = v5.7.1 = v6.0.0 (195/195); these share 184/195 with v2.0, and the differences come from start refinement.
- The SF releases share 195 with v2.0 and add one ORF (SPCC4G3.03.1:1038-1352), which reflects the training-set differences.
- devel shares 159 with v2.0.

This small, stranded, complete-transcript subset has few multi-ORF transcripts, so it does not separate the
v4/v5/v6 selection rules. The full benchmark (especially the unstranded and fragment inputs) is needed for that.
One smoke observation for devel: 28 more ORFs on `rand_*` sequences (52 to 86 codons; 300-transcript training only),
which is worth checking against the full-size decoy results.

Runtime: Predict in v4.1 to v6.0 has a fixed cost of about 75 s on this input, almost all of it PWM start-refinement training
(`train_start_PWM.pl` feature scoring and ROC). On full inputs, expect the ORF-scoring steps to dominate as usual.

## 5. Dependencies, workarounds and caveats for benchmarking

1. **DB_File (Perl) is missing from env/td.** It is needed by every version up to v5.2.0, including all SF releases and v2.0 to v5.0.0 (TiedHash /
   Gene_obj_indexer). Failure: `Can't locate DB_File.pm` in `index_gff3_files_by_isoform.pl`. **Workaround (nothing
   installed):** `work/history/compat/` contains two symlinks into the existing mamba package cache.
   `compat/perl5` points to `env/mamba/pkgs/.../perl-db_file-1.858-pl5321hb9d3cd8_0/lib/perl5/5.32/site_perl`, and `compat/lib/libdb-6.2.so`
   points to `.../libdb-6.2.32-h9c3ff4c_0/lib/libdb-6.2.so`. Prepend them to PERL5LIB and LD_LIBRARY_PATH (the `COMPAT` line above).
   They were built for Perl 5.32, matching env/td. Not needed for v5.3.0 and later.
2. **URI::Escape is not in env/td.** It comes only from `~/perl5/lib/perl5` via the user's login `PERL5LIB`, and it is needed
   by every version before v5.7.0. A runner that resets PERL5LIB, as I did once by mistake, makes Gene_obj.pm fail to load.
   Keep the existing PERL5LIB and prepend to it. For a non-interactive or other-user run, add `~/perl5/lib/perl5` explicitly.
3. **cd-hit-est** is needed by r20131117 through v3.0.1, and is optional in r20131110. env/td has CD-HIT 4.8.1, while the releases bundled 4.5.4/4.6.1.
   The newer version may shift the training set slightly. **ffindex/ParaFly/hmmer** are only for the old built-in Pfam
   search and are not needed when no Pfam search is requested. The SF archives contain prebuilt ELF binaries (ffindex in `3rd_party/`).
   These were **not copied or run**: `sf_runs/` holds only the Perl scripts, and `util/bin` is empty, so `cd-hit-est` resolves to env/td.
4. **The R package seqLogo is missing** (v4.0 to v6.0 start-refinement logos). Its absence is non-fatal: only the logo PDFs are skipped and the predictions are unaffected.
5. **Genetic code flags differ by version.**
   - v5.1.0 and later: `-G Universal` (NCBI names; `Tetrahymena` for the ciliate), as `run_tool.sh` already passes.
   - v2.0 to v5.0.x: `universal` is the default, but an explicit **`-G universal` crashes Predict** (`gff3_file_to_proteins.pl`: "code universal is not recognized"), so **omit -G for universal**. Pass `-G Tetrahymena` to LongOrfs (and to Predict for v3.0.0 and later).
   - **v2.0 Predict has no `-G`**, so its `.pep` is always translated with the universal code. GFF3 coordinates are still correct, because the ORFs come from LongOrfs.
   - The 2012 and 2014 single-step scripts take `-G` directly.
6. **Output locations differ.**
   - r2012: `best_candidates.eclipsed_orfs_removed.gff3` in the current directory. Model IDs are `m.N` and gene IDs are not prefixed by the transcript; seqid is column 1, so the converter is fine.
   - r2013 and r2014: `<input>.transdecoder.gff3`, plus `transdecoder.tmp.<pid>/`.
   - v2.0 to v5.7.1: `<input>.transdecoder.gff3` plus `<input>.transdecoder_dir/` in the current directory.
   - v6.0.0 and devel: the phase scripts are under `util/` (the top-level `TransDecoder` is the wrapper).
   - **`convert_predictions.py transdecoder_gff3` parsed every version's GFF3 without error**, and the row counts equal the mRNA counts. That includes the URI-encoded Name attributes of versions before v5.7 and the CDS-only gene spans of 2012.
7. **Behaviour quirks to report in the history text** (both from code reading):
   - From r20131110 through v4.1.0, `remove_eclipsed_ORFs.pl` is a no-op, because gene spans are whole transcripts. In a v2.0 test on all 233 candidate ORFs, including 26 nested within another ORF on the same transcript, it removed 0.
   - In v5.0.0 to v5.3.0, the default output is effectively a single ORF per transcript.
8. **Strand handling.** All versions support `-S`. The 6-frame rule (up to v4.1) lets an antisense frame veto a forward ORF even under `-S`, because scoring is done on the CDS sequence in all 6 frames. From v5.0 only same-strand frames compete.
9. Set `PERL_HASH_SEED=0`, as the benchmark already does. Several old scripts iterate over hashes, and in v5.0 to v6.0 the
   ranking ties are broken by input order.
10. Worktrees created here are detached at their tags (`git worktree add --detach`): `TD_hist_v2.0`, `TD_hist_v3.0.1`,
   `TD_hist_v4.1.0`, `TD_hist_v5.0.0` and `TD_hist_TransDecoder-v5.5.0`. TD_base (v6.0.0+1), TD_v571 and TD_fix were reused unchanged.
   No repository was modified or committed to, and nothing was installed.
