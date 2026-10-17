# TransDecoder v6.0.0 — algorithm reconstruction (from code)

Source: github.com/TransDecoder/TransDecoder @ 66d4712 (v6.0.0, 2026-03-20). Paths are relative to the repo root.
Built by reading the code on 2026-10-17. Use this as raw material for the Methods section.

## 0. Pipeline layout
- `TransDecoder` wrapper (new in v6) runs these phases:
  - 0: GTF+genome → cDNA (genome mode)
  - 1: `util/TransDecoder.LongOrfs`
  - 1.5: optional DIAMOND/BLASTP + hmmsearch
  - 2: `util/TransDecoder.Predict`
  - 3: genome propagation
- Checkpointing: `PerlLib/Pipeliner.pm`.
- Model inspiration: GeneID (Parra et al. 2000, Genome Res 10:511), cited in `notes`.

## 1. LongOrfs (`util/TransDecoder.LongOrfs`, `PerlLib/Longest_orf.pm`)

**Background base composition** (`util/compute_base_probs.pl`):
- Counts are pooled over all transcripts.
- Without `-S`, counts are complement-symmetrized, so P(A)=P(T) and P(C)=P(G).
- Output: `base_freqs.dat`.

**ORF enumeration**:
- Scans the + strand, plus the reverse complement unless `-S`.
- Stops come from the active genetic code.
- Starts are literal ATG only, for every genetic code.
- 3′-partial pseudo-stops are placed at L, L-1, L-2. 5′-partial pseudo-starts are placed at 0, 1, 2 (unless a stop is there).
- Each stop is paired with its most-upstream valid start, so there is one (longest) ORF per stop.
- `--complete_orfs_only` turns off both partial modes.

**Minimum length**: `-m 100` codons. This count includes the stop codon (or a phantom codon for 3′ partials).

**ORF types**:
- complete
- 5prime_partial
- 3prime_partial
- internal

**Genetic codes** (`PerlLib/Nuc_translator.pm`):
- Universal.
- Ciliate-type variants: Tetrahymena, Ciliate, Dasycladacean, Hexamita, Acetabularia, Candida, Euplotid, SR1_Gracilibacteria, Pachysolen, Mesodinium, Peritrich.
- 11+ mitochondrial codes.
- The genetic code affects stops and translation only.

## 2. Predict (`util/TransDecoder.Predict`)

### 2a. Training set
1. Take the 10·T longest CDS. Skip any with CDS longer than 5000 nt.
2. Redundancy removal (`exclude_similar_proteins.pl`; CD-HIT was dropped in v4.1):
   - Proteins are processed longest first.
   - Query words are the non-overlapping (tiled) 5-mers of the protein. Retained proteins are indexed by all overlapping 5-mers.
   - Skip a protein if its distinct tiled 5-mers are fewer than 30% of len/5 (low complexity).
   - Proxy identity = the maximum, over retained proteins, of shared tiled 5-mers ÷ the query's distinct tiled 5-mers.
   - Keep the protein if proxy identity ≤ 80%.
3. Take the top T. The default is T=500.

Training ORFs may be partial and from either strand.

### 2b. Coding model (`util/seq_n_baseprobs_to_loglikelihood_vals.pl`, `util/score_CDS_likelihood_all_6_frames.pl`)
- 5th-order, 3-periodic inhomogeneous Markov chain, compared against a 0th-order background q.
- Transition probability:
  P(x_i | x_{i-k..i-1}, f) = (N(x_{i-k..i}, f) + 1) / (N(prefix, (f-1) mod 3) + 4),  with f = i mod 3
- Per-position log-odds: ℓ = ln[ P / q(x_i) ]
- ORF score: S(s) = Σ_i ℓ(x_{i-m..i}, i mod 3), with m = min(i, 5)
- The terminal stop codon is excluded. These stops are hard-coded universal ones (TAA/TAG/TGA).
- Six scores are computed: S1..S3 for the forward frame offsets 0/1/2, and S4..S6 for the reverse complement.

### 2c. Dynamic long-ORF threshold
- GC% = 2·P(C)·100.
- Lookup table, giving ORF length in nt at the 0.999 quantile of random-sequence ORF lengths:

  | GC% | 25 | 30 | 35 | 40 | 45 | 50 | 55 | 60 | 65 | 70 | 75 | 80 |
  |---|---|---|---|---|---|---|---|---|---|---|---|---|
  | ORF length (nt) | 465 | 510 | 555 | 590 | 645 | 749 | 797 | 927 | 1086 | 1358 | 1743 | 2422 |

  Above 80% GC the threshold is disabled.
- **The simulation code is NOT in any repo. It must be regenerated or described.**
- The usage text says "1% FDR", but the 0.999 quantile is what the code uses.

### 2d. Selection (`util/select_best_ORFs_per_transcript.pl`)
- An ORF passes if any of these holds:
  - homology_count > 0
  - CDS length ≥ the long-ORF threshold
  - S1 > 0 AND S1 > max(S2, S3). Only forward frames are compared, since v5.0.0.
- Ranking: homology count (descending), then S1, then length.
- Greedy overlap removal:
  - Reject an ORF if it overlaps an accepted ORF by more than 10% of either ORF's span.
  - **This check is strand-agnostic.**
  - Homology-supported ORFs are always kept.
- `--single_best_only`: rank by homology, then length (S1 is ignored), and keep the top ORF.
- BLAST outfmt6 / Pfam domtblout: every line counts. No E-value filter is applied inside Predict.

### 2e. Start refinement (on by default)

**PWM** (`util/PWM/*`, `PerlLib/PWM.pm`):
- Window: 20 nt upstream + ATG + 10 nt downstream.
- Positive features: training-ORF starts.
- Negative features: downstream ATGs in any frame.
- Pseudocount: 0.1.

**Sub-window selection**:
- 5 rounds of 75/25 cross-validation.
- All 200 (u,d) sub-windows are scored with length-normalized log-odds (+ vs −).
- 10-threshold ROC; AUC computed per (u,d).

**Noise depletion**:
- Start from 30% of positives and swap in better-scoring features.
- Drop members with score ≤ 0.
- Re-run the sub-window selection on the enhanced positive set.

**Application** (`util/start_codon_refinement.pl`):
- Applies only to ORFs whose first codon is not ATG.
- Candidates: in-frame ATGs within max(30 nt, 15% of ORF length) downstream.
- Use the best-AUC (u,d) and the best-F1 threshold. Pick the highest-scoring candidate that is at or above the threshold.
- Relabel the ORF type: 5′partial → complete, internal → 3′partial.

**Seeds**: `srand(1)` in PWM.pm and `srand(1234)` elsewhere.

### 2f. Outputs
- `.transdecoder.{gff3, bed, pep, cds}`.
- The GFF3 includes UTRs.
- Name attributes carry the score and the BLAST/Pfam hits.

## 3. Genome propagation
- `gtf_genome_to_cdna_fasta.pl`, `gtf_to_alignment_gff3.pl`, `cdna_alignment_orf_to_genome_orf.pl`.
- An ORF on the − strand of a multi-exon transcript is dropped.
- An ORF on the − strand of a single-exon transcript flips the transcript strand (v5.6.0).
- Inputs supported: PASA, GMAP, Cufflinks, StringTie.

## 4. v6 wrapper
- `--blast_search_pep`: DIAMOND (default) or blastp, run with `-k 1` and evalue 1e-5.
- `--pfam-search-db`: hmmsearch with **no E-value/GA cutoff**.

## 5. Doc/code discrepancies and possible bugs (decide what to fix before the paper)
1. The wiki says the score is compared against "other 5 frames"; the code compares against frames 2–3 only.
2. The wiki says an "encapsulated ORF → longer one" is kept; the code uses 10% overlap with a homology > score > length ranking, regardless of strand.
3. The wiki file names are stale (e.g., `.scores.selected`).
4. "1% FDR" vs the 0.999 quantile, and the simulation is missing.
5. `max_prot_length_for_training=5000` is applied in nt (≈1666 aa).
6. **Start refinement:**
   - The threshold is learned on length-normalized scores but applied to unnormalized sums.
   - The scoring window is shifted 1 nt upstream relative to training.
   - **These are likely bugs.**
7. `start_codon_refinement.pl`: the usage text says `--adj_pct` but the option is `atg_pct`, and the default ordering for `--workdir` is wrong.
8. Euplotes vs Euplotid naming.
9. Universal stops are hard-coded in Markov train/score.
10. With `-S`, GC% = 2·P(C) is approximate.
11. No homology E-value thresholds are applied in Predict, and the wrapper's hmmsearch has no cutoff, so any reported domain counts.
12. Unused options: `--cpu`, `--search_pfam`, `--reuse`, `--debug`.
13. In `cdna_alignment_orf_to_genome_orf.pl:88`, the code splits on `::` but LongOrfs uses `~~`.
14. Stale `__testing` fixtures.
15. The Docker wiki uses old paths.
16. Many unused utilities remain in `util/`.
17. `--single_best_only` ignores the Markov score.

## 6. Follow-up findings (2026-10-17)

**Long-ORF threshold table origin**
- First appears in TransDecoder commit `33247f5` (2017-08-24). The full quantile table (0.95, 0.99, 0.999, 0.9999) is pasted as a comment.
- No script that generated it was found in the git history of TransDecoder, TransDecoder_benchmarking_data, PyTransDecoder or lrgTests.
- It is also absent from ~/GITHUB/HAASUTIL, ~/utilities and ~/bin.
- Plan: regenerate it with a documented, seeded simulation, and confirm the result reproduces the 2017 numbers.

**PyTransDecoder start refinement**
- Ported on 2026-03-15.
- It reproduces both Perl bugs:
  - Training scores are normalized by (u+d) at `predict.py:1174`.
  - Application uses the raw sum over the window shifted by −1 at `predict.py:1354-1364`.
- Fix both implementations together.
