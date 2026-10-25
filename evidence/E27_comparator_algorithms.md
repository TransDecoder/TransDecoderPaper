# E27. How the benchmarked ORF/CDS predictors work (source-code and literature review, 2026-10-25)

## Question
What statistical model, training regime, ORF-selection rule, start model, partial-CDS handling,
genetic-code support, strand handling, noncoding control and external evidence does each benchmarked
method use? Which of these features plausibly explain the benchmark patterns in
`notebooks/02_benchmark.html`? For the Background paragraph and a supplementary table
(`data/benchmark/comparator_features.tsv`).

## Sources
- **Installed code** (the versions that were benchmarked; paths from `benchmark_v2/scripts/run_tool.sh`).
  `$P = /broad/hptmp/bhaas/TRANSDECODER_PAPER`.
  - CPAT 3.0.5: `$P/env/cpat/lib/python3.12/site-packages/cpmodule/{__init__,find_orfs,FrameKmer,fickett,utils}.py`; prebuilt models `$P/tools/models/cpat/*_logitModel.RData` (formulas and coefficients printed with the env's `Rscript`).
  - CPC2 1.0.1: `$P/env/cpc2/bin/CPC2.py`, `$P/env/cpc2/data/cpc2.{model,range}` (libsvm).
  - Borf 1.2: `$P/env/borf/lib/python3.9/site-packages/borf/{borf,get_orfs}.py`.
  - RNAsamba 0.2.5: `$P/env/rnasamba/lib/python3.7/site-packages/rnasamba/{cli.py,core/*.py}`. The ORF finder and k-mer counter are compiled Rust (`core/orf.cpython-37m-*.so`, `kmer...so`; built from `src/orf.rs`, `src/kmer.rs`, which are not installed). Its CLI prints version 0.2.4 although the package metadata says 0.2.5.
  - CodAn 1.1: `$P/tools/CodAn/bin/{codan.py,predict}`; models `$P/tools/CodAn/models/<CLADE>_{full,partial}/` (ToPS GHMM files: `ghmm/model/ghmm_intronless.model`, `ghmm.<0-9>/model/*.model`, `cnf/*.cfg`, `cnf/meta.cnf`).
  - GeneMarkS-T (`gmst.pl` 5.1, GeneMark.hmm `gmhmmp` 3.20): `~/software/genemarks_t/dist/{gmst.pl,par_1.default,README.GeneMarkS-T}`, `gmhmmp` usage text. (Licence key not read.)
  - Prodigal 2.6.3: binary only (`$P/env/td/bin/prodigal`); usage text read. **Source not inspected**; algorithm from the paper.
  - ANGEL 3.0: `$P/tools/ANGEL/{README.md,src/DumbORF.py,src/SmartORF.py,src/ORFscores.py}`.
  - TransDecoder v7 / PyTD: `manuscript/sections/methods_algorithm.tex` (MA:line below), which documents the release code.
- **Papers** (`$P/literature_pdfs/`, text extracted with `env/latex/bin/pdftotext` into the session scratchpad):
  Hyatt2010_Prodigal, Tang2015_GeneMarkS-T, Nachtigall2021_CodAn, Wang2013_CPAT, Kang2017_CPC2,
  Camargo2020_RNAsamba, Signal2021_Borf (bioRxiv 10.1101/2021.04.12.439551). ANGLE (Shimizu 2006): not read (abstract only, E16).
- **Benchmark numbers**: `data/benchmark/summary.interim_20261025.tsv` (the file the notebook reads) and
  per-transcript tallies from `$P/work/bench_v3/runs/<org>/<variant>/<tool>/preds.tsv` against
  `$P/data/bench_v3/inputs/<org>/<variant>.truth.tsv` (scripts `by_type.py`, `start_codons.py` in the
  session scratchpad; same matching rule as `score_predictions.py`: stop within 2 nt, exact = start
  also within 2 nt). Medians over the 7 organisms other than *T. thermophila* unless stated.

Labels: **[code]** verified in installed source/model files; **[paper]** from the method's paper or
README only; **[run]** observed in our benchmark outputs; **[?]** uncertain/not verified.

## Per-method record

### TransDecoder v7 / PyTD (this paper)
1. Purpose: CDS in eukaryotic transcript assemblies (de novo or genome-guided).
2. Model: three-periodic fifth-order Markov chain (log-odds vs zeroth-order background), six-frame scores; acceptance by frame contrast C = S1 - max(S2,S3) (MA:65-82, 88-98).
3. Training: unsupervised on the input: the T = 500 longest non-redundant ORFs (MA:47-63). Presets from a related transcriptome for small inputs (MA:191-214).
4. Selection: candidates pass homology, a GC-calibrated length L*(GC), or q <= 0.05 against dinucleotide-shuffled decoy ORFs; ranked by homology, S1, length; one ORF per transcript by default, several non-overlapping with `--retain_multiple_orfs` (MA:88-188).
5. Start: ORFs begin at the most upstream in-frame ATG after the previous stop or run off the 5' end; open 5' ends are trimmed to the first ATG when a coding-vs-UTR score with a Kozak log-odds term favours it (MA:23-30, 237-302).
6. Partial: 5'-, 3'- and internal-partial ORFs reported (MA:26-31). Codes: 24 named codes (`-G`), starts ATG only (MA:19-21).
7. Strand: both, or forward only with `-S` (MA:18-19).
8. Noncoding control: target-decoy q-value per stratum (ORF length x transcript length x complete/partial) (MA:115-150).
9. External evidence: optional DIAMOND/BLASTP and Pfam/HMMER hits retain and rank ORFs (MA:169-176).

### Prodigal 2.6.3 (Hyatt et al. 2010)
1. Purpose: protein-coding genes in prokaryotic genomes and metagenomes [paper L18].
2. Model: in-frame hexamer log-ratio table C(w) = log(G(w)/B(w)) (gene training set vs whole-sequence background), "corresponds to a 5th-order Markov model" [paper L237-260]; start score S = 4.25(R + T + 0.4U) + C combining RBS motif (R), start codon type (T), upstream composition (U) and coding (C) [paper L641-650]; genes < 250 bp penalised [paper ~L672-685].
3. Training: unsupervised on the input; training genes from GC frame-plot bias and a first dynamic-programming pass; then hexamer, start-type and RBS statistics [paper L113-262]. Single mode was used; `-n` forces a full motif scan instead of the Shine-Dalgarno trainer [code: usage text; `run_tool.sh`].
4. Selection: dynamic programming over all start/stop nodes ("tiling path"), overlaps up to 60 bp same strand / 200 bp opposite strand, genes with negative score removed [paper L169-232, 701-712]. **Multiple genes per sequence** are reported [paper; run: 8.7% of complete transcripts had >1 prediction].
5. Start: ATG, GTG, TTG nodes scored by RBS/start-type model; within 15 bp the start score alone decides [paper L183-185, 712-716]. **[run]** On standard-code complete inputs every non-edge Prodigal start was ATG (GTG/TTG never chosen); 4-39% (human 38.5%) of predictions on complete transcripts run off the 5' edge with no start codon.
6. Partial: start and stop nodes added at sequence edges, so genes may run off either end unless `-c` (closed ends) [paper L186-189; code: usage]. Codes: any NCBI table via `-g` [code: usage]; run with `-g 6` for *T. thermophila*.
7. Strand: always both strands [paper]; for stranded inputs the benchmark kept + strand calls only (`run_tool.sh`, `--plus-only`).
8. Noncoding control: none beyond the score > 0 rule and the length penalty [paper]; designed for gene-dense genomes.
9. External evidence: none.

### GeneMarkS-T (gmst.pl 5.1; Tang et al. 2015)
1. Purpose: CDS in eukaryotic (default) or prokaryotic (`--prok`) transcripts [code: gmst.pl:200-201; paper].
2. Model: hidden semi-Markov model (GeneMark.hmm); coding = three-periodic Markov chain of order 4 (`--order 4`), non-coding order 2 (gmst.pl:78-79, `--ORDM $order --order_non $order_non` in `train`); gamma (coding) and exponential (non-coding) length distributions [paper L134-137].
3. Training: unsupervised and iterative on the input: initial prediction with the MetaGeneMark heuristic model, then rebuild model / re-predict until 99% identity or 10 iterations (gmst.pl:84-85, `train` sub); transcripts split into up to 3 G+C clusters, one model each (`--bins 0` = automatic; gmst.pl `cluster`). Inputs < 50 kb skip training and use the heuristic model (par_1.default:11; gmst.pl:414-430). Training genes >= 450 nt, partial genes excluded (par_1.default:51,57).
4. Selection: Viterbi parse of each transcript; `--filter 1` (default) passes `-b` = "output one best prediction per sequence" (highest log-odds) (gmst.pl:91, 557; gmhmmp usage). Multiple CDS possible with `--filter 0`.
5. Start: ATG only for code 1 (par_1.default:30-32); a 12-nt positional start-context (Kozak) model learned in the iterations, fixed position, A of ATG at position 7 (gmst.pl:83,87; paper L125-130), used via `gmhmmp -r`.
6. Partial: incomplete genes at sequence ends predicted by default (`gmhmmp -n` turns this off; `--shape partial` default). Codes: **1, 4 and 11 only** (gmst.pl:111, 185, 283). No ciliate option.
7. Strand: both by default; `--strand direct` removes reverse-strand states (strand-specific HSMM) [code gmst.pl:89,178; paper L130-134].
8. Noncoding control: implicit; the Viterbi parse may contain no coding region (non-coding states with order-2 model). No explicit threshold.
9. External evidence: none.

### CodAn 1.1 (Nachtigall et al. 2021)
1. Purpose: CDS (and UTRs) in eukaryotic transcripts, full-length or partial, strand-specific or strand-blind [paper L18, 118-124].
2. Model: generalized HMM (ToPS) with states begin(5'UTR) - start - CDS - stop - end(3'UTR), exactly one CDS per transcript (`ghmm/model/ghmm_intronless.model`; paper L123-124). CDS emission: three-phase inhomogeneous Markov chain with BIC-selected variable context up to order 4 (`cnf/cds.cfg`: InterpolatedPhasedMarkovChain, orders 0-4; `ghmm.5/model/cds.model`: p0,p1,p2, contexts <= 4 nt); UTR: variable-length Markov chain up to order 4 (`intergenic.cfg`); explicit CDS length distribution (`cds_length`, smoothed histogram). Ten GC-specific submodels between 32% and 72% G+C; each transcript uses the band closest to its own G+C (`cnf/meta.cnf`: isochore_min 0.32, max 0.72, nband 10; `bin/predict` `gc_content`, `get_closest_ghmm_id`).
3. Training: **pre-trained**, supervised on RefSeq/WormBase full-length transcripts of species mixes for four clades (VERT, INV, PLANTS, FUNGI); no retraining by users ("contact the authors") [paper L130-161; E16: 5 of our 8 organisms are training species].
4. Selection: Viterbi decoding (`tops-viterbi_decoding` in `bin/predict`); one CDS per transcript by architecture.
5. Start: start state = fixed ATG followed by a 4-nt weight-array "initial pattern" (`ghmm/model/start_composed_motifless.model`, `atg.model`); no upstream (Kozak) motif in the model used, although an upstream start WAM is trained (`start.cfg`, unused by `ghmm_intronless.model`) [code].
6. Partial: separate `_partial` model whose initial state may be CDS in phase 0/1/2 (5'-truncated) (`diff` of full vs partial `ghmm_intronless.model`: only the initial probabilities differ); the full model starts only in the 5'UTR state. [run] Full model: no prediction on complete transcripts lacked an ATG start. Codes: **no option**; stop codons TAA/TAG/TGA are built into the trained models (`cnf/stop.cfg`, `cnf/cds.cfg` RemoveSequence list); the CLI has only -t -m -s -c -o -b -H (codan.py:420-428).
7. Strand: `-s plus|minus|both`; with both, each strand is decoded separately and the higher Viterbi score wins (codan.py:171-184). Note: the help text says "[default=plus]" but the code default is "both" (codan.py:424); the benchmark set `-s` explicitly.
8. Noncoding control: implicit in the Viterbi parse (an all-UTR path is allowed; no score threshold). [run] Median 95.6% of complete coding transcripts received a CDS; random negatives (complete.wRand) received one in 0.01% (P. falciparum) to 51% (C. reinhardtii) of cases, median 5.2% (human 5.7%, D. melanogaster 17%).
9. External evidence: optional `-b` BLASTX of predicted ORFs against a protein DB; it only labels hits, it does not change the CDS (codan.py `_blast_search_`).

### CPAT 3.0.5 (Wang et al. 2013)
1. Purpose: coding-potential classification of transcripts (lncRNA vs mRNA) [paper L7, 26-31].
2. Model: logistic regression. CPAT 3 features: transcript length, ORF length, Fickett TESTCODE score of the ORF, in-frame hexamer log-likelihood ratio (coding vs noncoding hexamer frequencies) averaged over hexamers (utils.py:29 `glm(Label ~ mRNA + ORF + Fickett + Hexamer)`; __init__.py:246-248; FrameKmer.py:82-105). The 2013 paper used ORF coverage instead of transcript length [paper L124-139]. Prebuilt human model coefficients: intercept -7.83, mRNA -5.5e-5, ORF +0.0112/nt, Fickett +3.45, Hexamer +6.07 (Rscript on `Human_logitModel.RData`); mouse and fly are similar in form.
3. Training: **pre-trained** per species (human model: 10,000 RefSeq/CCDS coding + 10,000 GENCODE noncoding [paper L217-222]); prebuilt human, mouse, fly, zebrafish; users can train (`make_hexamer_tab`, `make_logitModel`).
4. Selection: ORF finder takes each ATG after the previous ORF's stop to the next stop (find_orfs.py:55-70); the 5 longest ORFs (`--top-orf 5`) are scored and the one with highest coding probability is reported (`--best-orf p`) (__init__.py:102-133, 340-356). One ORF per transcript.
5. Start: first ATG after the preceding in-frame ORF; no start model (`--start`, default ATG).
6. Partial: 3'-open ORFs allowed (ORF runs to the sequence end, find_orfs.py:59-66); 5'-open ORFs not (ATG required). Min ORF 75 nt (`--min-orf`). Codes: **configurable** start and stop codons (`--start`, `--stop`, __init__.py:79-94); run with `--stop TGA` for *T. thermophila*.
7. Strand: sense only by default; `--antisense` adds the reverse complement (__init__.py:70-78; find_orfs.py:96-105).
8. Noncoding control: coding probability cutoff (published 0.364 human, 0.44 mouse, 0.39 fly; `Human_cutoff.txt`), applied by our converter (`convert_predictions.py cpat`).
9. External evidence: none.

### CPC2 1.0.1 (Kang et al. 2017)
1. Purpose: coding-potential classification, "species-neutral" [paper L19, 36].
2. Model: RBF-kernel SVM (libsvm, `cpc2.model`: c_svc, rbf, gamma 0.125) on four features: Fickett score of the **whole transcript**, peptide length, isoelectric point and ORF integrity of the longest ORF (CPC2.py:274-297; paper L46-53).
3. Training: **pre-trained** on human only: 17,984 coding and 10,452 noncoding transcripts [paper L55-57, 180-181]. No user training.
4. Selection: the single longest ATG-initiated ORF (ties: upstream), reported if the SVM label is coding (CPC2.py:52-141, 333-341).
5. Start: first ATG of the longest ORF; no start model.
6. Partial: 3'-open ORFs allowed but flagged (integrity -1, a classifier feature) (CPC2.py:99-104); no 5'-open ORFs. No minimum ORF length. Codes: **no option**; stops hard-coded (CPC2.py:121), translation by Biopython default table 1 (CPC2.py:244-245); CLI has only -i, -o, -r, --ORF (CPC2.py:26-29).
7. Strand: forward; `-r` also scans the reverse complement, which wins only if its longest ORF is strictly longer (the forward `self.longest` carries over) (CPC2.py:121-141, 111-118).
8. Noncoding control: SVM class label (libsvm `-b 1`; label = higher-probability class) (CPC2.py:319-341).
9. External evidence: none.

### RNAsamba 0.2.5 (Camargo et al. 2020)
1. Purpose: coding-potential classification of RNA sequences, including partial transcripts [paper L24-27].
2. Model: neural network with two branches: whole-sequence branch (nucleotide embedding + two IGLOO layers on up to 3,000 nt) and longest-ORF branch (protein IGLOO on up to 1,000 aa, ORF 2-4-mer frequencies, amino-acid frequencies, ORF length); an ORF-length-driven attention weight combines them; softmax output (core/model.py:209-339; core/inputs.py; paper L147-231).
3. Training: **pre-trained** on human sequences (four published human training sets, non-redundant) [paper L688-690]; separate full-length and partial-length weights (partial trained on complete and truncated human sequences [paper L346-349]); users can train (`rnasamba train`).
4. Selection: the longest ATG-initiated ORF over the three forward frames (stop or sequence end) [paper L170-178; code: compiled `orf.rs`, not inspected]. The protein is written for sequences classified coding (model.py:102-110). One ORF per transcript.
5. Start: first ATG of the longest ORF; Kozak-like signal may be learned implicitly by the network [paper], no explicit start choice.
6. Partial: 3'-open ORFs allowed; no 5'-open ORF (no ATG = no ORF branch input) [paper L173-178]. Codes: **no option** (CLI: output, fasta, weights, -p, -v; cli.py:51-75); stops TAG/TAA/TGA [paper L173]; ORF finder compiled.
7. Strand: forward strand only (no option) [code cli.py; paper].
8. Noncoding control: argmax of softmax (coding score >= 0.5) (model.py:85-89).
9. External evidence: none.

### Borf 1.2 (Signal & Kahlke 2021, bioRxiv; GitHub)
1. Purpose: ORF prediction in de novo assemblies, emphasising strandedness and 5'-incomplete transcripts [preprint abstract].
2. Model: **none**; length only (get_orfs.py). No coding statistic.
3. Training: none.
4. Selection: per frame the longest ATG-to-stop ORF (stop or end), then the longest over frames per sequence (get_orfs.py:46, 91-92, 255-308). Ties in length keep all tied ORFs [code: `transform(max) ==`]. `-a` reports all ORFs >= cutoff.
5. Start: first ATG of the longest ORF; extended to the transcript start (5'-incomplete) when the frame has no stop upstream of the ATG and the upstream segment exceeds 50 aa (`-u 50`) (get_orfs.py:378-428; preprint L446-453).
6. Partial: 5'- and 3'-incomplete ORFs labelled (get_orfs.py:640-667). Min ORF 100 aa (`-l`). Sequences with non-ACGT characters are skipped (get_orfs.py:187-207). Codes: **no option** (CLI borf.py:16-23); translation with scikit-bio `DNA.translate()` default (standard code) (get_orfs.py:194, 201) [default table 1 from scikit-bio docs].
7. Strand: forward by default; `-s` both strands (borf.py:18). The strandedness check described in the preprint was not found in 1.2's code [?].
8. Noncoding control: none other than the 100-aa minimum.
9. External evidence: none.

### ANGEL 3.0 (PacBio; reimplements the ANGLE classifier, Shimizu et al. 2006)
1. Purpose: robust ORF prediction in cDNA/Iso-Seq transcripts with residual errors (README). Deprecated since 2021 (README).
2. Model: "dumb" stage: longest ORF in any frame, scored with a self-trained three-periodic hexamer/pentamer log-odds (fifth-order) model from the top 500 non-redundant longest ORFs (CD-HIT) — TransDecoder's scheme (src/DumbORF.py:106-160). ANGEL stage: AdaBoost (50 stumps) on windowed amino-acid, di-amino-acid and codon-usage features (96-nt windows), followed by dynamic programming over frames (src/SmartORF.py:89-131, src/ORFscores.py:104-170, src/findPath.py).
3. Training: self-trained on the input via the dumb-ORF training set (README workflow; `run_tool.sh angel`).
4. Selection: `--output_mode best` keeps the longest of the ANGEL and dumb ORFs (SmartORF.py:139-238).
5. Start: furthest in-frame ATG without an intervening stop (SmartORF/ORFscores `extend_in_frame`) [code, not traced in detail].
6. Partial: 5'/3'-partial ORFs (DumbORF.py:30-56). Codes: **standard only** (ORFscores.py:100-101 hard-coded; Biopython standard table).
7. Strand: forward; `--use_rev_strand` adds reverse.
8. Noncoding control: classifier score plus minimum lengths (50 aa ANGEL, 100 aa dumb).
9. External evidence: none. Dropped from the benchmark (17-22 h per *A. thaliana* run; E07).

## Genetic-code capability (checked from source, for the n/a treatment of *T. thermophila*)
| Method | Code/stop option | Evidence |
|---|---|---|
| TransDecoder | yes, `-G` (24 codes) | MA:19-21 |
| Prodigal | yes, `-g` NCBI table | usage text |
| CPAT | yes, `--start`/`--stop` lists | cpmodule/__init__.py:79-94 |
| GeneMarkS-T | codes 1, 4, 11 only | gmst.pl:111, 185, 283 |
| CodAn | none; stops built into trained clade models | codan.py:420-428; cnf/stop.cfg, cnf/cds.cfg |
| CPC2 | none; hard-coded stops, Biopython table 1 | CPC2.py:26-29, 121, 244-245 |
| RNAsamba | none (CLI); stops fixed in compiled ORF finder | cli.py:51-75; paper L173 |
| Borf | none; scikit-bio standard translation | borf.py:16-23; get_orfs.py:194, 201 |
| ANGEL | none; hard-coded | ORFscores.py:100-101 |

**No tool marked n/a for *T. thermophila* has a user-settable genetic code.** Only code changes
(CPC2, Borf, ANGEL: one line; RNAsamba: Rust rebuild) or retraining (CodAn) would make them usable.
GeneMarkS-T code 4 (TGA = Trp) is the opposite of the ciliate code (TGA only stop).

## Claims supported (for Background / supplementary table)
1. **Common features.** Every method enumerates ORFs from in-frame start/stop codons, and all except Borf score coding sequence with codon-position statistics: periodic Markov chains (TransDecoder order 5, Prodigal hexamers, GeneMarkS-T order 4, CodAn up to order 4), in-frame hexamer log-likelihoods (CPAT, ANGEL), Fickett TESTCODE (CPAT, CPC2) or learned k-mer/sequence representations (RNAsamba).
2. **Discriminating features.**
   - Training: self-trained on the input (TransDecoder, Prodigal, GeneMarkS-T, ANGEL); pre-trained supervised models (CodAn by clade; CPAT, CPC2, RNAsamba on human or a few model species); none (Borf).
   - Model type: generative gene-structure models with a Viterbi/DP parse (Prodigal, GeneMarkS-T, CodAn) versus discriminative classifiers of one candidate ORF (CPAT logistic regression, CPC2 SVM, RNAsamba neural network) versus scoring and ranking all candidate ORFs (TransDecoder).
   - ORF choice: longest ATG ORF (CPC2, RNAsamba, Borf), best of the 5 longest (CPAT), parse (Prodigal, GeneMarkS-T, CodAn), rank by homology/score/length after a decoy test (TransDecoder). Multiple CDS per transcript: Prodigal by design; GeneMarkS-T and TransDecoder optional.
   - Start: first/most-upstream ATG (CPAT, CPC2, RNAsamba, Borf, TransDecoder before refinement); RBS/start-type model (Prodigal; ATG/GTG/TTG); learned positional Kozak model (GeneMarkS-T, TransDecoder's refinement term); ATG + downstream 4-nt pattern (CodAn).
   - Partial CDS: 5'-open ORFs only in TransDecoder, Prodigal, GeneMarkS-T, CodAn (partial model), Borf (> 50 aa rule) and ANGEL; the three classifiers need an ATG.
   - Noncoding control: decoy q-value (TransDecoder); probability/label cutoff (CPAT, CPC2, RNAsamba); implicit in the parse (GeneMarkS-T, CodAn); score > 0 with length penalty (Prodigal); length only (Borf).
   - Genetic code: table above.
   - External evidence: TransDecoder (selection); CodAn (annotation only).
3. **Explaining benchmark patterns** (medians over 7 organisms, interim summary; mechanism statements are plausible, not tested by ablation):
   - *Classifiers score well on complete transcripts at the exact level but poorly on fragments.* Complete exact F1: CPAT 0.945, CPC2 0.935, RNAsamba 0.961, Borf 0.923 (TransDecoder 0.916). On a complete MANE-like transcript the annotated CDS is usually the longest ORF from its first ATG, which these tools report. On fragments, exact F1 falls to 0.447, 0.183, 0.486, 0.538 (TransDecoder 0.825). Split by fragment type [run]: where the start is missing (5'-partial truth), exact sensitivity is 0.015 (CPAT), 0.015 (CPC2), 0.018 (RNAsamba), because they must start at an internal ATG; Borf reaches 0.262 through its 50-aa upstream rule; TransDecoder 0.744, CodAn 0.910, Prodigal 0.938. Where the stop is missing (3'-partial truth), CPC2 calls only 24.9% of transcripts (stop sensitivity 0.247), consistent with ORF integrity being an SVM feature; CPAT 0.792, RNAsamba 0.907.
   - *Prodigal: high sensitivity, low precision and poor starts.* Stop sensitivity is the highest on complete (0.981) and fragment (0.951) inputs, because genes may run off both ends and any ORF with a positive score is called. Precision with random negatives is 0.664 (fragments 0.703): there is no noncoding model beyond score > 0, and multiple genes per transcript are reported (8.7% of complete transcripts had > 1). Exact F1 on complete transcripts is 0.653: 4-39% of its predictions on complete transcripts extend through the 5' UTR to the transcript edge (human 38.5%), as expected from an RBS-based start model applied to eukaryotic transcripts; on stop-missing fragments its exact sensitivity is 0.661.
   - *GeneMarkS-T and CodAn starts.* GeneMarkS-T 13.6% (human) of complete-transcript predictions extend to the 5' edge; exact F1 0.783. CodAn's full model must begin at ATG (0% edge starts) and has exact F1 0.830.
   - *Length-driven methods lose precision on random sequence, most at high G+C.* Median precision with GC-matched random negatives: longest-ORF 0.823, Borf 0.838, RNAsamba 0.721, against TransDecoder 0.994, GeneMarkS-T 0.996, CPAT 0.992, CPC2 0.982. Per organism the drop is largest where random sequence holds long ORFs (C. reinhardtii 68% G+C: Borf 0.582, CPAT 0.629, CodAn 0.652, RNAsamba 0.695, longest ORF 0.567; TransDecoder 0.992) and smallest in P. falciparum (22% G+C, all >= 0.993 except Prodigal 0.866). CPAT's human model weights ORF length heavily (+0.0112 logit per nt), so long random ORFs pass; human precision 0.854.
   - *Pre-trained models outside their training range.* CodAn on P. falciparum (INV model; transcripts at 22% G+C lie below its lowest GC band, 32%): stop sensitivity 0.552 (complete) and 0.417 (fragments), versus 0.94-0.97 elsewhere; plausible cause is the band/clade mismatch [?]. RNAsamba (human weights) fragment stop sensitivity 0.587 in P. falciparum vs 0.78-0.84 elsewhere. CPAT with the human model: lowest complete exact F1 in C. reinhardtii (0.786) — but every method, including self-trained ones, is lowest there (TransDecoder 0.747), so this is not specific to pre-training.
   - *Genetic code.* Five methods (GeneMarkS-T, CodAn, CPC2, RNAsamba, Borf) have no ciliate-code option and are not applicable to T. thermophila; their T. thermophila runs are not results.
4. **Caveats.** Prodigal's source was not inspected; its "only ATG on table 1" behaviour is observed, not read in code. RNAsamba's ORF finder is compiled Rust; behaviour is from the paper. CodAn rows: 62/80 runs in the interim summary (E07); medians computed over available organisms (all 7 non-ciliate organisms present for complete, fragments and complete.wRand). The per-type tallies used a simplified re-implementation of the scorer (no one-credit-per-truth bookkeeping across duplicate predictions), so they are for mechanism, not for quoting as headline numbers.

## New references needed (not in `manuscript/references.bib`)
- Signal B, Kahlke T. Borf: improved ORF prediction in de-novo assembled transcriptome annotation. bioRxiv 2021. doi:10.1101/2021.04.12.439551.
- Fickett JW. Recognition of protein coding regions in DNA sequences. Nucleic Acids Res. 1982;10:5303-5318. doi:10.1093/nar/10.17.5303.
- Lukashin AV, Borodovsky M. GeneMark.hmm: new solutions for gene finding. Nucleic Acids Res. 1998;26:1107-1115. doi:10.1093/nar/26.4.1107.
- Besemer J, Lomsadze A, Borodovsky M. GeneMarkS: a self-training method for prediction of gene starts in microbial genomes. Nucleic Acids Res. 2001;29:2607-2618. doi:10.1093/nar/29.12.2607.
(Check each against the paper before citing, per E14.)

## Files
- `data/benchmark/comparator_features.tsv` (one row per method).
- Draft Background paragraph and supplementary table: session scratchpad `background_comparators.tex` (not in the repo).
- Analysis scripts (scratchpad, not in git): `scripts/by_type.py` (sensitivity by fragment type), `scripts/start_codons.py` (start codons and 5'-edge starts of predictions).
