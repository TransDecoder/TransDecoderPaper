# True CDSs that fail TransDecoder's coding criteria

Brian's question (2026-10-21): real long ORFs sometimes fail the coding criteria. How often, and why?

Full tables, commands and checks are in `evidence/E12_long_orf_failures.md`; every number below comes from `results/*.tsv`.

## Setup
- **Data and code:** bench_v3 `complete` inputs (8 organisms); PyTD final defaults (`pytd_dev_v2` runs, code 48c870c); GM12878 StringTie/LRAA `final_abinitio`.
- **Reproducibility:** `markov.py` re-implements PyTD's 5th-order Markov scoring and reproduces `hexamer.scores` and `longest_orfs.cds.scores` exactly.
- **Rule 3 (coding test), no homology:** S1 > 0, and the frame contrast S1 − max(S2, S3) > 0. For CDS < 300 nt, also contrast ≥ 200·μ (the short-ORF guard).
- **Rule 2 (long-ORF rule):** ORF length ≥ L\*(GC). It is off for TGA-only codes (*T. thermophila*).

## How often (shown)
- **Overall:** 4.9% of 135,429 matched true CDSs fail rule 3.
  - 2,347 fail through the short-ORF guard, 3,091 because S1 ≤ 0, and 1,189 because a shifted frame scores ≥ S1.
- **ORFs ≥ 300 nt:** 2.38% fail.
  - By organism: *T. thermophila* 6.2%, fly 2.8%, *C. reinhardtii* 2.1%, mouse 1.7%; human, *S. pombe*, *A. thaliana* ≤ 0.9%; *P. falciparum* 0.1%.
- **By length:** <300 nt 67%, 300–599 nt 9.4%, 600–999 nt 2.3%, 1–2 kb 0.6%, ≥2 kb 0.3–0.4%.
- **Long examples:** 398 true CDSs ≥ 1 kb fail, e.g. Muc14A (48.7 kb), kl-2, kl-3, Nipped-A, MUC4, FLG, HRNR, FLG2, COL4A4, ELN.
- **The long-ORF rule:** it rescues 447 failures (397 end up predicted). Above L\* every failing CDS is rescued, except in *T. thermophila*, where the rule is off and 76 failing CDSs ≥ 1 kb are lost. Most failures are shorter than L\*: 6,180 true CDSs (4.6%) are still lost.
- **Predictions kept by length alone:** 0–1.4% of all predictions; 398 of 402 are true on complete inputs.
  - On negatives they make as many or more false calls as true ones: human lncRNA 47 true vs 99 false, human random 40 vs 82, *C. reinhardtii* random 8 vs 58.
  - In GM12878 StringTie: 116 distinct CDSs; 11 MANE-exact, 35 GENCODE-compatible, 46 with no GENCODE overlap.
- **GM12878 MANE assembled but not predicted exactly** (651 StringTie, 412 LRAA):
  - 62–68% have the correct stop but a different start (a start-codon problem).
  - 17–30% fail rule 3, mostly through the short-ORF guard.
  - Only 1.4–2.1% of all assembled MANE CDSs fail rule 3.

## Why (shown)
1. **Composition heterogeneity within a genome: a minority codon-usage class.** One Markov model and one background cannot fit every gene.
   - Failures concentrate in a GC3 extreme: fly lowest GC3 quintile 13.5% vs 0.04–0.3% in the others; fly chr4 40/80, chrY 8/22; *C. reinhardtii* lowest quintile 8.4%; *A. thaliana* highest quintile 2.6%.
   - Highly expressed, codon-biased *T. thermophila* ribosomal proteins fail at 21/60.
   - **Decisive test:** rescoring with a model trained on 500 CDSs from the same GC3 quintile rescues 25% of failures (fly 69%, mouse 52%, human 36%; 17 of the 21 Tetrahymena ribosomal proteins).
   - Models trained on all annotated CDSs, 500 random ones or the 500 longest give failure rates of 2.08–2.49%, vs 2.38% for TD's own model. TD's training set is 99.8–100% true CDSs.
   - So the training set is not the problem; a single model is.
2. **Low-complexity and repetitive proteins: the main cause in mammals.** A shifted frame scores at least as well as the true frame (human 134/162, mouse 318/358).
   - 78% of human and mouse failures are low-complexity or repetitive, vs 13% of controls.
   - Families: KRTAP (human 52/72, mouse 78/92), LCE, SRSF, histones, KRAB-ZNF arrays, FLG/HRNR, mucins, collagen/elastin.
   - Failing human CDSs: Cys 14% vs 2.5%; SEG-like low-complexity fraction 0.42 vs 0.11.
3. **Unconserved, partly dubious gene models.** Mostly *T. thermophila* "hypothetical proteins", plus *S. pombe*, *C. reinhardtii* and mouse Gm genes.
   - 97% of Tetrahymena failures have no Swiss-Prot or Pfam hit (controls 58%).
   - Human failure rate is 14.7% without a Swiss-Prot hit vs 0.4% with one.
   - No truth-trained model rescues most of them. A quarter of Tetrahymena failures have a stronger ORF in another frame, which may be an annotation frame error (speculative).

Exclusive attribution of the 3,114 failures ≥ 300 nt: low complexity 32%, GC3/codon class 33%, stronger out-of-frame ORF 6%, no homology 26%, unexplained 2.4%.

## Not supported
- **Organellar genes:** only 7 in the truth sets, and all pass.
- **Long UTRs:** they affect ranking against other ORFs, not the coding test.
- **Extension into the 5′ UTR:** rescoring only the annotated CDS gives 5.16% vs 4.89%.
- **Horizontal transfer:** not testable here; mechanism 1 is how such genes would fail.

## Implications (speculative, untested)
- The long-ORF rule rescues the long members of classes 1–2. It cannot reach most failures, which are shorter than L\*, and it is off in Tetrahymena.
- Candidate improvements:
  - a mixture of GC3-class Markov models, or a GC3-conditioned background;
  - homology evidence for repeat proteins.
- **Dropping S1 > 0 is not a fix** (exploratory): it would add only 0.06–2% of true CDSs, while transcripts on random sequence with a candidate would rise 10–40×.
- Reduced codon selection on fly chromosomes 4 and Y is a known literature theme; not tested here.

## Caveats
- **Truth sets** are annotation-derived, and unconserved models may not be real proteins.
- **Homology flags** come from the older `pytd_homology` runs (−m 100), so there are none for CDSs < 100 aa.
- **Low-complexity measures** are proxies.
- **Attribution** depends on the order the flags are applied.
- **Gene split:** the truth-trained models split genes two ways by gene, not by family.
- **Negative-set counts** are taken before single-best selection, so they are upper bounds.
- **Run commits:** the `pytd_dev_v2` runs predate the 48c870c commit; see E10.
