# Inventory of TransDecoder org repos (surveyed 2026-10-17)

| Repo | Last push | Contents |
|---|---|---|
| TransDecoder | 2026-03 (v6.0.0) | Perl source, sample_data, `__testing` |
| PyTransDecoder | 2026-03 | Python port v0.1.0. Ported with Claude Sonnet 4.5 under bhaas's guidance, Jan–Mar 2026 |
| TransDecoder_benchmarking_data | 2018-10 | 2017–18 accuracy benchmark (see below) |
| TransDecoder.lrgTests | 2018-04 | 3 large genome-propagation runtime examples (Cufflinks de novo, mm9 ref, Cuffmerge). No truth sets |
| transdecoder.github.io | 2017 | Redirects to the wiki |
| cdhit | 2017 | Vendored mirror; unused since v4.1 |

Local clones are in `/broad/hptmp/bhaas/TRANSDECODER_PAPER/repos/`.

## TransDecoder_benchmarking_data (2017–2018)

**Truth sets**
- Built from RefSeq, one transcript per gene: the longest CDS exactly contained in an mRNA, with UTRs trimmed to ≤ 250 bp.
- Only full-length CDS that start with M, end with a stop, and have no internal stops are kept.

| Organism | Transcripts | GC% |
|---|---|---|
| A. thaliana (TAIR10) | 27,328 | 42 |
| D. melanogaster (r6) | 13,539 | 51 |
| M. musculus (GRCm38.p5) | 21,825 | 52 |
| S. pombe (ASM294v2) | 5,118 | 38 |

**Fragments**
- Each transcript is cut once inside the CDS, at a position drawn from N(0.5, 0.1).
- The kept piece is 5′ (start present, stop lost) or 3′ (stop present, start lost).

**Negatives**
- Each transcript has a GC-matched random sequence of the same length.
- Any prediction on a random sequence counts as an FP.

**Tools compared**
- GeneMarkS-T (± strand-specific)
- Prodigal 2.6.3 (± strand filtering)
- TransDecoder (circa v5.x, ± SS, ± refine-starts, ± single-best)
- LongOrfs (all ORFs)
- ANGEL and CodAn are NOT included.

**Metrics**
- TP if the stop end is within 2 bp of the truth (3′ mode), or both ends are within 2 bp (strict 3′+5′ mode).
- The ROC sweeps the minimum ORF length from 750 to 300 bp.
- The axis labelled "specificity" is actually precision.
- AUCs are all about 0.99, so they discriminate poorly between methods.

**Selected AUCs (3′ mode)**

| Dataset | AUCs |
|---|---|
| Mouse complete | GMST 0.997, Prodigal 0.995, TD 0.995 |
| Mouse fragments (TDvary) | TD 0.986, GMST 0.983, Prodigal 0.970, longorfs 0.863 |
| Dmel complete | TD 0.994, Prodigal 0.993, GMST 0.989 |

**`__rpart/`**
- Decision-tree analysis of ORF-selection features (length, 6-frame scores).
- The generic rule "S1 > 0 & max + length rescue" beat the per-organism trees.
- **This is the empirical justification for the v5 selection rule.** The PDFs are in `complete/*/__custom_filters/`.

**Repro debt**
- The scripts are Python 2.
- Paths to Broad-specific tools are hard-coded.
- The GeneMarkS-T license is needed.
- The deleted `longest_orfs.*` files must be regenerated.
- Some analysis PDFs are recoverable via `git show da05542^:<path>`.

## PyTransDecoder
- About 3.9k lines of Python. Pure Python plus biopython and tqdm.
- Single-command driver: `pyTransdecoder` (diamond/hmmsearch built in).
- Perl parity, measured before Python got PWM start refinement:

  | Dataset | Result |
  |---|---|
  | Trinity sample, LongOrfs | 845/845 identical |
  | Trinity sample, Predict | Perl 679 vs Python 682 |
  | PASA, Cufflinks | 100% |

- The undocumented Perl strand-agnostic overlap rule was discovered during the port.
- **Parity and speed have not been re-measured since start refinement was added (Mar 2026).** The "10× faster" claim is not a fair comparison.
- Tests: 14 unit tests.
