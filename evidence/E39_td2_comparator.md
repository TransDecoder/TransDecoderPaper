# E39. TD2 as a comparator method (2026-10-31)

Decision (Brian Haas, 2026-10-31): include TD2 in the manuscript, cited and compared, reversing the earlier private-only rule. **Co-authorship caveat:** Brian Haas, sole author of this paper, is a co-author of TD2 (preprint CRediT: software, validation, review). The manuscript discloses this (Declarations). Runs are ab initio only, like the other comparators; TD2's homology options were not run (Brian, 2026-10-31).

## Claims supported
Table 1 (TD2 row); Fig. 2; Fig. 4d; Supplementary Fig. S4 (ranks); Supplementary Tables S1, S3, S4, S5; Supplementary Figs S5-S7; Results "Accuracy across eight eukaryotes" (ranks, noncoding, starts, fragments, balanced types, ciliate code, simulated assemblies), short-ORF lncRNA paragraph; Background comparator paragraph; Discussion, "A self-trained coding model"; Methods, "Compared methods".

## Software
- `td2` 1.1.0, latest on PyPI (`$P/env/td2`; console scripts `TD2.LongOrfs`, `TD2.Predict`), PSAURON 1.1.3, recorded per run in `td2_version.txt`.
- Preprint: Mao, Ji, Haas, Salzberg, Sommer, bioRxiv 2025, doi:10.1101/2025.04.13.648579 (paperclip bio_ab7dbef7775b; authors and title checked on Crossref). PSAURON: Sommer, Zimin, Salzberg, NAR Genom Bioinform 7(1):lqae189, 2025 (Crossref). PSAURON is a neural network: a temporal convolutional network (six layers, 383,457 parameters) over amino-acid sequence, trained on CDSs of 1,090 RefSeq plant and animal genomes (paper, PMC11704789 L17, L23, L26).
- Relation to TransDecoder: TD2 keeps TransDecoder's two-step structure and output format (`TD2/Predict.py` comments "TransDecoder style out-of-frame", "match TransDecoder output format"). When the preprint appeared (April 2025) the current TransDecoder release was v5.7.1 (2023-07-16); v6.0.0 came on 2026-03-20 and predicts identically (E23). Algorithm details: E27, TD2 section.

## Commands
```bash
cd $P/work/bench_v3
bash launch_td2_v1.sh    # copy: benchmark_v2/launch_v3/launch_td2_v1.sh; jobs td2_jobs_v1.tsv (120: input, outdir, org; largest first)
#   xargs -P 4 -> benchmark_v2/scripts/run_td2.sh <input.fa> <outdir> <org>
#   TD2.LongOrfs -t input.fa [-S] -G <NCBI table> -@ 2 ; TD2.Predict -t input.fa -G <NCBI table>
bash $T/pipeline/derive/10_td2_tables.sh      # score, collect; log runs/logs/td2_tables_v1.log
python $T/pipeline/derive/10b_td2_mechanisms.py
```
- `-S` on stranded inputs (names without `unstranded`). `-G 6` for *T. thermophila*, otherwise 1. Defaults otherwise: minimum 90 aa, PSAURON cutoff 0.5, dynamic long-ORF rule at FDR 0.10.
- One ORF per transcript: TD2.Predict's default ("Selecting single best ORF per transcript"); `run_td2.sh` fails if any transcript has more than one prediction. None did, including unstranded inputs.
- Threads: 2 for LongOrfs, `OMP_NUM_THREADS=MKL_NUM_THREADS=2` for PSAURON (it otherwise takes all cores). Other tools mostly used one thread, so TD2's times are not strictly comparable.
- `run_td2.sh` is a separate script so `run_tool.sh` was not edited. Fresh runs, no resume; `work/td2_private` was not reused.

## Inputs and outcome
120 of 120 finished with exit 0 (2026-10-31 12:01-15:05; 4 jobs, shared host, load 26-40):
84 main inputs (8 organisms incl. *T. thermophila* with table 6, complete, complete_fullutr, fragments, +random, +lncRNA, stranded and unstranded), 24 mixed t1-t3, 6 human and 6 *A. thaliana* simulated depths (d1-d32). Run directories: `$P/work/bench_v3/runs/<org>/<variant>/td2/`, `$P/work/simtrinity/runs/td2`, `$P/work/simtrinity/depth/<level>/runs/td2`, `$P/work/simtrinity_athal/depth/<level>/runs/td2`. No input failed.

## Outputs
`final_v7/all_scores.tsv` (+84 rows), `mixed_scores.tsv`, `partial_types.tsv`, `five_prime_by_distance.tsv`, `mixed_medians.tsv`, `table1.tsv`, `lnc_by_era.human.v7.tsv`, `td2_fragments_by_cds_length.tsv`, `td2_random_calls_by_rule.tsv`; `data/simtrinity/*.tsv`, `data/simtrinity/depth/`, `data/simtrinity_athal/depth/`; `data/supplementary_tables/benchmark_per_input.tsv` (1,315 rows); `data/benchmark/comparator_features.tsv`, `comparator_table.tex`.

## Results
Table 1 (median over 8 organisms; TD2 supports the ciliate code):

| method | complete stop / exact F1 | fragments stop / exact F1 | mixed stop / exact F1 | complete + random precision |
|---|---|---|---|---|
| TransDecoder | 0.980 / 0.912 | 0.930 / 0.822 | 0.970 / 0.843 | 0.994 |
| TD2 | 0.982 / 0.941 | 0.876 / 0.421 | 0.929 / 0.673 | 0.904 |

- Exact F1 on complete inputs is higher than TransDecoder in 7 of 8 organisms (*P. falciparum* 0.939 vs 0.940).
- Ranks (210 categories; 9 methods): TD2 first in 5, top three in 73, median rank 4, median gap 0.029 (second smallest). TransDecoder first 49 (was 50), top three 126 (was 134), median gap 0.015.
- Mixed types, exact F1: complete 0.940, 3' partial 0.644, 5' partial 0.182, internal 0.903.
- Fragments by type (median): 3' partial stop / exact sensitivity 0.791 / 0.544, 5' partial 0.785 / 0.194.
- CDSs of at least 90 codons, stop sensitivity on fragments: TD2 0.958-1.000 (median 0.987), TransDecoder 0.935-0.995 (0.976); TD2 >= TransDecoder in 6 of 8 organisms. 18-26% of fragment CDSs are below 90 codons (TD2 stop sensitivity on them 0.02-0.17, TransDecoder 0.26-0.67).
- Random sequences with a prediction (complete.wRand): 0.6% (*P. falciparum*) to 15.2% (*C. reinhardtii*), and 95.3% in *T. thermophila*. There 66% had no ORF above the PSAURON cutoff, so they passed on length: `length_at_fdr` in `TD2/Predict.py` assumes TAA/TAG/TGA (source comment "TODO alternate translation tables"). Elsewhere length-only calls were 5-25% (organisms with more than 100 predictions).
- lncRNAs with a prediction: human 15.2%, mouse 14.2% (TransDecoder 4.9%, 5.2%); same counts on complete.wLnc and fragments.wLnc, because TD2 trains nothing on the input. GENCODE era: 15.5% in v28, 15.1% after (ratio 0.97).
- Run time on complete inputs: median 389 s (82-796 s), 2 threads; ratio to PyTD 1 thread over 84 inputs: median 1.03 (0.54-1.58). Shared host; indicative.

Simulated assemblies, d1 / d2 / d4 / d8 / d16 / d32:

| | stop F1 | exact F1 | negatives with a prediction (%) |
|---|---|---|---|
| human TD2 | 0.899 / 0.907 / 0.911 / 0.904 / 0.884 / 0.857 | 0.822 / 0.825 / 0.819 / 0.789 / 0.753 / 0.723 | 13.8 / 13.5 / 13.2 / 13.0 / 13.4 / 14.6 |
| human TransDecoder | 0.916 / 0.920 / 0.924 / 0.924 / 0.929 / 0.925 | 0.844 / 0.846 / 0.846 / 0.832 / 0.825 / 0.806 | 8.1 / 8.9 / 10.0 / 12.0 / 14.3 / 16.6 |
| *A. thaliana* TD2 | 0.959 / 0.953 / 0.938 / 0.906 / 0.863 / 0.817 | 0.907 / 0.886 / 0.850 / 0.795 / 0.740 / 0.698 | 6.4 / 5.1 / 4.4 / 3.4 / 1.3 / 1.4 |
| *A. thaliana* TransDecoder | 0.972 / 0.970 / 0.973 / 0.968 / 0.969 / 0.971 | 0.882 / 0.873 / 0.859 / 0.838 / 0.822 / 0.828 | 6.0 / 6.1 / 4.6 / 4.7 / 2.6 / 4.0 |

*A. thaliana* d1 and d2: TD2 has the highest exact F1 of the ab initio methods, which changed the manuscript's sentence that TransDecoder led down to 1/8. Human: TransDecoder remains highest at every depth. Homology mode (TransDecoder; 0.9138 at *A. thaliana* d1) is above TD2 at every depth, but TD2 was run without homology, so this is not a like-for-like comparison and the text says so.
By type with negatives counted, human: complete 0.881-0.923 (TransDecoder 0.869-0.894), 3' partial 0.420-0.537, 5' partial 0.177-0.237, internal 0.747-0.821.

## Checks
- Non-td2 rows of every regenerated table are byte-identical to the pre-E39 tables (`10_td2_tables.sh` compares them; previous tables in `$P/_trash/e39_pre_20261031_150546/`); tables without a tool column are identical.
- `by_length.v7.tsv`, `lnc_call_lengths.human.tsv`, `sens_by_length.v7.tsv`, `lnc_riboseq_orfs.v7.tsv` unchanged.
- TD2 medians agree with the private runs of 2026-10-23 (`work/td2_private/RESULTS.md`: complete 0.982 / 0.941, fragments 0.876 / 0.421, precision 0.904), which were not reused.
- Notebook 02 and 04 hold new `check_text` rows for the TD2 values in the text.

## Caveats
- The author is a TD2 co-author (above).
- Mechanism statements for fragments are from code plus the two tallies; no ablated TD2 was run.
- PSAURON's training set includes all RefSeq plant and animal genomes of August 2023, so human, mouse, *D. melanogaster* and *A. thaliana* are training species (1,000,000 of 33.9 million CDSs were sampled); *P. falciparum*, *S. pombe* and *T. thermophila* are not plants or animals. Whether *C. reinhardtii* was included was not checked.
- PSAURON translates with a built-in standard table and drops `*`, so TAA/TAG (Gln in ciliates) are removed before scoring on *T. thermophila*.
- TD2's `--precise` mode was not run for the paper. The private precise runs are not used.
- Runs were 2-threaded and 4 at a time, unlike the 6-16 single-thread lanes of the main benchmark.
- Preprint, not peer-reviewed. E16's row on it stays as "not independent".
