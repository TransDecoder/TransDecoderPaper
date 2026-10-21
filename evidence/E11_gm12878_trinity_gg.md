# E11. Trinity genome-guided assembly of GM12878 (and upstream GM12878 data)

## Claims supported
Manuscript: Methods, application data (GM12878 transcript assemblies).

1. The Trinity assembly is Trinity 2.15.2 genome-guided assembly of the STAR alignments of ENCSR000AED rep 1 (stranded RF).
   - It gave 428,960 transcripts in 319,068 Trinity "genes", from 77,695 read partitions, none failed.
2. The Trinity transcripts were aligned to GRCh38 with minimap2 `splice:hq`: 428,223 (99.83%) have a primary alignment.
3. The same STAR alignments feed StringTie (E03) and the expression estimates (E05). The PacBio long reads feed LRAA (E03).

## Data
Downloaded 2026-10-19 from the ENCODE portal with `$T/genome_annotation/scripts/data_prep/download.sh <outdir> <url> <md5>` (curl with resume, then an md5 check). The selection rationale and ENCODE metadata are in `$P/data/gm12878/README.md`.

| Role | Experiment / library | File | Reads | md5 (verified, `logs/dl_*.log`) |
|---|---|---|---|---|
| Short reads R1 | ENCSR000AED rep 1 (CSHL; polyA+; stranded, read 1 antisense; HiSeq 2000 PE 101; library ENCLB037ZZZ, GEO GSM2072351) | ENCFF001REK.fastq.gz | 97,548,052 | 0055876c86d1909b7b7ff0e5497bd92b |
| Short reads R2 | same | ENCFF001REJ.fastq.gz | 97,548,052 | b484a76ccf390e100d60f373a9c54dd6 |
| Long reads | ENCSR838WFC rep 2 (UCI; PacBio Sequel II Iso-Seq FLNC; library ENCLB735WVC, GEO GSM5331162) | ENCFF694DIE.fastq.gz | 2,538,701 | 38c1eaea4283c60b6d8f0b3025d540cd |

- **Local paths:** `$P/data/gm12878/reads/{short,long}/`.
- **Download URLs** were passed as arguments and not logged. The ENCODE file download URL has the form `https://www.encodeproject.org/files/<ENCFF>/@@download/<ENCFF>.fastq.gz`.
- **Genome:** GENCODE release 50 `GRCh38.primary_assembly.genome.fa.gz` (md5 da1a11258be075cfa7af718162c894e7).
- **Annotation:** `gencode.v50.primary_assembly.annotation.gtf.gz` (md5 289b91e5e95e8b0450d223246f10a12e); both are in `$P/data/gm12878/ref/`.
- **STAR alignment** (`star_index.sh`, `star_align_stringtie.sh`; STAR 2-pass, index with GENCODE v50 junctions, overhang 100):
  - 97,548,052 input pairs, 83.31% uniquely mapped (`GM12878_ENCSR000AED_rep1.Log.final.out`);
  - output: `$P/data/gm12878/star_stringtie/GM12878_ENCSR000AED_rep1.STAR.bam`.

## Software
| Tool | Version (from the binary) | Environment |
|---|---|---|
| Trinity | Trinity-v2.15.2 (bioconda build pl5321h077b44d_6) | `$P/env/trinity215` |
| salmon | 1.10.3 | `$P/env/trinity215` |
| bowtie2 / jellyfish / samtools | 2.5.5 / 2.3.1 / 1.23 | `$P/env/trinity215` |
| STAR / minimap2 / samtools | 2.7.11b / 2.31-r1302 / 1.24 | `$P/env/rnaseq_asm` |

Environment:
```bash
env/bin/micromamba create -y -p env/trinity215 -c conda-forge -c bioconda --override-channels trinity=2.15.2 salmon=1.10.3
```
- Pin salmon 1.10.3. Unpinned salmon 2.x breaks Trinity's `salmon_runner.pl`.
- An earlier env, `$P/env/trinity`, was broken by a salmon install (Trinity was downgraded to 2.1.1; `$P/env/trinity.BROKEN_README.txt`).
- Trinity's Perl needs `PERL5LIB` etc. unset (`trinity_gg.sh` does this).

**ParaFly replacement.**
- **Problem:** the bioconda ParaFly in Trinity 2.15.2 runs commands serially, whatever `-CPU` says.
  - The binary links libgomp and calls `omp_set_num_threads`, but has no `GOMP_parallel*` calls. So its `#pragma omp parallel` loop was compiled without `-fopenmp`; per the issue below, the recipe adds `-fopenmp` to `LDFLAGS` only.
  - Reported as bioconda-recipes issue #70548 (2026-10-21), which also covers the salmon 2.x pin.
- **Fix:** a bash/xargs drop-in, committed as `$T/genome_annotation/scripts/data_prep/ParaFly_xargs.sh`, installed in place of the binary:
```bash
B=$P/env/trinity215/bin/trinity-plugins/BIN
mv $B/ParaFly $B/ParaFly.bioconda_serial
cp $T/genome_annotation/scripts/data_prep/ParaFly_xargs.sh $B/ParaFly && chmod +x $B/ParaFly
```
- **Behaviour:** the drop-in keeps ParaFly's interface and resume semantics:
  - options `-c`, `-CPU`, `-shuffle`, `-v`, `-failed_cmds`, `-max_retry`;
  - successful commands are appended to `<cmds>.completed` and skipped on rerun;
  - failures are written to `FailedCommands`, with a non-zero exit;
  - commands run with `xargs -P <CPU>`.

## Commands
`$T/genome_annotation/scripts/data_prep/trinity_gg.sh` (identical to `$P/data/gm12878/bin/trinity_gg.sh`). It was run with `CPU=16`, as `Trinity.timing` and the minimap2 log record:
```bash
Trinity --genome_guided_bam $D/star_stringtie/GM12878_ENCSR000AED_rep1.STAR.bam \
  --genome_guided_max_intron 100000 --SS_lib_type RF \
  --max_memory 50G --CPU 16 --output $D/trinity/trinity_gg_out
ln -sf $D/trinity/trinity_gg_out/Trinity-GG.fasta $D/trinity/Trinity.fasta
minimap2 -t 16 -ax splice:hq -uf --secondary=no $D/ref/GRCh38.primary_assembly.genome.fa $D/trinity/Trinity.fasta \
  | samtools sort -@ 4 -m 2G -o $D/trinity/Trinity.minimap2.bam -
samtools index $D/trinity/Trinity.minimap2.bam
samtools flagstat $D/trinity/Trinity.minimap2.bam > $D/trinity/Trinity.minimap2.flagstat
```
- `D=$P/data/gm12878`.
- `-uf` aligns transcripts on their given (RF-derived) strand only.
- Trinity's defaults otherwise apply (genome-guided read normalisation to max coverage 200, `norm_200`).

Timeline (2026-10-20/21):
- **Earlier invocation:** partitioning was done by an invocation started before the ParaFly fix. It wrote the normalised BAM at 13:40, the partitions at 17:18 and `trinity_GG.cmds` at 17:45.
- **Resume:** after the fix, the script was rerun at 18:11. Trinity resumed from its checkpoints, and ParaFly ran 77,441 of 77,695 partition commands (254 had already completed). The run finished at 04:59:42.
- **Aggregation and alignment:** aggregation took until 05:00:50; minimap2 took 164 s wall-clock; `TRINITY_ROUTE_DONE` was printed by 05:04.

## Outputs
Not in git:
- `$P/data/gm12878/trinity/trinity_gg_out/Trinity-GG.fasta` (= `Trinity.fasta`) and `.gene_trans_map`;
- `$P/data/gm12878/trinity/Trinity.minimap2.bam` (+ `.bai`, `.flagstat`), `minimap2.log`;
- logs: `$P/data/gm12878/logs/trinity_gg.log`, `trinity_gg_out/Trinity.timing`, `trinity_GG.cmds{,.completed}`, `FailedCommands`.

## Results
| Metric | Value | Source |
|---|---|---|
| Partition commands | 77,695 | `trinity_GG.cmds` (lines) |
| Completed / failed | 77,695 / 0 | `trinity_GG.cmds.completed` (lines); `FailedCommands` (0 bytes); log "All commands completed successfully" |
| Transcripts | 428,960 | `grep -c '^>' Trinity-GG.fasta` |
| Trinity "genes" | 319,068 | distinct genes in `.gene_trans_map` |
| Primary alignments mapped | 428,223 (99.83%) | `Trinity.minimap2.flagstat` |
| Supplementary alignments | 6,588 | flagstat |

Downstream (E03):
- 97,088 transcript ORFs were predicted with the final defaults, and 92,464 projected to the genome. The 4,611 CDSs extending outside the aligned part and the 13 unmapped transcripts were skipped.

## Checks
- A timing test of 4 × `sleep 2` with `-CPU 4` took 8.0 s with `ParaFly.bioconda_serial` (also with `OMP_NUM_THREADS` unset) and 2.3 s with the drop-in.
- All data-prep scripts in `$T/genome_annotation/scripts/data_prep/` are identical to those run from `$P/data/gm12878/bin/`.
- The number of mapped transcripts equals the assembled-transcript count used by `eval_genome_cds.py` (428,223).

## Caveats
- **`logs/trinity_gg.log` covers only the resumed run from 18:11.** The earlier invocation's log (read normalisation and partitioning) was overwritten, and its exact start time and ParaFly version are not recorded. Its outputs were reused through Trinity's checkpoints (`*.ok` files).
- **Trinity de novo** (`trinity_denovo.sh`) was started but abandoned: it was projected to take weeks on this machine. Its partial output is in `$P/data/gm12878/trinity_denovo_killed/`. It is not used.
- **This is a single short-read replicate,** so Trinity GG produces many fragments and low-expression isoforms. That is reflected in the partial-ORF counts in E03.
