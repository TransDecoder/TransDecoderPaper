#!/bin/bash
# GM12878 Ribo-seq (GSE65912): adapter trim -> rRNA/tRNA removal -> STAR (end-to-end, unique) -> merged BAM.
set -euo pipefail
P=/broad/hptmp/bhaas/TRANSDECODER_PAPER; R=$P/data/gm12878/riboseq_cenik2015
export PATH=$P/env/riboseq/bin:$PATH; STAR=$P/env/rnaseq_asm/bin/STAR
cd $R; mkdir -p trim star logs
until grep -q DL_DONE download.log; do grep -q -i error download.log && exit 1; sleep 60; done
for r in SRR1803147 SRR1803148; do
  [ -s trim/$r.norrna.fq.gz ] || {
    cutadapt -j 8 -a AGATCGGAAGAGCACACGTCT -O 3 --discard-untrimmed --trim-n -m 25 -M 34 \
      -o trim/$r.trim.fq.gz fastq/$r.fastq.gz > logs/$r.cutadapt.log
    bowtie2 -p 12 --very-sensitive -L 15 -x rrna/rrna_trna -U trim/$r.trim.fq.gz --un-gz trim/$r.norrna.fq.gz \
      -S /dev/null 2> logs/$r.rrna_bowtie2.log; }
  [ -s star/$r.Aligned.sortedByCoord.out.bam ] || \
  $STAR --runThreadN 12 --genomeDir $P/data/gm12878/star_stringtie/star_index --readFilesIn trim/$r.norrna.fq.gz \
    --readFilesCommand zcat --outFileNamePrefix star/$r. --outSAMtype BAM SortedByCoordinate \
    --alignEndsType EndToEnd --outFilterMismatchNmax 2 --outFilterMultimapNmax 1 --outSAMattributes NH NM \
    --limitBAMsortRAM 20000000000 > logs/$r.star.log 2>&1
  echo "ok $r"
done
samtools merge -f -@ 8 footprints.genome.bam star/SRR1803147.Aligned.sortedByCoord.out.bam star/SRR1803148.Aligned.sortedByCoord.out.bam
samtools index footprints.genome.bam
echo RIBO_PROC_DONE
