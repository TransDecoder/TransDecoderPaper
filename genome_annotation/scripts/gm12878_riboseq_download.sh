#!/bin/bash
# GM12878 Ribo-seq, Cenik et al. 2015 Genome Res (GSE65912): CHX footprints, RNase A+T1.
# .sra from the NCBI AWS open-data bucket (ENA FTP ran at ~0.15 MB/s), then fasterq-dump.
set -euo pipefail
P=/broad/hptmp/bhaas/TRANSDECODER_PAPER; export PATH=$P/env/sratools/bin:$PATH
cd $(dirname $(readlink -f $0)); mkdir -p sra fastq tmp
for r in SRR1803147 SRR1803148; do
  curl -s -S --retry 5 -C - -o sra/$r.sra https://sra-pub-run-odp.s3.amazonaws.com/sra/$r/$r
  vdb-validate sra/$r.sra 2>&1 | tail -1
  fasterq-dump -e 8 -t tmp -O fastq sra/$r.sra
  pigz -p 8 fastq/$r.fastq
  echo "ok $r $(( $(zcat fastq/$r.fastq.gz | wc -l) / 4 )) reads"
done
echo DL_DONE
