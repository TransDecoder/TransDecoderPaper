#!/bin/bash
# E06 release rerun: run the TE screen recipe of PyTD docs/homology_te_screen.md (PyTD d3892b5) on the
# GM12878 homology series (Swiss-Prot minus human, BLAST only, vs ab initio; PyTD a5b24c3 runs) and
# summarise what it flags. Output: data/gm12878/gm12878_te_screen_v7h.tsv
# Needs bedtools ($P/env/bedtools) and the UCSC hg38 RepeatMasker table ($P/data/repeats/ucsc_hg38_rmsk).
# usage (from the paper repo root): bash pipeline/derive/08_te_screen_check_v7h.sh
set -euo pipefail
P=/broad/hptmp/bhaas/TRANSDECODER_PAPER; T=$(pwd); W=$P/work/gm12878_td; X=$P/work/gm12878_td/te_screen_v7h
export PATH=$P/env/bedtools/bin:$P/env/td/bin:$PATH; export LC_ALL=C
mkdir -p $X; cd $X
zcat $P/data/homology/uniprot_sprot.fasta.gz | awk '/^>/ {id = substr($1, 2); $1 = ""; sub(/^ /, ""); print id"\t"$0}' |
  sort -k1,1 > sprot_names.tsv
zcat $P/data/repeats/ucsc_hg38_rmsk/rmsk.txt.gz |
  awk -F'\t' -v OFS='\t' '$12 ~ /^(SINE|LINE|LTR|DNA|RC|Retroposon)/ {print $6, $7, $8}' |
  sort -k1,1 -k2,2n | bedtools merge > te.bed
ids() { grep '^>' "$1" | cut -c2- | cut -d' ' -f1 | sort -u; }
out=$T/data/gm12878/gm12878_te_screen_v7h.tsv
echo -e "assembly\thomology_run\tpredicted\tabinitio\thomology_only\tflag_by_name\tflag_by_overlap\tflag_either\thomology_only_genome_mapped\tbest_hit_UTY_PANTR\tUTY_PANTR_genome_mapped\tUTY_PANTR_flag_by_overlap" > $out
for a in lraa stringtie trinity; do
  H=$W/$a/hseries_v7h; R=blast_9606_nopfam; mkdir -p $a; cd $a
  ids $H/$R/transcripts.fa.transdecoder.pep > homology.ids
  ids $H/abinitio/transcripts.fa.transdecoder.pep > abinitio.ids
  comm -23 homology.ids abinitio.ids > homology_only.ids
  awk -F'\t' '!seen[$1]++ {print $1"\t"$2}' $H/blastp.minus_9606.outfmt6 | sort -k1,1 > best_hit.tsv
  join -t $'\t' homology_only.ids best_hit.tsv | sort -t $'\t' -k2,2 |
    join -t $'\t' -1 2 -2 1 - ../sprot_names.tsv |
    awk -F'\t' -v OFS='\t' '{print $2, $1, $3}' | sort -k1,1 > homology_only.best_hit.tsv
  grep -iE 'LINE-1|transpos|reverse transcriptase|endogenous retrovirus|retrovirus-related' \
    homology_only.best_hit.tsv | cut -f1 > te_by_name.ids || true
  awk -F'\t' -v OFS='\t' '$3 == "CDS" { match($9, /Parent=[^;]+/); print $1, $4 - 1, $5, substr($9, RSTART + 7, RLENGTH - 7) }' \
    $H/$R/transcripts.fa.transdecoder.genome.gff3 | sort -k1,1 -k2,2n > cds.bed
  bedtools coverage -a cds.bed -b ../te.bed 2> bedtools.warnings |
    awk -F'\t' -v OFS='\t' '{cov[$4] += $6; len[$4] += $7} END {for (o in len) print o, cov[o] / len[o]}' |
    sort -k1,1 > cds_te_fraction.tsv
  join -t $'\t' homology_only.ids cds_te_fraction.tsv > homology_only.te_fraction.tsv
  awk -F'\t' '$2 >= 0.5 {print $1}' homology_only.te_fraction.tsv > te_by_overlap.ids
  sort -u te_by_name.ids te_by_overlap.ids > te_flagged.ids
  awk -F'\t' '$2 ~ /UTY_PANTR/ {print $1}' homology_only.best_hit.tsv | sort > uty.ids
  echo -e "$a\t$R\t$(wc -l < homology.ids)\t$(wc -l < abinitio.ids)\t$(wc -l < homology_only.ids)\t$(wc -l < te_by_name.ids)\t$(wc -l < te_by_overlap.ids)\t$(wc -l < te_flagged.ids)\t$(wc -l < homology_only.te_fraction.tsv)\t$(wc -l < uty.ids)\t$(join uty.ids cds_te_fraction.tsv | wc -l)\t$(comm -12 uty.ids te_by_overlap.ids | wc -l)" >> $out
  cd ..
done
cat $out
