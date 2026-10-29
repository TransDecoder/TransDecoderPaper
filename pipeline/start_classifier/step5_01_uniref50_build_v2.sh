#!/bin/bash
# v2: XML parser step5_01_uniref50_xml_v2.awk (field split on quotes, ~2x faster than v1; v1 run stopped during pass 1)
# E33 step 5: build the UniRef50 2026_03 DIAMOND db with taxonomy, streaming (no FASTA/XML on disk).
#  pass 1: XML -> representative taxid per cluster (+ which benchmark organisms are cluster members)
#  pass 2: FASTA -> headers stripped to the cluster id, piped into diamond makedb with
#          --taxonmap (representative taxid) --taxonnodes/--taxonnames (NCBI taxdump);
#          per-sequence length, first residue, n= and LCA TaxID recorded on the way (gz).
set -euo pipefail
P=/broad/hptmp/bhaas/TRANSDECODER_PAPER; T=$P/TransDecoderPaper; D=$T/pipeline/start_classifier
B=$P/data/start_classifier/step5; DB=$B/db; TX=$B/taxonomy
export PATH=$P/env/td/bin:$PATH
URL=https://ftp.uniprot.org/pub/databases/uniprot/uniref/uniref50
mkdir -p $DB
date; echo "pass 0: own taxsets"
python3 -I $D/step5_00_own_taxsets_v1.py $TX/nodes.dmp > $DB/own_taxsets.tsv
if [ ! -s $DB/uniref50_xml_rep.tsv.gz ]; then
  date; echo "pass 1: XML"
  curl -sS --fail $URL/uniref50.xml.gz | pigz -dc -p 2 \
    | LC_ALL=C grep -F -e '<entry id=' -e '<representativeMember>' -e '<member>' -e '"NCBI taxonomy"' -e '"common taxon ID"' -e '"member count"' \
    | LC_ALL=C awk -F"\"" -f $D/step5_01_uniref50_xml_v2.awk $DB/own_taxsets.tsv - | gzip -c > $DB/uniref50_xml_rep.tsv.gz.tmp
  mv $DB/uniref50_xml_rep.tsv.gz.tmp $DB/uniref50_xml_rep.tsv.gz
fi
date; echo "pass 1 rows: $(zcat $DB/uniref50_xml_rep.tsv.gz | wc -l)"
# taxonmap: representative taxid; if the XML lacks it, the cluster's common taxid
zcat $DB/uniref50_xml_rep.tsv.gz | awk 'BEGIN{OFS="\t"; print "accession","accession.version","taxid","gi"} {t=($2=="NA"?$3:$2); print $1,$1,t,0}' | gzip -c > $DB/uniref50_rep_taxonmap.gz
date; echo "pass 2: FASTA -> diamond makedb"
curl -sS --fail $URL/uniref50.fasta.gz | pigz -dc -p 2 \
  | LC_ALL=C awk '
      function flush() { if (id != "") print id "\t" len "\t" first "\t" n "\t" lca | "gzip -c > '"$DB"'/uniref50_seqinfo.tsv.gz.tmp" }
      /^>/ { flush(); id = substr($1, 2); n = "NA"; lca = "NA"
             for (i = 2; i <= NF; i++) { if ($i ~ /^n=/) n = substr($i, 3); else if ($i ~ /^TaxID=/) lca = substr($i, 7) }
             len = 0; first = ""; print ">" id; next }
      { if (first == "") first = substr($0, 1, 1); len += length($0); print }
      END { flush() }' \
  | nice -n 19 diamond makedb --in /dev/stdin -d $DB/uniref50_2026_03 --threads 8 \
      --taxonmap $DB/uniref50_rep_taxonmap.gz --taxonnodes $TX/nodes.dmp --taxonnames $TX/names.dmp
mv $DB/uniref50_seqinfo.tsv.gz.tmp $DB/uniref50_seqinfo.tsv.gz
date; echo "seqinfo rows: $(zcat $DB/uniref50_seqinfo.tsv.gz | wc -l)"
diamond dbinfo -d $DB/uniref50_2026_03.dmnd
ls -la $DB; du -sh $DB
echo BUILD_DONE
