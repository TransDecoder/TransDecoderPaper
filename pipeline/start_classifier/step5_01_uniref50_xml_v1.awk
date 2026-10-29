# E33 step 5: one row per UniRef50 cluster from the (grep-filtered) UniRef50 XML stream.
# Input 1: own_taxsets.tsv (taxid, organism taxid); input 2: filtered XML lines.
# Output: cluster id, representative's NCBI taxid, common (LCA) taxid, member count,
#         benchmark organisms (taxid incl. descendants) present among rep + members.
function val(   s) { s = $0; sub(/.*value="/, "", s); sub(/".*/, "", s); return s }
function out() { print id "\t" (rep == "" ? "NA" : rep) "\t" ct "\t" n "\t" (flags == "" ? "-" : flags) }
FNR == NR { if (FNR > 1) own[$1] = $2; next }
/^<entry id=/ {
  if (id != "") out()
  s = $0; sub(/^<entry id="/, "", s); sub(/".*/, "", s); id = s
  rep = ""; ct = ""; n = ""; state = 0; flags = ""; split("", seen); next
}
/"member count"/ { n = val(); next }
/"common taxon ID"/ { ct = val(); next }
/^<representativeMember>/ { state = 1; next }
/^<member>/ { state = 2; next }
/"NCBI taxonomy"/ {
  v = val()
  if (state == 1 && rep == "") rep = v
  if (v in own) { o = own[v]; if (!(o in seen)) { seen[o] = 1; flags = flags (flags == "" ? "" : ",") o } }
}
END { if (id != "") out() }
