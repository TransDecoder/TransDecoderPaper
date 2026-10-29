# E33 step 5 (v2 = v1 without regex substitutions; run with -F'"'): one row per UniRef50 cluster from
# the grep-filtered UniRef50 XML stream.
# Input 1: own_taxsets.tsv (taxid, organism taxid; read with the same -F, so lines are split by hand);
# input 2: filtered XML lines.
# Output: cluster id, representative's NCBI taxid, common (LCA) taxid, member count,
#         benchmark organisms (taxid incl. descendants) present among rep + members.
function out() { print id "\t" (rep == "" ? "NA" : rep) "\t" ct "\t" n "\t" (flags == "" ? "-" : flags) }
FNR == NR { if (FNR > 1) { split($0, z, "\t"); own[z[1]] = z[2] } next }
$2 == "NCBI taxonomy" {
  v = $4
  if (state == 1 && rep == "") rep = v
  if (v in own) { o = own[v]; if (!(o in seen)) { seen[o] = 1; flags = flags (flags == "" ? "" : ",") o } }
  next
}
$1 == "<member>" { state = 2; next }
$1 == "<entry id=" {
  if (id != "") out()
  id = $2; rep = ""; ct = ""; n = ""; state = 0; flags = ""; split("", seen); next
}
$2 == "member count" { n = $4; next }
$2 == "common taxon ID" { ct = $4; next }
$1 == "<representativeMember>" { state = 1; next }
END { if (id != "") out() }
