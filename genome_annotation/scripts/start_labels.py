"""Label 5'-open ORFs (trimmed in the trim_all run) by MANE truth and dump the start-rule features.
usage: start_labels.py <run_dir (trim_all)> <out.tsv>
label: utr_first_atg  - trimmed CDS == MANE CDS (transcript holds the MANE start; trimming is right)
       truncated      - MANE start not in the transcript; trimmed CDS is an in-frame piece of MANE (keep open)
       upstream_atg   - trimmed CDS starts upstream of MANE start, same stop (first ATG is upstream in-frame)
       start_downstream_in_tx - MANE start in transcript but trimmed start is downstream (should not happen)
"""
import re, sys
sys.path.insert(0, "/broad/hptmp/bhaas/TRANSDECODER_PAPER/TransDecoderPaper/genome_annotation/scripts")
from eval_genome_cds import read_ref_gtf, read_pred_gff3, Index, cds_overlap, compatible
run, out = sys.argv[1], sys.argv[2]
feat, cds5 = {}, {}
for l in open(f"{run}/transcripts.fa.transdecoder.gff3"):
    f = l.split("\t")
    if len(f) == 9 and f[2] == "CDS":
        cds5[re.search(r"Parent=([^;\s]+)", f[8]).group(1)] = int(f[3]) - 1   # all ORFs on + (run with -S)
for l in open(f"{run}/transcripts.fa.transdecoder.gff3"):
    f = l.split("\t")
    if len(f) == 9 and f[2] == "mRNA" and "start_revised=true" in f[8]:
        a = dict(x.split("=", 1) for x in f[8].strip().split(";") if "=" in x)
        feat[a["ID"]] = (float(a["upstream_frame_contrast"]), float(a["kozak_score"]),
                         None)
preds = read_pred_gff3(f"{run}/transcripts.fa.transdecoder.genome.gff3")
exons = {}
for l in open(f"{run}/transcripts.fa.transdecoder.genome.gff3"):
    f = l.split("\t")
    if len(f) == 9 and f[2] == "exon":
        exons.setdefault(re.search(r"Parent=([^;\s]+)", f[8]).group(1), []).append((int(f[3]), int(f[4])))
mane, _ = read_ref_gtf("/broad/hptmp/bhaas/TRANSDECODER_PAPER/data/annot/human/MANE.GRCh38.v1.5.ensembl_genomic.gtf.gz", None)
mane = {t: r for t, r in mane.items() if "MANE_Select" in r["tags"]}
idx = Index(mane)
def ends(x):
    s = x["segs"]
    return (s[0][0], s[-1][1]) if x["strand"] == "+" else (s[-1][1], s[0][0])  # (start, stop)
def in_exons(pos, ex): return any(a <= pos <= b for a, b in ex)
with open(out, "w") as fh:
    fh.write("model_id\tlabel\tupstream_frame_contrast\tkozak_score\tupstream_nt\tcds_len\ttx_upstream_nt\n")
    for mid, (con, koz, tx_up) in feat.items():
        tx_up = cds5.get(mid)
        p = preds.get(mid)
        if p is None: continue
        pst, psp = ends(p)
        hits = [mane[k] for k in idx.overlapping(p) if ends(mane[k])[1] == psp]
        if not hits: continue
        r = hits[0]; rst, _ = ends(r)
        ex = exons[mid]
        tx5 = min(a for a, b in ex) if p["strand"] == "+" else max(b for a, b in ex)
        up = abs(pst - tx5)   # genomic distance (approx; introns included)
        if p["segs"] == r["segs"]: lab = "utr_first_atg"
        elif compatible(p, r) and not in_exons(rst, ex): lab = "truncated"
        elif (pst < rst) == (p["strand"] == "+") and pst != rst: lab = "upstream_atg"
        elif in_exons(rst, ex): lab = "start_downstream_in_tx"
        else: lab = "other"
        L = sum(b - a + 1 for a, b in p["segs"])
        fh.write(f"{mid}\t{lab}\t{con}\t{koz}\t{up}\t{L}\t{tx_up}\n")
