"""Add per-base frame contrasts (upstream segment, downstream ORF, downstream window) to start labels.
usage: start_features.py <run_dir> <labels.tsv> <out.tsv>"""
import csv, re, sys
from pytransdecoder.predict import TransDecoderPredict
run, lab_f, out = sys.argv[1:4]
P = TransDecoderPredict.__new__(TransDecoderPredict)
P.stop_codons = {"TAA", "TAG", "TGA"}
d = f"{run}/transcripts.fa.transdecoder_dir"
scores = P._load_hexamer_scores(f"{d}/hexamer.scores")
seqs, name = {}, None
for l in open(f"{run}/transcripts.fa"):
    if l.startswith(">"): name = l[1:].split()[0]; seqs[name] = []
    else: seqs[name].append(l.strip().upper())
seqs = {k: "".join(v) for k, v in seqs.items()}
cds = {}
for l in open(f"{run}/transcripts.fa.transdecoder.gff3"):
    f = l.split("\t")
    if len(f) == 9 and f[2] == "CDS":
        cds[re.search(r"Parent=([^;\s]+)", f[8]).group(1)] = (f[0], int(f[3]), int(f[4]))
fc = lambda s: P._frame_contrast(s, scores) / max(len(s), 1)
rows = list(csv.DictReader(open(lab_f), delimiter="\t"))
with open(out, "w") as fh:
    cols = list(rows[0].keys()) + ["up_pb", "orf_pb", "down_win_pb"]
    fh.write("\t".join(cols) + "\n")
    for r in rows:
        tid, s, e = cds[r["model_id"]]
        seq = seqs[tid]; a = s - 1; off = a % 3
        up, orf = seq[off:a], seq[a:e]
        win = orf[:max(len(up), 90) // 3 * 3]
        r.update(up_pb=f"{fc(up):.4f}", orf_pb=f"{fc(orf):.4f}", down_win_pb=f"{fc(win):.4f}")
        fh.write("\t".join(str(r[c]) for c in cols) + "\n")
