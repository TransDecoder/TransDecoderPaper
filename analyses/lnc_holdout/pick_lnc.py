"""One lncRNA per gene (longest transcript, >= 300 nt) whose exons overlap no GENCODE CDS on the same
strand; sequences renamed lnc_<transcript_id>.
usage: pick_lnc.py <lnc.gtf> <lnc.fa> <annotation gtf|gff3> <out.fa>"""
import re, sys
from collections import defaultdict
lnc_gtf, lnc_fa, annot, out = sys.argv[1:5]
BIN = 100000
cds = defaultdict(list)
op = open
for l in op(annot):
    f = l.split("\t")
    if len(f) == 9 and f[2] == "CDS":
        a, b = int(f[3]), int(f[4])
        for k in range(a // BIN, b // BIN + 1):
            cds[(f[0], f[6], k)].append((a, b))
ex, gene = defaultdict(list), {}
for l in open(lnc_gtf):
    f = l.split("\t")
    tid = re.search(r'transcript_id "([^"]+)"', f[8]) or re.search(r"transcript_id=([^;]+)", f[8])
    gid = re.search(r'gene_id "([^"]+)"', f[8]) or re.search(r"gene_id=([^;]+)", f[8])
    ex[tid.group(1)].append((f[0], f[6], int(f[3]), int(f[4]))); gene[tid.group(1)] = gid.group(1)
def hits(t):
    for c, s, a, b in ex[t]:
        for k in range(a // BIN, b // BIN + 1):
            if any(x <= b and a <= y for x, y in cds.get((c, s, k), ())):
                return True
    return False
seqs, name = {}, None
for l in open(lnc_fa):
    if l[0] == ">": name = l[1:].split()[0]; seqs[name] = []
    else: seqs[name].append(l.strip())
seqs = {k: "".join(v) for k, v in seqs.items()}
best = {}
for t, s in seqs.items():
    if len(s) >= 300 and t in gene and not hits(t):
        g = gene[t]
        if g not in best or len(s) > len(seqs[best[g]]): best[g] = t
with open(out, "w") as fh:
    for g, t in sorted(best.items()):
        fh.write(f">lnc_{t}\n{seqs[t]}\n")
print(len(best), "lncRNAs", file=sys.stderr)
