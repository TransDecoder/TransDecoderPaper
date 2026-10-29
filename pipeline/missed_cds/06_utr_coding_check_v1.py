#!/usr/bin/env python3
"""E34: are the 5' UTR segments that score as coding actually coding in another isoform of the same gene?

usage: python 06_utr_coding_check_v1.py <manifest.tsv> <big_dir> <out_tsv>
For truth complete / 3' partial CDSs whose candidate ORF is 5'-open, whose true start is the first in-frame ATG
and whose stop was predicted (start5open.tsv.gz, outcome exact or kept_open), the segment u (open end to the ATG)
is translated in the ORF frame. u is "annotated coding" when its translation shares an identical 10-aa word with
any annotated protein of the same gene (kept.tsv/kept.fa of E30/E32: the simulated coding transcripts), i.e. the
same sequence is translated in this frame in another isoform (an N-terminal extension). Segments shorter than
30 nt (< 10 aa) are reported as too short to test.
"""
import sys

import pandas as pd

man_path, big, out_tsv = sys.argv[1:4]
P = "/broad/hptmp/bhaas/TRANSDECODER_PAPER"
REF = {"human": f"{P}/data/simtrinity/ref", "athal": f"{P}/data/simtrinity_athal/ref"}
BASES = "TCAG"
AA = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
CODON = {a + b + c: AA[16 * i + 4 * j + k] for i, a in enumerate(BASES) for j, b in enumerate(BASES)
         for k, c in enumerate(BASES)}
ALPHA = 0.4


def translate(s):
    return "".join(CODON.get(s[i:i + 3], "X") for i in range(0, len(s) - 2, 3))


def stream_fasta(path, keep=None):
    name, buf = None, []
    for line in open(path):
        if line.startswith(">"):
            if name and (keep is None or name in keep):
                yield name, "".join(buf).upper()
            name, buf = line[1:].split()[0], []
        else:
            buf.append(line.strip())
    if name and (keep is None or name in keep):
        yield name, "".join(buf).upper()


man = pd.read_csv(man_path, sep="\t")
rows = []
for sp, M in man.groupby("species"):
    kept = pd.read_csv(f"{REF[sp]}/kept.tsv", sep="\t")
    kept = kept[kept["class"] == "coding"].set_index("transcript_id")
    work = {}
    genes = set()
    for _, m in M.iterrows():
        S = pd.read_csv(f"{big}/{m.dataset}/start5open.tsv.gz", sep="\t")
        S = S[S.true_is_first_atg.fillna(False).astype(bool) & S.outcome.isin(["exact", "kept_open"])]
        t = pd.read_csv(f"{m.labels_prefix}.truth.tsv", sep="\t", usecols=["transcript_id", "source_transcripts"])
        t["src"] = t.source_transcripts.str.split(",").str[0]
        S = S.merge(t[["transcript_id", "src"]], on="transcript_id")
        S["gene"] = S.src.map(kept.gene_id)
        if m.condition == "B":
            H = pd.read_csv(f"{big}/{m.dataset}/homology_decisions.tsv.gz", sep="\t", dtype=str).set_index("orf_id")
            S["source"] = S.cand_id.map(H.source)
        else:
            S["source"] = "ab_initio"
        work[m.dataset] = (m, S)
        genes |= set(S.gene.dropna())
    need = set(kept.index[kept.gene_id.isin(genes)])
    words = {}
    for tid, seq in stream_fasta(f"{REF[sp]}/kept.fa", need):
        r = kept.loc[tid]
        p = translate(seq[int(r.cds_start) - 1:int(r.cds_end)])
        w = words.setdefault(r.gene_id, set())
        w.update(p[i:i + 10] for i in range(len(p) - 9))
    for ds, (m, S) in work.items():
        need_tx = set(S.transcript_id)
        tx = dict(stream_fasta(f"{m.run_dir}/input.fa", need_tx))
        lo = {}
        for line in open(f"{m.run_dir}/input.fa.transdecoder_dir/longest_orfs.gff3"):
            x = line.split("\t")
            if len(x) == 9 and x[2] == "mRNA":
                lo[x[8].split("ID=")[1].split(";")[0]] = int(x[3])
        res = []
        for r in S.itertuples():
            u = tx[r.transcript_id][lo[r.cand_id] - 1:lo[r.cand_id] - 1 + int(r.atg_first)]
            pu = translate(u)
            if len(pu) < 10:
                res.append("too_short")
                continue
            w = words.get(r.gene, set())
            res.append("annotated_coding" if any(pu[i:i + 10] in w for i in range(len(pu) - 9)) else "not_coding")
        S["u_status"] = res
        S["scores_coding"] = S.coding_u >= ALPHA * S.atg_first
        S.to_csv(f"{big}/{ds}/utr_coding_check.tsv.gz", sep="\t", index=False)
        for typ in ["complete", "3prime_partial"]:
            for oc in ["kept_open", "exact"]:
                for sc in [True, False]:
                    g = S[(S.type == typ) & (S.outcome == oc) & (S.scores_coding == sc)]
                    if oc == "kept_open":
                        g = g[g.source != "homology"]
                    rows.append(dict(dataset=ds, type=typ, outcome=oc, u_scores_coding=sc, n=len(g),
                                     annotated_coding=int((g.u_status == "annotated_coding").sum()),
                                     not_coding=int((g.u_status == "not_coding").sum()),
                                     too_short=int((g.u_status == "too_short").sum())))
        print(ds, "done", file=sys.stderr)
D = pd.DataFrame(rows)
D["annotated_coding_frac_of_testable"] = (D.annotated_coding / (D.annotated_coding + D.not_coding)).round(3)
D.to_csv(out_tsv, sep="\t", index=False)
print("written", out_tsv)
