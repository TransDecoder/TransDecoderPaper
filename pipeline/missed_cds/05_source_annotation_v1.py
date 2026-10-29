#!/usr/bin/env python3
"""E34: what kind of annotated CDS are the class (b) misses? (source-transcript annotation)

usage: python 05_source_annotation_v1.py <manifest.tsv> <big_dir> <out_tsv>
For each truth CDS (first source transcript): its gene's longest annotated protein among the simulated coding
transcripts (kept.tsv / kept.fa of E30 or E32), and
  rel_main      main = the CDS protein is (a substring of) the gene's longest protein, or that protein itself
                shared_frame = shares >= 1 identical 10-aa word with it (same frame, alternative exons/ends)
                novel = no shared 10-aa word (a different reading frame or unrelated sequence)
  frac_of_main  CDS length / longest protein length of the gene
  tags          (human) GENCODE tags of the source transcript: MANE_Select, Ensembl_canonical, basic
Writes one summary row per dataset x truth type x group (missed_b / found) and per-CDS rows to <big_dir>/<ds>/source_annot.tsv.gz.
"""
import collections
import gzip
import sys

import pandas as pd

man_path, big, out_tsv = sys.argv[1:4]
P = "/broad/hptmp/bhaas/TRANSDECODER_PAPER"
REF = {"human": f"{P}/data/simtrinity/ref", "athal": f"{P}/data/simtrinity_athal/ref"}
BASES = "TCAG"
AA = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
CODON = {a + b + c: AA[16 * i + 4 * j + k] for i, a in enumerate(BASES) for j, b in enumerate(BASES)
         for k, c in enumerate(BASES)}


def translate(s):
    return "".join(CODON.get(s[i:i + 3], "X") for i in range(0, len(s) - 2, 3)).rstrip("*")


def stream_fasta(path, keep):
    name, buf = None, []
    for line in open(path):
        if line.startswith(">"):
            if name in keep:
                yield name, "".join(buf).upper()
            name, buf = line[1:].split()[0], []
        else:
            buf.append(line.strip())
    if name in keep:
        yield name, "".join(buf).upper()


man = pd.read_csv(man_path, sep="\t")
cache = {}
rows = []
for sp, M in man.groupby("species"):
    kept = pd.read_csv(f"{REF[sp]}/kept.tsv", sep="\t")
    kept = kept[kept["class"] == "coding"].set_index("transcript_id")
    tags = {}
    if sp == "human":
        g = pd.read_csv(f"{REF[sp]}/gencode_v50.tx.tsv.gz", sep="\t", usecols=["transcript_id", "tags"])
        tags = dict(zip(g.transcript_id, g.tags.fillna(".")))
    # genes needed
    truths = {}
    genes = set()
    for _, m in M.iterrows():
        t = pd.read_csv(f"{m.labels_prefix}.truth.tsv", sep="\t", usecols=["transcript_id", "source_transcripts"])
        t["src"] = t.source_transcripts.str.split(",").str[0]
        truths[m.dataset] = t
        genes |= set(kept.gene_id.reindex(t.src).dropna())
    need = set(kept.index[kept.gene_id.isin(genes)])
    prot = {}
    for tid, seq in stream_fasta(f"{REF[sp]}/kept.fa", need):
        r = kept.loc[tid]
        prot[tid] = translate(seq[int(r.cds_start) - 1:int(r.cds_end)])
    main = {}
    for tid, p in prot.items():
        gid = kept.at[tid, "gene_id"]
        if gid not in main or len(p) > len(main[gid]):
            main[gid] = p
    print(sp, len(prot), "proteins,", len(main), "genes", file=sys.stderr)
    for _, m in M.iterrows():
        D = pd.read_csv(f"{big}/{m.dataset}/per_cds.tsv.gz", sep="\t",
                        usecols=["transcript_id", "type", "status", "miss_class", "n_codons"])
        D = D.merge(truths[m.dataset][["transcript_id", "src"]], on="transcript_id")
        out = []
        for r in D.itertuples():
            p = prot.get(r.src)
            gid = kept.gene_id.get(r.src)
            mp = main.get(gid)
            if p is None or mp is None:
                out.append((".", float("nan"), "."))
                continue
            if p in mp:
                rel = "main"
            else:
                words = {p[i:i + 10] for i in range(0, max(1, len(p) - 9))}
                rel = "shared_frame" if any(mp[i:i + 10] in words for i in range(0, max(1, len(mp) - 9))) else "novel"
            tg = tags.get(r.src, ".")
            tclass = ("MANE_Select" if "MANE_Select" in tg else "Ensembl_canonical" if "Ensembl_canonical" in tg
                      else "basic" if "basic" in tg else "other") if sp == "human" else "."
            out.append((rel, round(len(p) / len(mp), 3), tclass))
        D["rel_main"], D["frac_of_main"], D["tag_class"] = zip(*out)
        D.to_csv(f"{big}/{m.dataset}/source_annot.tsv.gz", sep="\t", index=False)
        D["group"] = D.status.map(lambda s: "missed" if s == "missed" else "found")
        D.loc[(D.group == "missed") & (D.miss_class == "b"), "group"] = "missed_b"
        for typ in ["complete", "5prime_partial", "3prime_partial", "internal", "all"]:
            d = D if typ == "all" else D[D.type == typ]
            for grp in ["missed_b", "found"]:
                g = d[d.group == grp]
                if not len(g):
                    continue
                r = dict(dataset=m.dataset, type=typ, group=grp, n=len(g),
                         frac_of_main_median=round(g.frac_of_main.median(), 3),
                         **{f"rel_{k}": round((g.rel_main == k).mean(), 3) for k in ["main", "shared_frame", "novel"]})
                if sp == "human":
                    r.update({f"tag_{k}": round((g.tag_class == k).mean(), 3)
                              for k in ["MANE_Select", "Ensembl_canonical", "basic", "other"]})
                # the same, within 50-99 codons (length-matched)
                h = g[g.n_codons < 100]
                r["n_50_99"] = len(h)
                r["rel_novel_50_99"] = round((h.rel_main == "novel").mean(), 3) if len(h) else float("nan")
                r["rel_main_50_99"] = round((h.rel_main == "main").mean(), 3) if len(h) else float("nan")
                rows.append(r)
pd.DataFrame(rows).to_csv(out_tsv, sep="\t", index=False)
print("written", out_tsv)
