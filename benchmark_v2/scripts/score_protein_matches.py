#!/usr/bin/env python3
"""Score predicted proteins against a reference proteome, as in the TD2 preprint (Mao et al. 2025).

Each predicted CDS is translated (standard code, terminal stop removed) and classed as
  complete       identical over its full length to a reference protein
  match30        not complete, but shares at least one exact 30-residue window with a reference protein
  nonmatching    neither (mostly false positives, though some may be unannotated proteins)
'true positives' in the preprint = complete + match30. Also reports distinct protein sequences and
distinct reference proteins recovered complete, since assemblies contain redundant transcripts.

usage: score_protein_matches.py <ref_proteins.fa[.gz]> <label=transcripts.fa:preds.tsv> ... > out.tsv
"""
import gzip
import sys

import numpy as np

K = 30
AA = "ACDEFGHIKLMNPQRSTVWYUXBZOJ"
CODE = np.zeros(256, dtype=np.uint64)
for i, a in enumerate(AA):
    CODE[ord(a)] = i + 1  # 0 is the separator
BASE = np.uint64(1000003)
RC = str.maketrans("ACGTacgt", "TGCAtgca")
BASES = "TCAG"
AMINO = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
CODON = {a + b + c: AMINO[16 * i + 4 * j + k]
         for i, a in enumerate(BASES) for j, b in enumerate(BASES) for k, c in enumerate(BASES)}


def read_fasta(path):
    op = gzip.open if path.endswith(".gz") else open
    name, buf = None, []
    with op(path, "rt") as fh:
        for line in fh:
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(buf)
                name, buf = line[1:].split()[0], []
            else:
                buf.append(line.strip())
    if name is not None:
        yield name, "".join(buf)


def kmer_hashes(seqs):
    """uint64 hashes of every K-residue window of each sequence (windows never span sequences)."""
    joined = "\0".join(seqs)
    c = CODE[np.frombuffer(joined.encode(), dtype=np.uint8)]
    n = len(c) - K + 1
    if n <= 0:
        return np.zeros(0, np.uint64), np.zeros(0, np.int64)
    h = np.zeros(n, dtype=np.uint64)
    with np.errstate(over="ignore"):
        for k in range(K):
            h = h * BASE + c[k:k + n]
    sep = np.concatenate([[0], np.cumsum(c == 0)])
    ok = (sep[K:K + n] - sep[0:n]) == 0  # no separator inside the window
    owner = np.cumsum(c == 0)[:n]  # sequence index of each window start
    return h[ok], owner[ok]


def translate(cds):
    p = "".join(CODON.get(cds[i:i + 3], "X") for i in range(0, len(cds) - 2, 3))
    return p[:-1] if p.endswith("*") else p


def main():
    ref = sorted({s.rstrip("*").upper() for _, s in read_fasta(sys.argv[1])})
    ref_set = set(ref)
    ref_h = np.unique(kmer_hashes(ref)[0])
    print(f"reference: {len(ref)} distinct proteins, {len(ref_h)} distinct {K}-mers", file=sys.stderr)
    print("label\tn_pred\tcomplete\tmatch30\ttrue_positives\tnonmatching\tdistinct_pred_proteins\t"
          "distinct_tp_proteins\tdistinct_ref_complete")
    for arg in sys.argv[2:]:
        label, rest = arg.split("=", 1)
        fa, pf = rest.split(":")
        tx = dict(read_fasta(fa))
        prots = []
        with open(pf) as fh:
            for line in fh:
                f = line.rstrip("\n").split("\t")
                if f[0] == "transcript_id" or len(f) < 4:
                    continue
                s, e = int(f[1]), int(f[2])
                cds = tx[f[0]][s - 1:e].upper()
                if f[3] == "-":
                    cds = cds.translate(RC)[::-1]
                prots.append(translate(cds))
        complete = [p in ref_set for p in prots]
        h, owner = kmer_hashes(prots)
        hit = np.zeros(len(prots), bool)
        if len(h):
            i = np.searchsorted(ref_h, h)
            i[i == len(ref_h)] = 0
            hit[np.unique(owner[ref_h[i] == h])] = True
        n_c = sum(complete)
        n_m = int(sum(1 for c, x in zip(complete, hit) if x and not c))
        tp = n_c + n_m
        tp_prot = {p for p, c, x in zip(prots, complete, hit) if c or x}
        print(f"{label}\t{len(prots)}\t{n_c}\t{n_m}\t{tp}\t{len(prots) - tp}\t{len(set(prots))}\t"
              f"{len(tp_prot)}\t{len({p for p, c in zip(prots, complete) if c})}")
        sys.stdout.flush()


if __name__ == "__main__":
    main()
