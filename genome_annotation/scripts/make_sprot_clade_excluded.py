#!/usr/bin/env python3
"""Swiss-Prot FASTA with every entry from a taxonomic clade removed (emulating a novel-organism target).

usage: make_sprot_clade_excluded.py uniprot_sprot.fasta.gz nodes.dmp <clade_taxid> > out.fasta

An entry is removed when its OX= taxid has <clade_taxid> in its NCBI lineage (inclusive). Entries with
taxids missing from nodes.dmp are kept and counted. Counts go to stderr.
"""
import gzip
import re
import sys


def main():
    fasta, nodes, clade = sys.argv[1], sys.argv[2], sys.argv[3]
    parent = {}
    for line in open(nodes):
        f = line.split("\t|\t", 2)
        parent[f[0]] = f[1]

    memo = {}

    def in_clade(t):
        path = []
        while t not in memo:
            if t == clade:
                memo[t] = True
                break
            p = parent.get(t)
            if p is None or p == t:
                memo[t] = False
                break
            path.append(t)
            t = p
        for x in path:
            memo[x] = memo[t]
        return memo[t]

    ox = re.compile(r" OX=(\d+)")
    kept = removed = unknown = 0
    keep = False
    out = sys.stdout
    for line in gzip.open(fasta, "rt"):
        if line.startswith(">"):
            t = ox.search(line).group(1)
            if t not in parent:
                unknown += 1
            keep = not in_clade(t)
            kept += keep
            removed += not keep
        if keep:
            out.write(line)
    print(f"clade {clade}: kept {kept}, removed {removed}, taxid not in nodes.dmp {unknown}", file=sys.stderr)


if __name__ == "__main__":
    main()
