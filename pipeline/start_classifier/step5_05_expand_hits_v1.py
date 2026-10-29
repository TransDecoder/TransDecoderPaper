#!/usr/bin/env python3
"""E33 step 5: expand the per-organism UniRef50 search outputs (unique-peptide query ids) to per-dataset hit files.

usage: python step5_05_expand_hits_v1.py <step5_dir>
(organisms whose searches have not finished are skipped; rerun when all are done)

Reads <step5_dir>/searches/uniref50/{k25,k1}.<org>.tsv.gz (columns: qseqid sseqid pident length mismatch gapopen
qstart qend sstart send evalue bitscore qlen slen staxids), <step5_dir>/queries/{cand,k1}_map.tsv.gz,
<step5_dir>/db/uniref50_seqinfo.tsv.gz (id, length, first residue, n=, LCA TaxID) and
<step5_dir>/db/uniref50_xml_rep.tsv.gz (id, representative taxid, common taxid, member count, benchmark
organisms among members).
Writes, per dataset, <step5_dir>/hits/uniref50/<dataset>/:
  k25.tsv.gz   candidates, -k 25 hits;  k1.tsv.gz  candidates, -k 1 hits
  (query = ORF id; DIAMOND's line order; + s_first, s_n, s_lca, s_rep_taxid, s_own)
  blastp.outfmt6  (full-mode datasets only) all LongOrfs, -k 1, the 12 standard columns, for PyTD Predict
Also <step5_dir>/hits/uniref50/subject_check.tsv: per organism, hits whose representative taxid is in the
excluded set (must be 0).
"""
import collections
import gzip
import os
import sys

B = sys.argv[1]
S = f"{B}/searches/uniref50"
O = f"{B}/hits/uniref50"
ORG_TAXID = {"human": "9606", "mouse": "10090", "athal": "3702", "dmel": "7227", "spombe": "284812",
             "creinhardtii": "3055", "pfalciparum": "36329", "tthermophila": "312017"}

own = collections.defaultdict(set)
with open(f"{B}/db/own_taxsets.tsv") as fh:
    fh.readline()
    for line in fh:
        t, o = line.split()
        own[o].add(t)

raw = {}
subjects = set()
for kind in ("k25", "k1"):
    for org in ORG_TAXID:
        p = f"{S}/{kind}.{org}.tsv.gz"
        if not os.path.exists(f"{S}/{kind}.{org}.done"):
            print("not finished, skipped:", p, file=sys.stderr)
            continue
        by_q = collections.defaultdict(list)
        with gzip.open(p, "rt") as fh:
            for line in fh:
                x = line.rstrip("\n").split("\t")
                by_q[x[0]].append(x)
                subjects.add(x[1])
        raw[kind, org] = by_q
        print(kind, org, len(by_q), "queries with hits", file=sys.stderr, flush=True)

info = {}
with gzip.open(f"{B}/db/uniref50_seqinfo.tsv.gz", "rt") as fh:
    for line in fh:
        i = line.find("\t")
        if line[:i] in subjects:
            x = line.rstrip("\n").split("\t")
            info[x[0]] = [x[2], x[3], x[4], "NA", "-"]
with gzip.open(f"{B}/db/uniref50_xml_rep.tsv.gz", "rt") as fh:
    for line in fh:
        i = line.find("\t")
        if line[:i] in info:
            x = line.rstrip("\n").split("\t")
            info[x[0]][3] = x[1]
            info[x[0]][4] = x[4]
miss = subjects - set(info)
print(len(subjects), "subjects,", len(miss), "without info", file=sys.stderr)
assert not miss

check = []
for (kind, org), by_q in raw.items():
    bad = sum(1 for lst in by_q.values() for x in lst if info[x[1]][3] in own[ORG_TAXID[org]])
    badtax = sum(1 for lst in by_q.values() for x in lst if x[14] in own[ORG_TAXID[org]])
    n = sum(len(v) for v in by_q.values())
    check.append(f"{kind}\t{org}\t{n}\t{bad}\t{badtax}\n")
os.makedirs(O, exist_ok=True)
with open(f"{O}/subject_check.tsv", "w") as fh:
    fh.write("search\torg\tn_hit_lines\tn_rep_taxid_in_excluded_taxa\tn_staxids_in_excluded_taxa\n")
    fh.writelines(check)


def maps(name):
    m = collections.defaultdict(list)
    with gzip.open(f"{B}/queries/{name}_map.tsv.gz", "rt") as fh:
        fh.readline()
        for line in fh:
            ds, org, oid, q = line.rstrip("\n").split("\t")
            m[ds].append((org, oid, q))
    return m


def write(path, rows, plain12=False):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    op = open if plain12 else gzip.open
    with op(path, "wt") as fh:
        for oid, x in rows:
            if plain12:
                fh.write("\t".join([oid] + x[1:12]) + "\n")
            else:
                fh.write("\t".join([oid] + x[1:] + info[x[1]]) + "\n")


for kind, name in (("k25", "cand"), ("k1", "k1")):
    for ds, lst in maps(name).items():
        if any((kind, org) not in raw for org, _, _ in lst):
            continue
        rows = [(oid, x) for org, oid, q in lst for x in raw[kind, org].get(q, [])]
        if kind == "k25":
            write(f"{O}/{ds}/k25.tsv.gz", rows)
        elif ds.endswith("#cand"):
            write(f"{O}/{ds[:-5]}/k1.tsv.gz", rows)
        else:
            write(f"{O}/{ds}/blastp.outfmt6", rows, plain12=True)
        print(kind, ds, len(lst), "queries", len({r[0] for r in rows}), "with hits", file=sys.stderr, flush=True)
