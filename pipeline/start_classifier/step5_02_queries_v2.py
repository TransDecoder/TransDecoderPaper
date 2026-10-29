#!/usr/bin/env python3
"""E33 step 5: query sets for the UniRef50 searches (identical peptides searched once per organism).

v2: one query file per organism (each search excludes that organism's taxon), i.e.
cand.<org>.pep / k1.<org>.pep; maps carry the organism.

usage: python step5_02_queries_v2.py <out_dir>

Candidate datasets (start-decision isolate; 21): the 5'-open candidate peptides (longest-ORF peptides of the
01_features.py candidates) of sim.human.{d1,d8,d32}, sim.athal.{d1,d16} and bench.<org>.{complete,fragments}
for the 8 organisms. Step-4 cand.pep files are reused; sim.athal.d16 is new (features in
$P/data/start_classifier/step5/features/).
Full-mode datasets (all LongOrfs peptides): sim.athal.{d1,d16}, sim.human.d8 (the a5b24c3 LongOrfs used by
the athal depth series / hstart runs; checked identical to the release runs' LongOrfs), and
bench.{human,athal,dmel,tthermophila}.{complete,fragments,complete.wRand,fragments.wRand} plus
bench.human.{complete,fragments}.wLnc (LongOrfs of the pytd_v7_homology runs).

Writes
  cand.<org>.pep + cand_map.tsv.gz (dataset, org, orf_id, qid): -k 25 search
  k1.<org>.pep   + k1_map.tsv.gz   : -k 1 search = all full-mode LongOrfs + all candidate peptides
  datasets.tsv                 : dataset, kind (cand/full), n peptides, pep source
"""
import gzip
import hashlib
import os
import sys

P = "/broad/hptmp/bhaas/TRANSDECODER_PAPER"
S4 = f"{P}/data/start_classifier/step4/searches"
R = f"{P}/work/bench_v3/runs"
ORGS = ["human", "mouse", "dmel", "athal", "spombe", "tthermophila", "pfalciparum", "creinhardtii"]
out = sys.argv[1]
os.makedirs(out, exist_ok=True)


def read_fasta(path, keep=None):
    seqs, name, buf = {}, None, []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if name is not None and (keep is None or name in keep):
                    seqs[name] = "".join(buf)
                name, buf = line[1:].split()[0], []
            else:
                buf.append(line.strip())
    if name is not None and (keep is None or name in keep):
        seqs[name] = "".join(buf)
    return seqs


cand = {}
for ds in ["sim.human.d1", "sim.human.d8", "sim.human.d32", "sim.athal.d1"] + \
        [f"bench.{o}.{v}" for o in ORGS for v in ("complete", "fragments")]:
    cand[ds] = (f"{S4}/{ds}/cand.pep", None)
# sim.athal.d16: candidate ids from the new feature table, peptides from the release run's LongOrfs
ids = set()
with gzip.open(f"{P}/data/start_classifier/step5/features/sim.athal.d16.tsv.gz", "rt") as fh:
    fh.readline()
    for line in fh:
        ids.add(line.split("\t")[2])
cand["sim.athal.d16"] = (f"{P}/work/simtrinity_athal/depth/d16/runs/pytd_v7/input.fa.transdecoder_dir/longest_orfs.pep", ids)

full = {
    "sim.athal.d1": f"{P}/work/start_classifier/athal_depth/longorfs/d1/input.fa.transdecoder_dir/longest_orfs.pep",
    "sim.athal.d16": f"{P}/work/start_classifier/athal_depth/longorfs/d16/input.fa.transdecoder_dir/longest_orfs.pep",
    "sim.human.d8": f"{P}/work/start_classifier/hstart_runs/sim.human.d8/on/input.fa.transdecoder_dir/longest_orfs.pep",
}
for o in ("human", "athal", "dmel", "tthermophila"):
    for v in ("complete", "fragments", "complete.wRand", "fragments.wRand") + \
            (("complete.wLnc", "fragments.wLnc") if o == "human" else ()):
        full[f"bench.{o}.{v}"] = f"{R}/{o}/{v}/pytd_v7_homology/input.fa.transdecoder_dir/longest_orfs.pep"


def qid(seq):
    return "q" + hashlib.md5(seq.encode()).hexdigest()[:16]


meta = []


def org_of(ds):
    return ds.split(".")[1]


def build(name, sets):
    uniq, rows = {}, []  # uniq[org][qid] = seq
    for ds, (path, keep) in sets.items():
        org = org_of(ds)
        s = read_fasta(path, keep)
        if keep is not None:
            assert len(s) == len(keep), (ds, len(s), len(keep))
        u = uniq.setdefault(org, {})
        for oid, seq in s.items():
            q = qid(seq)
            if q in u:
                assert u[q] == seq
            u[q] = seq
            rows.append(f"{ds}\t{org}\t{oid}\t{q}\n")
        meta.append(f"{ds}\t{name}\t{org}\t{len(s)}\t{path}\n")
    for org, u in uniq.items():
        with open(f"{out}/{name}.{org}.pep", "w") as fh:
            for q, seq in u.items():
                fh.write(f">{q}\n{seq}\n")
        print(name, org, len(u), "unique,", sum(map(len, u.values())), "aa", file=sys.stderr)
    with gzip.open(f"{out}/{name}_map.tsv.gz", "wt") as fh:
        fh.write("dataset\torg\torf_id\tqid\n")
        fh.writelines(rows)
    print(name, len(rows), "peptides", file=sys.stderr)


build("cand", cand)
k1 = {ds: (p, None) for ds, p in full.items()}
for ds, v in cand.items():
    k1[f"{ds}#cand"] = v  # candidates of every dataset (full-mode ones are a subset of their LongOrfs)
build("k1", k1)
with open(f"{out}/datasets.tsv", "w") as fh:
    fh.write("dataset\tquery_set\torg\tn_peptides\tpep_source\n")
    fh.writelines(meta)
