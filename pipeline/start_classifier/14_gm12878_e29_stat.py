#!/usr/bin/env python3
"""E33 adoption: E29 5'-end statistic of GM12878 predictions against the transcript-level MANE truth.

usage: python 14_gm12878_e29_stat.py <truth.tsv> <label>=<preds.tsv> ...

truth.tsv: 03_utr_gm12878_truth.py output (transcripts carrying a MANE Select stop with the reference
structure to their 5' end; '+' = transcript orientation; 5' partials start at the first whole codon).
preds.tsv: convert_predictions.py output of the run's transdecoder.gff3.
Per reference type: references; references whose stop is predicted (+-2 nt, same strand); of those,
exact 5' end (+-2 nt); for 5' partials, trimmed (5' end downstream of the transcript start: the
E29 error); for complete references, 5' end upstream of the true start (left open through it or
trimmed to an upstream ATG) and downstream of it.
"""
import collections
import sys

TOL = 2


def ends(s, e, strand):
    return (s, e) if strand == "+" else (e, s)


truth = {}
with open(sys.argv[1]) as fh:
    hdr = fh.readline().rstrip("\n").split("\t")
    for line in fh:
        r = dict(zip(hdr, line.rstrip("\n").split("\t")))
        truth[r["transcript_id"]] = (int(r["cds_start"]), int(r["cds_end"]), r["strand"], r["type"])
print("label\tref_type\tn_ref\tn_stop\texact\tfive_prime_downstream\tfive_prime_upstream\texact_given_stop")
for arg in sys.argv[2:]:
    label, pf = arg.split("=", 1)
    preds = collections.defaultdict(list)
    with open(pf) as fh:
        fh.readline()
        for line in fh:
            x = line.rstrip("\n").split("\t")
            preds[x[0]].append((int(x[1]), int(x[2]), x[3]))
    c = collections.defaultdict(collections.Counter)
    for tid, (s, e, strand, typ) in truth.items():
        k = c[typ]
        k["n"] += 1
        t5, t3 = ends(s, e, strand)
        hit = None
        for p in preds.get(tid, []):
            if p[2] != strand:
                continue
            p5, p3 = ends(*p)
            if abs(p3 - t3) <= TOL:
                if hit is None or abs(p5 - t5) < abs(hit - t5):
                    hit = p5
        if hit is None:
            continue
        k["stop"] += 1
        # distance along the transcript in the CDS direction (positive = downstream)
        dd = (hit - t5) if strand == "+" else (t5 - hit)
        if abs(dd) <= TOL:
            k["exact"] += 1
        elif dd > 0:
            k["down"] += 1
        else:
            k["up"] += 1
    for typ in sorted(c):
        k = c[typ]
        print(f"{label}\t{typ}\t{k['n']}\t{k['stop']}\t{k['exact']}\t{k['down']}\t{k['up']}\t"
              f"{k['exact'] / k['stop']:.4f}")
