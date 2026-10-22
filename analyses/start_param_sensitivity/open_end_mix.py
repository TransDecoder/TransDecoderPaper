#!/usr/bin/env python3
"""Composition of the ORFs that start refinement acts on, per benchmark input.

usage: open_end_mix.py [--grid $P/work/start_param_grid] > results/open_end_mix.tsv

Uses the --no-refine-starts predictions (norefine config). A prediction is 5'-open when its
first codon is not ATG (stranded inputs: + strand only). Open ORFs that match a truth CDS
at the stop are classed as
  truncated     truth is 5' partial and the prediction already starts at the truth 5' end
                (keeping the ORF open is correct)
  utr_extended  the prediction starts upstream of the truth start (trimming to the right ATG
                is needed)
  other         the remainder (e.g. truth start not reachable)
The truncated fraction among (truncated + utr_extended) is the quantity to compare with
real assemblies (GM12878 start labels, record E09: utr_first_atg vs truncated).
"""
import argparse
import collections

P = "/broad/hptmp/bhaas/TRANSDECODER_PAPER"
ORGS = "athal creinhardtii dmel human mouse pfalciparum spombe tthermophila".split()


def fasta(path):
    seqs, name, buf = {}, None, []
    for line in open(path):
        if line.startswith(">"):
            if name:
                seqs[name] = "".join(buf).upper()
            name, buf = line[1:].split()[0], []
        else:
            buf.append(line.strip())
    if name:
        seqs[name] = "".join(buf).upper()
    return seqs


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--grid", default=f"{P}/work/start_param_grid")
    a = ap.parse_args()
    print("org\tinput\tn_pred\topen_5prime\ttruncated\tutr_extended\tother\ttruncated_frac")
    for inp in ("complete", "fragments"):
        tot = collections.Counter()
        for org in ORGS:
            seqs = fasta(f"{P}/data/bench_v3/inputs/{org}/{inp}.fa")
            truth = {}
            with open(f"{P}/data/bench_v3/inputs/{org}/{inp}.truth.tsv") as fh:
                hdr = fh.readline().rstrip("\n").split("\t")
                for line in fh:
                    r = dict(zip(hdr, line.rstrip("\n").split("\t")))
                    truth[r["transcript_id"]] = (int(r["cds_start"]), int(r["cds_end"]), r["strand"], r.get("type", ""))
            c = collections.Counter()
            for line in open(f"{a.grid}/runs/{org}/{inp}/norefine/preds.tsv"):
                if line.startswith("transcript_id"):
                    continue
                tid, s, e, strand = line.split("\t")[:4]
                s, e = int(s), int(e)
                c["n_pred"] += 1
                if strand != "+" or seqs[tid][s - 1:s + 2] == "ATG":
                    continue
                c["open"] += 1
                t = truth.get(tid)
                if not t or t[2] != "+" or abs(e - t[1]) > 2:
                    c["other"] += 1
                elif abs(s - t[0]) <= 2 and t[3] == "5prime_partial":
                    c["truncated"] += 1
                elif s < t[0] - 2:
                    c["utr_extended"] += 1
                else:
                    c["other"] += 1
            tot.update(c)
            frac = c["truncated"] / max(1, c["truncated"] + c["utr_extended"])
            print(f"{org}\t{inp}\t{c['n_pred']}\t{c['open']}\t{c['truncated']}\t{c['utr_extended']}\t{c['other']}\t{frac:.3f}")
        frac = tot["truncated"] / max(1, tot["truncated"] + tot["utr_extended"])
        print(f"ALL\t{inp}\t{tot['n_pred']}\t{tot['open']}\t{tot['truncated']}\t{tot['utr_extended']}\t{tot['other']}\t{frac:.3f}")


if __name__ == "__main__":
    main()
