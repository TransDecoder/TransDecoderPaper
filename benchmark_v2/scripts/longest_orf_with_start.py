#!/usr/bin/env python3
"""Longest-ORF baseline that requires a start codon.

Reads a LongOrfs longest_orfs.cds. An ORF that is open at its 5' end (type 5prime_partial or
internal) is trimmed to its first in-frame ATG, and dropped if it has none; ORFs of type complete
or 3prime_partial already begin with ATG. The longest of the resulting ORFs is kept per transcript
and written in the preds.tsv format of convert_predictions.py (transcript_id, start, end, strand,
score; start < end). No minimum length is re-applied after trimming.

usage: longest_orf_with_start.py longest_orfs.cds > longest_orf_atg.preds.tsv
"""
import re
import sys

HDR = re.compile(r"^>(\S+)\s+type:(\S+)(?:\s.*)?\s(\S+):(\d+)-(\d+)\(([+-])\)\s*$")


def records(path):
    hdr, seq = None, []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if hdr:
                    yield hdr, "".join(seq)
                hdr, seq = line, []
            else:
                seq.append(line)
    if hdr:
        yield hdr, "".join(seq)


def main():
    best, order = {}, []
    for hdr, seq in records(sys.argv[1]):
        m = HDR.match(hdr)
        if not m:
            sys.exit(f"unparsed header: {hdr}")
        _, otype, tid, a, b, strand = m.groups()
        a, b = int(a), int(b)
        off = 0
        if otype in ("5prime_partial", "internal"):
            s = seq.upper()
            off = next((i for i in range(0, len(s) - 2, 3) if s[i:i + 3] == "ATG"), None)
            if off is None:
                continue
        if strand == "+":
            lo, hi = a + off, b
        else:
            lo, hi = b, a - off
        if tid not in best:
            order.append(tid)
        if tid not in best or hi - lo > best[tid][1] - best[tid][0]:
            best[tid] = (lo, hi, strand)
    print("transcript_id\tstart\tend\tstrand\tscore")
    for tid in order:
        lo, hi, strand = best[tid]
        print(f"{tid}\t{lo}\t{hi}\t{strand}\t.")


if __name__ == "__main__":
    main()
