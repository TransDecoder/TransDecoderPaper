#!/usr/bin/env python3
"""Pool raw ORF lengths (gc<TAB>length) from simulation chunks and report quantiles per GC level.

usage: pool_chunk_quantiles.py <rawdir> [--quantiles 0.95,0.99,0.999,0.9999] [--boot N]

With --boot, also reports a bootstrap 95% interval for the 0.999 quantile.
"""
import argparse
import glob
import os
import random
import sys
from collections import defaultdict

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from simulate_random_orf_lengths import quantile  # noqa: E402


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("rawdir")
    ap.add_argument("--quantiles", default="0.95,0.99,0.999,0.9999")
    ap.add_argument("--boot", type=int, default=0)
    args = ap.parse_args()
    qs = [float(q) for q in args.quantiles.split(",")]

    lens = defaultdict(list)
    for path in sorted(glob.glob(os.path.join(args.rawdir, "*.tsv"))):
        with open(path) as fh:
            for line in fh:
                gc, x = line.split()
                lens[int(gc)].append(int(x))

    hdr = "gc\tn_orfs\t" + "\t".join(f"q{q}" for q in qs)
    if args.boot:
        hdr += "\tq0.999_lo95\tq0.999_hi95"
    print(hdr)
    rng = random.Random(1)
    for gc in sorted(lens):
        v = sorted(lens[gc])
        row = [str(gc), str(len(v))] + [f"{quantile(v, q):.1f}" for q in qs]
        if args.boot:
            b = sorted(quantile(sorted(rng.choices(v, k=len(v))), 0.999) for _ in range(args.boot))
            row += [f"{quantile(b, 0.025):.1f}", f"{quantile(b, 0.975):.1f}"]
        print("\t".join(row))


if __name__ == "__main__":
    sys.exit(main())
