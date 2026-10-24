#!/usr/bin/env python3
"""Reference transcriptomes for the bundled presets: per organism, all complete coding transcripts plus
their coding fragments (IDs suffixed .frag), plus annotated lncRNAs where the benchmark has them (human).
Built from the full annotation, so bundled presets must not be evaluated on the benchmark inputs (the
held-out evaluation of presets is E21).

usage: make_bundled_preset_refs.py <bench inputs dir> <out dir>
"""
import os
import sys

CLADES = {"human": "mammal", "dmel": "insect", "athal": "plant", "spombe": "fungus",
          "pfalciparum": "apicomplexan", "tthermophila": "ciliate", "creinhardtii": "green_alga"}


def read_fasta(p):
    name, buf = None, []
    for line in open(p):
        if line.startswith(">"):
            if name:
                yield name, "".join(buf)
            name, buf = line[1:].split()[0], []
        else:
            buf.append(line.strip())
    if name:
        yield name, "".join(buf)


def main():
    src, out = sys.argv[1], sys.argv[2]
    os.makedirs(out, exist_ok=True)
    for org, clade in CLADES.items():
        recs = list(read_fasta(f"{src}/{org}/complete.fa"))
        n_c = len(recs)
        recs += [(k + ".frag", s) for k, s in read_fasta(f"{src}/{org}/fragments.fa")]
        n_f = len(recs) - n_c
        lnc = f"{src}/{org}/complete.wLnc.fa"
        if os.path.exists(lnc):
            recs += [(k, s) for k, s in read_fasta(lnc) if k.startswith("lnc_")]
        with open(f"{out}/{clade}.fa", "w") as fh:
            for k, s in recs:
                fh.write(f">{k}\n" + "\n".join(s[i:i + 60] for i in range(0, len(s), 60)) + "\n")
        print(f"{clade}\t{org}\tcomplete={n_c}\tfragments={n_f}\tlnc={len(recs) - n_c - n_f}")


if __name__ == "__main__":
    main()
