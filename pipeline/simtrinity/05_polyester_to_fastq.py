#!/usr/bin/env python3
"""Convert one Polyester chunk (sample_01_{1,2}.fasta.gz, strand_specific=TRUE) to RF (dUTP) FASTQ.

usage: 05_polyester_to_fastq.py <chunk_dir> <chunk_tag> <out_prefix>

Polyester (strand_specific = TRUE) writes mate 1 as the first 101 nt of the fragment in transcript
sense and mate 2 as the reverse complement of its last 101 nt (get_reads(); no reverse_complement()
call). The real library is RF (read 1 antisense), so:  R1 = Polyester mate 2,  R2 = Polyester mate 1.

Read names: <chunk_tag>_<n>:<transcript_id>:<frag_start>-<frag_end> (1-based transcript coordinates
of the fragment, from Polyester's mate1 start and mate2 end; '1-T' = the whole transcript, when it
is shorter than the drawn fragment length). The name keeps the true source transcript.
No '/' or spaces, so STAR (readNameSeparator '/') keeps the whole name. Base qualities are a constant
'I' (Polyester writes FASTA; the error model has already been applied to the bases).
Writes <out_prefix>_R1.fastq.gz, <out_prefix>_R2.fastq.gz (via pigz -p 2) and prints the pair count.
"""
import gzip
import subprocess
import sys


def fasta_records(path):
    with gzip.open(path, "rt") as fh:
        name, seq = None, []
        for line in fh:
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(seq)
                name, seq = line[1:].rstrip("\n"), []
            else:
                seq.append(line.strip())
        if name is not None:
            yield name, "".join(seq)


def parse_name(name):
    # read123/ENST...;mate1:a-b;mate2:c-d   or   read123/ENST...;mate1Start:1;mate2Start:1
    num, rest = name.split("/", 1)
    f = rest.split(";")
    tid = f[0]
    if f[1].startswith("mate1:"):
        s = int(f[1][6:].split("-")[0])
        e = int(f[2][6:].split("-")[1])
    else:
        s, e = 1, "T"  # transcript shorter than the drawn fragment: the whole transcript is the fragment
    return num[4:], tid, s, e


def main():
    d, tag, out = sys.argv[1:4]
    outs = []
    for r in ("R1", "R2"):
        fh = open(f"{out}_{r}.fastq.gz", "wb")
        p = subprocess.Popen(["pigz", "-p", "2", "-c"], stdin=subprocess.PIPE, stdout=fh)
        outs.append((p, fh))
    w1, w2 = outs[0][0].stdin, outs[1][0].stdin
    n = 0
    for (n1, s1), (n2, s2) in zip(fasta_records(f"{d}/sample_01_1.fasta.gz"),
                                  fasta_records(f"{d}/sample_01_2.fasta.gz")):
        if n1 != n2:
            sys.exit(f"mate names differ: {n1} {n2}")
        num, tid, s, e = parse_name(n1)
        rn = f"{tag}_{num}:{tid}:{s}-{e}"
        # R1 = Polyester mate 2 (antisense), R2 = Polyester mate 1 (sense)
        w1.write(f"@{rn}\n{s2}\n+\n{'I' * len(s2)}\n".encode())
        w2.write(f"@{rn}\n{s1}\n+\n{'I' * len(s1)}\n".encode())
        n += 1
    for p, fh in outs:
        p.stdin.close()
        p.wait()
        fh.close()
    print(n)


if __name__ == "__main__":
    main()
