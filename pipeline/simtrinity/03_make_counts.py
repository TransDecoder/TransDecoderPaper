#!/usr/bin/env python3
"""Read counts and fragment-length parameters for the Polyester simulation (evidence E30).

usage: 03_make_counts.py <salmon_quant_dir> <kept.tsv> <kept.fa> <outdir> <n_chunks>

- counts: round(NumReads) of each kept transcript (coding with complete CDS, lncRNA) from salmon run
  on the full GENCODE v50 transcriptome; transcripts with 0 reads are not simulated.
- fragment length: mean and sd of salmon's fragment-length distribution, as salmon reports them in
  aux_info/meta_info.json (frag_length_mean, frag_length_sd). aux_info/fld.gz (gzipped int32 counts
  for lengths 0..1000; salmon writes a histogram of 10,000 draws from the distribution) is copied to
  fld.tsv, and its mean/sd are reported as fld_sample_mean/sd for comparison.
- chunks: transcripts with reads are dealt to <n_chunks> chunks in decreasing count order (snake
  order), so the chunks have similar read totals; each chunk gets chunk_<k>.fa and
  chunk_<k>.counts.tsv (same order).
Writes <outdir>/counts.tsv (transcript_id class numreads tpm count), fld.tsv, sim_params.tsv and prints
a summary.
"""
import array
import gzip
import json
import math
import os
import sys


def main():
    qdir, kept_tsv, kept_fa, outdir, nch = sys.argv[1:6]
    nch = int(nch)
    os.makedirs(outdir, exist_ok=True)
    q = {}
    with open(os.path.join(qdir, "quant.sf")) as fh:
        fh.readline()
        for line in fh:
            f = line.rstrip("\n").split("\t")
            q[f[0]] = (float(f[3]), float(f[4]))
    kept = []
    with open(kept_tsv) as fh:
        fh.readline()
        for line in fh:
            f = line.rstrip("\n").split("\t")
            kept.append((f[0], f[1]))
    meta = json.load(open(os.path.join(qdir, "aux_info", "meta_info.json")))
    raw = gzip.open(os.path.join(qdir, "aux_info", "fld.gz")).read()
    fld = array.array("i")
    fld.frombytes(raw)
    tot = sum(fld)
    mean = sum(i * c for i, c in enumerate(fld)) / tot
    sd = math.sqrt(sum(c * (i - mean) ** 2 for i, c in enumerate(fld)) / tot)
    with open(os.path.join(outdir, "fld.tsv"), "w") as fh:
        fh.write("length\tcount\n")
        for i, c in enumerate(fld):
            fh.write(f"{i}\t{c}\n")

    rows = []
    missing = 0
    with open(os.path.join(outdir, "counts.tsv"), "w") as fh:
        fh.write("transcript_id\tclass\tnumreads\ttpm\tcount\n")
        for tid, cls in kept:
            if tid not in q:
                missing += 1
                continue
            tpm, nr = q[tid]
            c = int(round(nr))
            fh.write(f"{tid}\t{cls}\t{nr}\t{tpm}\t{c}\n")
            if c > 0:
                rows.append((c, tid))
    rows.sort(key=lambda x: (-x[0], x[1]))
    chunks = [[] for _ in range(nch)]
    for i, r in enumerate(rows):
        k = i % (2 * nch)
        chunks[k if k < nch else 2 * nch - 1 - k].append(r)
    want = {tid: k for k, ch in enumerate(chunks) for c, tid in ch}
    cnt = {tid: c for c, tid in rows}
    fhs = [open(os.path.join(outdir, f"chunk_{k}.fa"), "w") for k in range(nch)]
    order = [[] for _ in range(nch)]
    name = None
    with open(kept_fa) as fh:
        for line in fh:
            if line.startswith(">"):
                name = line[1:].split()[0]
                k = want.get(name)
                if k is not None:
                    order[k].append(name)
            if name in want:
                fhs[want[name]].write(line)
    for k in range(nch):
        fhs[k].close()
        with open(os.path.join(outdir, f"chunk_{k}.counts.tsv"), "w") as fh:
            for tid in order[k]:
                fh.write(f"{tid}\t{cnt[tid]}\n")
    kd = dict(kept)
    summ = {
        "salmon_num_processed": meta.get("num_processed"),
        "salmon_num_mapped": meta.get("num_mapped"),
        "salmon_percent_mapped": meta.get("percent_mapped"),
        "salmon_library_types": ",".join(meta.get("library_types", [])),
        "kept_transcripts": len(kept),
        "kept_missing_in_quant": missing,
        "kept_with_reads": len(rows),
        "coding_with_reads": sum(1 for c, t in rows if kd[t] == "coding"),
        "lncRNA_with_reads": sum(1 for c, t in rows if kd[t] == "lncRNA"),
        "pairs_total": sum(c for c, t in rows),
        "pairs_coding": sum(c for c, t in rows if kd[t] == "coding"),
        "pairs_lncRNA": sum(c for c, t in rows if kd[t] == "lncRNA"),
        "numreads_all_gencode": round(sum(v[1] for v in q.values())),
        "fld_mean": round(meta["frag_length_mean"], 2),
        "fld_sd": round(meta["frag_length_sd"], 2),
        "fld_sample_mean": round(mean, 2),
        "fld_sample_sd": round(sd, 2),
        "fld_sample_n": tot,
    }
    for k in range(nch):
        summ[f"chunk_{k}_pairs"] = sum(c for c, t in chunks[k])
    with open(os.path.join(outdir, "sim_params.tsv"), "w") as fh:
        for k, v in summ.items():
            fh.write(f"{k}\t{v}\n")
            print(k, v, sep="\t")


if __name__ == "__main__":
    main()
