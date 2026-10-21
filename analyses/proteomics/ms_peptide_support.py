#!/usr/bin/env python3
"""MS peptide support for predicted ORFs, from published MaxQuant search output.

usage: ms_peptide_support.py --peptides peptides.txt [--peptides ...]
                             --ref ref1.fa[.gz] [--ref ...]
                             --pred <name>=<transdecoder.pep>,<pred_status.tsv> [...]
                             -o out_dir [--prefix ms_peptide_support]

Peptides: MaxQuant peptides.txt rows (already at the authors' 1% PSM/protein FDR). Rows with
Reverse == '+' or Potential contaminant == '+' are dropped. Sequences are normalised I->L and
deduplicated; per normalised peptide the MS/MS counts are summed, the best Score / lowest PEP is kept,
and the experiments (protease_fragmentation) with evidence are listed.

Mapping: exact substring matching (Aho-Corasick, pyahocorasick) of normalised peptides against
normalised predicted proteins (terminal '*' stripped) and against the reference proteome.
A peptide is "novel-specific" if it does not occur in any --ref sequence.

Classes (from pred_status.tsv, distinct representatives only):
  exact if gencode_exact=1; else compatible if gencode_compatible=1; else novel if
  no_gencode_cds_overlap=1; else other_overlap. short = aa < 100.

Tryptic status of a peptide is evaluated in the context of the predicted ORF (Trypsin/P rules,
no proline rule): N-side OK if the preceding residue is K/R or the peptide starts at protein
position 1 (or 2 after an initiator M); C-side OK if the peptide ends in K/R or at the protein end.
full = both sides OK, semi = one, non = neither. Reported for the first occurrence.

Outputs (in out_dir):
  <prefix>.peptide_counts.tsv
  <prefix>.summary.tsv
  <prefix>.<name>.model_hits.tsv        models with >= 1 peptide
  <prefix>.novel_specific_hits.tsv      one row per (model, novel-specific peptide), any class
  <prefix>.novel_with_peptides.tsv      novel-class models with >= 1 peptide: the reference protein
                                        sharing the most peptides (gene from GENCODE-style '|' headers
                                        or UniProt GN=), and whether the whole ORF occurs in a
                                        reference protein (I/L-normalised substring)
"""
import argparse
import csv
import gzip
import os
import sys
from collections import defaultdict

import ahocorasick

csv.field_size_limit(sys.maxsize)
CLASSES = ["exact", "compatible", "other_overlap", "novel"]


def norm(s):
    return s.upper().replace("I", "L")


def opener(p):
    return gzip.open(p, "rt") if p.endswith(".gz") else open(p)


def ref_label(header):
    """short label for a reference FASTA header: GENCODE 'ENSP|ENST|ENSG|..|..|tx_name|gene|len'
    -> 'gene:ENSP'; UniProt 'sp|ACC|ID ... GN=gene' -> 'gene:ACC'."""
    f = header.split()[0].split("|")
    if len(f) >= 7 and f[0].startswith("ENSP"):
        return f"{f[6]}:{f[0]}"
    gn = [w[3:] for w in header.split() if w.startswith("GN=")]
    acc = f[1] if len(f) >= 3 else f[0]
    return f"{gn[0] if gn else '?'}:{acc}"


def read_fasta(path, full_header=False):
    name, chunks = None, []
    with opener(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(chunks)
                hdr = line[1:].strip()
                name, chunks = (hdr if full_header else (hdr.split()[0] if hdr else "")), []
            else:
                chunks.append(line.strip())
    if name is not None:
        yield name, "".join(chunks)


def read_peptides(paths):
    pep = {}  # norm seq -> dict
    counts = defaultdict(int)
    for path in paths:
        with open(path, newline="") as fh:
            rd = csv.reader(fh, delimiter="\t")
            hdr = next(rd)
            ix = {h: i for i, h in enumerate(hdr)}
            exp_cols = [(h[len("Experiment "):], i) for h, i in ix.items() if h.startswith("Experiment ")]
            for r in rd:
                counts["rows"] += 1
                if r[ix["Reverse"]] == "+":
                    counts["reverse"] += 1
                    continue
                if r[ix["Potential contaminant"]] == "+":
                    counts["contaminant"] += 1
                    continue
                counts["kept_rows"] += 1
                seq = r[ix["Sequence"]]
                k = norm(seq)
                d = pep.get(k)
                if d is None:
                    d = pep[k] = {"seqs": set(), "msms": 0, "score": 0.0, "pep": 1.0, "exps": set(),
                                  "mutated": False, "proteins": set()}
                d["seqs"].add(seq)
                d["msms"] += int(r[ix["MS/MS Count"]] or 0)
                try:
                    d["score"] = max(d["score"], float(r[ix["Score"]]))
                    d["pep"] = min(d["pep"], float(r[ix["PEP"]]))
                except ValueError:
                    pass
                for e, i in exp_cols:
                    if r[i] not in ("", "0"):
                        d["exps"].add(e)
                if "Mutated" in ix and r[ix["Mutated"]] == "Yes":
                    d["mutated"] = True
                d["proteins"].update(p for p in r[ix["Proteins"]].split(";") if p)
    counts["distinct_IL_normalised"] = len(pep)
    return pep, counts


def build_automaton(keys):
    A = ahocorasick.Automaton()
    for k in keys:
        A.add_word(k, k)
    A.make_automaton()
    return A


def scan(A, seq):
    """yield (peptide, start0) for every occurrence."""
    for end, k in A.iter(seq):
        yield k, end - len(k) + 1


def tryptic(prot, start, pep_len):
    end = start + pep_len  # exclusive
    n_ok = start == 0 or (start == 1 and prot[0] == "M") or prot[start - 1] in "KR"
    c_ok = end == len(prot) or prot[end - 1] in "KR"
    return {2: "full", 1: "semi", 0: "non"}[n_ok + c_ok]


def read_status(path):
    st = {}
    with open(path) as fh:
        rd = csv.DictReader(fh, delimiter="\t")
        for r in rd:
            if r["gencode_exact"] == "1":
                c = "exact"
            elif r["gencode_compatible"] == "1":
                c = "compatible"
            elif r["no_gencode_cds_overlap"] == "1":
                c = "novel"
            else:
                c = "other_overlap"
            aa = int(r["aa"])
            st[r["model_id"]] = (c, aa, "short" if aa < 100 else "long", r["orf_type"])
    return st


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--peptides", action="append", required=True)
    ap.add_argument("--ref", action="append", required=True)
    ap.add_argument("--pred", action="append", required=True, help="name=pep_fasta,pred_status.tsv")
    ap.add_argument("-o", "--out-dir", required=True)
    ap.add_argument("--prefix", default="ms_peptide_support")
    a = ap.parse_args()
    os.makedirs(a.out_dir, exist_ok=True)
    P = os.path.join(a.out_dir, a.prefix)

    pep, counts = read_peptides(a.peptides)
    print(f"[peptides] {dict(counts)}", file=sys.stderr)
    A = build_automaton(pep)
    print("[automaton] built", file=sys.stderr)

    in_ref = set()
    n_ref = 0
    for path in a.ref:
        for _, s in read_fasta(path):
            n_ref += 1
            for k, _ in scan(A, norm(s.rstrip("*"))):
                in_ref.add(k)
    counts["ref_sequences"] = n_ref
    counts["peptides_in_ref"] = len(in_ref)
    counts["peptides_not_in_ref"] = len(pep) - len(in_ref)
    counts["peptides_not_in_ref_mutated"] = sum(1 for k, d in pep.items() if k not in in_ref and d["mutated"])
    print(f"[ref] {n_ref} sequences; {len(in_ref)} peptides in ref", file=sys.stderr)

    summary = []
    novel_rows = []
    novel_models = []  # (assembly, model_id, aa, normalised seq, peptides)
    for spec in a.pred:
        name, files = spec.split("=", 1)
        pep_fa, status_tsv = files.split(",")
        st = read_status(status_tsv)
        hits = {}  # model -> {pep: (start, tryptic)}
        any_pep_matched = set()
        for mid, s in read_fasta(pep_fa):
            if mid not in st:
                continue
            prot = norm(s.rstrip("*"))
            h = {}
            for k, start in scan(A, prot):
                if k not in h:
                    h[k] = (start, tryptic(prot, start, len(k)))
            if h:
                hits[mid] = h
                any_pep_matched.update(h)
                if st[mid][0] == "novel":
                    novel_models.append((name, mid, st[mid][1], prot, set(h)))
        counts[f"{name}:distinct_models"] = len(st)
        counts[f"{name}:peptides_matched"] = len(any_pep_matched)

        agg = defaultdict(lambda: [0, 0, 0, 0, 0])  # n, >=1, >=2, >=1 novel-specific, >=2 novel-specific
        for mid, (c, aa, sz, _) in st.items():
            h = hits.get(mid, {})
            nspec = sum(1 for k in h if k not in in_ref)
            for key in ((c, sz), (c, "all")):
                v = agg[key]
                v[0] += 1
                v[1] += len(h) >= 1
                v[2] += len(h) >= 2
                v[3] += nspec >= 1
                v[4] += nspec >= 2
        for c in CLASSES:
            for sz in ("short", "long", "all"):
                n, g1, g2, s1, s2 = agg[(c, sz)]
                summary.append([name, c, sz, n, g1, g2, f"{g1 / n:.4f}" if n else "NA",
                                f"{g2 / n:.4f}" if n else "NA", s1, s2])

        with open(f"{P}.{name}.model_hits.tsv", "w") as out:
            out.write("model_id\tclass\taa\tsize\torf_type\tn_peptides\tn_novel_specific\t"
                      "total_msms\tpeptides\n")
            for mid in sorted(hits, key=lambda m: (CLASSES.index(st[m][0]), m)):
                c, aa, sz, ot = st[mid]
                h = hits[mid]
                ns = [k for k in h if k not in in_ref]
                out.write(f"{mid}\t{c}\t{aa}\t{sz}\t{ot}\t{len(h)}\t{len(ns)}\t"
                          f"{sum(pep[k]['msms'] for k in h)}\t{','.join(sorted(h))}\n")
                for k in ns:
                    d = pep[k]
                    novel_rows.append([name, mid, c, aa, ot, k, "|".join(sorted(d["seqs"])), h[k][0] + 1,
                                       h[k][1], d["msms"], f"{d['score']:.2f}", f"{d['pep']:.3g}",
                                       "yes" if d["mutated"] else "no", ";".join(sorted(d["exps"])),
                                       ";".join(sorted(d["proteins"]))[:300]])

    # diagnostics for novel-class models with peptides
    need_peps = set().union(*(m[4] for m in novel_models)) if novel_models else set()
    A2 = build_automaton(need_peps) if need_peps else None
    orf_seqs = {m[3] for m in novel_models}
    A3 = build_automaton(orf_seqs) if orf_seqs else None
    pep2ref = defaultdict(set)
    orf_in_ref = defaultdict(set)
    if A2 is not None:
        for path in a.ref:
            for hdr, s in read_fasta(path, full_header=True):
                rs = norm(s.rstrip("*"))
                lab = None
                for k, _ in scan(A2, rs):
                    lab = lab or ref_label(hdr)
                    pep2ref[k].add(lab)
                for k, _ in scan(A3, rs):
                    orf_in_ref[k].add(lab or ref_label(hdr))
    with open(f"{P}.novel_with_peptides.tsv", "w") as out:
        out.write("assembly\tmodel_id\taa\tn_peptides\tn_novel_specific\ttotal_msms\t"
                  "best_ref_protein\tn_peptides_shared_with_best\tbest_ref_genes\torf_contained_in_ref\n")
        for name, mid, aa, seq, peps in novel_models:
            cnt = defaultdict(int)
            for k in peps:
                for lab in pep2ref.get(k, ()):
                    cnt[lab] += 1
            if cnt:
                mx = max(cnt.values())
                best = sorted(l for l, v in cnt.items() if v == mx)
                genes = sorted({l.split(":")[0] for l in best})
                best_s, genes_s = best[0], ",".join(genes)
            else:
                mx, best_s, genes_s = 0, "NA", "NA"
            contained = sorted({l.split(":")[0] for l in orf_in_ref.get(seq, ())})
            out.write(f"{name}\t{mid}\t{aa}\t{len(peps)}\t{sum(1 for k in peps if k not in in_ref)}\t"
                      f"{sum(pep[k]['msms'] for k in peps)}\t{best_s}\t{mx}\t{genes_s}\t"
                      f"{','.join(contained) if contained else 'no'}\n")

    with open(f"{P}.summary.tsv", "w") as out:
        out.write("assembly\tclass\tsize\tn\tn_ge1_peptide\tn_ge2_peptides\tfrac_ge1\tfrac_ge2\t"
                  "n_ge1_novel_specific\tn_ge2_novel_specific\n")
        for r in summary:
            out.write("\t".join(map(str, r)) + "\n")
    with open(f"{P}.novel_specific_hits.tsv", "w") as out:
        out.write("assembly\tmodel_id\tclass\taa\torf_type\tpeptide_IL\tpeptide_seqs\tstart_in_orf\t"
                  "tryptic_in_orf\tmsms_count\tbest_score\tbest_PEP\tmaxquant_mutated\texperiments\t"
                  "maxquant_proteins\n")
        for r in novel_rows:
            out.write("\t".join(map(str, r)) + "\n")
    with open(f"{P}.peptide_counts.tsv", "w") as out:
        for k, v in counts.items():
            out.write(f"{k}\t{v}\n")
    print(f"[done] {P}.*", file=sys.stderr)


if __name__ == "__main__":
    main()
