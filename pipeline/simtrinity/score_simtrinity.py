#!/usr/bin/env python3
"""Score the ORF predictors on the simulated GM12878 Trinity assembly (evidence E30).

usage: score_simtrinity.py <input.fa> <labels_prefix> <runs_dir> <outdir>
  <labels_prefix>.{truth.tsv,labels.tsv,negatives.txt,excluded.txt} from label_assembly.py
  <runs_dir>/<tool>/preds.tsv (run_tool.sh); the longest-ORF baseline is pytd_v7/longest_orf.preds.tsv

Scoring convention
- Scored set = truth transcripts (complete / 5' / 3' / internal) + negatives (lncRNA-only contigs, not in
  the truth file). Predictions on the excluded transcripts (unresolved, ambiguous, short_*) are removed
  before scoring, so they count neither as true nor as false positives; they are reported separately.
- Transcript IDs in preds.tsv are cut at the first whitespace (GeneMarkS-T writes the whole header).
- The remaining predictions are scored with benchmark_v2/scripts/score_predictions.py (stop within 2 nt;
  exact = start also within 2 nt; one credit per truth CDS; every other prediction, including any on a
  negative, is a false positive and counts in precision).
Writes to <outdir>:
  scores_overall.tsv     tool + score_predictions.py metrics + runtime_sec + n_pred_excluded
  scores_by_type.tsv     per tool x truth type: n, stop / exact sensitivity
  scores_by_expr.tsv     per tool x TPM bin of the source transcript(s): n, stop / exact sensitivity
  negatives.tsv          per tool: negatives with >= 1 prediction (fraction), predictions on negatives
  excluded_preds.tsv     per tool x excluded label (unresolved split by reason): transcripts, transcripts with a prediction, predictions
  pred_types.tsv         per tool: predicted ORFs by their own ends (complete / 5prime_partial /
                         3prime_partial / internal), on the scored set and on all transcripts
"""
import collections
import os
import subprocess
import sys

T = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..")
sys.path.insert(0, os.path.join(T, "benchmark_v2", "scripts"))
from score_predictions import TOL, ends, read_preds  # noqa: E402

TOOLS = ["pytd_v7", "perl_td_v7", "perl_td_v600", "prodigal", "gmst", "cpat", "cpc2", "rnasamba", "borf",
         "longest_orf"]
TYPES = ["complete", "5prime_partial", "3prime_partial", "internal"]
EXPR = [(0, 1, "<1"), (1, 5, "1-5"), (5, 20, "5-20"), (20, 100, "20-100"), (100, float("inf"), ">100")]
STOPS = {"TAA", "TAG", "TGA"}
COMP = str.maketrans("ACGTN", "TGCAN")


def read_fasta(path):
    seqs, name, buf = {}, None, []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if name:
                    seqs[name] = "".join(buf).upper()
                name, buf = line[1:].split()[0], []
            else:
                buf.append(line.strip())
    if name:
        seqs[name] = "".join(buf).upper()
    return seqs


def read_tsv(path):
    with open(path) as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        return [dict(zip(hdr, line.rstrip("\n").split("\t"))) for line in fh]


def expr_bin(tpm):
    return next(lab for lo, hi, lab in EXPR if lo <= tpm < hi)


def best_match(t, plist):
    best = None
    t5, t3 = ends(int(t["cds_start"]), int(t["cds_end"]), t["strand"])
    for p in plist:
        if p[2] != t["strand"]:
            continue
        p5, p3 = ends(*p)
        if abs(p3 - t3) <= TOL:
            lvl = "exact" if abs(p5 - t5) <= TOL else "stop"
            if best != "exact":
                best = lvl
    return best


def orf_type(seq, s, e, strand):
    c = seq[s - 1:e]
    if strand == "-":
        c = c.translate(COMP)[::-1]
    st, sp = c[:3] == "ATG", c[-3:] in STOPS
    return {(True, True): "complete", (False, True): "5prime_partial",
            (True, False): "3prime_partial", (False, False): "internal"}[(st, sp)]


def main():
    fa, lp, runs, outdir = sys.argv[1:5]
    os.makedirs(outdir, exist_ok=True)
    seqs = read_fasta(fa)
    truth = read_tsv(lp + ".truth.tsv")
    negs = {l.strip() for l in open(lp + ".negatives.txt") if l.strip()}
    excl = {}
    for line in open(lp + ".excluded.txt"):
        tid, lab, reason = line.rstrip("\n").split("\t")
        excl[tid] = f"unresolved:{reason}" if lab == "unresolved" else lab
    tset = {t["transcript_id"] for t in truth}
    assert not (tset & negs) and not (tset & set(excl)) and not (negs & set(excl))
    assert len(tset) + len(negs) + len(excl) == len(seqs), "labels do not cover the input FASTA"
    scorer = os.path.join(T, "benchmark_v2", "scripts", "score_predictions.py")
    tmp = os.path.join(outdir, "_filtered_preds")
    os.makedirs(tmp, exist_ok=True)

    overall, by_type, by_expr, neg_rows, exc_rows, type_rows = [], [], [], [], [], []
    hdr = None
    excl_n = collections.Counter(excl.values())
    for tool in TOOLS:
        pf = os.path.join(runs, "pytd_v7", "longest_orf.preds.tsv") if tool == "longest_orf" \
            else os.path.join(runs, tool, "preds.tsv")
        if not os.path.exists(pf):
            continue
        rtf = os.path.join(runs, tool, "runtime_sec.txt")
        rt = open(rtf).read().strip() if tool != "longest_orf" and os.path.exists(rtf) else "NA"
        # transcript IDs are cut at the first whitespace (GeneMarkS-T reports the whole FASTA header line)
        norm = os.path.join(tmp, tool + ".all_preds.tsv")
        with open(pf) as fh, open(norm, "w") as out:
            for line in fh:
                if not line.strip() or line.startswith(("#", "transcript_id\t")):
                    continue
                x = line.split("\t", 1)
                out.write(x[0].split()[0] + "\t" + x[1])
        preds = read_preds(norm)
        kept_lines, n_exc = [], 0
        with open(norm) as fh:
            for line in fh:
                if line.split("\t", 1)[0] in excl:
                    n_exc += 1
                else:
                    kept_lines.append(line)
        ff = os.path.join(tmp, tool + ".preds.tsv")
        with open(ff, "w") as fh:
            fh.writelines(kept_lines)
        out = subprocess.run([sys.executable, scorer, "--truth", lp + ".truth.tsv", "--preds", ff, "--label", tool],
                             capture_output=True, text=True, check=True).stdout.strip().split("\n")
        hdr = "tool\t" + out[0].split("\t", 1)[1] + "\truntime_sec\tn_pred_excluded"
        overall.append(f"{tool}\t" + out[1].split("\t", 1)[1] + f"\t{rt}\t{n_exc}")

        acc_t = collections.defaultdict(lambda: [0, 0, 0])
        acc_e = collections.defaultdict(lambda: [0, 0, 0])
        for t in truth:
            m = best_match(t, preds.get(t["transcript_id"], []))
            for acc, key in ((acc_t, t["type"]), (acc_e, expr_bin(float(t["tpm"])))):
                acc[key][0] += 1
                acc[key][1] += m is not None
                acc[key][2] += m == "exact"
        for key in TYPES:
            n, s, x = acc_t[key]
            by_type.append(f"{tool}\t{key}\t{n}\t{s}\t{x}\t{s / n:.4f}\t{x / n:.4f}")
        for _, _, key in EXPR:
            n, s, x = acc_e[key]
            if n:
                by_expr.append(f"{tool}\t{key}\t{n}\t{s}\t{x}\t{s / n:.4f}\t{x / n:.4f}")
        nneg_hit = sum(1 for tid in negs if tid in preds)
        npred_neg = sum(len(preds[tid]) for tid in negs if tid in preds)
        neg_rows.append(f"{tool}\t{len(negs)}\t{nneg_hit}\t{nneg_hit / len(negs):.4f}\t{npred_neg}")
        for lab in sorted(excl_n):
            ids = [tid for tid, l in excl.items() if l == lab]
            hit = sum(1 for tid in ids if tid in preds)
            npr = sum(len(preds[tid]) for tid in ids if tid in preds)
            exc_rows.append(f"{tool}\t{lab}\t{len(ids)}\t{hit}\t{hit / len(ids):.4f}\t{npr}")
        for scope in ("scored", "all"):
            c = collections.Counter()
            for tid, pl in preds.items():
                if scope == "scored" and tid in excl:
                    continue
                for s, e, strand in pl:
                    c[orf_type(seqs[tid], s, e, strand)] += 1
            n = sum(c.values())
            type_rows.append(f"{tool}\t{scope}\t{n}\t" + "\t".join(
                f"{c[k]}\t{c[k] / n:.4f}" if n else "0\tNA" for k in TYPES))

    def write(name, header, rows):
        with open(os.path.join(outdir, name), "w") as fh:
            fh.write(header + "\n" + "\n".join(rows) + "\n")

    write("scores_overall.tsv", hdr, overall)
    sh = "n\ttp_stop\ttp_exact\tsensitivity_stop\tsensitivity_exact"
    write("scores_by_type.tsv", "tool\ttype\t" + sh, by_type)
    write("scores_by_expr.tsv", "tool\ttpm_bin\t" + sh, by_expr)
    write("negatives.tsv", "tool\tn_negatives\tnegatives_with_pred\tfrac_negatives_with_pred\tn_pred_on_negatives",
          neg_rows)
    write("excluded_preds.tsv", "tool\tlabel\tn_transcripts\twith_pred\tfrac_with_pred\tn_pred", exc_rows)
    write("pred_types.tsv", "tool\tscope\tn_pred\t" + "\t".join(f"{k}\t{k}_frac" for k in TYPES), type_rows)


if __name__ == "__main__":
    sys.exit(main())
