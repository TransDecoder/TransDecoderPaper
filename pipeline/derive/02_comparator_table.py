#!/usr/bin/env python3
"""Write the detailed supplementary table of the compared methods (LaTeX longtable).

Features come from data/benchmark/comparator_features.tsv (evidence E27: installed code and papers);
benchmark settings from methods_benchmark.tex; benchmark numbers from data/benchmark/final_v7/table1.tsv
and partial_types.tsv (medians over organisms; 7 for tools without the ciliate code).
usage: 02_comparator_table.py > data/benchmark/comparator_table.tex   (run from the repo root)
"""
import csv
import re
import statistics

FEATURES = [("original_domain", "Designed for"), ("coding_model", "Coding model"), ("training", "Training"),
            ("orf_choice", "CDS choice"), ("multiple_cds_per_transcript", "Several CDSs per sequence"),
            ("start_selection", "Start codon"), ("partial_cds", "Partial CDSs"),
            ("genetic_code_option", "Genetic code"), ("strand", "Strand"),
            ("noncoding_control", "Noncoding control"), ("external_evidence", "External evidence")]
KEY = {"TransDecoder": "pytd_v7", "Prodigal": "prodigal", "GeneMarkS-T": "gmst", "CodAn": "codan", "CPAT": "cpat",
       "CPC2": "cpc2", "RNAsamba": "rnasamba", "Borf": "borf"}
SETTINGS = {
    "TransDecoder": "Release defaults (Perl v7.0.0 and PyTD give identical predictions); ciliate code for T. thermophila.",
    "Prodigal": "Single mode, full motif scan (-n), the organism's translation table.",
    "GeneMarkS-T": "Self-trained on each input; not applicable to T. thermophila.",
    "CodAn": "Pre-trained model of the closest clade; full model for complete inputs, partial model for fragments; not applicable to T. thermophila.",
    "CPAT": "Prebuilt human, mouse or fly model with its published cutoff; human model for other organisms; TGA as the only stop for T. thermophila.",
    "CPC2": "ORF reported for transcripts classified as coding (--ORF); not applicable to T. thermophila.",
    "RNAsamba": "Human-trained full-length or partial weights; stranded inputs only; not applicable to T. thermophila.",
    "Borf": "Default settings; not applicable to T. thermophila.",
    "ANGEL": "Dropped: each of three A. thaliana runs took 17--22 h.",
}
NOTES = {
    "TransDecoder": "Most often first or in the top three of the method ranks; high precision with random sequence.",
    "Prodigal": "Highest stop-level sensitivity, but reports ORFs in many random sequences (score > 0 suffices). Runs partial genes to the sequence ends, so best exact sensitivity on 5' partials and poor start placement on complete transcripts.",
    "GeneMarkS-T": "Self-trained like TransDecoder; highest precision with random sequence, lower fragment sensitivity and exact F1.",
    "CodAn": "Strong on fragments with its partial model; clade model fits AT-rich P. falciparum poorly (complete stop F1 0.710); slowest method (median 5.8 h per input).",
    "CPAT": "High exact F1 on complete transcripts; needs an ATG, so almost no exact 5'-partial starts; precision with random sequence falls on GC-rich C. reinhardtii (0.629).",
    "CPC2": "ORF integrity is a feature, so few calls on 3' partials (stop-level sensitivity 0.247).",
    "RNAsamba": "Highest F1 on complete transcripts; reports ORFs in many random sequences.",
    "Borf": "Longest ATG ORF of at least 100 aa without a coding score.",
}


CITE = {"Prodigal": "Hyatt2010", "GeneMarkS-T": "Tang2015", "CodAn": "Nachtigall2021", "CPAT": "Wang2013",
        "CPC2": "Kang2017", "RNAsamba": "Camargo2020", "Borf": "Signal2021", "ANGEL": "Shimizu2006"}


def tex(s):
    s = re.sub(r"(?<![0-9])--(?=[a-z])", "-{}-", s)
    for a, b in [("\\", "\\textbackslash{}"), ("&", "\\&"), ("%", "\\%"), ("_", "\\_"), ("#", "\\#"),
                 (">=", "$\\ge$"), ("<=", "$\\le$"), (">", "$>$"), ("<", "$<$"), ("~", "$\\sim$"),
                 ("5'", "5$'$"), ("3'", "3$'$"), ("T. thermophila", "\\textit{T.~thermophila}"),
                 ("A. thaliana", "\\textit{A.~thaliana}"), ("P. falciparum", "\\textit{P.~falciparum}"),
                 ("C. reinhardtii", "\\textit{C.~reinhardtii}")]:
        s = s.replace(a, b)
    return s


def main():
    feats = list(csv.DictReader(open("data/benchmark/comparator_features.tsv"), delimiter="\t"))
    t1 = {r["tool"]: r for r in csv.DictReader(open("data/benchmark/final_v7/table1.tsv"), delimiter="\t")}
    pt = {}
    for r in csv.DictReader(open("data/benchmark/final_v7/partial_types.tsv"), delimiter="\t"):
        if r["variant"] != "fragments" or (r["org"] == "tthermophila" and r["tool"] in ("gmst", "codan", "cpc2", "rnasamba", "borf")):
            continue
        pt.setdefault((r["tool"], r["type"]), []).append(float(r["sensitivity_exact"]))
    out = ["\\setlength{\\LTcapwidth}{\\textwidth}", "\\begin{longtable}{@{}>{\\raggedright}p{3.2cm}>{\\raggedright\\arraybackslash}p{12.3cm}@{}}",
           "\\caption{\\textbf{The compared methods in detail.} Features from the installed code and the "
           "methods' papers (evidence E27; ANGEL was not benchmarked). Benchmark values are medians over the "
           "organisms of complete stop / exact F1, fragment stop / exact F1, stop-level precision on complete "
           "transcripts with random sequence added, and exact sensitivity on 3$'$- and 5$'$-partial fragments "
           "(Table~1 and Supplementary Table~S8).}\\label{tab:comparators_detail}\\\\",
           "\\toprule", "Feature & Description \\\\", "\\midrule", "\\endfirsthead",
           "\\toprule", "Feature & Description \\\\", "\\midrule", "\\endhead", "\\bottomrule", "\\endfoot"]
    for r in feats:
        m = r[list(r)[0]]
        cite = f" \\citep{{{CITE[m]}}}" if m in CITE else ""
        out.append(f"\\multicolumn{{2}}{{@{{}}l}}{{\\textbf{{{tex(m)}}}{cite} ({tex(r['version_run'])})}} \\\\*")
        for k, lab in FEATURES:
            out.append(f"{lab} & {tex(r[k])} \\\\*")
        out.append(f"Run here & {tex(SETTINGS[m])} \\\\*" if m in KEY else f"Run here & {tex(SETTINGS[m])} \\\\")
        if m in KEY:
            k = KEY[m]
            v = t1[k]
            f = lambda c: f"{float(v[c]):.3f}"
            e3 = statistics.median(pt[(k, "3prime_partial")])
            e5 = statistics.median(pt[(k, "5prime_partial")])
            out.append(f"Benchmark & complete {f('complete_stop_f1')} / {f('complete_exact_f1')}; fragments "
                       f"{f('fragments_stop_f1')} / {f('fragments_exact_f1')}; precision +random "
                       f"{f('complete_rand_precision')}; exact 3$'$ / 5$'$ partial {e3:.3f} / {e5:.3f} \\\\*")
            out.append(f"Compared & {tex(NOTES[m])} \\\\")
        out.append("\\addlinespace[6pt]" if r is not feats[-1] else "")
    out.append("\\end{longtable}")
    print("\n".join(x for x in out if x))


if __name__ == "__main__":
    main()
