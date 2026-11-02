#!/usr/bin/env python3
"""Would choosing a later in-frame ATG by start context help? Release code only (E41).

For each complete-type truth CDS of the stranded complete input, the candidate starts are the in-frame
ATGs from the first one after the upstream in-frame stop (or the transcript start) up to the annotated
stop. Each is scored with the released Kozak model trained on that input (start_refinement.kozak.tsv of
the pinned release run, bits). The first ATG is what the released rule uses when it trims.

Rule tested: move from the first ATG to the highest-scoring later ATG within W codons downstream of the
first ATG when that ATG scores at least M bits above the first ATG.
  gain: annotated start is a later ATG, and the rule moves to it
  loss: annotated start is the first ATG, and the rule moves away from it
  later_wrong: annotated start is a later ATG, but the rule moves to a different later ATG (still wrong)
Also reported: how often the annotated start is a later ATG, and, among first-ATG CDSs, how often some
later ATG within W codons simply scores higher than the first (M = 0).

usage: later_atg_kozak_v1.py > data/start_refinement/later_atg_kozak.tsv
"""
import sys

sys.path.insert(0, "/broad/hptmp/bhaas/TRANSDECODER_PAPER/work/trees/PY_rel_d3892b5")
from pytransdecoder.core import kozak  # noqa: E402

P = "/broad/hptmp/bhaas/TRANSDECODER_PAPER"
RUNS, INPUTS = f"{P}/work/release_pin/bench", f"{P}/data/bench_v3/inputs"
ORGS = "athal creinhardtii dmel human mouse pfalciparum spombe tthermophila".split()
STOPS = {"tthermophila": {"TGA"}}
COMP = {"A": "T", "C": "G", "G": "C", "T": "A"}
WINDOWS = (10, 50, 100000)
MARGINS = (0, 1, 2, 3, 4)


def rc(s):
    return "".join(COMP.get(c, "N") for c in reversed(s))


def fasta(path):
    d, k, b = {}, None, []
    for line in open(path):
        if line[0] == ">":
            if k:
                d[k] = "".join(b)
            k, b = line[1:].split()[0], []
        else:
            b.append(line.strip().upper())
    if k:
        d[k] = "".join(b)
    return d


print("org\twindow_codons\tmargin_bits\tn_cds\tannot_later\tfirst_outscored\tgain\tloss\tlater_wrong\tnet")
for org in ORGS:
    stops = STOPS.get(org, {"TAA", "TAG", "TGA"})
    model = kozak.read_model(f"{RUNS}/{org}/complete/rel/input.fa.transdecoder_dir/start_refinement.kozak.tsv")
    tx = fasta(f"{INPUTS}/{org}/complete.fa")
    cases = []  # (index of annotated ATG among candidates, [(offset_codons_from_first, score)])
    for line in open(f"{INPUTS}/{org}/complete.truth.tsv"):
        f = line.split("\t")
        if f[0] == "transcript_id" or f[7] != "complete" or f[0] not in tx:
            continue
        s, a, b = tx[f[0]], int(f[1]), int(f[2])
        if f[3] == "-":
            s = rc(s)
            a, b = len(s) - b + 1, len(s) - a + 1
        i = a - 1
        if s[i:i + 3] != "ATG":
            continue
        j = i
        while j - 3 >= 0 and s[j - 3:j] not in stops:
            j -= 3
        atgs = [k for k in range(j, b - 3, 3) if s[k:k + 3] == "ATG"]
        if i not in atgs:
            continue
        cases.append((atgs.index(i), [((k - atgs[0]) // 3, kozak.score(model, s, k)) for k in atgs]))
    n = len(cases)
    later = sum(c[0] > 0 for c in cases)
    for w in WINDOWS:
        for m in MARGINS:
            gain = loss = wrong = out = 0
            for idx, cand in cases:
                first = cand[0][1]
                rest = [(sc, t) for t, (off, sc) in enumerate(cand) if t > 0 and off <= w]
                if idx == 0 and rest and max(rest)[0] > first:
                    out += 1
                if not rest:
                    continue
                best, t = max(rest)
                if best < first + m or (m == 0 and best <= first):
                    continue
                if idx == 0:
                    loss += 1
                elif t == idx:
                    gain += 1
                else:
                    wrong += 1
            wl = "all" if w >= 100000 else str(w)
            print(f"{org}\t{wl}\t{m}\t{n}\t{later}\t{out}\t{gain}\t{loss}\t{wrong}\t{gain - loss}")
