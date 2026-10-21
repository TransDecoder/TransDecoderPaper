import sys, math, collections
sys.argv = ["x", "data/bench_v3/inputs", "work/bench_v3/runs", sys.argv[1], "complete.wRand", "pytd_dev_single"]
src = open("work/analysis/kozak/sim_kozak.py").read().split("results = []")[0]
exec(src)
U, D = 6, 2
pos, neg = [], []
for tid, st, s, a, b in train:
    pos.append(window(s, a, U, D)); neg.extend(window(s, i, U, D) for i in inframe_atgs(s, a, b))
P, N = build(pos, U + D), build(neg, U + D)
labels = [f"{k}" for k in range(-U, 0)] + [f"+{k}" for k in range(4, 4 + D)]
print("self-trained PWM log2(P/N):  (positions relative to A of ATG = +1)")
for k, lab in enumerate(labels):
    print(f"  {lab:>3} " + " ".join(f"{b}:{math.log2(P[k][b]/N[k][b]):+.2f}" for b in "ACGT"))
# truth-trained on annotated-UTR genes for comparison
tr = {}
for l in open(f"data/bench_v3/inputs/{org}/{variant}.truth.tsv"):
    x = l.rstrip("\n").split("\t")
    if x[0] != "transcript_id" and x[3] == "+" and x[8] == "annotated": tr[x[0]] = (int(x[1]) - 1, int(x[2]))
tpos, tneg = [], []
for tid, (a, b) in tr.items():
    s = seqs[tid]
    if a >= U: tpos.append(window(s, a, U, D)); tneg.extend(window(s, i, U, D) for i in inframe_atgs(s, a, b))
TP, TN = build(tpos, U + D), build(tneg, U + D)
def llr_with(Pm, Nm, s, i): return sum(math.log(Pm[k][c] / Nm[k][c]) for k, c in enumerate(window(s, i, U, D)) if c and c in "ACGT")
# decision cases: ATG-start, stop-correct predictions; fixable = truth is a later in-frame ATG
fix, ok = [], []
for tid, st, s, a, b in preds:
    t = truth_oriented(tid, st, len(s))
    if t is None or abs(t[1] - b) > 2 or s[a:a+3] != "ATG": continue
    atgs = inframe_atgs(s, a, b)
    if abs(t[0] - a) <= 2:
        if atgs: ok.append((s, a, atgs[0]))
    elif t[0] in atgs:
        fix.append((s, a, t[0]))
def auc(xs, ys):
    import bisect
    ys = sorted(ys); n = 0
    for x in xs: n += bisect.bisect_left(ys, x) + 0.5 * (bisect.bisect_right(ys, x) - bisect.bisect_left(ys, x))
    return n / (len(xs) * len(ys))
for name, (Pm, Nm) in [("self", (P, N)), ("truth", (TP, TN))]:
    dfix = [llr_with(Pm, Nm, s, j) - llr_with(Pm, Nm, s, a) for s, a, j in fix]
    dok = [llr_with(Pm, Nm, s, j) - llr_with(Pm, Nm, s, a) for s, a, j in ok]
    print(f"{name}-trained: n_fixable={len(fix)} n_correct={len(ok)}  AUC(delta: fixable vs correct)={auc(dfix, dok):.3f}  "
          f"median delta fixable={sorted(dfix)[len(dfix)//2]:.2f} correct={sorted(dok)[len(dok)//2]:.2f}")
