import sys, re, collections
best = {}
for l in open(sys.argv[1]):
    x = l.rstrip("\n").split("\t")
    if len(x) != 9 or x[2] != "CDS": continue
    pass
# scores live on mRNA lines; CDS shares coords in PyTD final gff3
for l in open(sys.argv[1]):
    x = l.rstrip("\n").split("\t")
    if len(x) != 9 or x[2] != "mRNA": continue
    sc = float(re.search(r"score=([-\d.]+)", x[8]).group(1)); L = int(x[4]) - int(x[3])
    k = (sc, L)
    if x[0] not in best or k > best[x[0]][0]:
        best[x[0]] = (k, x)
print("transcript_id\tstart\tend\tstrand\tscore")
for t, (k, x) in best.items():
    print(f"{t}\t{x[3]}\t{x[4]}\t{x[6]}\t{k[0]}")
