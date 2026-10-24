"""Held-out preset test (E21): split each organism's coding transcripts by an ID hash into halves A/B.
A (+ hashed-in lncRNAs for human/mouse) = reference transcriptome for building the preset.
B = source of small test sets (coding + random negatives 1:1, complete and fragments), sizes 50-5000."""
import hashlib, os, random, sys
src, out = sys.argv[1], sys.argv[2]
def read_fasta(p):
    name, buf = None, []
    for line in open(p):
        if line.startswith(">"):
            if name: yield name, "".join(buf)
            name, buf = line[1:].split()[0], []
        else: buf.append(line.strip())
    if name: yield name, "".join(buf)
def half(k):
    base = k[5:] if k.startswith(("rand_", "lnc_")) else k
    return "A" if int(hashlib.md5(base.encode()).hexdigest(), 16) % 2 == 0 else "B"
def write(path, recs):
    with open(path, "w") as fh:
        for k, s in recs:
            fh.write(f">{k}\n" + "\n".join(s[i:i + 60] for i in range(0, len(s), 60)) + "\n")
os.makedirs(f"{out}/ref", exist_ok=True)
for org in ["human", "mouse", "athal", "dmel"]:
    ref = [(k, s) for k, s in read_fasta(f"{src}/{org}/complete.fa") if half(k) == "A"]
    lnc = f"{src}/{org}/complete.wLnc.fa"
    if os.path.exists(lnc):
        ref += [(k, s) for k, s in read_fasta(lnc) if k.startswith("lnc_") and half(k) == "A"]
    write(f"{out}/ref/{org}.A.fa", ref)
    print(org, "reference A:", len(ref), file=sys.stderr)
    for v in ["complete.wRand", "fragments.wRand"]:
        hdr, *rows = open(f"{src}/{org}/{v}.truth.tsv").read().splitlines()
        truth = {r.split("\t")[0]: r for r in rows}
        recs = [(k, s) for k, s in read_fasta(f"{src}/{org}/{v}.fa") if half(k) == "B"]
        coding = [r for r in recs if r[0] in truth]
        rand = [r for r in recs if r[0].startswith("rand_")]
        os.makedirs(f"{out}/tests/{org}", exist_ok=True)
        for n, reps in [(50, 3), (200, 3), (1000, 3), (5000, 1)]:
            for rep in range(1, reps + 1):
                rng = random.Random(f"{org}:{v}:{n}:{rep}")
                pick = rng.sample(coding, n // 2) + rng.sample(rand, n // 2)
                name = f"{v}.B{n}.r{rep}"
                write(f"{out}/tests/{org}/{name}.fa", pick)
                with open(f"{out}/tests/{org}/{name}.truth.tsv", "w") as fh:
                    fh.write(hdr + "\n" + "".join(truth[k] + "\n" for k, _ in pick if k in truth))
