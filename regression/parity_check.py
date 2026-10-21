import sys,os
R=sys.argv[1]
def fasta(p):
    d={};k=None
    for l in open(p):
        l=l.rstrip("\n")
        if l.startswith(">"): k=l[1:].split()[0]; d[k]=[]
        elif k: d[k].append(l.strip())
    return {k:"".join(v) for k,v in d.items()}
def hexs(p):
    d={}
    for l in open(p):
        f=l.split()
        if len(f)>=2:
            try: d[f[0]]=float(f[-1])
            except: pass
    return d
def gffids(p,ft="mRNA"):
    s=set()
    for l in open(p):
        f=l.rstrip("\n").split("\t")
        if len(f)==9 and f[2]==ft:
            s.add((f[0],f[3],f[4],f[6]))
    return s
for ds in ["trinity","pasa","stringtie","supertx","cufflinks"]:
    a=f"{R}/{ds}/perl/{ds}.fa.transdecoder_dir"; b=f"{R}/{ds}/py/{ds}.fa.transdecoder_dir"
    pa,pb=fasta(a+"/longest_orfs.pep"),fasta(b+"/longest_orfs.pep")
    ca,cb=fasta(a+"/longest_orfs.cds"),fasta(b+"/longest_orfs.cds")
    ha,hb=hexs(a+"/hexamer.scores"),hexs(b+"/hexamer.scores")
    hk=set(ha)&set(hb); hmax=max(abs(ha[k]-hb[k]) for k in hk) if hk else None
    ba=gffids(a+"/longest_orfs.cds.best_candidates.gff3"); bb=gffids(b+"/best_candidates.gff3")
    pep_ids_same = list(pa)==list(pb)
    print(ds, f"ORFs perl/py {len(pa)}/{len(pb)}", "ids+order same" if pep_ids_same else "ids/order DIFF",
          "pep seqs same" if pa==pb else f"pep diff {sum(pa.get(k)!=pb.get(k) for k in set(pa)|set(pb))}",
          "cds same(case-insens)" if {k:v.upper() for k,v in ca.items()}=={k:v.upper() for k,v in cb.items()} else "cds DIFF",
          "cds case differs" if ca!=cb else "cds case same",
          f"hexamer keys {len(ha)}/{len(hb)} shared {len(hk)} max|d| {hmax}",
          f"preRefine mRNA {len(ba)}/{len(bb)} same" if ba==bb else f"preRefine DIFF {len(ba^bb)}")
