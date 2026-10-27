#!/usr/bin/env python3
"""Why are simulated Trinity transcripts structurally incompatible? (evidence E30)
For each transcript labelled unresolved/incompatible by label_assembly.py, against all overlapping simulated
transcripts on its strand: novel_intron (an intron in no overlapping simulated transcript),
terminal_diff_le10nt (best match differs by <= 10 exonic bases), introns_ok_exon_diff (introns agree
within the overlap but an exon end differs by > 10 nt), known_introns_new_combination (all introns are
annotated but no single transcript has the combination).
usage: 10_incompatible_breakdown.py <labels_prefix> <ref_dir> <Trinity.minimap2.bam>  > table.tsv
"""
import sys, collections, gzip
import os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import label_assembly as L
LP, REF, BAM = sys.argv[1:4]
lab={}
for l in open(LP + '.labels.tsv'):
    f=l.split('\t')
    if f[2]=='incompatible': lab[f[0]]=1
kept={l.split('\t')[0]:l.split('\t')[1] for l in open(REF + '/kept.tsv')}
idx=L.Index()
with gzip.open(REF + '/gencode_v50.tx.tsv.gz', 'rt') as fh:
    fh.readline()
    for line in fh:
        f=line.rstrip('\n').split('\t')
        if f[0] in kept:
            ex=L.segs(f[8]); r=dict(tid=f[0],exons=ex,s=ex[0][0],e=ex[-1][1])
            idx.add(f[6],f[7],r['s'],r['e'],r)
aln,supp=L.read_bam(BAM)
def introns(ex): return [(ex[i][1],ex[i+1][0]) for i in range(len(ex)-1)]
cat=collections.Counter()
for t in lab:
    a=aln[t]; T=a['exons']; ts,te=T[0][0],T[-1][1]
    ti=introns(T)
    best=None
    allintr=set()
    for r in idx.query(a['chrom'],a['strand'],ts,te):
        if r['s']>te or r['e']<ts: continue
        ri=set(introns(r['exons'])); allintr|=ri
        lo,hi=max(ts,r['s']),min(te,r['e'])
        ct,cr=L.clip(T,lo,hi),L.clip(r['exons'],lo,hi)
        a1=set(p for s,e in ct for p in range(s,e+1)); b1=set(p for s,e in cr for p in range(s,e+1))
        d=len(a1^b1)
        intr_ok=all(i in ri for i in ti if lo<=i[0] and i[1]<=hi)
        if best is None or d<best[0]: best=(d,intr_ok,r['tid'])
    if best is None: cat['no_cand']+=1; continue
    d,intr_ok,_=best
    novel=[i for i in ti if i not in allintr]
    if novel: k='novel_intron(%s)'%('1' if len(novel)==1 else '2+')
    elif d<=10: k='terminal_diff_le10nt'
    elif intr_ok: k='introns_ok_exon_diff>10nt(retained/extended_end)'
    else: k='known_introns_new_combination'
    cat[k]+=1
print('category\tn')
for k,v in cat.most_common(): print(k,v,sep='\t')
