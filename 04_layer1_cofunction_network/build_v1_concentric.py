#!/usr/bin/env python3
#!/usr/bin/env python3
"""Build V1-style nodes/edges from top-300 QQMan hits for concentric annotation plot."""
import pandas as pd, re, urllib.parse, sys
from math import log10
from collections import OrderedDict

QQMAN, ANNO, GFF, GO_TAB, OUT = sys.argv[1:6]

# Parse QQMan top-300
df = pd.read_csv(QQMAN, sep='\t', header=None,
    names=['o','c','p','pv'], dtype={'o':str,'c':str,'p':int,'pv':float})
RE = re.compile(r'LOC\d+|TCONS_\d+')
df['gid'] = df['o'].apply(lambda x: RE.search(x).group(0) if RE.search(x) else x)
df = df.loc[df.groupby('gid')['pv'].idxmin()]
df['nlp'] = -(df['pv'] + 1e-300).apply(log10)
df = df.sort_values('pv').head(300)

# eggNOG
em = OrderedDict()
with open(ANNO) as f:
    for line in f:
        if '##' in line: continue
        h = line.strip().lstrip('#').split('\t'); break
    cm = {n:i for i,n in enumerate(h)}
    for line in f:
        p = line.rstrip('\n').split('\t')
        if len(p)>=8: em[p[0]] = p

# GFF fallback
gff = {}
with open(GFF) as f:
    for line in f:
        if line.startswith('#') or not line.strip(): continue
        p = line.rstrip('\n').split('\t')
        if len(p)<9 or p[2]!='gene': continue
        a = dict(kv.split('=',1) for kv in [x.strip() for x in p[8].split(';') if '=' in x])
        gid = a.get('ID','').replace('gene-','')
        if not gid.startswith('LOC'): gid = a.get('Name',gid)
        gff[gid] = a.get('description',a.get('product',''))

# GO
go_map = {}
with open(GO_TAB) as f:
    f.readline()
    for line in f:
        p = line.rstrip().split('\t')
        if len(p)>=2: go_map[p[0]]=p[1]

ROOT = {'GO:0003674','GO:0008150','GO:0005575'}

COG = {c:n for c,n in [('A','RNA processing'),('B','Chromatin'),('C','Energy'),('D','Cell cycle'),
    ('E','Amino acid'),('F','Nucleotide'),('G','Carbohydrate'),('H','Coenzyme'),
    ('I','Lipid'),('J','Translation'),('K','Transcription'),('L','Replication/repair'),
    ('M','Cell wall'),('N','Cell motility'),('O','PTM/chaperones'),('P','Inorganic ion'),
    ('Q','Secondary metabolites'),('R','General'),('S','Unknown'),
    ('T','Signal transduction'),('U','Trafficking'),('V','Defense'),
    ('W','Extracellular'),('X','Mobilome'),('Y','Nuclear'),('Z','Cytoskeleton')]}

ph = 'chlorophyll_a+b'
vn, ve, seen = [], [], set()

def add_node(nid, tp, dp):
    if nid not in seen:
        seen.add(nid)
        vn.append({'id':nid,'type':tp,'display_name':dp,'significance':'','description':'',
            'source':tp,'chr':'','pos':'','pvalue':'','neg_log10_p':'',
            'COG_category':'','COG_full_name':'','gene_biotype':'','Preferred_name':''})

add_node(ph, 'phenotype', 'Chlorophyll a+b')

for _,r in df.iterrows():
    gid, pv, nlp = r['gid'], r['pv'], r['nlp']
    rec = em.get(gid)
    desc, cog_raw, gos, pfams = '','','',''
    if rec:
        desc = urllib.parse.unquote(rec[cm['Description']]) if cm['Description']<len(rec) and rec[cm['Description']] not in ('-','') else ''
        cog_raw = rec[cm['COG_category']] if cm['COG_category']<len(rec) and rec[cm['COG_category']]!='-' else ''
        gos = rec[cm['GOs']] if cm['GOs']<len(rec) and rec[cm['GOs']]!='-' else ''
        pfams = rec[cm['PFAMs']] if cm['PFAMs']<len(rec) and rec[cm['PFAMs']]!='-' else ''
    elif gid.startswith('LOC') and gid in gff: desc = gff[gid]
    cog_full = '; '.join(f'[{ch}] {COG.get(ch,ch)}' for ch in cog_raw if ch in COG) if cog_raw else ''

    add_node(gid, 'transcript' if gid.startswith('TCONS') else 'gene', gid)
    # update with annotation
    for nd in vn:
        if nd['id']==gid:
            nd.update(significance=f'{nlp:.2f}', description=desc,
                chr=r['c'], pos=r['p'], pvalue=f'{pv:.4e}', neg_log10_p=nlp,
                COG_category=cog_raw, COG_full_name=cog_full)
    ve.append({'source':ph,'target':gid,'weight':nlp,'type':'association','pvalue':pv})

    # Annotation edges
    if desc and len(desc)>5:
        add_node(f'DESC:{desc[:100]}','description',desc[:100])
        ve.append({'source':gid,'target':f'DESC:{desc[:100]}','weight':1,'type':'annotation','pvalue':''})
    if cog_full:
        add_node(f'COG:{cog_full[:100]}','COG_category',cog_full[:100])
        ve.append({'source':gid,'target':f'COG:{cog_full[:100]}','weight':1,'type':'annotation','pvalue':''})
    for pf in pfams.split(',')[:3]:
        pf=pf.strip()
        if pf:
            add_node(f'PFAM:{pf}','PFAM_domain',pf)
            ve.append({'source':gid,'target':f'PFAM:{pf}','weight':1,'type':'annotation','pvalue':''})
    for go_id in gos.split(',')[:3]:
        go_id=go_id.strip()
        if go_id in ROOT: continue
        gd=go_map.get(go_id,go_id)
        add_node(f'GO:{gd[:80]}','GO_term',gd[:80])
        ve.append({'source':gid,'target':f'GO:{gd[:80]}','weight':1,'type':'annotation','pvalue':''})

pd.DataFrame(vn).to_csv(f'{OUT}/v1_nodes.tsv', sep='\t', index=False)
pd.DataFrame(ve).to_csv(f'{OUT}/v1_edges.tsv', sep='\t', index=False)
print(f'V1 nodes: {len(vn)}  edges: {len(ve)}')
