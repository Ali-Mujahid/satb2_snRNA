"""Build public evidence tables from original study files. Arguments in README."""
import argparse,csv,gzip,hashlib,json,math,sqlite3,statistics
from collections import Counter,defaultdict
from pathlib import Path
from lxml import html
import openpyxl

parser=argparse.ArgumentParser()
for name in ('workspace','sources','output'): parser.add_argument('--'+name,type=Path,required=True)
args=parser.parse_args()
src=args.sources; out=args.output; out.mkdir(exist_ok=True,parents=True)
ws=args.workspace; prior=ws/'SATB2_CA1_EPHYS_PLANNING_20261125'
def read(p):
    with (gzip.open(p,'rt') if p.suffix=='.gz' else p.open()) as f: return list(csv.DictReader(f,delimiter='\t'))
def tab(name,rows,fields):
    p=out/name
    with (gzip.open(p,'wt') if p.suffix=='.gz' else p.open('w')) as f:
        w=csv.DictWriter(f,fields,delimiter='\t',lineterminator='\n',extrasaction='ignore')
        w.writeheader();w.writerows(rows)
def num(v):
    try:
        v=float(v);return v if math.isfinite(v) else None
    except (ValueError,TypeError):return None
def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()
def showep(x):
    return {'CA1_Superficial':'Superficial','CA1_Deep':'Deep'}.get(x,x.replace('IN_','IN '))
config=json.loads((ws/'SATB2_CA1_IN_EXPANDED_20261102/workflow/config/analysis.json').read_text())
pairs=read(prior/'strict_DEG_ligand_receptor_pairs.tsv')
assert len(pairs)==14 and all('_' not in v for p in pairs for v in (p['ligand'],p['receptor']))
pairgenes=sorted({p[k] for p in pairs for k in ('ligand','receptor')})
panel='Grid2ip Cacna1g Stac Dlgap4 Dagla Foxp4 Grik4 Zfp385c Mir9-3hg Abr Ryr1 Cacna2d2 Arhgap31 Ablim2 Hpcal1 Col5a1 Myo10 Arap1 Plxdc1 Pde9a 9530059O14Rik Ppp1r16b Plekhd1 Clmn Bcl11a Sh3pxd2a Gm20635 Cntnap5c Col18a1 Srgap1'.split()
assert len(set(panel))==30
endpoints=['CA1_Superficial','CA1_Deep']+['IN_C'+str(i) for i in range(1,21)]
degs=[];cm=[];stats=[];tested={};reference_ids=defaultdict(set); endpoint_files={}
for pop in endpoints:
    cfg=config['endpoints'][pop]
    if cfg['status']!='TESTED':
        tested[pop]=None;continue
    if pop=='CA1_Superficial':p=prior/'superficial_C2_C3_DESeq2.tsv'
    else:
        p=Path(cfg['stable_audit_path']); assert sha(p)==cfg['sha256']
    endpoint_files[pop]=sha(p)
    rr=read(p)
    sym=lambda r:r.get('gene_symbol',r.get('source_gene_symbol',''))
    for r in rr:
        if sym(r) not in ('NA',''):reference_ids[sym(r)].add(r['gene_id'])
    freq=Counter(map(sym,rr))
    valid={sym(r):r for r in rr if num(r['padj']) is not None and sym(r) not in ('NA','') and freq[sym(r)]==1}
    called={g for g,r in valid.items() if float(r['padj'])<.05}
    for r in rr:
        if num(r['padj']) is not None and float(r['padj'])<.05:
            degs.append(dict(endpoint=pop,gene=sym(r),gene_id=r['gene_id'],log2FC=float(r['log2FoldChange']),FDR=float(r['padj'])))
    for gene in panel:
        r=valid.get(gene); matched=[x for x in rr if sym(x)==gene]
        status='DEG' if gene in called else 'tested' if r else 'FDR unavailable or ambiguous' if matched else 'outside fitted gene list'
        cm.append(dict(endpoint=pop,gene=gene,reference_figure_order=panel.index(gene)+1,status=status,log2FC=float(r['log2FoldChange']) if r else None,FDR=float(r['padj']) if r else None))
    N=len(valid);K=len(set(panel)&valid.keys());D=len(called);hits=set(panel)&called;A=len(hits)
    # Hypergeometric upper tail with a panel of at most 30 genes.
    prob=sum(math.comb(D,i)*math.comb(N-D,K-i) for i in range(A,min(K,D)+1) if K-i<=N-D)/math.comb(N,K)
    stats.append(dict(endpoint=pop,tested_unique_symbols=N,DEG_unique_symbols=D,panel_tested=K,panel_DEG=A,up=sum(float(valid[g]['log2FoldChange'])>0 for g in hits),down=sum(float(valid[g]['log2FoldChange'])<0 for g in hits),p=prob,genes=', '.join(sorted(hits))))
    tested[pop]=sum(num(r['padj']) is not None and float(r['padj'])<.05 for r in rr)
for pop in endpoints:
    if tested[pop] is None:
        for gene in panel:cm.append(dict(endpoint=pop,gene=gene,reference_figure_order=panel.index(gene)+1,status='no eligible model',log2FC=None,FDR=None))
ordered=sorted(stats,key=lambda x:x['p']);q=1.
for i in range(len(ordered)-1,-1,-1):
    q=min(q,ordered[i]['p']*len(ordered)/(i+1));ordered[i]['BH_FDR']=q
genes={r['gene'] for r in degs if r['gene'] not in ('','NA')}|set(pairgenes)
byupper=defaultdict(set)
for g in reference_ids:byupper[g.upper()].add(g)
synd=[]; byfrac=defaultdict(list); areas=['SO','SP','SR','SLM']
w=openpyxl.load_workbook(src/'41467_2025_63119_MOESM4_ESM.xlsx',data_only=True,read_only=True)
for sheet,frac in [('strata syn FCs','FASS'),('strata tissue FCs','Tissue')]:
    rows=iter(w[sheet].values)
    for row in rows:
        if row[0]=='Protein.ID':keys=row;break
    for row in rows:
        r=dict(zip(keys,row));names={v.strip().upper() for v in str(r['Genes']).split(';') if v.strip()}
        for key in names:
            mapped=byupper.get(key,set()) & genes
            for g in sorted(mapped):
                unambig=len(names)==1 and len(byupper[key])==1 and len(reference_ids[g])==1
                for area in areas:
                    fc=num(r['log2fc_'+area]);qval=num(r['qvalue_'+area])
                    item=dict(gene=g,fraction=frac,stratum=area,protein_group=r['Protein.ID'],reported_symbols=r['Genes'],unique_peptides=r['unique_peptides'],unique_gene=unambig,log2FC=fc,q=qval,author_significant=(r['sig_'+area] is True) if fc is not None else None)
                    synd.append(item)
                    if unambig and fc is not None:byfrac[g,frac,area].append(item)
w.close()
synplot=[]
for g in sorted(genes):
    for fr in ['FASS','Tissue']:
        for a in areas:
            rr=byfrac.get((g,fr,a),[])
            if rr:
                synplot.append(dict(gene=g,fraction=fr,stratum=a,log2FC=statistics.median(r['log2FC'] for r in rr),groups=len(rr),significant_groups=sum(r['author_significant'] is True for r in rr)))
db=sqlite3.connect(next(src.glob('*.db.sqlite')).as_uri()+'?mode=ro',uri=True)
db.row_factory=sqlite3.Row
synrows=[dict(r) for r in db.execute('SELECT DISTINCT MouseName AS gene, HumanName AS human_symbol, Localisation AS compartment, SpeciesTaxID AS species_taxid, BrainRegion AS region, PaperPMID AS pmid, Paper AS study FROM FullGeneFullPaperFullRegion') if r['gene'] in genes]
native=read(src.parent/'native_synaptome_R_lookup.tsv')
assert {(r['MouseName'],r['Localisation'],r['BrainRegion'],str(r['PaperPMID']),str(r['SpeciesTaxID'])) for r in native}=={(r['gene'],r['compartment'],r['region'],str(r['pmid']),str(r['species_taxid'])) for r in synrows if r['gene'] in pairgenes}
ppi={}
for r in db.execute('SELECT a.MouseName AS g1,b.MouseName AS g2,p.type,p.method,pp.PMID AS pmid FROM PPI p JOIN Gene a ON a.ID=p.A JOIN Gene b ON b.ID=p.B LEFT JOIN PaperPPI pp ON pp.PPID=p.ID'):
    if r['g1'] in pairgenes and r['g2'] in pairgenes:
        k=tuple(sorted([r['g1'],r['g2']]));ppi.setdefault(k,[]).append(dict(r))
idnames=dict(db.execute('select ID,MouseName from Gene'))
native_edges={tuple(sorted([idnames[int(r['A'])],idnames[int(r['B'])]])) for r in read(src.parent/'native_synaptome_R_ppi.tsv')}
assert native_edges==set(ppi)
db.close()
def htmlrows(t):
    return [[el.text_content().strip() for el in tr.xpath('./td|./th')] for tr in t.xpath('.//tr')]
g2=[];mgi={};id_map={i:g for g in reference_ids if len(reference_ids[g])==1 for i in reference_ids[g]}
bay=html.parse(str(src/'human-mouse-psp.html'))
st,=bay.xpath('//table[./thead//th[normalize-space(.)="Detected Mouse"]]')
assert [e.text_content() for e in st.xpath('./thead//th')][1:4]==['Detected Mouse','Detected Human (PMID:21170055)','Appr Gene (m)']
def member(g,collection,source,value):
    if g in genes:g2.append(dict(gene=g,collection=collection,reference_table=source,status=value))
for r in htmlrows(st):
    g=id_map.get(r[7],r[3])
    if r[4].startswith('MGI:'):mgi[r[4]]=g
    for col,coll in [(1,'Mouse cortical PSD'),(2,'Human cortical PSD ortholog')]:
        if r[col]=='YES':member(g,coll,'Bayes et al. 2012 Table S1','listed')
        elif r[col]=='AMBIGUOUS':member(g,coll,'Bayes et al. 2012 Table S1','ambiguous')
        else:member(g,coll,'Bayes et al. 2012 Table S1','unmarked')
for t in bay.xpath('//table'):
    if 'mGluR5 Complex' not in [i.text_content() for i in t.xpath('./thead//th')]:continue
    for r in htmlrows(t):
        for idx,coll in [(2,'PSD95 complex'),(3,'NR / MASC'),(4,'mGluR5 complex')]:
            member(mgi.get(r[0],''),coll,'Bayes et al. 2012 Table S2d','listed' if r[idx]=='YES' else 'unmarked')
tap=html.parse(str(src/'tap-psd-95.html'))
for t in tap.xpath('//table'):
    rows=htmlrows(t)
    if len(rows)<3 or 'PSD-95 Core Complex' not in rows[1]:continue
    assert rows[1][6:9]==['MASC/NRC','AMPA','mGluR5']
    for r in rows[2:]:
        g=mgi.get(r[2],r[0])
        for i,coll in [(5,'PSD95 core'),(6,'NR / MASC'),(7,'AMPA complex'),(8,'mGluR5 complex')]:
            member(g,coll,'Fernandez et al. 2009 Table S1','listed' if r[i].upper()=='X' else 'unmarked')
group_members=defaultdict(set)
for r in g2:group_members[r['gene'],r['collection']].add(r['status'])
gcols=['Mouse cortical PSD','Human cortical PSD ortholog','PSD95 core','PSD95 complex','NR / MASC','AMPA complex','mGluR5 complex']
def gstate(g,c):
    s=group_members.get((g,c),set());return next((v for v in ['listed','ambiguous','unmarked'] if v in s),'no eligible row')
by_gene=defaultdict(list)
for r in synrows:by_gene[r['gene']].append(r)
def areas_for(g,fr):return {a for a in areas if byfrac.get((g,fr,a))}
for r in degs:
    g=r['gene']
    r['SynDive_FASS_strata']=','.join(a for a in areas if a in areas_for(g,'FASS'))
    r['SynDive_tissue_strata']=','.join(a for a in areas if a in areas_for(g,'Tissue'))
    r['Synaptome_mouse_hippocampus_PMIDs']=','.join(map(str,sorted({x['pmid'] for x in by_gene[g] if x['region']=='Hippocampus' and x['species_taxid']==10090})))
    r['G2C_lists']=' | '.join(c for c in gcols if gstate(g,c)=='listed')
    r['Connectome_seq_30_marker_panel']=g in panel
for r in pairs:
    g,h=r['ligand'],r['receptor']
    r['FASS_common_strata']=','.join(a for a in areas if a in areas_for(g,'FASS')&areas_for(h,'FASS'))
    r['mouse_hippocampus_both']=all(any(x['species_taxid']==10090 and x['region']=='Hippocampus' for x in by_gene[z]) for z in (g,h))
    r['G2C_shared_lists']=' | '.join(c for c in gcols if gstate(g,c)==gstate(h,c)=='listed')
    infos=ppi.get(tuple(sorted([g,h])),[])
    r['PPI_pmids']=','.join(str(x) for x in sorted({i['pmid'] for i in infos if i['pmid']}))
    r['PPI_types']=','.join(sorted({i['type'] for i in infos}))
    r['connectome_marker_ends']=','.join(z for z in (g,h) if z in panel)
    for k in ('magnitude_rank','specificity_rank','ligand_log2FC','ligand_FDR','receptor_log2FC','receptor_FDR'):r[k]=float(r[k])
ranks={(r['source'],r['target'],r['ligand'],r['receptor']):r for r in pairs}
edges=[]
for key in sorted({tuple(sorted([r['ligand'],r['receptor']])) for r in pairs}):edges.extend(ppi.get(key,[]))
for n,rows,fields in [
    ('lr_pair_evidence.tsv',pairs,list(pairs[0])),('all_DEGs_reference_lookup.tsv.gz',degs,list(degs[0])),
    ('SynDive_protein_groups.tsv.gz',synd,list(synd[0])),('SynDive_gene_strata.tsv.gz',synplot,list(synplot[0])),
    ('Synaptome_study_hits.tsv.gz',synrows,list(synrows[0])),('Synaptome_pair_PPI.tsv',edges,['g1','g2','type','method','pmid']),
    ('G2C_memberships.tsv.gz',g2,list(g2[0])),('Connectome_seq_30_genes_by_endpoint.tsv',cm,list(cm[0])),('Connectome_seq_overlap_tests.tsv',stats,list(stats[0]))]:tab(n,rows,fields)
data=dict(pairs=pairs,genes=pairgenes,strata=areas,collections=gcols,markers=panel,status=tested,
          syn=[x for x in synplot if x['gene'] in pairgenes],
          db=[x for x in synrows if x['gene'] in pairgenes],
          g2=[dict(gene=g,collection=c,status=gstate(g,c)) for g in pairgenes for c in gcols],
          con=cm,stats=stats)
(out/'data.json').write_text(json.dumps(data,separators=(',',':'),allow_nan=False))
files=['connect_fig15.jpg','41467_2025_63119_MOESM4_ESM.xlsx','human-mouse-psp.html','tap-psd-95.html','synaptome.db.tar.gz',next(src.glob('*.db.sqlite')).name]
(out/'provenance.json').write_text(json.dumps({'sources_sha256':{n:sha(src/n) for n in files},'strict_input_sha256':sha(prior/'strict_DEG_ligand_receptor_pairs.tsv'),'endpoint_source_sha256':endpoint_files,'excluded':'CA1 C8','IN_C18':'no eligible original mouse-level model','DEG_criterion':'per-endpoint adjusted P < 0.05','reference_panel':'Connectome-seq Extended Data Fig. 9b, top 30, published 2026','synaptome':'AH107282 supplied by synaptome.db 0.99.17','G2C':'G2C programme publication tables Bayes S1/S2d, Fernandez S1'},indent=2))
print(tested,'\n',len(degs),'DEG rows.',len(pairgenes),'pair genes.',len(edges),'PPI assay records.',len(set(r['gene'] for r in cm if r['status']=='DEG')),'panel genes in DEGs.')
print([(s['endpoint'],s['panel_DEG'],round(s['BH_FDR'],4)) for s in stats if s['panel_DEG']]);print([(r['ligand'],r['receptor'],r['FASS_common_strata'],r['PPI_pmids']) for r in pairs])
