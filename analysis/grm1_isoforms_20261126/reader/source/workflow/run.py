import argparse,base64,csv,gzip,hashlib,itertools,json,shlex,subprocess,sys
from pathlib import Path
ROOT=Path(__file__).resolve().parent

def save_json(path,obj):path.write_text(json.dumps(obj,indent=2,allow_nan=False)+'\n')
def tsv(path,rows):
 assert rows
 with path.open('w') as f:
  w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rows)
def ssh_task(cfg,mode,payload):
 code=(ROOT/'core.py').read_text()+'\n'+(ROOT/'remote_counter.py').read_text()
 b=base64.b64encode(json.dumps(payload,separators=(',',':')).encode()).decode()
 code+='\ncfg=json.loads(base64.b64decode('+repr(b)+'))\nprint(json.dumps('+mode+'(cfg),separators=(",",":"),allow_nan=False))\n'
 p=subprocess.run(['/usr/bin/ssh','-o','BatchMode=yes','-o','ConnectTimeout=15',cfg['ssh_alias'],'/usr/bin/python3 -'],input=code,text=True,capture_output=True)
 if p.returncode:raise RuntimeError(p.stderr[-5000:]+'\n'+p.stdout[-1000:])
 return json.loads(p.stdout)
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--config',required=True);args=ap.parse_args()
 cfg=json.loads(Path(args.config).read_text());out=Path(cfg['results_dir']);out.mkdir(parents=True,exist_ok=True)
 assert not (out/'COMPLETE.json').exists(),'Existing completed analysis cannot be overwritten'
 ref=ssh_task(cfg,'read_reference',cfg);save_json(out/'reference.json',ref);ann=ref['annotation']
 (out/'Grm1.reference.gtf').write_text(ref['gtf_excerpt'])
 codon=dict(zip(map(''.join,itertools.product('TCAG',repeat=3)),'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'))
 for t in ann['transcripts']:
  s=t.pop('coding_sequence');phase=(min(t['CDS_phases']) if ann['gene']['strand']=='+' else max(t['CDS_phases']))[1] if t['CDS_phases'] else 0;s=s[phase:];protein=''.join(codon.get(s[i:i+3],'X') for i in range(0,len(s)-2,3));t['protein_length']=len(protein);t['C_terminal_30aa']=protein[-30:];t['protein_sha256']=hashlib.sha256(protein.encode()).hexdigest()
  assert '*' not in protein and 'X' not in protein,'CDS translation is incomplete or out of frame'
  t['protein_family']='noncoding' if not t['CDS'] else 'alpha' if protein.endswith('DYKQSSSTL') and len(protein)==1199 else 'beta' if protein.endswith('KKRQPEFSPSSQCPSAHVQL') and len(protein)==906 else 'unclassified'
  assert t['protein_family']!='unclassified','Unverified protein mapping'
  t['tail_assignment']='Noncoding RNA; no annotated CDS' if t['protein_family']=='noncoding' else 'mGluR1 '+t['protein_family']
 assert {t['name']:t['protein_family'] for t in ann['transcripts']}=={'Grm1-201':'alpha','Grm1-202':'beta','Grm1-203':'beta','Grm1-204':'noncoding','Grm1-205':'noncoding','Grm1-206':'noncoding'}
 assert len({t['protein_sha256'] for t in ann['transcripts'] if t['protein_family']=='beta'})==1
 save_json(out/'annotation.json',ann)
 tsv(out/'transcripts.tsv',[{k:v for k,v in t.items() if k not in ['exons','CDS','UTR','junctions','tags','CDS_phases']}|{'exon_count':len(t['exons']),'CDS_start_incomplete':'cds_start_NF' in t['tags'],'tags':';'.join(t['tags'])} for t in ann['transcripts']])
 features=[dict(f,transcripts=';'.join(f['transcripts'])) for f in ann['features']]
 tsv(out/'features.tsv',features)
 rosters={};genotypes={}
 for r in csv.DictReader(Path(cfg['roster']).open(),delimiter='\t'):
  rosters.setdefault(r['sample'],{})[r['barcode']]={'groups':r['groups'].split(';'),'total_gene_UMIs':int(float(r['total_gene_UMIs']))};genotypes[r['sample']]=r['genotype']
 samples=sorted(cfg['samples'],key=lambda s:(s.startswith('cKO'),int(s.split('_')[1])))
 results=[]
 for sample in samples:
  p=out/(sample+'.json')
  if p.exists():raise RuntimeError('Refusing to replace sample checkpoint '+str(p))
  r=ssh_task(cfg,'count_sample',dict(cfg,sample=sample,sample_paths=cfg['samples'][sample],annotation=ann,roster=rosters[sample],genotype=genotypes[sample]));save_json(p,r);results.append(r)
  print(sample,'completed;',r['diagnostics'].get('accepted_records',0),'accepted Grm1 records',flush=True)
 novel={f[0]:(f[1],f[2]) for r in results for f in r['novel_junctions']}
 for fid,(s,e) in sorted(novel.items()):features.append(dict(id=fid,kind='novel_junction',start=s,end=e,transcripts='',detail='Not present in matching GTF'))
 tsv(out/'features.tsv',features)
 rows=[];coverage=[]
 for r in results:
  counts={(x['group'],x['metric'],x['feature']):x['count'] for x in r['counts']}
  for group,meta in r['groups'].items():
   common=dict(sample=r['sample'],genotype=r['genotype'],population=group,nuclei=meta['nuclei'],total_gene_UMIs=meta['total_gene_UMIs'],Grm1_locus_UMIs=counts.get((group,'gene','gene'),0))
   coverage.append(common)
   # Explicit zeros require observed valid population coverage and completed extraction.
   keys={('gene','gene'),('exonic_gene','gene'),('spliced_gene','gene'),('intronic_evidence','gene'),('records','gene')}
   keys|={(metric,t['name']) for metric in ['compatible_upper_bound','diagnostic_feature','splice_supported_unique'] for t in ann['transcripts']}
   keys|={(metric,fam) for metric in ['family_compatible_upper_bound','family_diagnostic_feature','family_splice_supported'] for fam in ['alpha','beta']}
   keys|={(metric,f['id']) for metric in ['feature_UMIs','feature_records'] for f in features if not f['id'].startswith('novel_junction')}
   keys|={(metric,feature) for (g,metric,feature) in counts if g==group and (metric in ['composition','family_composition','equivalence_class'] or feature.startswith('novel_junction'))}
   for metric,feature in sorted(keys):
    n=counts.get((group,metric,feature),0);den=common['Grm1_locus_UMIs']
    rows.append(dict(common,metric=metric,feature=feature,count=n,per_1000_nuclei=n/meta['nuclei']*1000,per_million_library_UMIs=n/meta['total_gene_UMIs']*1e6,per_1000_Grm1_locus_UMIs=n/den*1000 if den else None))
 tsv(out/'all_counts_long.tsv',rows);tsv(out/'population_mouse_coverage.tsv',coverage)
 tsv(out/'isoform_evidence_per_mouse.tsv',[r for r in rows if r['metric'] in ['compatible_upper_bound','diagnostic_feature','splice_supported_unique']])
 tsv(out/'alpha_beta_evidence_per_mouse.tsv',[r for r in rows if r['metric'].startswith('family_')])
 fdict={f['id']:f for f in features}
 for kind in ['exon','junction','novel_junction','CDS','UTR','exonic_bin']:
  subset=[dict(r,start=fdict[r['feature']]['start'],end=fdict[r['feature']]['end'],transcripts=fdict[r['feature']]['transcripts']) for r in rows if r['metric']=='feature_UMIs' and r['feature'] in fdict and fdict[r['feature']]['kind']==kind]
  if subset:tsv(out/(kind+'_UMI_counts_per_mouse.tsv'),subset)
 tsv(out/'gene_UMI_counts_per_mouse.tsv',[r for r in rows if r['metric'] in ['gene','exonic_gene','spliced_gene','intronic_evidence']])
 summary=[]
 for group in sorted({r['population'] for r in coverage}):
  for t in ann['transcripts']:
   rs=[r for r in rows if r['population']==group and r['metric']=='splice_supported_unique' and r['feature']==t['name']]
   ds=[r for r in rows if r['population']==group and r['metric']=='diagnostic_feature' and r['feature']==t['name']]
   summary.append(dict(population=group,transcript=t['name'],biotype=t['biotype'],tail=t['tail_assignment'],strict_UMIs=sum(r['count'] for r in rs),mice_with_population=len(rs),mice_any_strict_support=sum(r['count']>0 for r in rs),WT_mice_any=sum(r['count']>0 and r['genotype']=='WT' for r in rs),cKO_mice_any=sum(r['count']>0 and r['genotype']=='cKO' for r in rs),mice_at_least_3_strict_UMIs=sum(r['count']>=3 for r in rs),replicated_support=sum(r['count']>=3 for r in rs)>=2,diagnostic_feature_UMIs=sum(r['count'] for r in ds)))
 tsv(out/'isoform_detection_summary.tsv',summary)
 families=[]
 for group in sorted({r['population'] for r in coverage}):
  for family in ['alpha','beta']:
   rs=[r for r in rows if r['population']==group and r['metric']=='family_splice_supported' and r['feature']==family]
   ds=[r for r in rows if r['population']==group and r['metric']=='family_diagnostic_feature' and r['feature']==family]
   families.append(dict(population=group,protein_family=family,strict_UMIs=sum(r['count'] for r in rs),WT_UMIs=sum(r['count'] for r in rs if r['genotype']=='WT'),cKO_UMIs=sum(r['count'] for r in rs if r['genotype']=='cKO'),mice_with_population=len(rs),WT_mice_any=sum(r['count']>0 and r['genotype']=='WT' for r in rs),cKO_mice_any=sum(r['count']>0 and r['genotype']=='cKO' for r in rs),replicated_support=sum(r['count']>=3 for r in rs)>=2,diagnostic_feature_UMIs=sum(r['count'] for r in ds)))
 tsv(out/'alpha_beta_detection_summary.tsv',families)
 original=Path(cfg['original_gene_counts']);orig=list(csv.DictReader(original.open(),delimiter='\t'))
 mapping={}
 benchmark=[]
 for r in coverage:
  gene=next(x['raw_UMI_count'] for x in orig if x['sample']==r['sample'] and x['population']==mapping.get(r['population'],r['population']))
  benchmark.append(dict(r,original_matrix_Grm1_UMIs=int(gene),BAM_minus_matrix=r['Grm1_locus_UMIs']-int(gene)))
 tsv(out/'original_matrix_vs_BAM_gene_counts.tsv',benchmark)
 manifest=dict(status='complete',mice=samples,n_mice=len(samples),n_WT=sum(genotypes[s]=='WT' for s in samples),n_cKO=sum(genotypes[s]=='cKO' for s in samples),roster_sha256=hashlib.sha256(Path(cfg['roster']).read_bytes()).hexdigest(),gtf_sha256=ref['gtf_sha256'],scripts={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in ROOT.glob('*.py')},notes=['Transcript evidence is not complete isoform reconstruction.','Compatible upper bounds overlap and cannot be summed.','Exon and junction feature counts overlap and cannot be summed.','All counts are mouse-level corrected CB/UB support.','cKO_933 uses the recovered original Cell Ranger 10.0 BAM; eleven others use 10.1.'])
 save_json(out/'COMPLETE.json',manifest)
if __name__=='__main__':main()
