#!/usr/bin/env python3
"""Save source-bound human-readable results and full-gene mouse display inputs.
Engineering export only: no RNA extraction, aggregation, normalization fitting,
DE fitting, embedding fitting or alteration of the completed native workflow.
"""
from __future__ import annotations
import argparse,base64,collections,csv,gzip,hashlib,html,json,math,os,re,shutil,subprocess,time
from pathlib import Path
SEVEN=['gene','baseMean','log2FoldChange','lfcSE','stat','pvalue','padj']
POP_NAMES={'allcells':'All nuclei','Excitatory':'Excitatory','Inhibitory':'Inhibitory','CA1_ProS':'CA1-ProS candidates'}
SHA_CACHE={}
def sha(p):
 p=Path(p).resolve();st=p.stat();key=(str(p),st.st_size,st.st_mtime_ns)
 if key not in SHA_CACHE:
  h=hashlib.sha256()
  with p.open('rb') as f:
   for x in iter(lambda:f.read(8*1024*1024),b''):h.update(x)
  SHA_CACHE[key]=h.hexdigest()
 return SHA_CACHE[key]
def read_tsv(p):
 with Path(p).open(encoding='utf-8',newline='') as f:return list(csv.DictReader(f,delimiter='\t'))
def write_tsv(p,rows,fields=None):
 p=Path(p);p.parent.mkdir(parents=True,exist_ok=True)
 if fields is None:fields=list(rows[0]) if rows else []
 with p.open('w',encoding='utf-8',newline='') as f:
  w=csv.DictWriter(f,fieldnames=fields,delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rows)
def dump(p,data):Path(p).write_text(json.dumps(data,ensure_ascii=True,separators=(',',':'),allow_nan=False)+'\n',encoding='utf-8')
def truth(x):return x is True or str(x).lower() in ['true','1']
def scalar(x):
 if x is None or str(x) in ['', 'NA','nan','None']:return None
 v=float(x)
 if not math.isfinite(v):return None
 return v
class Evidence:
 def __init__(self):self.inputs={};self.exports=[]
 def bind(self,p,expected=None):
  p=Path(p).resolve();digest=sha(p)
  if expected and digest!=expected:raise ValueError('Source SHA mismatch: '+str(p))
  self.inputs[str(p)]=digest;return p
 def clone(self,source,dest):
  source=self.bind(source);dest=Path(dest);dest.parent.mkdir(parents=True,exist_ok=True)
  if dest.exists():raise ValueError('Refusing existing export: '+str(dest))
  # APFS clone: independent inode, copy-on-write storage; never shared-inode chmod.
  cp=subprocess.run(['/bin/cp','-c',str(source),str(dest)],capture_output=True)
  mode='APFS_COPY_ON_WRITE_CLONE'
  if cp.returncode:
   if dest.exists():dest.unlink()
   shutil.copy2(source,dest);mode='VERIFIED_BYTE_COPY'
  if sha(dest)!=sha(source) or dest.stat().st_ino==source.stat().st_ino:raise ValueError('Export is not source-isolated and byte-identical')
  self.exports.append({'source_path':str(source),'saved_path':str(dest),'sha256':sha(source),'bytes':source.stat().st_size,'copy_method':mode})
  return dest

def payload(e,base,registry,key):
 meta=registry[key];p=e.bind(Path(base)/meta['path'],meta['sha256']);t=p.read_text(encoding='utf-8')
 entry=json.loads(t.split(',entry=',1)[1].split(';(w.SATB2_SHARDS',1)[0]);raw=gzip.decompress(base64.b64decode(entry['data']))
 if hashlib.sha256(raw).hexdigest()!=meta['payload_sha256']:raise ValueError('Catalog payload hash mismatch')
 return json.loads(raw)
def slug(s):return re.sub('[^A-Za-z0-9._-]+','_',str(s)).strip('_')[:100] or 'unlabeled'
def source_table_counts(e,p):
 p=e.bind(p)
 with p.open(encoding='utf-8',newline='') as f:
  r=csv.reader(f,delimiter='\t')
  if next(r)!=SEVEN:raise ValueError('Nonstandard seven-column results')
  n=0
  for row in r:
   if len(row)!=7:raise ValueError('Malformed result row')
   n+=1
 return n

def write_count_shard(e,counts,expected,genes,samples,bundle):
 counts=e.bind(counts,expected);rows=[];sums=[0]*12
 with gzip.open(counts,'rt',encoding='utf-8',newline='') as f:
  r=csv.reader(f,delimiter='\t');hdr=next(r)
  if hdr!=['gene_id','gene_symbol']+samples:raise ValueError('Full gene source sample order differs')
  n=0
  for n,row in enumerate(r,1):
   i=n-1
   if i>=len(genes) or row[0]!=genes[i]['gene_id'] or len(row)!=14:raise ValueError('Count gene axis mismatch')
   symbol=row[1] if row[1] not in ['', 'NA'] else None
   if symbol!=genes[i]['gene_symbol']:raise ValueError('Count gene-symbol axis mismatch')
   a=[int(v) for v in row[2:]]
   if any(v<0 or v>2**53-1 or str(v)!=tok for v,tok in zip(a,row[2:])):raise ValueError('Source counts are not exact nonnegative safe integers')
   for j,v in enumerate(a):sums[j]+=v
   if any(a):rows.append([i,a])
  if n!=len(genes):raise ValueError('Source full gene universe incomplete')
 data={'gene_axis':'global','n_genes':len(genes),'sample_order':samples,'rows':rows,'omitted_rows':'EXPLICIT_ALL_ZERO','source_counts_sha256':expected,'all12_integer_column_totals':sums}
 raw=json.dumps(data,separators=(',',':'),allow_nan=False).encode();packed=gzip.compress(raw,compresslevel=6,mtime=0);ident='mouse_counts__'+expected
 entry={'encoding':'gzip-base64','data':base64.b64encode(packed).decode('ascii')}
 script=';(function(w){var id='+json.dumps(ident)+',entry='+json.dumps(entry,separators=(',',':'))+';(w.SATB2_SHARDS=w.SATB2_SHARDS||{})[id]=entry;if(typeof w.registerShard===\'function\'){w.registerShard(id,entry.data);}})(window);\n'
 p=bundle/'data/mouse_shards'/(ident+'.js');p.parent.mkdir(parents=True,exist_ok=True)
 if p.exists():raise ValueError('Refusing existing mouse shard')
 p.write_text(script,encoding='utf-8')
 return ident,{'path':str(p.relative_to(bundle)),'sha256':sha(p),'encoding':'gzip-base64','payload_sha256':hashlib.sha256(raw).hexdigest(),'source_counts_sha256':expected,'source_counts_path':str(counts),'uncompressed_bytes':len(raw),'gzip_bytes':len(packed),'n_nonzero_genes':len(rows),'n_genes':len(genes),'all12_integer_column_totals':sums}

def build(args):
 cfg_path=Path(args.config).resolve();cfg=json.load(cfg_path.open());study=Path(cfg['study_root']);native=Path(cfg['results_root']);saved=Path(args.outdir).resolve();bundle=Path(args.mouse_bundle_root).resolve();old=study/'reports/offline_dashboard_v1';e=Evidence();e.bind(cfg_path)
 if saved.exists() and any(saved.iterdir()):raise ValueError('Refusing nonempty saved export directory')
 if (bundle/'data/mouse_catalog.js').exists() or (bundle/'MOUSE_CATALOG_COMPLETE.json').exists():raise ValueError('Refusing existing mouse catalog')
 saved.mkdir(parents=True,exist_ok=True);bundle.mkdir(parents=True,exist_ok=True)
 cat_path=e.bind(old/'DE_CATALOG.json');cat=json.load(cat_path.open());cat_receipt=e.bind(old/'DE_CATALOG_COMPLETE.json');cr=json.load(cat_receipt.open())
 if cat['status']!='COMPLETE' or cr['status']!='COMPLETE' or cat['config_sha256']!=sha(cfg_path):raise ValueError('Final exact DE catalog is not closed for current config')
 for branch in ['original_full','paired_original','paired_cellbender']:
  cp=e.bind(native/'deseq2'/branch/'COMPLETE.json');rp=json.load(cp.open())
  if rp['status']!='EXECUTED_REVIEW_REQUIRED':raise ValueError('Native DE producer incomplete')
 acceptance=json.load(e.bind(native/'count_validation/COMPLETE.json').open());comparison=json.load(e.bind(native/'comparison/COMPLETE.json').open())
 if not acceptance['exact_integer_counts_PASS'] or comparison['status']!='EXECUTED_REVIEW_REQUIRED':raise ValueError('Native count acceptance/comparison is not complete')
 gene_path=e.bind(native/'populations/original_full/GENE_ID_SYMBOL_MANIFEST.tsv');genes=read_tsv(gene_path)
 if len(genes)!=32285 or len({x['gene_id'] for x in genes})!=len(genes):raise ValueError('Unexpected original fullgene axis')
 for x in genes:x['gene_symbol']=x['gene_symbol'] if x['gene_symbol'] not in ['', 'NA'] else None
 fits=[]
 for ds in cat['datasets']:
  for f in payload(e,old,cat['shards'],ds['fits_shard'])['fits']:
   # Exact prevalidated fit ID and native metadata; coordinate lists unnecessary here.
   f.pop('members',None);fits.append(f)
 if len(fits)!=1060 or len({f['id'] for f in fits})!=1060:raise ValueError('Incomplete exact fit registry')
 first=next(x for x in fits if x['analysis_kind']=='total' and x['branch']=='original_full' and x['endpoint']=='PRIMARY');samples=[x['sample'] for x in first['coverage']];conditions={x['sample']:x['condition'] for x in first['coverage']}
 if len(samples)!=12 or len(set(samples))!=12:raise ValueError('Need all12 independent mouse axis')
 idx={}
 for branch in ['original_full','paired_original','paired_cellbender']:
  ip=e.bind(native/'populations'/branch/'POPULATION_INDEX.tsv');idx[branch]={x['population_key']:x for x in read_tsv(ip)}
 result_rows=[];mouse_fits={};mouse_shards={};counts_registry=[];count_to_shard={};allcoverage=[];source_count_exports={}
 for pos,f in enumerate(fits):
  if pos%50==0:print(f'Exporting endpoints {pos}/{len(fits)}; distinct count arrays {len(count_to_shard)}',flush=True)
  row=idx[f['branch']][f['population_key']];countsha=row['counts_sha256']
  if countsha!=f['counts_source_sha256']:raise ValueError('Fit/count registry binding differs')
  counts=e.bind(row['counts_file'],countsha)
  if countsha not in count_to_shard:
   ident,meta=write_count_shard(e,counts,countsha,genes,samples,bundle);count_to_shard[countsha]=ident;mouse_shards[ident]=meta
   saved_count=e.clone(counts,saved/'pseudobulk_source_counts'/(countsha+'.tsv.gz'));source_count_exports[countsha]=saved_count
   counts_registry.append({'counts_shard':ident,'counts_sha256':countsha,'n_genes':len(genes),'n_nonzero_genes':meta['n_nonzero_genes'],'source_path':str(counts),'saved_path':str(saved_count.relative_to(saved))})
  cov={x['sample']:x for x in f['coverage']}
  if set(cov)!=set(samples) or any(cov[s]['condition']!=conditions[s] for s in samples):raise ValueError('Coverage full mouse identity/genotype differs')
  audit=e.bind(f['native_source_paths']['endpoint_audit']);ep=audit.parent.parent;actual=read_tsv(audit)
  if len(actual)!=1 or actual[0]['status']!=f['status'] or actual[0]['n_WT']!=str(f['n_WT']):raise ValueError('Endpoint audit/catalog mismatch')
  sf_path=ep/'AUDIT/FITTED_SIZE_FACTORS.tsv';sf={};sf_sha=None
  eligible=[truth(cov[s]['included_in_model']) for s in samples]
  if sf_path.exists():
   e.bind(sf_path);sf_sha=sha(sf_path)
   for sr in read_tsv(sf_path):
    val=scalar(sr['fitted_size_factor'])
    if sr['sample'] in sf or val is None or val<=0:raise ValueError('Invalid recorded size factors')
    sf[sr['sample']]=val
   if set(sf)!={s for s,b in zip(samples,eligible) if b}:raise ValueError('Recorded SF does not exactly match model mouse roster')
  elif f['status']=='TESTED_REVIEW_REQUIRED':raise ValueError('Tested model lacks recorded fitted SF')
  mouse_fits[f['id']]={'dataset_id':f['dataset_id'],'tool':f['tool'],'branch':f['branch'],'population':f['population'],'population_key':f['population_key'],'endpoint':f['endpoint'],'analysis_kind':f['analysis_kind'],'category_family':f['category_family'],'analysis_label':f['analysis_label'],'status':f['status'],'reason':f['reason'],'counts_shard':count_to_shard[countsha],'sample_order':samples,'size_factors':[sf.get(s) for s in samples],'eligible':eligible,'n_nuclei':[int(cov[s]['n_nuclei']) for s in samples],'counts_sha256':countsha,'membership_sha256':f['membership_sha256'],'n_WT':f['n_WT'],'n_cKO':f['n_cKO'],'residual_df':f['residual_df'],'coverage_flag':f['coverage_flag'],'warning_count':f['warning_count'],'source_endpoint_audit_path':str(audit),'source_endpoint_audit_sha256':sha(audit),'source_size_factors_path':str(sf_path) if sf_sha else None,'source_size_factors_sha256':sf_sha,'normalization_display':'source integer count / source fitted size factor; no new fitting','excluded_normalized_value':None}
  segment='primary' if f['endpoint']=='PRIMARY' else 'sensitivity'
  dest=saved/'DESeq2'/segment/f['branch']/f['tool']/f['population']/f['category_family']/(f['population_key']+'__'+slug(f['analysis_label']))/f['endpoint']
  dest.mkdir(parents=True,exist_ok=True)
  base={k:f[k] for k in ['id','dataset_id','tool','branch','population','population_key','endpoint','analysis_kind','category_family','analysis_label','status','reason','n_nuclei','n_WT','n_cKO','residual_df','coverage_flag','warning_count','n_endpoint_genes','n_pvalue','n_padj','DEG_FDR05','membership_sha256','execution_origin']};base['population_display_name']=POP_NAMES[f['population']];base['physical_layer_validated']=False;base['counts_sha256']=countsha;base['source_output_directory']=str(ep);base['saved_directory']=str(dest.relative_to(saved));base['all_results']='';base['all_results_sha256']='';base['deg_fdr05']='';base['deg_fdr05_sha256']='';base['stable_id_audit']='';base['stable_id_audit_sha256']=''
  if f['status']=='TESTED_REVIEW_REQUIRED':
   n=source_table_counts(e,ep/'ALL_RESULTS.txt');nd=source_table_counts(e,ep/'DEG_FDR05.txt')
   if n!=f['n_endpoint_genes'] or nd!=f['DEG_FDR05']:raise ValueError('Saved result row count disagrees with native endpoint audit')
   for src,name,pkey in [(ep/'ALL_RESULTS.txt','ALL_RESULTS_7_COLUMNS.tsv','all_results'),(ep/'DEG_FDR05.txt','DEG_FDR05_7_COLUMNS.tsv','deg_fdr05'),(ep/'AUDIT/ALL_RESULTS_STABLE_ID_AUDIT.tsv','ALL_RESULTS_STABLE_ID_AUDIT.tsv','stable_id_audit')]:
    out=e.clone(src,dest/name);base[pkey]=str(out.relative_to(saved));base[pkey+'_sha256']=sha(src)
  for name in ['EXECUTION_AUDIT.tsv','MICE_ALL12_AUDIT.tsv','MODEL_SETTINGS.tsv','FITTED_SIZE_FACTORS.tsv']:
   src=ep/'AUDIT'/name
   if src.exists():e.clone(src,dest/'AUDIT'/name)
  base['coverage_audit']=str((dest/'AUDIT/MICE_ALL12_AUDIT.tsv').relative_to(saved));base['execution_audit']=str((dest/'AUDIT/EXECUTION_AUDIT.tsv').relative_to(saved));result_rows.append(base)
  dump(dest/'STATUS_AND_SCOPE.json',{**base,'contrast':'cKO / WT','BH_scope':'within this endpoint only','all12_coverage':[cov[s] for s in samples],'model_settings':f['model_settings'],'source_count_file':str(counts),'source_recorded_size_factors_sha256':sf_sha,'expression_programs_are_anatomical_layers':False})
  for s in samples:allcoverage.append({'fit_id':f['id'],'tool':f['tool'],'branch':f['branch'],'population':f['population'],'population_key':f['population_key'],'endpoint':f['endpoint'],'sample':s,'condition':conditions[s],'n_nuclei':cov[s]['n_nuclei'],'included_in_model':eligible[samples.index(s)],'recorded_size_factor':sf.get(s),'inclusion_or_exclusion_reason':cov[s]['inclusion_or_exclusion_reason']})
 primary=[r for r in result_rows if r['endpoint']=='PRIMARY'];sens=[r for r in result_rows if r['endpoint']!='PRIMARY']
 write_tsv(saved/'PRIMARY_RESULTS_INDEX.tsv',primary);write_tsv(saved/'SENSITIVITY_RESULTS_INDEX.tsv',sens);write_tsv(saved/'ALL12_MOUSE_COVERAGE_INDEX.tsv',allcoverage);write_tsv(saved/'PSEUDOBULK_SOURCE_COUNTS_INDEX.tsv',counts_registry)
 # Existing plot books and completed diagnostics, copied by source-isolated clones.
 figrows=[];asset_path=e.bind(old/'data/asset_map.json');assets=json.load(asset_path.open());dedup={}
 for a in assets:
  src=Path(a['source_path'])
  if src.suffix.lower()!='.pdf':continue
  src=e.bind(src,a['sha256']);rel=Path(a['bundle_path'])
  target=saved/'figures'/rel.relative_to('assets/pdf');target=e.clone(src,target)
  figrows.append({'title':src.stem.replace('_',' '),'scope':('original frozen main inference UMAP and marker book' if '/main/' in a['bundle_path'] else 'secondary representation or diagnostic figure; inspect caption'),'saved_pdf':str(target.relative_to(saved)),'source_pdf':str(src),'sha256':a['sha256']})
 write_tsv(saved/'FIGURE_INDEX.tsv',figrows)
 # Full gene lookup and portable mouse catalog, with exact fit IDs compatible with DE selector.
 catalog={'schema_version':'satb2-mouse-offline-1','status':'COMPLETE','genes':genes,'samples':[{'sample':s,'condition':conditions[s]} for s in samples],'shards':mouse_shards,'fits':mouse_fits,'default_contrast':{'numerator':'cKO','denominator':'WT'},'config_sha256':sha(cfg_path),'de_catalog_sha256':sha(cat_path),'gene_axis_source_sha256':sha(gene_path),'n_genes':len(genes),'scope':'Existing all-gene integer mouse aggregates and recorded endpoint size factors only. No nucleus extraction, new aggregation or fitting. Unsupported/excluded normalized values remain null. All-zero source rows are explicit zeros.','statistical_inference_for_all_genes':False,'inference_scope':'Only a produced DE endpoint result supports inference; full-gene display availability does not add a fitted result.'}
 dump(bundle/'data/MOUSE_CATALOG.json',catalog);(bundle/'data/mouse_catalog.js').write_text('window.SATB2_MOUSE_CATALOG='+json.dumps(catalog,ensure_ascii=True,separators=(',',':'),allow_nan=False)+';\n',encoding='utf-8')
 write_tsv(saved/'GENE_ID_SYMBOL_MANIFEST.tsv',[{'gene_id':g['gene_id'],'gene_symbol':g['gene_symbol'] or 'NA'} for g in genes])
 write_tsv(saved/'SOURCE_FILE_EXPORT_INDEX.tsv',e.exports)
 build_html(saved,primary,sens,figrows)
 (saved/'README.txt').write_text('DESeq2 cKO/WT saved results\n\nOpen index.html. PRIMARY_RESULTS_INDEX.tsv lists all948 primary endpoints, including unsupported clusters. Full and BH FDR<0.05 files retain the seven standard source columns and gene symbols; stable IDs are separate audits. No extra fold-change filter. All tables/PDFs are unchanged, verified source-isolated copies or APFS clones.\n\nCA1-ProS candidates and layer-associated expression programs are the actual available definitions. They are not validated pure CA1 or physical sublayers. Fine labels remain provisional. BH is within each endpoint and overlapping categories are not independent discoveries. IN uses all available nuclei; mice contributing fewer than5 nuclei are exploratory. All12 mouse coverage and unsupported reasons are retained.\n\nPseudobulk mouse display inputs include32285 full source genes per array, all12 mice, and recorded fit size factors only. Missing normalized values stay unavailable. Full-gene display does not imply the gene was tested.\n',encoding='utf-8')
 outputs={str(p.relative_to(saved)):sha(p) for p in saved.rglob('*') if p.is_file()};mouse_outputs={str(p.relative_to(bundle)):sha(p) for p in [bundle/'data/MOUSE_CATALOG.json',bundle/'data/mouse_catalog.js']};mouse_outputs.update({m['path']:m['sha256'] for m in mouse_shards.values()})
 receipt={'status':'COMPLETE','stage':'saved_primary_results_and_mouse_display_inputs','config_sha256':sha(cfg_path),'source_code_sha256':sha(__file__),'input_sha256':e.inputs,'output_sha256':outputs,'n_primary_endpoints':len(primary),'n_tested_primary':sum(r['status']=='TESTED_REVIEW_REQUIRED' for r in primary),'n_unsupported_primary':sum(r['status']!='TESTED_REVIEW_REQUIRED' for r in primary),'n_sensitivity_endpoints':len(sens),'n_figure_pdfs':len(figrows),'n_mouse_endpoints':len(mouse_fits),'n_distinct_fullgene_count_arrays':len(mouse_shards),'n_genes':len(genes),'source_isolated_copy_methods':dict(collections.Counter(r['copy_method'] for r in e.exports)),'scientific_fitting_performed':False,'raw_nucleus_count_aggregation_performed':False,'native_outputs_modified':False,'pure_CA1_population_available':False,'physical_sublayer_assignment_validated':False,'mouse_catalog_root':str(bundle),'mouse_catalog_output_sha256':mouse_outputs}
 dump(saved/'COMPLETE.json',receipt);dump(bundle/'MOUSE_CATALOG_COMPLETE.json',{k:v for k,v in receipt.items() if k!='output_sha256'}|{'output_sha256':mouse_outputs})
 print(json.dumps({k:v for k,v in receipt.items() if k.startswith('n_') or k=='source_isolated_copy_methods'},indent=2),flush=True)

def build_html(saved,primary,sens,figures):
 def esc(x):return html.escape(str(x))
 def link(path,label):return '<a href="'+esc(path)+'">'+esc(label)+'</a>' if path else '<span class="na">not fitted</span>'
 tab=[]
 for r in primary:
  key=' '.join(str(r[k]) for k in ['tool','branch','population','analysis_kind','category_family','analysis_label','status','coverage_flag'])
  tab.append('<tr data-search="'+esc(key.lower())+'"><td>'+esc(r['tool'])+'</td><td>'+esc(r['branch'])+'</td><td>'+esc(r['population_display_name'])+'</td><td>'+esc(r['analysis_kind'])+'</td><td>'+esc(r['analysis_label'])+'</td><td>'+esc(r['status'].replace('_REVIEW_REQUIRED',''))+'</td><td>'+esc(r['n_nuclei'])+'</td><td>'+esc(str(r['n_WT'])+'/'+str(r['n_cKO']))+'</td><td>'+esc(r['DEG_FDR05'])+'</td><td>'+link(r['deg_fdr05'],'DEG FDR05')+' · '+link(r['all_results'],'Full 7 columns')+'</td><td>'+link(r['coverage_audit'],'12-mouse coverage')+'<br>'+esc(r['coverage_flag'])+' '+esc(r['reason'])+'</td></tr>')
 figs=''.join('<li>'+link(r['saved_pdf'],r['title'])+' <small>'+esc(r['scope'])+'</small></li>' for r in figures)
 text='''<!doctype html><html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>Saved DESeq2 results and figures</title><style>body{font:16px system-ui,sans-serif;margin:28px;color:#182c40;background:#fbfcfd}main{max-width:1700px;margin:auto}h1{margin-bottom:8px}p{max-width:1000px;line-height:1.5}a{color:#005fa8}input{padding:12px;width:min(850px,90%);border:1px solid #aebdcb;border-radius:5px;font-size:16px}table{border-collapse:collapse;width:100%;font-size:13px;background:white;margin-top:18px}th,td{text-align:left;border-bottom:1px solid #dde5ec;padding:9px;vertical-align:top}th{position:sticky;top:0;background:#e5eef6}.note{background:#edf5fa;border-left:4px solid #337ba5;padding:12px}.na{color:#66788a}small{color:#66788a}li{margin:12px 0}details{margin-top:28px}.tablebox{overflow:auto}code{overflow-wrap:anywhere}</style></head><body><main><h1>Saved DESeq2 results and figures</h1><p>cKO versus WT · all12 independent mice · Scanpy and Seurat original Harmony0.4 populations · three count branches.</p><p class="note">CA1/ProS candidates are the actual available CA1-related population. Deep- and superficial-associated expression programs are not validated physical sublayers. Gene symbols are retained in the seven-column downloads; stable gene IDs and full mouse coverage remain in the audits.</p><p>Every PRIMARY endpoint is listed, including unsupported clusters. <strong>FDR05</strong> means BH adjusted P&lt;0.05 within that endpoint, with no added fold-change cutoff. Missing P values/FDR are unavailable. Low-nucleus inhibitory models are exploratory; overlapping cluster/category endpoints are not independent discoveries.</p><p><a href="PRIMARY_RESULTS_INDEX.tsv">PRIMARY index</a> · <a href="SENSITIVITY_RESULTS_INDEX.tsv">Partition/leave-one-mouse-out index</a> · <a href="ALL12_MOUSE_COVERAGE_INDEX.tsv">All12 mouse coverage</a> · <a href="FIGURE_INDEX.tsv">Figure index</a> · <a href="PSEUDOBULK_SOURCE_COUNTS_INDEX.tsv">Full-gene source mouse count index</a> · <a href="README.txt">Scope and file notes</a></p><h2>Per-total, per-cluster and reviewed category DEG lists</h2><p><input id="search" placeholder="Filter: CA1, superficial, Inhibitory, total, cluster, paired_cellbender…"><span id="shown"></span></p><div class="tablebox"><table><thead><tr><th>Tool</th><th>Count branch</th><th>Population</th><th>Definition</th><th>Cluster / category</th><th>Model</th><th>Nuclei</th><th>WT/cKO mice</th><th>BH05 genes</th><th>Saved gene-symbol tables</th><th>Coverage / caveat</th></tr></thead><tbody>'''+''.join(tab)+'''</tbody></table></div><details open><summary><strong>Saved UMAP, representative/given-marker and effect figure PDFs</strong></summary><p>The original-frozen books show the actual inference memberships on existing coordinates. Secondary fresh-branch books are descriptive representation sensitivity. Marker pages use existing native full-cohort summaries, and retain the supplied/canonical panel scope.</p><ul>'''+figs+'''</ul></details><p><small>All files are verified unchanged source-isolated exports. No scientific fitting or nucleus-level aggregation was performed to build this page. COMPLETE.json and SOURCE_FILE_EXPORT_INDEX.tsv record the source and export SHA256 hashes.</small></p></main><script>const q=document.querySelector('#search'),rs=[...document.querySelectorAll('tbody tr')];function filter(){let n=0;const s=q.value.trim().toLowerCase();for(const r of rs){r.hidden=!r.dataset.search.includes(s);if(!r.hidden)n++}document.querySelector('#shown').textContent=' '+n+' / '+rs.length+' endpoints'}q.addEventListener('input',filter);filter();</script></body></html>'''
 (saved/'index.html').write_text(text,encoding='utf-8')

if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--config',required=True);ap.add_argument('--outdir',required=True);ap.add_argument('--mouse-bundle-root',required=True);a=ap.parse_args();build(a)
