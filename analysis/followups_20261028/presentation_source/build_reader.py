#!/usr/bin/env python3
"""Package a reviewed, portable reader of native results; never infer or rescore.

Only final seals, unchanged native assets and explicitly approved interpretation
are used. The JSON embedded in HTML is a navigation index, not a replacement
scientific result. No run directory, GitHub repository or release is modified.
"""
from __future__ import annotations
import argparse, csv, gzip, hashlib, json, shutil, tempfile
from pathlib import Path
from datetime import datetime, timezone
from urllib.parse import quote, urlparse

HERE = Path(__file__).resolve().parent
SEALS = {
 'composition': ('satb2-composition-final-validation/v1', {'COMPLETE','COMPLETE_WITH_INCONCLUSIVE_OR_UNTESTABLE_SCOPES'}, 'outputs'),
 'cellchat': ('satb2-cellchat-final-validation/v1', {'COMPLETE'}, 'files'),
 'monocle3': ('satb2-monocle3-final-validation/v1', {'COMPLETE'}, 'outputs'),
}
TITLES={'composition':'1 · Cell-type proportions','cellchat':'2 · CellChat communication','monocle3':'3 · Monocle3 state graphs','velocity':'4 · RNA velocity'}
RAW_AXES={'nucleus_metadata.tsv','input_cells.tsv','input_genes.tsv','feature_mapping.tsv','frozen_candidates_before_after_keep.tsv.gz'}

def sha(path):
 h=hashlib.sha256()
 with Path(path).open('rb') as stream:
  for chunk in iter(lambda:stream.read(1024*1024),b''):h.update(chunk)
 return h.hexdigest()

def receipt(path):
 p=Path(path).resolve(strict=True);return {'path':str(p),'bytes':p.stat().st_size,'sha256':sha(p)}

def verify(record, root=None):
 p=Path(record['path'])
 if not p.is_absolute():
  if root is None:raise ValueError('Relative receipt requires native root')
  p=root/p
 p=p.resolve(strict=True)
 if root is not None and not p.is_relative_to(root.resolve()):raise ValueError(f'Receipt escapes native root: {p}')
 if not p.is_file() or p.stat().st_size!=record.get('bytes',record.get('size_bytes')) or sha(p)!=record['sha256']:
  raise ValueError(f'Native/source receipt mismatch: {p}')
 return p

def native_seal(kind, spec):
 root=Path(spec['native_root']).resolve(strict=True);seal_path=verify(spec['terminal_seal'])
 if not seal_path.is_relative_to(root):raise ValueError('Terminal seal outside its native root')
 seal=json.loads(seal_path.read_text());schema,states,key=SEALS[kind]
 if seal.get('schema')!=schema or seal.get('status') not in states:raise ValueError(f'{kind} native task is not in an accepted terminal state')
 # Fail closed on synthetic fixture roots or declared synthetic analysis IDs.
 if 'synthetic' in str(root).lower() or 'synthetic' in str(seal.get('analysis_id','')).lower():
  raise ValueError('Synthetic results cannot be published as study results')
 records=seal.get(key,[])
 if not records:raise ValueError(f'{kind} seal has no native output receipts')
 paths={}
 for rec in records:
  p=verify(rec,root);rel=p.relative_to(root)
  if str(rel) in paths:raise ValueError('Duplicate native receipt')
  paths[str(rel)]=p
 return root,seal_path,seal,paths

def descriptor(kind, rel, source, root, native_paths=None):
 parts=rel.parts;d={'task':kind,'scope':'All native scopes','cohort':'All saved cohorts','method':'Native','direction':'All directions','plot':'Source file'}
 if kind=='composition' and len(parts)>3 and parts[0]=='fits':
  d.update(scope=parts[1],cohort=parts[2],method='scCODA')
 elif kind=='composition':d.update(method='scCODA',scope='Cross-scope source tables')
 elif kind=='cellchat':
  d['method']=parts[1] if parts[0]=='figures' and len(parts)>2 else parts[0] if parts[0] in ('triMean','truncatedMean10') else 'Shared native sources'
  name=rel.name
  for direction in ('CA1_to_IN','IN_to_CA1','bidirectional'):
   if name.startswith(direction):d['direction']=direction
  if '_WT_' in name:d['cohort']='WT'
  elif '_cKO_' in name:d['cohort']='cKO'
  else:d['cohort']='Both genotypes / saved inference scopes'
  # Read only a manifest already verified as a sealed native output.  The mouse
  # boxplot filename has no direction prefix, but its native view is bidirectional.
  manifest_relative=(rel.parent/'PLOT_MANIFEST.tsv').as_posix()
  if source.suffix.lower() in ('.png','.pdf') and native_paths is not None and manifest_relative in native_paths:
   with native_paths[manifest_relative].open(newline='') as stream:
    rows=[row for row in csv.DictReader(stream,delimiter='\t') if row.get('base')==source.stem]
   if len(rows)>1:raise ValueError('Duplicate native plot identity in sealed plot manifest')
   if rows:
    d['direction']=rows[0]['direction']
    d['cohort']=rows[0]['condition'] if rows[0]['condition'] in ('WT','cKO') else 'Both genotypes / saved inference scopes'
  d['scope']='CA1 programs ↔ IN C1/C6/C8'
 elif kind=='monocle3':
  d.update(method='Monocle3',cohort='All 12 mice; native geometry')
  if 'models' in parts:
   i=parts.index('models');d['scope']=parts[i+1]
  elif len(parts)>1 and parts[0]=='monocle3':d['scope']=parts[1]
 name=source.stem
 d['plot']=name.replace('_',' ')
 d['relative_native_path']=rel.as_posix()
 d['status']='Native file'
 # Surface original per-model statuses without recomputing calls or diagnostics.
 if kind=='composition' and len(parts)>3 and parts[0]=='fits':
  status=root/'fits'/parts[1]/parts[2]/'RESULT_STATUS.json'
  if status.exists():d['status']=json.loads(status.read_text()).get('status','Unavailable')
 elif kind=='monocle3' and 'models' in parts:
  directory=root.joinpath(*parts[:parts.index('models')+2])
  status=directory/'SCOPE_VALIDATION.json'
  if status.exists():d['status']=json.loads(status.read_text()).get('status','Unavailable')
 return d

def choose_asset(kind, rel, path, max_table):
 if any(part.startswith('chain_') for part in rel.parts):return False,'Native chain objects/diagnostics retained in full archive'
 if path.suffix.lower() in ('.png','.pdf'):return True,None
 if path.name in RAW_AXES:return False,'Full nucleus/feature axes retained in native archive'
 if path.suffix.lower() in ('.tsv','.csv','.json','.md','.txt') and path.stat().st_size<=max_table:return True,None
 return False,'Large native models, raw axes or large tables retained in full archive'

def table_preview(path, budget):
 if path.suffix.lower() not in ('.tsv','.csv') or budget<=0:return None
 delimiter='\t' if path.suffix.lower()=='.tsv' else ','
 with path.open(newline='') as stream:
  reader=csv.reader(stream,delimiter=delimiter);rows=[];size=0
  for row in reader:
   size+=sum(len(x) for x in row)
   if size>budget or len(rows)>=61:break
   rows.append(row)
 return rows

def layout_override(spec, native_root, scientific_seal, native_paths):
 """Verify a separately sealed display-only repair. Never replace native data."""
 override=spec.get('presentation_override')
 if not override:return None
 if override.get('approved_for_display') is not True or not override.get('visual_review_receipts'):
  raise ValueError('Layout override requires explicit display approval and independent visual-review receipt')
 for record in override['visual_review_receipts']:verify(record)
 root=Path(override['presentation_root']).resolve(strict=True)
 if root==native_root or root.is_relative_to(native_root) or native_root.is_relative_to(root):
  raise ValueError('Presentation repair must remain separate from original scientific results')
 seal_path=verify(override['terminal_seal']);seal=json.loads(seal_path.read_text())
 if not seal_path.is_relative_to(root) or seal.get('schema')!='satb2-cellchat-layout-only-final/v1' or seal.get('status') not in {'COMPLETE','COMPLETE_PENDING_VISUAL_REVIEW'}:
  raise ValueError('Presentation repair lacks an accepted actual terminal seal')
 for flag in ('scientific_fit_or_rescore','scientific_aggregation_formulas_changed','colors_and_numerical_scales_changed','original_native_outputs_modified','new_statistical_results'):
  if seal.get(flag) is not False:raise ValueError(f'Layout-only preservation flag absent/violated: {flag}')
 original=seal.get('native_final_validation',{})
 if original.get('sha256')!=sha(scientific_seal) or verify(original)!=scientific_seal:
  raise ValueError('Layout repair belongs to a different original scientific seal')
 outputs={}
 for record in seal.get('outputs',[]):
  path=verify(record,root)
  if str(path) in outputs:raise ValueError('Duplicate repaired output receipt')
  outputs[str(path)]=path
 mapping_path=root/'NATIVE_TO_PRESENTATION_PLOT_MAPPING.json'
 if str(mapping_path) not in outputs:raise ValueError('Native-to-layout mapping is not sealed')
 mapping=json.loads(mapping_path.read_text())
 if not isinstance(mapping,list) or not mapping:raise ValueError('Empty or invalid layout mapping')
 checked=[];seen_original=set();seen_corrected=set()
 for item in mapping:
  original_path=verify(item['native_original'],native_root);corrected=verify(item['corrected_presentation'],root)
  relative=original_path.relative_to(native_root).as_posix()
  if relative not in native_paths or corrected.as_posix() not in outputs:
   raise ValueError('Mapped original/corrected plot is absent from its sealed outputs')
  if corrected.relative_to(root).as_posix()!=relative or original_path.suffix not in ('.pdf','.png') or corrected.suffix!=original_path.suffix:
   raise ValueError('Repaired figure identity/relative path changed')
  if relative in seen_original or str(corrected) in seen_corrected:raise ValueError('Duplicate native/corrected plot mapping')
  seen_original.add(relative);seen_corrected.add(str(corrected));checked.append((relative,original_path,corrected))
 expected={relative for relative,path in native_paths.items() if path.suffix.lower() in ('.png','.pdf')}
 if seen_original!=expected:raise ValueError('Layout mapping does not cover every sealed original PNG/PDF')
 return {'root':root,'seal_path':seal_path,'mapping_path':mapping_path,'checked':checked,'reviews':override['visual_review_receipts']}

def approved_configuration(path):
 cfg=json.loads(Path(path).read_text())
 if cfg.get('schema')!='satb2-followups-portable-presentation/v1' or cfg.get('approved_for_portable_presentation') is not True:
  raise ValueError('Root-approved reviewed interpretation configuration is required')
 if set(cfg.get('tasks',{}))!=set(TITLES):raise ValueError('All four task sections are mandatory')
 for kind,spec in cfg['tasks'].items():
  if not spec.get('findings') or not spec.get('limitations') or not spec.get('review_receipts'):
   raise ValueError(f'{kind} requires approved findings, limitations and reviews')
  for rec in spec['review_receipts']:verify(rec)
  if spec.get('archive_url'):
   parsed=urlparse(spec['archive_url'])
   if parsed.scheme!='https' or not parsed.netloc or parsed.username or parsed.password:
    raise ValueError('Archive download links must be ordinary HTTPS URLs')
 if cfg['tasks']['velocity'].get('status')!='INPUT_BLOCKED':raise ValueError('Velocity must remain INPUT_BLOCKED under this assessment')
 return cfg

def protected_result_roots(cfg):
 """Protect scientific outputs and their conventional complete run trees."""
 roots=set()
 for spec in cfg.get('tasks',{}).values():
  candidates=[spec.get('native_root'),spec.get('presentation_override',{}).get('presentation_root')]
  for value in candidates:
   if not value:continue
   root=Path(value).resolve();roots.add(root)
   for parent in (root,*root.parents):
    if parent.name=='results':
     roots.add(parent.parent);break
 return roots

def reject_protected_destination(destination,roots):
 destination=Path(destination).resolve()
 if any(destination.is_relative_to(root) for root in roots):
  raise ValueError('Presentation destination is inside a protected native/layout run tree')

def lossless_gzip(source,target):
 """Deterministic download transport of original bytes; never reformat a table."""
 target.parent.mkdir(parents=True,exist_ok=True)
 with Path(source).open('rb') as raw,target.open('wb') as encoded:
  with gzip.GzipFile(filename='',mode='wb',fileobj=encoded,compresslevel=9,mtime=0) as transport:
   shutil.copyfileobj(raw,transport,1024*1024)
 decoded=hashlib.sha256()
 with gzip.open(target,'rb') as stream:
  for chunk in iter(lambda:stream.read(1024*1024),b''):decoded.update(chunk)
 if decoded.hexdigest()!=sha(source):raise ValueError('Lossless native table transport verification failed')
 return target

def build(config_path,destination):
 cfg=approved_configuration(config_path);protected=protected_result_roots(cfg)
 reject_protected_destination(destination,protected)
 reject_protected_destination(tempfile.gettempdir(),protected)
 # Compression staging is outside every run/native tree and always removed.
 with tempfile.TemporaryDirectory(prefix='satb2_reader_transport_') as scratch:
  return _build(config_path,destination,Path(scratch))

def _build(config_path,destination,scratch):
 cfg=approved_configuration(config_path);destination=Path(destination).resolve()
 reject_protected_destination(destination,protected_result_roots(cfg))
 reject_protected_destination(scratch,protected_result_roots(cfg))
 if destination.exists():raise ValueError('Portable destination must be fresh; never overwrite a reader')
 # Validate all task seals/reviews before writing a single presentation file.
 native={k:native_seal(k,cfg['tasks'][k]) for k in SEALS}
 overrides={kind:layout_override(cfg['tasks'][kind],root,seal_path,paths) for kind,(root,seal_path,seal,paths) in native.items()}
 if any(overrides[k] is not None for k in ('composition','monocle3')):
  raise ValueError('This verified layout override is explicitly scoped to CellChat')
 velocity=cfg['tasks']['velocity'];assessment=verify(velocity['assessment_receipt'])
 velocity_native=json.loads(assessment.read_text())
 if velocity_native.get('task_status')!='INPUT_BLOCKED_NOT_EXECUTED' or velocity_native.get('scientific_jobs_started') is not False:
  raise ValueError('Velocity assessment no longer matches the declared blocked state')
 total_limit=int(cfg.get('maximum_package_bytes',45_000_000));max_table=int(cfg.get('maximum_table_bytes',2_000_000))
 preview_budget=int(cfg.get('maximum_embedded_table_preview_bytes',750_000))
 tasks=[];copies=[];omitted=[];compressed=[];planned_bytes=0
 for kind in TITLES:
  spec=cfg['tasks'][kind]
  task={'id':kind,'title':TITLES[kind],'headline':spec['headline'],'findings':spec['findings'],'limitations':spec['limitations'],
    'assets':[],'archive_url':spec.get('archive_url'),'archive_note':spec.get('archive_note','Large native models remain in the separately preserved native results.'),'reviews':[]}
  if kind in native:
   root,seal_path,seal,paths=native[kind];task['native_status']=seal['status'];task['native_schema']=seal['schema']
   # Native statuses/scopes remain verbatim for the visible provenance panel.
   task['native_summary']={k:seal[k] for k in ('n_scopes','n_models','n_models_inconclusive','scope_results','biological_temporal_root','genotype_DE_performed') if k in seal}
   for relative,path in paths.items():
    rel=Path(relative);include,reason=choose_asset(kind,rel,path,max_table)
    if kind=='monocle3' and path.name=='graph_test_all.tsv' and path.stat().st_size>max_table:
     # Only this native gene-result table is approved for compressed transport.
     encoded=lossless_gzip(path,scratch/kind/Path(relative+'.gz'))
     target=Path('native')/kind/Path(relative+'.gz')
     original={'path':relative,'bytes':path.stat().st_size,'sha256':sha(path)}
     compressed.append({'task':kind,'original_native':original,'transport':{'path':target.as_posix(),'bytes':encoded.stat().st_size,'sha256':sha(encoded)},'decoded_sha256':original['sha256'],'encoding':'gzip; mtime=0; empty filename; original bytes unchanged'})
     task['assets'].append({**descriptor(kind,rel,path,root,paths),'href':quote(target.as_posix(),safe='/'),'bytes':encoded.stat().st_size,'sha256':sha(encoded),'extension':'.gz','rendering':'Native table/source','transport':'Lossless gzip of complete original table','original_native_receipt':original})
     copies.append((encoded,target));planned_bytes+=encoded.stat().st_size
    elif include:
     target=Path('native')/kind/rel;descriptor_row=descriptor(kind,rel,path,root,paths)
     entry={**descriptor_row,'href':quote(target.as_posix(),safe='/'),'bytes':path.stat().st_size,'sha256':sha(path),'extension':path.suffix.lower(),
       'rendering':('Archival native (known display issues)' if kind=='cellchat' else 'Original native') if path.suffix.lower() in ('.png','.pdf') else 'Native table/source'}
     preview=table_preview(path, min(preview_budget,60_000))
     if preview:entry['preview']=preview;preview_budget-=sum(len(v) for row in preview for v in row)
     task['assets'].append(entry);copies.append((path,target));planned_bytes+=path.stat().st_size
    else:omitted.append({'task':kind,'native_path':rel.as_posix(),'bytes':path.stat().st_size,'sha256':sha(path),'reason':reason})
   copies.append((seal_path,Path('provenance')/kind/'FINAL_VALIDATION.json'));planned_bytes+=seal_path.stat().st_size
   task['seal_href']=f'provenance/{kind}/FINAL_VALIDATION.json'
   if kind=='cellchat':
    task['original_layout_warning']='Archival native figures preserve known display issues: some titles or long labels are clipped. Original mouse-boxplot dots include vertical jitter, so their plotted y position does not equal the saved raw score. Readable layout fixes display only and uses exact raw point heights when a reviewed repair is available.'
   override=overrides[kind]
   if override is not None:
    task['default_rendering']='Readable layout'
    task['layout_provenance']={'note':'Only presentation margins/headers and vertical jitter=0 were corrected. Original scientific values, aggregation, colors, scales and all native downloads remain unchanged.',
      'seal_href':f'provenance/{kind}/layout_FINAL_VALIDATION.json','mapping_href':f'provenance/{kind}/NATIVE_TO_PRESENTATION_PLOT_MAPPING.json'}
    for relative,original,corrected in override['checked']:
     target=Path('presentation')/kind/Path(relative)
     task['assets'].append({**descriptor(kind,Path(relative),original,root,paths),'href':quote(target.as_posix(),safe='/'),
       'bytes':corrected.stat().st_size,'sha256':sha(corrected),'extension':corrected.suffix.lower(),'rendering':'Readable layout',
       'original_href':quote((Path('native')/kind/Path(relative)).as_posix(),safe='/')})
     copies.append((corrected,target));planned_bytes+=corrected.stat().st_size
    for source,target in [(override['seal_path'],Path('provenance')/kind/'layout_FINAL_VALIDATION.json'),
                         (override['mapping_path'],Path('provenance')/kind/'NATIVE_TO_PRESENTATION_PLOT_MAPPING.json')]:
     copies.append((source,target));planned_bytes+=source.stat().st_size
    for i,record in enumerate(override['reviews'],1):
     source=verify(record);target=Path('reviews')/kind/f'layout-{i:02d}_{source.name}'
     copies.append((source,target));planned_bytes+=source.stat().st_size
     task['reviews'].append({'label':'Layout visual review: '+source.name,'href':quote(target.as_posix(),safe='/'),'sha256':record['sha256']})
  else:
   task['native_status']='INPUT_BLOCKED';task['native_schema']=velocity_native['schema'];task['native_summary']={'task_status':velocity_native['task_status'],'scientific_jobs_started':False}
   copies.append((assessment,Path('provenance/velocity/inventory_evidence.json')));planned_bytes+=assessment.stat().st_size
   task['seal_href']='provenance/velocity/inventory_evidence.json'
  for i,rec in enumerate(spec['review_receipts'],1):
   path=verify(rec);target=Path('reviews')/kind/f'{i:02d}_{path.name}'
   copies.append((path,target));planned_bytes+=path.stat().st_size
   task['reviews'].append({'label':path.name,'href':quote(target.as_posix(),safe='/'),'sha256':rec['sha256']})
  tasks.append(task)
 if planned_bytes>total_limit:raise ValueError(f'Portable assets require {planned_bytes:,} bytes, exceeding explicit limit {total_limit:,}; no assets silently dropped')
 # All paths came from verified, contained native receipts; no source reformatting.
 destination.mkdir(parents=True)
 copied=[]
 for source,target in copies:
  target=destination/target;target.parent.mkdir(parents=True,exist_ok=True)
  if target.exists():raise ValueError('Duplicate portable destination asset')
  shutil.copyfile(source,target)
  if sha(source)!=sha(target):raise ValueError('Asset copy hash mismatch')
  copied.append({'path':target.relative_to(destination).as_posix(),'bytes':target.stat().st_size,'sha256':sha(target)})
 nav={'title':cfg.get('title','SATB2 · Four follow-up analyses'),'created_utc':datetime.now(timezone.utc).isoformat(),'tasks':tasks,
  'reader_role':'Read-only navigation of preserved native results and approved interpretation. No computation or reanalysis occurs in this page.',
  'omitted_native_artifacts':omitted}
 template=(HERE/'template.html').read_text();payload=json.dumps(nav,ensure_ascii=False,separators=(',',':')).replace('<','\\u003c')
 (destination/'index.html').write_text(template.replace('/*__REVIEWED_NATIVE_NAVIGATION__*/',payload))
 (destination/'NATIVE_ARCHIVE_ONLY_FILES.json').write_text(json.dumps(omitted,indent=2)+'\n')
 (destination/'README.txt').write_text('Open index.html in a modern browser. Keep the entire folder together. No server or network is required for plots, native PDF/table links or embedded table previews. Selectors change only the displayed native assets. Original models/raw axes remain in separately preserved native archives.\n')
 (destination/'APPROVED_PRESENTATION_CONFIGURATION.json').write_text(Path(config_path).read_text())
 for p in [destination/'index.html',destination/'NATIVE_ARCHIVE_ONLY_FILES.json',destination/'README.txt',destination/'APPROVED_PRESENTATION_CONFIGURATION.json']:
  copied.append({'path':p.relative_to(destination).as_posix(),'bytes':p.stat().st_size,'sha256':sha(p)})
 package_bytes=sum(p.stat().st_size for p in destination.rglob('*') if p.is_file())
 if package_bytes>total_limit:raise ValueError('HTML/preview overhead exceeds package limit; preserved fresh draft requires explicit review, never delete native files')
 result={'schema':'satb2-followups-portable-reader-manifest/v1','status':'BUILT_PENDING_RENDERED_REVIEW','configuration':receipt(config_path),
  'package_bytes_before_manifest':package_bytes,'maximum_package_bytes':total_limit,'files':copied,
  'native_seals':{k:receipt(v[1]) for k,v in native.items()},'source_builder':receipt(__file__),
  'no_fit_rescore_or_native_mutation':True,'native_archive_only_file_count':len(omitted),'lossless_compressed_native_tables':compressed}
 (destination/'READER_MANIFEST.json').write_text(json.dumps(result,indent=2)+'\n')
 if package_bytes+(destination/'READER_MANIFEST.json').stat().st_size>total_limit:
  raise ValueError('Final receipt overhead exceeds explicit reader limit; no native files changed or silently omitted')
 print(json.dumps({'portable_directory':str(destination),'bytes':package_bytes,'files':len(copied),'status':result['status']}))
 return result

if __name__=='__main__':
 parser=argparse.ArgumentParser();parser.add_argument('--configuration',required=True,type=Path);parser.add_argument('--destination',required=True,type=Path)
 a=parser.parse_args();build(a.configuration,a.destination)
