#!/usr/bin/env python3
"""Lossless web transport for the completed offline SATB2 atlas.

Only unwraps existing base64 gzip bytes. No numerical reserialization, data
selection, RNA aggregation, fitting or change to source scientific payloads.
Frontend code is owned separately. The original offline atlas stays unchanged.
"""
from __future__ import annotations
import argparse,base64,collections,copy,csv,gzip,hashlib,json,re,time
from pathlib import Path

def sha_bytes(data):return hashlib.sha256(data).hexdigest()
def sha_file(path):
 h=hashlib.sha256()
 with Path(path).open('rb') as f:
  for b in iter(lambda:f.read(8*1024*1024),b''):h.update(b)
 return h.hexdigest()
def json_bytes(data):return (json.dumps(data,ensure_ascii=True,separators=(',',':'),allow_nan=False)+'\n').encode()
def write_json(path,data):
 path=Path(path);path.parent.mkdir(parents=True,exist_ok=True);temp=path.with_suffix(path.suffix+'.tmp');temp.write_bytes(json_bytes(data));temp.replace(path)
def write_tsv(path,rows,fields=None):
 path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
 with path.open('w',encoding='utf-8',newline='') as f:
  w=csv.DictWriter(f,fieldnames=fields or list(rows[0]),delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rows)
def write_identical(path,data):
 path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
 if path.exists():
  if path.read_bytes()!=data:raise ValueError('Refusing different binary shard at existing path: '+str(path))
 else:path.write_bytes(data)
 if path.read_bytes()!=data:raise ValueError('Written gzip bytes differ from source')
def load_json(path):
 with Path(path).open(encoding='utf-8') as f:return json.load(f)
def unwrap(source,descriptor):
 raw_script=Path(source).read_bytes();script_sha=sha_bytes(raw_script)
 if descriptor.get('sha256') and descriptor['sha256']!=script_sha:raise ValueError('Original source script SHA mismatch: '+str(source))
 strings=re.findall(rb'(?:\bdata\s*:|"data"\s*:)\s*("[A-Za-z0-9+/=]+")',raw_script)
 if len(strings)!=1:raise ValueError('Expected one native gzip-base64 field: '+str(source))
 packed=base64.b64decode(json.loads(strings[0]),validate=True)
 if packed[:2]!=b'\x1f\x8b':raise ValueError('Native shard does not contain gzip')
 if descriptor.get('gzip_bytes') is not None and len(packed)!=descriptor['gzip_bytes']:raise ValueError('Native gzip byte count differs')
 payload_raw=gzip.decompress(packed);payload_sha=sha_bytes(payload_raw)
 if descriptor.get('payload_sha256') and descriptor['payload_sha256']!=payload_sha:raise ValueError('Native payload SHA differs')
 if descriptor.get('uncompressed_bytes') is not None and len(payload_raw)!=descriptor['uncompressed_bytes']:raise ValueError('Native decoded byte count differs')
 payload=json.loads(payload_raw,parse_constant=lambda x:(_ for _ in ()).throw(ValueError('Nonstandard JSON number '+x)))
 is_download=isinstance(payload,dict) and isinstance(payload.get('text'),str)
 if is_download:
  if 'sha256' not in payload or sha_bytes(payload['text'].encode('utf-8'))!=payload['sha256']:raise ValueError('Exact native text download SHA differs')
 return packed,payload_raw,payload,script_sha,payload_sha,is_download

def build(args):
 source=Path(args.source_bundle).resolve();repo=Path(args.repo_root).resolve();site=repo/'site';staging=repo/'publication_build';staging.mkdir(parents=True,exist_ok=True)
 if (staging/'WEB_DATA_BUILD_COMPLETE.json').exists():raise ValueError('Refusing a second closed web transport build')
 stamp=time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime());raw_prefix=f'https://raw.githubusercontent.com/{args.owner}/{args.repository}/{args.release_tag}/data_downloads/'
 release_prefix=f'https://github.com/{args.owner}/{args.repository}/releases/download/{args.release_tag}/'
 inputs={'base_manifest':source/'data/manifest.json','de_catalog':source/'DE_CATALOG.json','mouse_catalog':source/'data/MOUSE_CATALOG.json','asset_map':source/'data/asset_map.json'}
 originals={k:load_json(v) for k,v in inputs.items()};original_sha={str(p.relative_to(source)):sha_file(p) for p in inputs.values()}
 for name in ['annotation_suggestions.js','ANNOTATION_REVIEW_PROVENANCE.json']:
  original_sha['data/'+name]=sha_file(source/'data'/name)
 if originals['de_catalog'].get('status')!='COMPLETE' or originals['mouse_catalog'].get('status')!='COMPLETE':raise ValueError('Source DE/mouse catalogs not complete')
 catalogs={k:copy.deepcopy(originals[k]) for k in ['base_manifest','de_catalog','mouse_catalog']};registry={};owners=collections.defaultdict(list)
 for name,cat in catalogs.items():
  for ident,desc in cat.get('shards',{}).items():
   if not re.fullmatch(r'[A-Za-z0-9_.-]+',ident):raise ValueError('Unsafe/noncanonical shard ID')
   if ident in registry and registry[ident]!=desc:raise ValueError('Different source descriptors share an ID')
   registry[ident]=desc;owners[ident].append(name)
 physical={str(p.resolve()) for d in ['data/shards','data/mouse_shards'] for p in (source/d).glob('*.js')}
 registered={str((source/d['path']).resolve()) for d in registry.values()}
 if physical!=registered:raise ValueError('Unregistered or missing lazy source scripts: '+str({'extra':len(physical-registered),'missing':len(registered-physical)}))
 if len(registry)!=3140:raise ValueError('Expected every3140 original lazy shard ID')
 rows=[];new_registry={};kind_stats=collections.defaultdict(lambda:{'count':0,'source_script_bytes':0,'gzip_bytes':0,'decoded_json_bytes':0});download_ids=[]
 for n,(ident,desc) in enumerate(registry.items(),1):
  source_path=source/desc['path'];packed,raw,payload,script_sha,payload_sha,is_download=unwrap(source_path,desc)
  # Native gzip is copied byte-for-byte; decoded JSON is never serialized again.
  webpath=('data_downloads/' if is_download else 'site/data/shards/')+ident+'.json.gz';destination=repo/webpath;write_identical(destination,packed)
  if gzip.decompress(destination.read_bytes())!=raw:raise ValueError('Strict decoded JSON byte comparison failed')
  meta={k:v for k,v in desc.items() if k not in ['path','sha256','encoding','source_counts_path']}
  meta.update({'path':raw_prefix+ident+'.json.gz' if is_download else 'data/shards/'+ident+'.json.gz','encoding':'gzip-bytes','sha256':sha_bytes(packed),'source_script_sha256':script_sha,'payload_sha256':payload_sha,'gzip_bytes':len(packed),'uncompressed_bytes':len(raw),'source_atlas_path':desc['path']})
  new_registry[ident]=meta
  if is_download:download_ids.append(ident)
  kind=ident.split('__')[0];s=kind_stats[kind];s['count']+=1;s['source_script_bytes']+=source_path.stat().st_size;s['gzip_bytes']+=len(packed);s['decoded_json_bytes']+=len(raw)
  rows.append({'shard_id':ident,'catalog_owners':';'.join(owners[ident]),'publication_role':'EXACT_SOURCE_TEXT_DOWNLOAD' if is_download else 'INTERACTIVE_DATA','source_atlas_path':desc['path'],'published_repo_path':webpath,'published_url_or_site_path':meta['path'],'source_script_sha256':script_sha,'gzip_sha256':meta['sha256'],'decoded_json_sha256':payload_sha,'source_script_bytes':source_path.stat().st_size,'gzip_bytes':len(packed),'decoded_json_bytes':len(raw),'strict_decoded_json_byte_comparison':'PASS'})
  if n%100==0:print(f'Unwrapped and verified {n}/{len(registry)} source shards',flush=True)
 if len(download_ids)!=2063:raise ValueError('Expected all2063 exact text download shards')
 for cat in catalogs.values():cat['shards']={ident:new_registry[ident] for ident in cat['shards']}
 # Every existing downloadable PDF/document receives one collision-free direct asset URL.
 public_assets=[];local_assets=[];asset_by_old={};used_names=set()
 for item in originals['asset_map']:
  old=item['bundle_path'];rel=Path(old)
  if not rel.parts or rel.parts[0]!='assets' or '..' in rel.parts:raise ValueError('Unsafe asset source path')
  asset_name='__'.join(rel.parts[1:])
  if asset_name in used_names:raise ValueError('Deterministic release filename collision')
  used_names.add(asset_name);p=source/old;digest=sha_file(p)
  if digest!=item['sha256']:raise ValueError('Native downloadable asset differs')
  url=release_prefix+asset_name
  record={'atlas_relative_path':old,'release_asset_name':asset_name,'url':url,'sha256':digest,'bytes':p.stat().st_size};public_assets.append(record);local_assets.append({'source_path':str(p),**record});asset_by_old[old]=record
 if len(public_assets)!=58:raise ValueError('Expected all58 original downloadable PDF/document assets')
 for d in catalogs['base_manifest'].get('downloads',[]):
  old=d['path']
  if old not in asset_by_old:raise ValueError('Source download has no exact published asset mapping: '+old)
  d['source_atlas_path']=old;d['path']=asset_by_old[old]['url'];d['sha256']=asset_by_old[old]['sha256']
 for d in catalogs['de_catalog'].get('review_metadata',{}).get('review_documents',[]):
  old=d.get('path') or d.get('bundle_path')
  if old not in asset_by_old:raise ValueError('Review document has no exact release mapping: '+str(old))
  d['source_atlas_path']=old;d['path']=asset_by_old[old]['url']
  if 'bundle_path' in d:d['bundle_path']=d['path']
  d['sha256']=asset_by_old[old]['sha256']
 for name,cat in catalogs.items():cat['web_transport']={'version':'atlas-web-gzip-1','scientific_payloads_changed':False,'native_gzip_bytes_preserved':True,'source_catalog_sha256':original_sha[str(inputs[name].relative_to(source))],'release_tag':args.release_tag,'external_exact_download_prefix':raw_prefix}
 # Keep the source linkage and bind the corresponding public transport catalog too.
 de=catalogs['de_catalog'];de_bytes=json_bytes(de);write_identical(site/'DE_CATALOG.json',de_bytes);write_identical(site/'data/DE_CATALOG.json',de_bytes)
 mouse=catalogs['mouse_catalog'];mouse['source_de_catalog_sha256']=mouse['de_catalog_sha256'];mouse['de_catalog_sha256']=sha_bytes(de_bytes)
 writes=[('base_manifest','manifest.json','manifest.js','SATB2_DATA'),('de_catalog',None,'de_catalog.js','SATB2_DE_CATALOG'),('mouse_catalog','MOUSE_CATALOG.json','mouse_catalog.js','SATB2_MOUSE_CATALOG')]
 for name,json_name,js_name,glob in writes:
  obj=catalogs[name]
  if json_name:write_identical(site/'data'/json_name,json_bytes(obj))
  write_identical(site/'data'/js_name,('window.'+glob+'='+json.dumps(obj,ensure_ascii=True,separators=(',',':'),allow_nan=False)+';\n').encode())
 for name in ['annotation_suggestions.js','ANNOTATION_REVIEW_PROVENANCE.json']:
  write_identical(site/'data'/name,(source/'data'/name).read_bytes())
 write_tsv(site/'data/PUBLIC_DOWNLOAD_ASSETS.tsv',public_assets);write_json(site/'data/PUBLIC_DOWNLOAD_ASSETS.json',public_assets)
 write_tsv(staging/'PUBLIC_DOWNLOAD_ASSET_INDEX.tsv',local_assets);write_tsv(staging/'SHARD_TRANSPORT_INDEX.tsv',rows)
 public_provenance={'status':'PASS','transport':'NATIVE_GZIP_BASE64_UNWRAPPED_TO_IDENTICAL_BINARY_GZIP','source_catalog_sha256':original_sha,'all_source_lazy_shards_preserved':True,'n_source_shards':len(rows),'n_interactive_shards':len(rows)-len(download_ids),'n_exact_text_download_shards':len(download_ids),'n_direct_release_assets':len(public_assets),'n_datasets':len(originals['de_catalog']['datasets']),'n_endpoints':len(originals['mouse_catalog']['fits']),'n_mouse_count_arrays':len(originals['mouse_catalog']['shards']),'n_full_source_genes':len(originals['mouse_catalog']['genes']),'numeric_payloads_reserialized':False,'strict_decoded_JSON_byte_comparison':'PASS_ALL3140','new_scientific_fitting_performed':False,'RNA_extraction_or_aggregation_performed':False,'source_atlas_modified':False,'release_tag':args.release_tag,'builder_code_sha256':sha_file(__file__)}
 write_json(site/'data/WEB_DATA_PROVENANCE.json',public_provenance)
 site_files=[p for p in site.rglob('*') if p.is_file()];site_bytes=sum(p.stat().st_size for p in site_files);data_bytes=sum((repo/r['published_repo_path']).stat().st_size for r in rows if r['publication_role']=='INTERACTIVE_DATA');download_bytes=sum((repo/r['published_repo_path']).stat().st_size for r in rows if r['publication_role']=='EXACT_SOURCE_TEXT_DOWNLOAD')
 if site_bytes>=args.max_site_bytes:raise ValueError(f'Current site budget {site_bytes} exceeds {args.max_site_bytes}')
 owned_static=['DE_CATALOG.json']+['data/'+n for n in ['manifest.json','manifest.js','DE_CATALOG.json','de_catalog.js','MOUSE_CATALOG.json','mouse_catalog.js','PUBLIC_DOWNLOAD_ASSETS.tsv','PUBLIC_DOWNLOAD_ASSETS.json','WEB_DATA_PROVENANCE.json','annotation_suggestions.js','ANNOTATION_REVIEW_PROVENANCE.json']]
 outputs={str((site/n).relative_to(repo)):sha_file(site/n) for n in owned_static};outputs.update({r['published_repo_path']:r['gzip_sha256'] for r in rows})
 receipt={**public_provenance,'completed_utc':stamp,'source_bundle_path':str(source),'repo_root':str(repo),'site_bytes_at_data_build':site_bytes,'interactive_gzip_bytes':data_bytes,'external_source_download_gzip_bytes':download_bytes,'site_budget_bytes':args.max_site_bytes,'remaining_site_budget_bytes':args.max_site_bytes-site_bytes,'kind_budgets':kind_stats,'original_catalog_input_sha256':original_sha,'output_sha256':outputs,'frontend_files_owned_separately':True,'remote_network_mutations_performed':False}
 write_json(staging/'WEB_DATA_BUILD_COMPLETE.json',receipt)
 print(json.dumps({k:receipt[k] for k in ['status','n_source_shards','n_interactive_shards','n_exact_text_download_shards','n_direct_release_assets','site_bytes_at_data_build','interactive_gzip_bytes','external_source_download_gzip_bytes','remaining_site_budget_bytes']},indent=2),flush=True)

if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--source-bundle',required=True);ap.add_argument('--repo-root',required=True);ap.add_argument('--owner',default='Ali-Mujahid');ap.add_argument('--repository',default='satb2_snRNA');ap.add_argument('--release-tag',default='atlas-v2.1.0');ap.add_argument('--max-site-bytes',type=int,default=950000000);build(ap.parse_args())
