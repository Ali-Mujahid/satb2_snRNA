#!/usr/bin/env python3
"""Copy byte-verified finished native fits into fresh export recovery; no fitting."""
from pathlib import Path
import argparse, hashlib, json, shutil

def need(ok, message):
    if not ok: raise ValueError(message)

def digest(path, algo='sha256'):
    with Path(path).open('rb') as f: return hashlib.file_digest(f, algo).hexdigest()

def load(path): return json.loads(Path(path).read_text())

def record(path):
    q=Path(path)
    need(q.is_file() and not q.is_symlink(), f'Expected regular file: {q}')
    return {'path':str(q.resolve()),'bytes':q.stat().st_size,'sha256':digest(q)}

def check(r):
    need(record(r['path'])==r, f'Changed bound source: {r["path"]}')

def child(root, rel):
    root=Path(root).resolve()
    q=Path(rel)
    need(not q.is_absolute() and q.parts and q.parts[0] in ['prepared','hdwgcna'] and '..' not in q.parts, f'Unsafe inventory path: {rel}')
    p=(root/q).resolve()
    need(p.is_relative_to(root), 'Inventory escapes source root')
    return p

def tree_files(root):
    out=set()
    for tree in ['prepared','hdwgcna']:
        q=root/tree
        need(q.is_dir() and not q.is_symlink(), f'Invalid tree: {q}')
        for f in q.rglob('*'):
            need(not f.is_symlink(), f'Refusing symlink: {f}')
            if f.is_file(): out.add(str(f.relative_to(root)))
    return out

def verify_rows(root, rows):
    need(len({r['path'] for r in rows})==len(rows), 'Duplicate inventory path')
    need(tree_files(root)=={r['path'] for r in rows}, 'Native file tree differs from frozen inventory')
    for r in rows:
        p=child(root,r['path'])
        need(p.stat().st_size==r['bytes'] and digest(p)==r['sha256'], f'Native bytes changed: {p}')

def run(config_path, results):
    config_path,results=Path(config_path).resolve(),Path(results).resolve()
    cfg_record=record(config_path);cfg=load(config_path); reuse=cfg['reuse_completed_native']
    need(reuse['schema']=='satb2-hdwgcna-copy-recovery/v1' and reuse['no_network_refits'] is True, 'Not a copy-only recovery')
    need(results==Path(cfg['run_root']).resolve()/'results', 'Wrong recovery results root')
    for k in ['native_config','inventory']: check(reuse[k])
    original=load(reuse['native_config']['path'])
    ignore={'run_root','workflow_version','reuse_completed_native'}
    need({k:v for k,v in cfg.items() if k not in ignore}=={k:v for k,v in original.items() if k not in ignore}, 'Scientific configuration changed; fits cannot be reused')
    source=Path(reuse['source_results']).resolve()
    need(source==Path(original['run_root']).resolve()/'results' and source!=results and not results.is_relative_to(source) and not source.is_relative_to(results), 'Invalid original/recovery root relationship')
    inv=load(reuse['inventory']['path'])
    need(inv['schema']=='satb2-hdwgcna-copy-inventory/v1' and inv['status']=='FROZEN_COMPLETE_FIT_COPY_SOURCE', 'Invalid completed-fit inventory')
    need(inv['native_config']==reuse['native_config'] and Path(inv['source_results']).resolve()==source, 'Inventory source identity differs')
    rows=inv['files'];verify_rows(source,rows)
    prep=load(source/'prepared/COMPLETE.json')
    need(prep.get('status')=='COMPLETE', 'Original preparation is incomplete')
    for rel,r in prep['files'].items():
        p=source/'prepared'/rel
        need(p.is_file() and p.stat().st_size==r['bytes'] and digest(p)==r['sha256'], 'Original preparation receipt mismatch')
    need(load(source/'prepared/source_hashes.json')['configuration']['expected_sha256']==reuse['native_config']['sha256'], 'Prepared scientific config differs')
    status=load(source/'hdwgcna/workflow_status.json')
    need(not status.get('primary_failures') and not status.get('optional_consensus_failures'), 'Native fit failures cannot be relabeled complete')
    for lineage in ['CA1','IN']:
        for network in ['pooled','consensus']:
            d=source/'hdwgcna'/lineage/network
            c=load(d/'completed.json')
            need(load(d/'status.json')['status']=='NETWORK_COMPLETE' and c['fingerprint']==load(d/'run_identity.json')['fingerprint'], 'Unfinished/mismatched native checkpoint')
            for rel,h in c['files'].items():
                p=(d/rel).resolve();need(p.is_relative_to(d) and p.is_file() and digest(p,'md5')==h, 'Native checkpoint MD5 mismatch')
    native_r=Path(__file__).resolve().parent/'run_hdwgcna.R'
    need(digest(native_r,'md5')==load(source/'hdwgcna/run_parameters.json')['script_md5'], 'Recovery changed the native fitting script')
    need(not (results/'REUSED_NATIVE_RECEIPT.json').exists(), 'Refusing to replace reuse receipt')
    for tree in ['prepared','hdwgcna']:
        q=results/tree
        need(not q.exists(), f'Refusing to reuse destination tree: {q}')
    results.mkdir(parents=True,exist_ok=True)
    for r in rows:
        src=child(source,r['path']);dst=child(results,r['path'])
        dst.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(src,dst)
    verify_rows(results,rows);verify_rows(source,rows)
    for k in ['native_config','inventory']: check(reuse[k])
    need(record(config_path)==cfg_record,'Recovery configuration changed during copy')
    out={'schema':'satb2-hdwgcna-reused-native-receipt/v1','status':'VERIFIED_EXACT_COPY','source_results':str(source),'native_config':reuse['native_config'],'inventory':reuse['inventory'],'recovery_config':cfg_record,'copied_files':rows,'scientific_fits_run':False,'all_four_completed_checkpoints_MD5_verified':True,'prepared_complete_SHA256_verified':True,'all_source_and_destination_SHA256_verified':True,'scientific_configuration_unchanged':True,'native_fitting_script_MD5_unchanged':True}
    with (results/'REUSED_NATIVE_RECEIPT.json').open('x') as f: json.dump(out,f,indent=2);f.write('\n')
    print(json.dumps({'status':out['status'],'files':len(rows),'bytes':sum(r['bytes'] for r in rows),'scientific_fits_run':False}))

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--config',required=True);p.add_argument('--results',required=True);a=p.parse_args();run(a.config,a.results)
