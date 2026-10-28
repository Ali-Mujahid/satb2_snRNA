from pathlib import Path
import tempfile,importlib.util,hashlib,json,copy
src=Path(__file__).resolve().parent/'copy_completed_native.py'
spec=importlib.util.spec_from_file_location('recovery',src);m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
def write(p,x):p.parent.mkdir(parents=True,exist_ok=True);p.write_text(json.dumps(x)+'\n')
def fixture(t):
 old=t/'old';res=old/'results';new=t/'new';oc=old/'config.json';write(oc,{'schema':'SYNTHETIC_ONLY','run_root':str(old),'workflow_version':'1.0.1','scientific_key':'UNCHANGED'})
 write(res/'prepared/source_hashes.json',{'configuration':{'expected_sha256':m.digest(oc)}})
 toy=res/'prepared/SYNTHETIC.txt';toy.write_text('SYNTHETIC; NO STUDY DATA\n');write(res/'prepared/COMPLETE.json',{'status':'COMPLETE','files':{'SYNTHETIC.txt':{'bytes':toy.stat().st_size,'sha256':m.digest(toy)}}})
 write(res/'hdwgcna/workflow_status.json',{'primary_failures':[],'optional_consensus_failures':[]})
 write(res/'hdwgcna/run_parameters.json',{'script_md5':m.digest(src.parent/'run_hdwgcna.R','md5')})
 for lin in ['CA1','IN']:
  for net in ['pooled','consensus']:
   d=res/'hdwgcna'/lin/net;write(d/'status.json',{'status':'NETWORK_COMPLETE'});write(d/'run_identity.json',{'fingerprint':'SYNTHETIC_'+lin+net})
   write(d/'completed.json',{'fingerprint':'SYNTHETIC_'+lin+net,'files':{f.name:m.digest(f,'md5') for f in d.iterdir() if f.is_file()}})
 rows=[{'path':r,'bytes':(res/r).stat().st_size,'sha256':m.digest(res/r)} for r in sorted(m.tree_files(res))]
 inv=t/'inventory.json';write(inv,{'schema':'satb2-hdwgcna-copy-inventory/v1','status':'FROZEN_COMPLETE_FIT_COPY_SOURCE','source_results':str(res),'native_config':m.record(oc),'files':rows})
 c=json.loads(oc.read_text());c.update(run_root=str(new),workflow_version='1.0.2',reuse_completed_native={'schema':'satb2-hdwgcna-copy-recovery/v1','source_results':str(res),'native_config':m.record(oc),'inventory':m.record(inv),'no_network_refits':True})
 cfg=t/'recovery.json';write(cfg,c);return cfg,new/'results',res,c
with tempfile.TemporaryDirectory(prefix='satb2_copy_SYNTHETIC_') as root:
 t=Path(root);cfg,out,res,c=fixture(t/'pass');m.run(cfg,out);a=m.load(out/'REUSED_NATIVE_RECEIPT.json');assert a['scientific_fits_run'] is False and a['status']=='VERIFIED_EXACT_COPY';m.verify_rows(res,a['copied_files']);m.verify_rows(out,a['copied_files'])
 for case in ['changed-science','changed-source','missing-source','symlink','nonempty-destination']:
  cfg,out,res,c=fixture(t/case)
  if case=='changed-science':c['scientific_key']='CHANGED';write(cfg,c)
  if case=='changed-source':(res/'prepared/SYNTHETIC.txt').write_text('CORRUPTED\n')
  if case=='missing-source':(res/'prepared/SYNTHETIC.txt').unlink()
  if case=='symlink':(res/'prepared/SYMLINK').symlink_to(res/'prepared/SYNTHETIC.txt')
  if case=='nonempty-destination':write(out/'prepared/KEEP.json',{'preserve':True})
  try:m.run(cfg,out)
  except ValueError:pass
  else:raise AssertionError(case+' was not rejected')
 print('PASS six synthetic copy/immutability gates; no study reads, transformation or fitting')
