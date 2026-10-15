"""Read-only provenance freezing for the native downstream revision."""
from pathlib import Path
import argparse, hashlib, json, shutil
import pandas as pd

def sha(p):
    h=hashlib.sha256()
    with open(p,'rb') as f:
        for block in iter(lambda:f.read(8*1024*1024),b''): h.update(block)
    return h.hexdigest()

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--config',required=True);ap.add_argument('--outdir',required=True);a=ap.parse_args()
    cfg=json.loads(Path(a.config).read_text());out=Path(a.outdir)
    if out.exists() and any(out.iterdir()):raise RuntimeError('Refusing output overwrite')
    out.mkdir(parents=True,exist_ok=True)
    old=Path(cfg['previous_run_root']);root=Path(cfg['study_root'])
    if root==old or old in root.parents:raise RuntimeError('Revision must be separate')
    rows=[]
    for path in cfg['frozen_input_paths']:
        p=Path(path).resolve(strict=True);rows.append({'path':str(p),'sha256':sha(p),'bytes':p.stat().st_size})
    pd.DataFrame(rows).to_csv(out/'INPUT_SHA256.tsv',sep='\t',index=False)
    shutil.copyfile(a.config,out/'MASTER_CONFIG.json')
    (out/'COMPLETE.json').write_text(json.dumps({'status':'COMPLETE','stage':'freeze_revision','inputs':len(rows),'config_sha256':sha(a.config),'previous_analysis_read_only':True},indent=2)+'\n')
