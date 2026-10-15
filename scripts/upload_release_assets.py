"""Resumable, source-SHA-bound GitHub CLI uploads; never print credentials."""
import argparse,concurrent.futures,hashlib,json,os,subprocess,time
from pathlib import Path

def digest(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for b in iter(lambda:f.read(1024*1024),b''):h.update(b)
    return h.hexdigest()

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--manifest',type=Path,required=True)
    ap.add_argument('--gh',required=True);ap.add_argument('--repo',default='Ali-Mujahid/satb2_snRNA')
    ap.add_argument('--tag',default='atlas-v2.1.0');ap.add_argument('--workers',type=int,default=3)
    ap.add_argument('--receipt',type=Path,required=True);a=ap.parse_args()
    env={**os.environ,'GH_PAGER':'cat','NO_COLOR':'1'}
    def gh(*args):return subprocess.run([a.gh,*args],capture_output=True,text=True,env=env)
    rows=json.loads(a.manifest.read_text());rows=rows.get('files',rows.get('assets',[]))
    release=json.loads(gh('api','repos/'+a.repo+'/releases/tags/'+a.tag).stdout or 'null')
    if not isinstance(release,dict) or 'id' not in release:
        ls=json.loads(gh('api','repos/'+a.repo+'/releases?per_page=100').stdout)
        release=next(x for x in ls if x['tag_name']==a.tag)
    rid=release['id'];observed=json.loads(gh('api',f'repos/{a.repo}/releases/{rid}/assets?per_page=100').stdout)
    existing={x['name']:x for x in observed};completed=[];errors=[];started=time.time()
    def job(row):
        p=Path(row['local_path']);name=row.get('asset_name',p.name)
        assert p.name==name and p.stat().st_size==row['bytes'] and row['bytes']<2*1024**3,name
        assert digest(p)==row['sha256'],name
        old=existing.get(name)
        if old:
            if old.get('digest')=='sha256:'+row['sha256'] and old['size']==row['bytes']:
                return {'asset_name':name,'status':'ALREADY_VERIFIED','bytes':row['bytes'],'sha256':row['sha256'],'id':old['id']}
            raise RuntimeError('Existing differently identified release asset: '+name)
        for attempt in range(3):
            r=gh('release','upload',a.tag,str(p),'--repo',a.repo)
            if r.returncode==0:break
            if attempt==2:raise RuntimeError(name+': '+r.stderr[-600:])
            time.sleep(5*(attempt+1))
        asset=json.loads(gh('api',f'repos/{a.repo}/releases/{rid}/assets?per_page=100').stdout)
        actual=next(x for x in asset if x['name']==name)
        assert actual['size']==row['bytes'] and actual.get('digest')=='sha256:'+row['sha256'],name
        return {'asset_name':name,'status':'UPLOADED_VERIFIED','bytes':row['bytes'],'sha256':row['sha256'],'id':actual['id']}
    a.receipt.parent.mkdir(parents=True,exist_ok=True)
    with concurrent.futures.ThreadPoolExecutor(max_workers=a.workers) as pool:
        jobs={pool.submit(job,row):row for row in rows}
        for f in concurrent.futures.as_completed(jobs):
            try:
                r=f.result();completed.append(r);print(json.dumps(r),flush=True)
            except Exception as e:
                r={'asset_name':jobs[f].get('asset_name'),'error':str(e)};errors.append(r);print(json.dumps(r),flush=True)
            receipt={'status':'COMPLETE' if len(completed)==len(rows) else 'UPLOADING','release_id':rid,'tag':a.tag,'completed':completed,'errors':errors,'expected_assets':len(rows),'elapsed_seconds':round(time.time()-started,1)}
            a.receipt.write_text(json.dumps(receipt,indent=2)+'\n')
    if errors:raise SystemExit('Some uploads failed; receipt retains verified completions.')

if __name__=='__main__':main()
