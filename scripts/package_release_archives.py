"""Repackage existing files only; never run or modify scientific analysis."""
from pathlib import Path
import argparse,hashlib,json,zipfile

def sha(path):
    h=hashlib.sha256()
    with path.open('rb') as f:
        for b in iter(lambda:f.read(1024*1024),b''):h.update(b)
    return h.hexdigest()

def main():
    a=argparse.ArgumentParser();a.add_argument('--source-zip',type=Path,required=True)
    a.add_argument('--repo',type=Path,required=True);a.add_argument('--out',type=Path,required=True)
    args=a.parse_args();args.out.mkdir(parents=True,exist_ok=True)
    runtime=args.out/'SATB2_atlas_v2_runtime.zip';saved=args.out/'SATB2_atlas_v2_saved_results.zip'
    assets=[]
    with zipfile.ZipFile(args.source_zip) as source,zipfile.ZipFile(runtime,'w',compression=zipfile.ZIP_DEFLATED,compresslevel=6,allowZip64=True) as r,zipfile.ZipFile(saved,'w',compression=zipfile.ZIP_DEFLATED,compresslevel=6,allowZip64=True) as s:
        for info in source.infolist():
            target=s if info.filename.startswith('offline_dashboard_v2/saved_results/') else r
            with source.open(info) as inp,target.open(info,'w',force_zip64=True) as out:
                for b in iter(lambda:inp.read(1024*1024),b''):out.write(b)
    for p,kind in [(runtime,'offline_runtime'),(saved,'offline_saved_results')]:
        assert p.stat().st_size<2*1024**3,(p,p.stat().st_size)
        assets.append({'asset_name':p.name,'local_path':str(p),'bytes':p.stat().st_size,'sha256':sha(p),'kind':kind})
    code=args.out/'SATB2_analysis_source_v2.zip'
    with zipfile.ZipFile(code,'w',compression=zipfile.ZIP_DEFLATED,compresslevel=6,allowZip64=True) as z:
        for sub in ['analysis','docs','scripts']:
            for p in sorted((args.repo/sub).rglob('*')):
                if p.is_file() and '__pycache__' not in p.parts and p.suffix!='.pyc':z.write(p,Path('satb2_analysis_source_v2')/p.relative_to(args.repo))
        for name in ['README.md','.github/workflows/pages.yml']:
            p=args.repo/name
            if p.exists():z.write(p,Path('satb2_analysis_source_v2')/name)
    assets.append({'asset_name':code.name,'local_path':str(code),'bytes':code.stat().st_size,'sha256':sha(code),'kind':'analysis_source'})
    manifest={'status':'PREPARED','source_portable_zip_sha256':sha(args.source_zip),'assets':assets}
    (args.out/'ARCHIVE_ASSETS.json').write_text(json.dumps(manifest,indent=2)+'\n')
    print(json.dumps(manifest))

if __name__=='__main__':main()
