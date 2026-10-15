"""Validate published transport identities without executing scientific code."""
from pathlib import Path
import gzip,hashlib,json,re,sys

ROOT=Path(__file__).resolve().parent.parent
SITE=ROOT/'site'
PREFIX='https://raw.githubusercontent.com/Ali-Mujahid/satb2_snRNA/atlas-v2.1.0/data_downloads/'

def main():
    maps=[json.loads((SITE/'data'/name).read_text()) for name in ['manifest.json','de_catalog.json','MOUSE_CATALOG.json']]
    count=0
    for registry in maps:
        for key,d in registry['shards'].items():
            path=d['path']
            if path.startswith(PREFIX):
                name=path[len(PREFIX):]
                assert re.fullmatch(r'[A-Za-z0-9_.-]+\.json\.gz',name),name
                p=ROOT/'data_downloads'/name
            else:
                assert path.startswith('data/shards/') and '..' not in path and '://' not in path,path
                p=SITE/path
            packed=p.read_bytes();assert hashlib.sha256(packed).hexdigest()==d['sha256'],key
            raw=gzip.decompress(packed);assert hashlib.sha256(raw).hexdigest()==d['payload_sha256'],key
            count+=1
    assert count==3140,count
    total=sum(p.stat().st_size for p in SITE.rglob('*') if p.is_file())
    assert total<950_000_000,total
    html=(SITE/'index.html').read_text();ids=re.findall(r'\bid="([^"]+)"',html)
    assert len(ids)==len(set(ids))
    for h in SITE.glob('*.html'):
        for ref in re.findall(r'(?:src|href)="([^"]+)"',h.read_text()):
            if ref.startswith(('https://','http://','#','data:','mailto:')):continue
            path=ref.split('#')[0].split('?')[0]
            if path:assert (h.parent/path).is_file(),(h.name,ref)
    assert len(maps[0]['datasets'])==24
    assert len(maps[2]['genes'])==32285 and len(maps[2]['fits'])==1060
    print(json.dumps({'status':'PASS','shards':count,'site_bytes':total,'n_datasets':24,'n_mouse_genes':32285,'n_mouse_endpoints':1060,'scientific_execution_performed':False}))

if __name__=='__main__':main()
