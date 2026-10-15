#!/usr/bin/env python3
"""Register hash-verified v3 release assets without changing historical URLs."""
import argparse, hashlib, json
from pathlib import Path

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--assets',type=Path,required=True);a=ap.parse_args()
    root=Path(__file__).resolve().parents[1];path=root/'site/data/published_downloads.js'
    data=json.loads(path.read_text().split('=',1)[1].strip().rstrip(';'))
    rows=json.loads(a.assets.read_text())['files'];prefix='https://github.com/Ali-Mujahid/satb2_snRNA/releases/download/atlas-v3.0.0/'
    new=[];figures=[]
    for row in rows:
        source=Path(row['local_path']);name=row['asset_name']
        assert source.name==name and source.stat().st_size==row['bytes']
        h=hashlib.sha256()
        with source.open('rb') as f:
            for b in iter(lambda:f.read(8*1024*1024),b''):h.update(b)
        assert h.hexdigest()==row['sha256'],name
        if name.endswith('.zip'):
            runtime='runtime' in name
            new.append({'label':('Recommended v3 offline runtime · screenshot panel and ordering' if runtime else 'V3 companion · unchanged native saved figures and results'),
                'url':prefix+name,'release_filename':name,'bytes':row['bytes'],'sha256':row['sha256'],
                'kind':'archive','recommended':True,
                'scope':'Extract with '+('SATB2_atlas_v3_saved_results.zip' if runtime else 'SATB2_atlas_v3_runtime.zip')+' into the same parent folder; common offline_dashboard_v3/ root. Earlier releases preserved.'})
        elif name.endswith('.pdf'):
            assert name.startswith(('v3_screenshot75_','v3_pattern75_'))
            tool='scanpy' if '_scanpy_' in name else 'seurat';fixed='screenshot75' in name
            title=('Exact screenshot order · 75 genes' if fixed else 'Cluster expression pattern · 75 genes · Scanpy reference')+' · '+tool
            caption=('Exact image gene order and fifteen five-gene blocks.' if fixed else 'One shared gene order defined by descriptive Scanpy cluster expression contrasts, used in both tools.')+' Inhibitory original-full cohort; source cluster rows descending. Dot area represents detected fraction; negative z dots are smaller and fainter. Blocks do not establish cross-dataset cluster identities.'
            new.append({'label':title,'url':prefix+name,'release_filename':name,'bytes':row['bytes'],'sha256':row['sha256'],'kind':'example_pdf','scope':caption})
            png='assets/gallery/'+source.stem+'.png';assert (root/'site'/png).is_file(),png
            figures.append({'title':title,'png_url':png,'pdf_url':prefix+name,'caption':caption,'tool':tool})
    names={r['release_filename'] for r in new}
    historical=[r for r in data['assets'] if r['release_filename'] not in names]
    for r in historical:
        if r.get('recommended'):
            r['recommended']=False
            if not r['label'].startswith('Legacy'):r['label']='Legacy v2.1 · '+r['label']
    data['assets']=new+historical;data['version']='atlas-v3.0.0';data['release_url']='https://github.com/Ali-Mujahid/satb2_snRNA/releases/tag/atlas-v3.0.0'
    data['figures']=figures+[r for r in data['figures'] if not r['pdf_url'].startswith(prefix)]
    path.write_text('window.SATB2_PUBLISHED_DOWNLOADS='+json.dumps(data,separators=(',',':'),ensure_ascii=False)+';\n')
    clean=[{k:v for k,v in row.items() if k!='local_path'} for row in rows]
    (root/'site/data/V3_DOWNLOAD_ASSETS.json').write_text(json.dumps({'version':'atlas-v3.0.0','assets':clean},indent=2)+'\n')
    print(json.dumps({'status':'PASS','added_assets':len(new),'all_assets':len(data['assets']),'gallery_figures':len(data['figures'])}))

if __name__=='__main__':main()
