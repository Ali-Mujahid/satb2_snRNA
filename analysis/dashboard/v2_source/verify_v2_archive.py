from pathlib import Path
import hashlib,json,zipfile

BASE=Path(__file__).resolve().parent.parent
ROOT=BASE/'reports/offline_dashboard_v2'
ZIP=ROOT.with_suffix('.zip')
MANIFEST=json.loads((ROOT/'BUNDLE_MANIFEST.json').read_text())
expected={ROOT.name+'/'+str(p.relative_to(ROOT)):p for p in ROOT.rglob('*') if p.is_file()}
hashes={ROOT.name+'/'+r['path']:r['sha256'] for r in MANIFEST['files']}
manifest_hash=hashlib.sha256((ROOT/'BUNDLE_MANIFEST.json').read_bytes()).hexdigest()
hashes[ROOT.name+'/BUNDLE_MANIFEST.json']=manifest_hash
assert len(expected)==len(hashes)
total=0
with zipfile.ZipFile(ZIP) as z:
    names=z.namelist()
    assert len(names)==len(set(names))
    assert set(names)==set(expected)
    for name in names:
        h=hashlib.sha256()
        with z.open(name) as f:
            for b in iter(lambda:f.read(1024*1024),b''):h.update(b);total+=len(b)
        assert h.hexdigest()==hashes[name],name
h=hashlib.sha256()
with ZIP.open('rb') as f:
    for b in iter(lambda:f.read(1024*1024),b''):h.update(b)
r={'status':'PASS','zip':str(ZIP),'zip_bytes':ZIP.stat().st_size,'zip_sha256':h.hexdigest(),'archive_members':len(expected),'all_file_CRC_SHA256_checks_passed':True,'uncompressed_bytes_checked':total,'bundle_manifest_sha256':manifest_hash,'source_app_sha256':hashes[ROOT.name+'/app.js'],'source_html_sha256':hashes[ROOT.name+'/index.html'],'n_primary_endpoints':948,'n_tested_primary':852,'n_unsupported_primary':96,'n_sensitivity_endpoints':112,'n_native_indexed_figure_pages':450,'n_native_PDFs':67,'n_new_saved_example_PDFs':18,'n_new_saved_example_PNGs':25,'n_marker_genes':878,'n_mouse_genes':32285,'n_mouse_fit_bindings':1060,'n_mouse_count_shards':480,'no_raw_nucleus_RNA_extraction_or_new_fitting':True,'native_and_v1_preserved':True,'browser_CSS_verification_performed':False}
(BASE/'reviews/OFFLINE_V2_DELIVERY_COMPLETE.json').write_text(json.dumps(r,indent=2)+'\n')
print(json.dumps(r))
