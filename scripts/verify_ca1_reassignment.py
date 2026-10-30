import hashlib,json,re
from pathlib import Path
p=Path('site/ca1_sublayers');prov=json.loads((p/'PROVENANCE.json').read_text());h=lambda f:hashlib.sha256((p/f).read_bytes()).hexdigest()
assert h('data.js')==prov['data_sha256'] and h('app.js')==prov['app_sha256']
html=(p/'index.html').read_text()
for name in ['data.js','app.js']:assert f'src="{name}?v={h(name)[:12]}"' in html
assert not any(t in html for t in ['RESULTS_ARCHIVE_URL','OFFLINE_ARCHIVE_URL','JSPDF_URL','REVIEW_URL'])
s=(p/'data.js').read_text();assert s.startswith('window.CA1_REVISION=') and s.strip().endswith(';');d=json.loads(s[len('window.CA1_REVISION='):].strip()[:-1]);assert len(d['coords'])==27703 and len(d['assignments'])==5 and all(x['strict70_assignment']=='Unassigned' for x in d['assignments'])
assert d['summary']['Deep']['nuclei_augmented']==11306 and d['summary']['Superficial']['nuclei_augmented']==16298
assert [d['summary'][k]['baseline'].get('up',0)+d['summary'][k]['baseline'].get('down',0) for k in ['Deep','Superficial']]==[742,1825]
assert 'ca1_sublayers/index.html' in Path('site/index.html').read_text()
print('PASS CA1 reassignment source-bound page, program counts and browser cache fingerprints')
