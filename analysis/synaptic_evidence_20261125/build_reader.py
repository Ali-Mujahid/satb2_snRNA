"""Build the offline and checksum-verified Pages reader from the same sources."""
import base64
import hashlib
import json
from pair_sets import compare_pairs
from binding_audit import annotate_pairs
import zipfile
from pathlib import Path

root = Path(__file__).resolve().parent
site = root.parents[1] / 'site'
source = root / 'reader'
release = 'atlas-v3.14.2-binding-audit-20261126'
remote = f'https://raw.githubusercontent.com/Ali-Mujahid/satb2_snRNA/{release}/analysis/{root.name}/offline.html'

def inline_script(path):
    code = path.read_text().replace('</script', r'<\/script')
    return '<script>\n' + code + '\n</script>'

page = (source / 'index.html').read_text()
page = page.replace('<link rel="stylesheet" href="evidence.css">', '<style>' + (source / 'evidence.css').read_text() + '</style>')
page = page.replace('<script src="../vendor/jspdf.umd.min.js"></script>', inline_script(site / 'vendor/jspdf.umd.min.js'))
payload = json.loads((root / 'results/data.json').read_text())
payload['pairs'], payload['bindingAudit'] = annotate_pairs(payload['pairs'], root / 'results')
payload['pairRoles'] = compare_pairs(payload, root / 'results')
data = json.dumps(payload, separators=(',', ':'), allow_nan=False).replace('<', r'\u003c')
page = page.replace('<script src="pair_views.js"></script>', inline_script(source / 'pair_views.js'))
page = page.replace('<script src="binding_views.js"></script>', inline_script(source / 'binding_views.js'))
page = page.replace('<script src="evidence.js"></script>', '<script type="application/json" id="sourceData">' + data + '</script>' + inline_script(source / 'evidence.js'))
page = page.replace('ATLAS_RELEASE_TAG', release)
(root / 'offline.html').write_text(page)
with zipfile.ZipFile(root / 'synaptic_evidence.zip', 'w', compression=zipfile.ZIP_DEFLATED) as archive:
    for file in [root / 'offline.html', *sorted((root / 'results').iterdir())]:
        archive.write(file, file.relative_to(root))
integrity = 'sha256-' + base64.b64encode(hashlib.sha256(page.encode()).digest()).decode()
shell = """<!doctype html><html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>Synaptic evidence · SATB2 atlas</title></head><body>
<main style="font:16px system-ui;max-width:68ch;margin:4rem auto;padding:0 1rem"><a href="../index.html">← SATB2 atlas</a><h1>Synaptic evidence in CA1</h1><p id="load" role="status">Loading the four reference analyses…</p><a href="REMOTE" download="SATB2_synaptic_evidence.html">Download the standalone explorer</a></main>
<script>
const url=location.pathname.includes('/site/synaptic_evidence/')?'../../analysis/ANALYSIS/offline.html':'REMOTE';
fetch(url,{integrity:'DIGEST',credentials:'omit'}).then(response=>{
  if(!response.ok)throw Error('HTTP '+response.status);
  return response.text();
}).then(html=>{document.open();document.write(html);document.close();}).catch(error=>{document.getElementById('load').textContent='Could not load the explorer. '+error.message+'. Download the standalone copy above to open it locally.';});
</script></body></html>
"""
(site / 'synaptic_evidence/index.html').write_text(shell.replace('REMOTE',remote).replace('ANALYSIS',root.name).replace('DIGEST',integrity))
