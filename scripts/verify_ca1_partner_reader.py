#!/usr/bin/env python3
"""Verify the additive CA1-side DEG / expressed IN partner reader."""
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SITE = ROOT / 'site'
READER = SITE / 'communication' / 'ca1_partner'

def main():
    receipt = json.loads((READER / 'READER_INTEGRITY.json').read_text())
    assert receipt['run_registry'] == '4276e1b6-bac1-4189-848f-9d161556ff4b'
    assert set(receipt['reader_files']) == {'index.html', 'style.css', 'app.js', 'model.js', 'data.js'}
    for name, expected in receipt['reader_files'].items():
        content = (READER / name).read_bytes()
        assert len(content) == expected['bytes'], name
        assert hashlib.sha256(content).hexdigest() == expected['sha256'], name
    prefix = 'window.SATB2_CA1_IN_DATA = '
    js = (READER / 'data.js').read_text()
    assert js.startswith(prefix) and js.endswith(';\n')
    data = json.loads(js[len(prefix):-2]); pairs = data['pairs']
    assert data['metadata']['state'] == 'reviewed_complete'
    assert (len(pairs),len(data['subunits']),len(data['mouse_expression'])) == (1364,2875,34500)
    truth = lambda v: str(v).lower() in ('true','1')
    assert sum(truth(r['passes_5pct']) for r in pairs) == 1364
    assert sum(truth(r['passes_10pct']) for r in pairs) == 937
    assert sum(truth(r['passes_5pct']) and truth(r['replicated_majority_5pct']) for r in pairs) == 1330
    assert sum(truth(r['passes_10pct']) and truth(r['replicated_majority_10pct']) for r in pairs) == 891
    ids = {r['record_id'] for r in pairs}
    assert len(ids) == len(pairs) and all(r['record_id'] in ids for r in data['subunits'] + data['mouse_expression'])
    assert any(r['padj'] == '' for r in data['subunits'])
    html = (READER / 'index.html').read_text()
    assert 'inhibitory cluster C8 has uncertain identity'.lower() in html.lower()
    assert 'interaction-level FDR' in html and 'Print / save plot as PDF' in html
    assert f'<script src="app.js?v={hashlib.sha256((READER / "app.js").read_bytes()).hexdigest()[:12]}"></script>' in html
    assert 'SATB2_CA1_LIANA_partner_portable_offline_20261030_v3.zip' in html
    assert 'SATB2_CA1_DE_IN_partner_reproducibility_20261030_v2.zip' in html
    assert 'communication/ca1_partner/index.html' in (SITE/'index.html').read_text()
    assert 'ca1_partner/index.html' in (SITE/'communication'/'index.html').read_text()
    print('CA1-side DEG / IN-expressed reader verified')

if __name__ == '__main__':
    main()
