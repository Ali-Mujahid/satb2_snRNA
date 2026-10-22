"""Read-only checks of saved communication transport and publication links (no fits)."""
from pathlib import Path
import base64
import gzip
import hashlib
import json
import math
import re

ROOT = Path(__file__).resolve().parent.parent
SITE = ROOT / 'site'
CCC = SITE / 'communication'
SEVEN = ['source', 'target', 'ligand_complex', 'receptor_complex', 'lr_logfc', 'lr_gmeans', 'magnitude_rank']
sha = lambda b: hashlib.sha256(b).hexdigest()


def main():
    manifest = json.loads((CCC / 'data/manifest.json').read_text())
    raw = (CCC / 'data/data.js').read_bytes()
    transport = manifest['transport']
    assert len(raw) == transport['script_bytes'] and sha(raw) == transport['script_sha256']
    encoded = re.search(rb'window.CCC_DATA_GZIP\s*=\s*"([A-Za-z0-9+/=]+)"', raw).group(1)
    compressed = base64.b64decode(encoded, validate=True)
    assert sha(compressed) == transport['gzip_sha256']
    unpacked = gzip.decompress(compressed)
    assert sha(unpacked) == transport['json_sha256']
    data = json.loads(unpacked)
    assert len(data['columns']) == len(set(data['columns'])) == 48
    assert data['provenance']['genotype_interaction_test'] is False
    assert data['provenance']['no_new_scientific_fits'] is True
    expected = {'IN_clusters|pooled': 20471, 'IN_clusters|WT': 19467, 'IN_clusters|cKO': 21680,
                'IN_total|pooled': 1108, 'IN_total|WT': 1040, 'IN_total|cKO': 1170}
    all_rows = 0
    for key, n in expected.items():
        dataset = data['datasets'][key]
        assert len(dataset['rows']) == n
        assert len(dataset['eligible_pairs']) == (80 if key.startswith('IN_clusters|') else 4)
        keys = set()
        for values in dataset['rows']:
            assert len(values) == len(data['columns'])
            row = dict(zip(data['columns'], values))
            k = tuple(row[c] for c in SEVEN[:4])
            assert k not in keys
            keys.add(k)
            assert -1e-9 <= row['lr_prop'] <= 1 + 1e-9
            assert all(v is None or not isinstance(v, float) or math.isfinite(v) for v in values)
            assert data['components'][row['source']][row['ligand_complex']]
            assert data['components'][row['target']][row['receptor_complex']]
        all_rows += n
    assert all_rows == 64936
    receipt = json.loads((CCC / 'data/TRANSPORT_RECEIPT.json').read_text())
    for item in receipt['files']:
        path = (ROOT / item['path']).resolve()
        assert path.is_relative_to(ROOT) and path.is_file(), path
        packed = path.read_bytes()
        assert sha(packed) == item['transport_sha256']
        assert len(packed) == item['transport_bytes']
        original = gzip.decompress(packed) if item['compression'] == 'gzip' else packed
        if item['native_bytes_preserved']:
            assert sha(original) == item['source_sha256']
            assert len(original) == item['source_bytes']
        else:
            # The download landing page points to reversible gzip for the large JSON.
            assert item['path'] == 'communication_data/full_result_tables/index.html'
            assert item['adaptation'] == 'Download index JSON link uses lossless gzip transport'
        if path.name.endswith('_rank_aggregate.tsv'):
            assert original.decode().splitlines()[0].split('\t') == SEVEN
    for page in [CCC / 'index.html', CCC / 'RESULTS.html', *list((CCC / 'data').glob('*/index.html'))]:
        html = page.read_text()
        ids = re.findall(r'\bid="([^"]+)"', html)
        assert len(ids) == len(set(ids)), page
        for reference in re.findall(r'(?:src|href)="([^"]+)"', html):
            if reference.startswith(('https://', 'http://', '#', 'data:', 'mailto:')):
                continue
            relative = reference.split('?')[0].split('#')[0]
            assert (page.parent / relative).is_file(), (page, reference)
    page = (CCC / 'index.html').read_text()
    assert 'communication/index.html' in (SITE / 'index.html').read_text()
    for key in ['palette', 'lowColor', 'midColor', 'highColor', 'dotRadius', 'edgeWidth', 'arrowSize',
                'curvature', 'arcWidth', 'arcGap', 'ribbonScale', 'nodeScale', 'nodeColor', 'edgeColor',
                'arcColor', 'sources', 'targets', 'compare', 'dotMetric', 'heatMetric']:
        assert 'id="' + key + '"' in page, key
    assert 'vendor/fflate.min.js' in page and 'vendor/jspdf.umd.min.js' in page
    total = sum(p.stat().st_size for p in SITE.rglob('*') if p.is_file())
    assert total < 950_000_000, total
    print(json.dumps({'status': 'PASS', 'native_rows': all_rows, 'transport_files': len(receipt['files']),
                      'six_fits': True, 'source_hashes_preserved': True,
                      'pages_total_bytes': total, 'no_scientific_fits': True}, indent=2))


if __name__ == '__main__':
    main()
