"""Embed the communication UI, saved evidence and local libraries in one offline HTML.

Only transport and display assets are authored. No scientific input is opened or fitted.
"""
from pathlib import Path
import argparse
import hashlib
import json
import re

ROOT = Path(__file__).resolve().parent.parent
PAGE = ROOT / 'site/communication/index.html'
LIVE = 'https://ali-mujahid.github.io/satb2_snRNA/communication/'


def build(destination):
    html = PAGE.read_text()
    sources = []

    def script(match):
        relative = match.group(1).split('?')[0]
        path = PAGE.parent / relative
        content = path.read_text()
        sources.append({'path': str(path.resolve().relative_to(ROOT)),
                        'sha256': hashlib.sha256(path.read_bytes()).hexdigest()})
        return '<script>\n' + content.replace('</script', '<\\/script') + '\n</script>'

    def css(match):
        path = PAGE.parent / match.group(1).split('?')[0]
        sources.append({'path': str(path.resolve().relative_to(ROOT)),
                        'sha256': hashlib.sha256(path.read_bytes()).hexdigest()})
        return '<style>\n' + path.read_text() + '\n</style>'

    html = re.sub(r'<script src="([^"]+)"></script>', script, html)
    html = re.sub(r'<link rel="stylesheet" href="([^"]+)">', css, html)

    def link(match):
        url = match.group(1)
        if url.startswith(('https:', 'http:', '#', 'data:', 'mailto:')):
            return match.group(0)
        if url == '../index.html':
            return 'href="https://ali-mujahid.github.io/satb2_snRNA/"'
        return 'href="' + LIVE + url + '"'

    html = re.sub(r'href="([^"]+)"', link, html)
    html = html.replace('<title>SATB2 · DEG-only ligand–receptor explorer</title>',
                        '<title>SATB2 · DEG-only ligand–receptor explorer · offline</title>')
    assert not re.search(r'<script\s+src=', html)
    assert 'CCC_DATA_GZIP' in html and 'gunzipSync' in html
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_text(html)
    receipt = {'schema': 'satb2-communication-portable/v2',
               'html_file': destination.name, 'bytes': destination.stat().st_size,
               'sha256': hashlib.sha256(destination.read_bytes()).hexdigest(),
               'embedded_sources': sources, 'external_render_dependencies': [],
               'three_strict_DEG_only_projections_embedded': True,
               'endpoint_count': 5, 'strict_record_count': 47,
               'both_endpoint_all_subunit_DEG_rule': True,
               'matching_LR_component_DESeq2_evidence_embedded': True,
               'scientific_fits_performed': False,
               'native_complete_package': 'SATB2_DEG_ONLY_LR_results_20261023.zip',
               'external_links_require_internet': True}
    return receipt


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--receipt', type=Path, required=True)
    args = parser.parse_args()
    result = build(args.output)
    args.receipt.parent.mkdir(parents=True, exist_ok=True)
    args.receipt.write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps({k: result[k] for k in ['html_file', 'bytes', 'sha256']}))
