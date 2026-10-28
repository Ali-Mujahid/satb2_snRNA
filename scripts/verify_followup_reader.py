"""Verify the additive online follow-up reader without executing scientific code."""
from pathlib import Path
from urllib.parse import unquote, urlparse
import hashlib
import json
import re

ROOT = Path(__file__).resolve().parent.parent
SITE = ROOT / 'site'
PAGE = SITE / 'followups'


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def local(relative):
    path = (PAGE / unquote(relative)).resolve(strict=True)
    assert path.is_relative_to(PAGE.resolve()) and path.is_file(), relative
    assert not (PAGE / unquote(relative)).is_symlink(), relative
    return path


def release_url(value):
    url = urlparse(value)
    assert url.scheme == 'https' and url.netloc == 'github.com'
    assert url.path.startswith('/Ali-Mujahid/satb2_snRNA/releases/download/followups-20261028-v1.0.0/')
    assert not url.query and not url.fragment and not url.username and not url.password


def main():
    manifest = json.loads((PAGE / 'ONLINE_TRANSPORT_MANIFEST.json').read_text())
    assert manifest['schema'] == 'satb2-followups-online-transport-manifest/v1'
    assert manifest['scientific_values_or_native_outputs_changed'] is False
    expected = set()
    for row in manifest['files']:
        p = local(row['path'])
        assert row['path'] not in expected
        expected.add(row['path'])
        assert p.stat().st_size == row['bytes'] and sha(p) == row['sha256'], row['path']
    actual = {p.relative_to(PAGE).as_posix() for p in PAGE.rglob('*') if p.is_file()}
    assert actual == expected | {'ONLINE_TRANSPORT_MANIFEST.json'}, 'Unsealed/missing online files'
    remote = {}
    for row in manifest['verified_release_assets']:
        release_url(row['url'])
        assert row['url'] not in remote
        assert row['remote_verified_bytes'] == row['bytes'] > 0
        assert row['remote_verified_sha256'] == row['sha256']
        assert re.fullmatch(r'[0-9a-f]{64}', row['sha256'])
        assert row['verification_method'] in {'github_server_digest', 'download_sha256'}
        remote[row['url']] = row
    nav = json.loads(re.search(r'const INDEX=(.*?);\nconst el=', (PAGE / 'index.html').read_text(), re.S)[1])
    assert [t['id'] for t in nav['tasks']] == ['composition', 'cellchat', 'monocle3', 'velocity']
    cfg = json.loads((PAGE / 'APPROVED_PRESENTATION_CONFIGURATION.json').read_text())
    assert cfg['approved_for_portable_presentation'] is True
    for task in nav['tasks']:
        assert task['headline'] == cfg['tasks'][task['id']]['headline']
        assert task['findings'] == cfg['tasks'][task['id']]['findings']
        assert task['limitations'] == cfg['tasks'][task['id']]['limitations']
        assert task['findings'] and task['limitations'] and task['reviews']
        seal = json.loads(local(task['seal_href']).read_text())
        assert task['native_schema'] == seal['schema']
        if task['id'] != 'velocity':
            assert task['native_status'] == seal['status']
        for review in task['reviews']:
            assert sha(local(review['href'])) == review['sha256']
        for asset in task['assets']:
            if asset['href'].startswith('https:'):
                r = remote[asset['href']]
                assert r['bytes'] == asset['bytes'] and r['sha256'] == asset['sha256']
                assert unquote(asset['portable_href']) == r['portable_path']
            else:
                p = local(asset['href'])
                assert p.stat().st_size == asset['bytes'] and sha(p) == asset['sha256']
        if task['id'] == 'cellchat':
            assert task['default_rendering'] == 'Readable layout'
            assert 'vertical jitter' in task['original_layout_warning']
        if task['id'] == 'velocity':
            assert task['native_status'] == 'INPUT_BLOCKED' and task['assets'] == []
            assert task['native_summary']['scientific_jobs_started'] is False
    ziprow = manifest['portable_zip']
    release_url(ziprow['url'])
    assert ziprow['bytes'] == ziprow['remote_verified_bytes']
    assert ziprow['sha256'] == ziprow['remote_verified_sha256']
    assert nav['delivery']['portable_zip_sha256'] == ziprow['sha256']
    assert nav['delivery']['portable_zip_url'] == ziprow['url']
    assert nav['delivery']['portable_zip_bytes'] == ziprow['bytes']
    page_size = sum(p.stat().st_size for p in PAGE.rglob('*') if p.is_file())
    site_size = sum(p.stat().st_size for p in SITE.rglob('*') if p.is_file())
    assert page_size <= manifest['maximum_online_bytes'] <= 10_000_000
    assert site_size <= manifest['maximum_site_bytes'] <= 950_000_000
    assert 'followups/index.html' in (SITE / 'index.html').read_text()
    print(json.dumps({'followups': 'PASS', 'online_files': len(expected),
                      'verified_remote_assets': len(remote), 'page_bytes': page_size, 'site_bytes': site_size}))


if __name__ == '__main__':
    main()
