#!/usr/bin/env python3
"""Stage a small online transport of a reviewed portable reader.

This does not upload, deploy, infer, reformat tables or change native results.
Release-hosted file bytes must already have been verified independently.  The
complete portable package remains the self-contained offline deliverable.
"""
from __future__ import annotations

import argparse
from copy import deepcopy
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import re
import shutil
import stat
import zipfile
from urllib.parse import quote, unquote, urlparse

from build_reader import HERE, receipt, sha, verify, protected_result_roots, reject_protected_destination


def relative_file(value):
    """A URL-encoded portable relative path, never a URL or filesystem escape."""
    if not isinstance(value, str) or urlparse(value).scheme or value.startswith('/'):
        raise ValueError('Expected a portable relative file path')
    decoded = unquote(value)
    path = Path(decoded)
    if not decoded or '\\' in decoded or path.is_absolute() or any(x in ('', '.', '..') for x in decoded.split('/')):
        raise ValueError('Unsafe portable relative path')
    return path.as_posix()


def read_navigation(root):
    match = re.search(r'const INDEX=(.*?);\nconst el=', (root / 'index.html').read_text(), flags=re.S)
    if not match:
        raise ValueError('Portable navigation payload is absent')
    return json.loads(match.group(1))


def portable_package(spec):
    root = Path(spec['portable_root']).resolve(strict=True)
    manifest_path = verify(spec['portable_manifest'])
    if manifest_path != root / 'READER_MANIFEST.json':
        raise ValueError('Portable manifest identity differs from the supplied root')
    manifest = json.loads(manifest_path.read_text())
    if manifest.get('schema') != 'satb2-followups-portable-reader-manifest/v1' or manifest.get('status') not in {
        'BUILT_PENDING_RENDERED_REVIEW', 'BUILT_REVIEWED'
    } or manifest.get('no_fit_rescore_or_native_mutation') is not True:
        raise ValueError('A preserved native portable reader is required')
    files = {}
    for record in manifest.get('files', []):
        path = verify(record, root)
        relative = path.relative_to(root).as_posix()
        if relative in files:
            raise ValueError('Duplicate portable file receipt')
        files[relative] = path
    if 'index.html' not in files:
        raise ValueError('Portable index is not sealed')
    files['READER_MANIFEST.json'] = manifest_path
    navigation = read_navigation(root)
    refs = {}
    for task in navigation['tasks']:
        for asset in task['assets']:
            relative = relative_file(asset['href'])
            if relative not in files or files[relative].stat().st_size != asset['bytes'] or sha(files[relative]) != asset['sha256']:
                raise ValueError('Portable navigation and sealed asset identity disagree')
            refs[relative] = asset
    return root, manifest_path, manifest, files, navigation, refs


def safe_release_url(url):
    parsed = urlparse(url)
    prefix = '/Ali-Mujahid/satb2_snRNA/releases/download/'
    if parsed.scheme != 'https' or parsed.netloc != 'github.com' or parsed.username or parsed.password or not parsed.path.startswith(prefix) or parsed.query or parsed.fragment:
        raise ValueError('Use the authorized repository\'s canonical HTTPS release download URL')
    tail = parsed.path[len(prefix):].split('/')
    if len(tail) != 2 or not all(tail) or any(unquote(x) in ('.', '..') or '/' in unquote(x) or '\\' in unquote(x) for x in tail):
        raise ValueError('Release download URL must identify one tag and one asset')
    return url


def verified_remote(record, local_path):
    safe_release_url(record['url'])
    digest = sha(local_path)
    size = local_path.stat().st_size
    if record.get('bytes') != size or record.get('sha256') != digest:
        raise ValueError('Release transport differs from the sealed portable file')
    if record.get('remote_verified_bytes') != size or record.get('remote_verified_sha256') != digest:
        raise ValueError('Server/download verification is missing or does not match local bytes')
    if record.get('verification_method') not in {'github_server_digest', 'download_sha256'}:
        raise ValueError('Remote verification must use server SHA256 or a downloaded SHA256')
    return record['url']


def verify_portable_zip(path, files, prefix=''):
    """Check complete offline contents without extracting or changing a file."""
    if prefix:
        prefix = relative_file(prefix).rstrip('/') + '/'
    seen = set()
    with zipfile.ZipFile(path) as archive:
        for entry in archive.infolist():
            if stat.S_ISLNK(entry.external_attr >> 16):
                raise ValueError('Portable ZIP cannot contain symbolic-link entries')
            name = relative_file(entry.filename.rstrip('/') if entry.is_dir() else entry.filename)
            if entry.is_dir() and prefix and name == prefix.rstrip('/'):
                continue
            if not name.startswith(prefix):
                raise ValueError('Portable ZIP prefix does not match its approved package')
            if entry.is_dir():
                continue
            relative = name[len(prefix):]
            if relative in seen or relative not in files:
                raise ValueError('Portable ZIP contains duplicate or unmanifested files')
            seen.add(relative)
            if entry.file_size != files[relative].stat().st_size:
                raise ValueError('Portable ZIP file size differs from the reviewed package')
            digest = hashlib.sha256()
            with archive.open(entry) as stream:
                for chunk in iter(lambda: stream.read(1024 * 1024), b''):
                    digest.update(chunk)
            if digest.hexdigest() != sha(files[relative]):
                raise ValueError('Portable ZIP file bytes differ from the reviewed package')
    if seen != set(files):
        raise ValueError('Portable ZIP omits reviewed package files')


def approved_config(path):
    cfg = json.loads(Path(path).read_text())
    if cfg.get('schema') != 'satb2-followups-online-transport/v1' or cfg.get('approved_for_online_presentation') is not True:
        raise ValueError('Root-approved online transport configuration is required')
    for key in ('portable_visual_review_receipts', 'remote_verification_receipts'):
        if not cfg.get(key):
            raise ValueError(f'{key} is required before online staging')
        for record in cfg[key]:
            verify(record)
    return cfg


def plan_release_assets(portable_root, tag, maximum_local_asset_bytes=20000):
    """A read-only upload proposal; never a verification or approval claim."""
    root = Path(portable_root).resolve(strict=True)
    spec = {'portable_root': str(root), 'portable_manifest': receipt(root / 'READER_MANIFEST.json')}
    _, manifest_path, _, files, _, refs = portable_package(spec)
    rows = []
    for relative, asset in sorted(refs.items()):
        path = files[relative]
        if asset['extension'] not in ('.png', '.pdf', '.gz') and path.stat().st_size <= maximum_local_asset_bytes:
            continue
        name = 'followups-' + hashlib.sha256(relative.encode()).hexdigest()[:12] + '-' + re.sub(r'[^A-Za-z0-9._-]', '_', path.name)
        url = f'https://github.com/Ali-Mujahid/satb2_snRNA/releases/download/{quote(tag, safe="")}/{name}'
        safe_release_url(url)
        rows.append({'portable_path': relative, 'local_path': str(path), 'asset_name': name, 'url': url,
                     'bytes': path.stat().st_size, 'sha256': sha(path), 'remote_verified_bytes': None,
                     'remote_verified_sha256': None, 'verification_method': None})
    return {'schema': 'satb2-followups-online-asset-proposal/v1', 'status': 'PROPOSAL_NOT_UPLOADED_NOT_VERIFIED',
            'portable_manifest': receipt(manifest_path), 'release_tag': tag, 'assets': rows,
            'remote_asset_bytes': sum(x['bytes'] for x in rows),
            'note': 'Upload exact bytes, verify GitHub bytes/digest or downloaded SHA256, and save those receipts before approving online staging.'}


def build(configuration, destination):
    cfg = approved_config(configuration)
    root, portable_manifest_path, manifest, files, original_nav, refs = portable_package(cfg)
    destination = Path(destination).resolve()
    source_cfg=json.loads(files['APPROVED_PRESENTATION_CONFIGURATION.json'].read_text())
    reject_protected_destination(destination,protected_result_roots(source_cfg))
    site = Path(cfg['existing_site_root']).resolve(strict=True)
    if destination.exists() or destination.is_relative_to(root) or root.is_relative_to(destination) or destination.is_relative_to(site) or site.is_relative_to(destination):
        raise ValueError('Online stage must be fresh and separate from portable package and existing site')
    remote = {}
    seen_urls = set()
    for record in cfg.get('verified_release_assets', []):
        relative = relative_file(record['portable_path'])
        if relative in remote or relative not in refs:
            raise ValueError('Duplicate or unused release transport mapping')
        url = verified_remote(record, files[relative])
        if url in seen_urls:
            raise ValueError('One release URL cannot represent different portable assets')
        seen_urls.add(url)
        remote[relative] = url
    if not remote:
        raise ValueError('An online transport requires verified external assets')
    zip_record = cfg['portable_zip']
    zip_path = verify(zip_record['local_receipt'])
    verify_portable_zip(zip_path, files, zip_record.get('archive_prefix', ''))
    zip_url = verified_remote(zip_record, zip_path)
    nav = deepcopy(original_nav)
    for task in nav['tasks']:
        for asset in task['assets']:
            original = relative_file(asset['href'])
            if original in remote:
                asset['portable_href'] = asset['href']
                asset['href'] = remote[original]
                asset['delivery'] = 'Verified GitHub release asset; same bytes as portable package'
            if asset.get('original_href'):
                relative = relative_file(asset['original_href'])
                if relative in remote:
                    asset['portable_original_href'] = asset['original_href']
                    asset['original_href'] = remote[relative]
    nav['delivery'] = {'mode': 'online', 'portable_zip_url': zip_url, 'portable_zip_bytes': zip_path.stat().st_size,
                       'portable_zip_sha256': sha(zip_path),
                       'note': 'This small online preview loads some original plots, PDFs and tables from verified GitHub release assets. Download and unpack the complete portable ZIP for offline operation.'}
    payload = json.dumps(nav, ensure_ascii=False, separators=(',', ':')).replace('<', '\\u003c')
    html = (HERE / 'template.html').read_text().replace('/*__REVIEWED_NATIVE_NAVIGATION__*/', payload)
    copies = [(path, Path(relative)) for relative, path in files.items() if relative not in remote and relative not in {'index.html', 'README.txt'}]
    planned = sum(p.stat().st_size for p, _ in copies) + len(html.encode())
    limit = int(cfg.get('maximum_online_bytes', 10_000_000))
    existing_bytes = sum(p.stat().st_size for p in site.rglob('*') if p.is_file())
    site_limit = int(cfg.get('maximum_site_bytes', 950_000_000))
    if planned > limit or existing_bytes + planned > site_limit:
        raise ValueError('Planned online transport exceeds explicit online or site limit; no files staged')
    destination.mkdir(parents=True)
    outputs = []
    for path, relative in copies:
        target = destination / relative
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(path, target)
        if sha(path) != sha(target):
            raise ValueError('Online small-file copy mismatch')
        outputs.append({'path': relative.as_posix(), 'bytes': target.stat().st_size, 'sha256': sha(target)})
    (destination / 'index.html').write_text(html)
    (destination / 'README.txt').write_text('Online preview: some plot/PDF/table bytes are delivered from verified GitHub release assets. Use the complete portable ZIP linked in the page for offline operation. No scientific values, filters or statistical results changed during online transport.\n')
    for name in ('index.html', 'README.txt'):
        p = destination / name
        outputs.append({'path': name, 'bytes': p.stat().st_size, 'sha256': sha(p)})
    result = {'schema': 'satb2-followups-online-transport-manifest/v1', 'status': 'STAGED_PENDING_LIVE_REVIEW',
              'created_utc': datetime.now(timezone.utc).isoformat(), 'source_portable_manifest': receipt(portable_manifest_path),
              'configuration': receipt(configuration), 'source_builder': receipt(__file__), 'files': outputs,
              'verified_release_assets': cfg['verified_release_assets'], 'portable_zip': zip_record,
              'scientific_values_or_native_outputs_changed': False, 'existing_site_bytes': existing_bytes,
              'maximum_site_bytes': site_limit, 'maximum_online_bytes': limit}
    (destination / 'ONLINE_TRANSPORT_MANIFEST.json').write_text(json.dumps(result, indent=2) + '\n')
    size = sum(p.stat().st_size for p in destination.rglob('*') if p.is_file())
    if size > limit or existing_bytes + size > site_limit:
        raise ValueError('Final online receipt overhead exceeds budget; preserve draft and review, do not deploy')
    return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    sub = parser.add_subparsers(dest='command', required=True)
    proposal = sub.add_parser('plan-assets')
    proposal.add_argument('--portable-root', required=True, type=Path)
    proposal.add_argument('--tag', required=True)
    proposal.add_argument('--maximum-local-asset-bytes', type=int, default=20000)
    stage = sub.add_parser('stage')
    stage.add_argument('--configuration', required=True, type=Path)
    stage.add_argument('--destination', required=True, type=Path)
    args = parser.parse_args()
    if args.command == 'plan-assets':
        print(json.dumps(plan_release_assets(args.portable_root, args.tag, args.maximum_local_asset_bytes), indent=2))
    else:
        result = build(args.configuration, args.destination)
        print(json.dumps({'status': result['status'], 'destination': str(args.destination.resolve()), 'files': len(result['files'])}))
