"""Upload exact reviewed follow-up files to an existing release; never overwrite assets.

The proposal is a transport manifest, not scientific output. GitHub server sizes
and SHA256 digests must match every local file. Supports more than 100 assets.
"""
import argparse
import concurrent.futures
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import threading


def sha(path):
    h = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--proposal', type=Path, required=True)
    ap.add_argument('--stage', type=Path, required=True)
    ap.add_argument('--receipt', type=Path, required=True)
    ap.add_argument('--gh', required=True)
    ap.add_argument('--tag', required=True)
    ap.add_argument('--workers', type=int, default=3)
    args = ap.parse_args()
    repo = 'Ali-Mujahid/satb2_snRNA'
    env = {**os.environ, 'GH_PAGER': 'cat', 'NO_COLOR': '1'}

    def gh(*command):
        result = subprocess.run([args.gh, *command], capture_output=True, text=True, env=env)
        if result.returncode:
            raise RuntimeError(result.stderr[-1000:])
        return result.stdout

    def assets(rid):
        pages = json.loads(gh('api', f'repos/{repo}/releases/{rid}/assets?per_page=100', '--paginate', '--slurp'))
        rows = [row for page in pages for row in page]
        if len({x['name'] for x in rows}) != len(rows):
            raise ValueError('Duplicate remote asset names')
        return {row['name']: row for row in rows}

    proposal = json.loads(args.proposal.read_text())
    rows = proposal['assets']
    assert rows and len({r['asset_name'] for r in rows}) == len(rows)
    releases = json.loads(gh('api', f'repos/{repo}/releases?per_page=100', '--paginate', '--slurp'))
    release = next(r for page in releases for r in page if r['tag_name'] == args.tag)
    rid = release['id']
    old = assets(rid)
    args.stage.mkdir(parents=True, exist_ok=False)
    args.receipt.parent.mkdir(parents=True, exist_ok=True)
    assert not args.receipt.exists(), 'Preserve prior receipt; use a new path'
    for row in rows:
        src = Path(row['local_path']).resolve(strict=True)
        name = row['asset_name']
        if Path(name).name != name or name in ('.', '..'):
            raise ValueError('Unsafe asset name')
        assert src.is_file() and src.stat().st_size == row['bytes'] < 1_900_000_000
        assert sha(src) == row['sha256'], str(src)
        dst = args.stage / name
        shutil.copyfile(src, dst)
        assert sha(dst) == row['sha256']
        if name in old:
            assert old[name]['state'] == 'uploaded'
            assert old[name]['size'] == row['bytes']
            assert old[name].get('digest') == 'sha256:' + row['sha256'], name
    completed, errors = [], []
    lock = threading.Lock()

    def upload(row):
        name = row['asset_name']
        if name not in old:
            gh('release', 'upload', args.tag, str(args.stage / name), '--repo', repo)
        actual = assets(rid)[name]
        assert actual['state'] == 'uploaded' and actual['size'] == row['bytes']
        assert actual.get('digest') == 'sha256:' + row['sha256'], name
        return {**row, 'remote_verified_bytes': actual['size'],
                'remote_verified_sha256': actual['digest'].removeprefix('sha256:'),
                'verification_method': 'github_server_digest', 'github_asset_id': actual['id'],
                'observed_browser_download_url': actual['browser_download_url']}

    with concurrent.futures.ThreadPoolExecutor(max_workers=max(1, min(args.workers, 5))) as pool:
        jobs = {pool.submit(upload, row): row for row in rows}
        for future in concurrent.futures.as_completed(jobs):
            try:
                completed.append(future.result())
            except Exception as error:
                errors.append({'asset_name': jobs[future]['asset_name'], 'error': str(error)})
            with lock:
                args.receipt.write_text(json.dumps({
                    'status': 'VERIFIED' if len(completed) == len(rows) else 'INCOMPLETE',
                    'repo': repo, 'release_id': rid, 'tag': args.tag,
                    'proposal_sha256': sha(args.proposal),
                    'expected_assets': len(rows), 'assets': sorted(completed, key=lambda r: r['asset_name']),
                    'errors': errors, 'no_overwrite': True,
                }, indent=2) + '\n')
            print(f'Verified {len(completed)}/{len(rows)} assets; errors {len(errors)}', flush=True)
    if errors:
        raise SystemExit('Upload incomplete; preserve receipt and inspect errors before continuing.')


if __name__ == '__main__':
    main()
