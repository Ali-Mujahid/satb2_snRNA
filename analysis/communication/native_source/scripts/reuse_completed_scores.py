"""Reuse hash-verified completed scores for rendering; never execute a fit."""
import argparse
import hashlib
import json
import shutil
from pathlib import Path

def sha(path):
    with Path(path).open('rb') as f:
        return hashlib.file_digest(f, 'sha256').hexdigest()

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--config', required=True)
    args = ap.parse_args()
    c = json.loads(Path(args.config).read_text())
    source = Path(c['completed_results_root'])
    out = Path(c['run_root']) / 'results'
    for item in c['completed_artifacts']:
        p = Path(item['path'])
        assert p.resolve().is_relative_to(source.resolve())
        assert p.stat().st_size == item['bytes'] and sha(p) == item['sha256'], p
    assert (source / 'communication/COMPLETE.json').is_file()
    for folder in ['communication', 'prepared', 'report']:
        assert not (out / folder).exists() or not any((out / folder).iterdir()), 'Refusing nonempty previous result tree'
        shutil.copytree(source / folder, out / folder, dirs_exist_ok=True,
                        ignore=shutil.ignore_patterns('*.h5ad') if folder == 'prepared' else None)
    record = {'status': 'HASH_VERIFIED_NATIVE_SCORE_REUSE_NO_REFITS',
              'source_results_root': str(source), 'verified_source_files': len(c['completed_artifacts']),
              'model_fits_executed': 0, 'native_scores_altered': False,
              'rendering_fix': 'Tolerate detection-fraction roundoff within 1e-9; clamp plotted fractions only'}
    (out / 'communication/REUSE_RECEIPT.json').write_text(json.dumps(record, indent=2))
    for item in c['completed_artifacts']:
        assert sha(item['path']) == item['sha256']
    print(json.dumps(record, indent=2))

if __name__ == '__main__':
    main()
