"""Check actual native outputs, projections, original hashes and portable assets."""
from __future__ import annotations
import argparse
import hashlib
import json
import shutil
from pathlib import Path

import numpy as np
import pandas as pd
from PIL import Image

KEY = ['source', 'target', 'ligand_complex', 'receptor_complex']
SEVEN = KEY + ['lr_logfc', 'lr_gmeans', 'magnitude_rank']


def sha(path):
    with Path(path).open('rb') as f:
        return hashlib.file_digest(f, 'sha256').hexdigest()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--config', required=True)
    ap.add_argument('--results', required=True)
    args = ap.parse_args()
    cfg = json.loads(Path(args.config).read_text())
    results = Path(args.results)
    score = results / 'communication'
    data = pd.read_csv(score / 'tables/interaction_evidence.tsv.gz', sep='\t')
    assert not data.empty and not data.duplicated(['scope', 'condition'] + KEY).any()
    for col in ['magnitude_rank', 'lr_gmeans', 'lr_logfc']:
        assert np.isfinite(data[col]).all(), col
    assert data['magnitude_rank'].between(0, 1).all()
    assert data['lr_gmeans'].ge(0).all()
    assert data['lr_prop'].ge(cfg['expr_prop'] - 1e-7).all()
    assert (data['deg_supported'] == (data['ligand_deg_any'] | data['receptor_deg_any'])).all()
    assert not (data['ligand_deg_any'] & ~data['source_de_tested']).any()
    assert not (data['receptor_deg_any'] & ~data['target_de_tested']).any()
    assert np.allclose(data['lr_logfc'], (data['ligand_logfc'] + data['receptor_logfc']) / 2,
                       rtol=3e-5, atol=3e-6)
    pairs = pd.read_csv(score / 'tables/condition_pair_eligibility.tsv', sep='\t')
    for scope, n in [('IN_clusters', 80), ('IN_total', 4)]:
        for condition in ['pooled', 'WT', 'cKO']:
            p = pairs[(pairs.scope == scope) & (pairs.condition == condition)]
            assert len(p) == n and not p.duplicated(['source', 'target']).any()
            assert ((p.source.str.startswith('CA1_') & p.target.str.startswith('IN_')) |
                    (p.source.str.startswith('IN_') & p.target.str.startswith('CA1_'))).all()
            f = data[(data.scope == scope) & (data.condition == condition)]
            for direction in ['CA1_to_IN', 'IN_to_CA1', 'bidirectional']:
                expected = f if direction == 'bidirectional' else f[f.direction == direction]
                for supported in [False, True]:
                    ex = expected[expected.deg_supported] if supported else expected
                    suffix = '_DESeq2_supported_rank_aggregate.tsv' if supported else '_rank_aggregate.tsv'
                    actual = pd.read_csv(score / 'tables' / scope / condition / (direction + suffix), sep='\t')
                    assert actual.columns.tolist() == SEVEN
                    assert len(actual) == len(ex)
                    ex = ex[SEVEN].sort_values(KEY).reset_index(drop=True)
                    actual = actual.sort_values(KEY).reset_index(drop=True)
                    assert actual[KEY].equals(ex[KEY])
                    assert np.allclose(actual[SEVEN[4:]], ex[SEVEN[4:]], rtol=1e-10, atol=1e-12)
            gm = pd.read_csv(score / 'native_scores' / f'{scope}_{condition}_geometric_mean.tsv.gz', sep='\t')
            assert np.allclose(gm.lr_gmeans, np.sqrt(gm.ligand_means * gm.receptor_means), rtol=3e-5, atol=3e-6)
    for condition in ['WT', 'cKO']:
        assert (data['sample_support_' + condition] <= data['sample_eligible_' + condition]).all()
        assert data['sample_eligible_' + condition].between(0, 6).all()
    differential = pd.read_csv(score / 'tables/WT_cKO_descriptive_score_comparison.tsv.gz', sep='\t')
    assert differential['genotype_interaction_pvalue'].isna().all()
    assert differential['genotype_interaction_test_status'].eq('NOT_TESTED_DESCRIPTIVE_COMPARISON_ONLY').all()
    report = results / 'report'
    manifest = json.loads((report / 'data/plot_manifest.json').read_text())
    assert len(manifest['summary']) == 36
    image_count = 0
    for record in manifest['figures']:
        p = report / record['prefix']
        for extension in ['.pdf', '.png', '.svg', '_SourceData.tsv']:
            path = Path(str(p) + extension)
            assert path.is_file() and path.stat().st_size > 0, path
            if extension == '.pdf':
                assert path.open('rb').read(5) == b'%PDF-'
            if extension == '.png':
                with Image.open(path) as im:
                    im.verify()
        image_count += 1
    assert image_count >= 144
    original_hashes = {}
    for name, source in cfg['inputs'].items():
        actual = sha(source['path'])
        assert actual == source['sha256'], name
        original_hashes[name] = actual
    assert sha(cfg['resource_file']) == cfg['resource_sha256']
    for endpoint in cfg['endpoints']:
        for field in ['all_results_path', 'deg_path', 'stable_audit_path']:
            if field + '_sha256' in endpoint:
                assert sha(endpoint[field]) == endpoint[field + '_sha256']
    # Self-contained report includes all raw seven-column projections and
    # untouched source DESeq2 exports, not just the display-prioritized subset.
    shutil.copytree(score / 'tables', report / 'data/full_result_tables')
    shutil.copytree(score / 'deseq2_source_exports', report / 'data/existing_DESeq2_exports')
    shutil.copytree(results / 'prepared', report / 'data/input_audits',
                    ignore=shutil.ignore_patterns('*.h5ad'))
    source_root = Path(__file__).resolve().parents[1]
    shutil.copytree(source_root, report / 'analysis_source',
                    ignore=shutil.ignore_patterns('__pycache__', '*.pyc'))
    shutil.copy2(args.config, report / 'analysis_source/executed_config.json')
    html = report / 'index.html'
    extra = '<p><a href="data/full_result_tables/RESULT_COUNTS.tsv" download>Complete result counts</a> · <a href="data/full_result_tables/" >All native rank tables and mouse support</a> · <a href="data/existing_DESeq2_exports/">Original DESeq2 tables used</a> · <a href="analysis_source/README.md">Analysis source</a></p>'
    # Explicit file index works in file:// browsers that disallow directory views.
    links = []
    for folder in [report / 'data/full_result_tables', report / 'data/existing_DESeq2_exports']:
        entries = sorted(p for p in folder.rglob('*') if p.is_file())
        index = '<!doctype html><meta charset="utf-8"><title>Result downloads</title><h1>Result downloads</h1><ul>'
        from html import escape
        index += ''.join('<li><a download href="' + escape(p.relative_to(folder).as_posix(), quote=True) + '">' + escape(p.relative_to(folder).as_posix()) + '</a></li>' for p in entries)
        index += '</ul>'
        (folder / 'index.html').write_text(index)
    extra = extra.replace('href="data/full_result_tables/"', 'href="data/full_result_tables/index.html"').replace('href="data/existing_DESeq2_exports/"', 'href="data/existing_DESeq2_exports/index.html"')
    html.write_text(html.read_text().replace('</footer>', extra + '</footer>'))
    receipt = {'status': 'VALIDATED_REVIEW_REQUIRED', 'nuclei': cfg['expected_nuclei'],
               'interaction_rows': len(data), 'figures_checked': image_count,
               'report_sha256': sha(html), 'source_hashes_unchanged': original_hashes,
               'checks': ['score ranges', 'native score formulae', 'unique exact LR keys',
                          '84 requested pairs in each condition', 'all seven-column projections',
                          'endpoint-specific DE support', 'eligible mouse denominators',
                          'untested genotype interaction P remains NA', 'all plot assets open',
                          'prior source content unchanged']}
    (results / 'FINAL_VALIDATION.json').write_text(json.dumps(receipt, indent=2))
    print(json.dumps(receipt, indent=2))


if __name__ == '__main__':
    main()
