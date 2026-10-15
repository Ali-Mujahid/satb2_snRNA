#!/usr/bin/env python3
"""Independent checks of viewer shards against existing produced tables only."""
from pathlib import Path
import argparse, base64, csv, gzip, hashlib, json, math, re
from collections import Counter

BASE = Path(__file__).resolve().parents[1]

def rows(path):
    with (gzip.open(path, 'rt') if str(path).endswith('.gz') else Path(path).open()) as f:
        yield from csv.DictReader(f, delimiter='\t')

def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for b in iter(lambda:f.read(1024*1024), b''): h.update(b)
    return h.hexdigest()

def number(v):
    if v is None or str(v) in ('', 'NA', 'NaN', 'nan', 'Inf', '-Inf'): return None
    n = float(v)
    return n if math.isfinite(n) else None

def close(a,b):
    if a is None or b is None: return a is b
    return math.isclose(float(a),float(b),rel_tol=2e-12,abs_tol=2e-10)

def load_shard(out, manifest, key):
    meta = manifest['shards'][key]
    path = out / meta['path']
    if sha(path) != meta['sha256']: raise AssertionError('Shard hash: '+key)
    text = path.read_text()
    match = re.search(r'data:("[A-Za-z0-9+/=]+")', text)
    if not match: raise AssertionError('Unsupported shard encoding: '+key)
    packed = base64.b64decode(json.loads(match[1]), validate=True)
    raw = gzip.decompress(packed)
    if len(raw) != meta['uncompressed_bytes'] or len(packed) != meta['gzip_bytes']:
        raise AssertionError('Shard length: '+key)
    return json.loads(raw)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--outdir', type=Path, default=BASE/'reports/offline_dashboard_v1')
    ap.add_argument('--config', type=Path, default=BASE/'config/analysis_v2_0_0.json')
    args = ap.parse_args(); out=args.outdir
    cfg=json.loads(args.config.read_text()); root=Path(cfg['results_root'])
    manifest=json.loads((out/'data/manifest.json').read_text())
    assert manifest['schema_version']=='satb2-offline-1'
    assert len(manifest['datasets'])==24 and len(manifest['genes'])==878
    assert len({g['gene_id'] for g in manifest['genes']})==878
    roster={x['sample']:x['condition'] for x in manifest['samples']}
    assert len(roster)==12 and Counter(roster.values())=={'WT':6,'cKO':6}
    all_genes={x['gene_id'] for x in manifest['genes']}
    for panel in manifest['gene_sets']: assert set(panel['gene_ids'])<=all_genes
    reviewed=next(x for x in manifest['gene_sets'] if x['id']=='provided_all')
    given=next(x for x in manifest['gene_sets'] if x['id']=='workbook_only')
    assert len(reviewed['gene_ids'])==139 and len(given['gene_ids'])==117
    acceptance=[]; paired_members={}; tool_members={}
    for dataset in manifest['datasets']:
        key=dataset['id']; tool,branch,pop=dataset['tool'],dataset['branch'],dataset['population']
        data=load_shard(out,manifest,dataset['coordinates_shard'])
        nuclei=data['nuclei']; ids=[x['nucleus_id'] for x in nuclei]
        assert ids==sorted(ids) and len(set(ids))==len(ids)==dataset['n_nuclei']
        assert {x['sample'] for x in nuclei}==set(roster)
        assert all(roster[x['sample']]==x['condition'] for x in nuclei)
        assert all(math.isfinite(x[k]) for x in nuclei for k in ['x_harmony','y_harmony','x_unintegrated','y_unintegrated'])
        clusters=Counter(x['cluster_id'] for x in nuclei)
        assert all(clusters.get(c['cluster_id'],0)==int(c['n_nuclei']) for c in data['clusters'])
        source=root/'plots'/tool/branch/pop/'original_frozen/rendered_nucleus_metadata.tsv.gz'
        native={x['nucleus_id']:x for x in rows(source) if x['keep_for_analysis'].lower()=='true'}
        assert set(ids)==set(native)
        labels=['sample','condition','cluster_id','reviewed_label','broad_label','label_confidence',
                'layer_program_label','expression_program_label','physical_layer_assignment']
        for x in nuclei:
            y=native[x['nucleus_id']]
            for k in labels:
                expected=None if y[k] in ('','NA') else y[k]
                assert x[k]==expected,(key,k,x[k],expected)
            for a,b in [('x_harmony','UMAP_harmony_1'),('y_harmony','UMAP_harmony_2'),
                        ('x_unintegrated','UMAP_unintegrated_1'),('y_unintegrated','UMAP_unintegrated_2')]:
                # Plot metadata retains native tool coordinate column naming.
                alt=b.replace('UMAP_harmony_','UMAPH_').replace('UMAP_unintegrated_','UMAPU_')
                value=y.get(b,y.get(alt))
                if value is not None: assert close(x[a],number(value)),(key,a)
        frozen=frozenset(ids)
        if (branch,pop) in tool_members: assert tool_members[(branch,pop)]==frozen
        else: tool_members[(branch,pop)]=frozen
        if branch.startswith('paired_'):
            if (tool,pop) in paired_members: assert paired_members[(tool,pop)]==frozen
            else: paired_members[(tool,pop)]=frozen
        assert (out/dataset['pdf_path']).is_file()
        markers=load_shard(out,manifest,dataset['markers_shard'])
        maps={scope:{(r['cluster_id'],r['gene_id']):r for r in markers[field]}
              for scope,field in [('pooled','rows'),('WT','WT_rows'),('cKO','cKO_rows')]}
        assert all(len(m)==len(data['clusters'])*878 for m in maps.values())
        numeric=['n_nuclei','sum_counts','n_detected','sum_log1p_normalized','mean_log1p_normalized',
                 'fraction_detected','mean_source_counts']
        checks=0
        for name,scope in [('gene_cluster_summary.tsv','pooled'),('gene_cluster_summary_WT.tsv','WT')]:
            for s in rows(root/'marker_summaries'/name):
                if (s['tool'],s['branch'],s['population'])!=(tool,branch,pop): continue
                r=maps[scope][(s['cluster_id'],s['gene_id'])]
                assert r['gene_symbol']==s['gene_symbol']
                for k in numeric+['scaled_expression_z']:
                    assert close(r[k],number(s[k])),(key,scope,k,s['cluster_id'],s['gene_id'])
                    checks+=1
        for s in rows(root/'marker_summaries/gene_cluster_condition_summary.tsv'):
            if (s['tool'],s['branch'],s['population'],s['condition_scope'])!=(tool,branch,pop,'cKO'):continue
            r=maps['cKO'][(s['cluster_id'],s['gene_id'])]
            for k in numeric:
                assert close(r[k],number(s[k])),(key,'cKO',k,s['cluster_id'],s['gene_id'])
                checks+=1
        assert all(r['scaled_expression_z'] is None for r in maps['cKO'].values())
        record={'dataset':key,'n_nuclei':len(nuclei),'n_clusters_audited':len(data['clusters']),
                'marker_values_compared':checks,'coordinates_labels_memberships_PASS':True,
                'produced_marker_values_PASS':True,'cKO_z_missing_preserved':True}
        acceptance.append(record); print(json.dumps(record),flush=True)
        del data,nuclei,native,markers,maps
    receipt={'status':'PASS','n_datasets':len(acceptance),'n_genes':878,
             'workbook_unique_gene_ids':117,'reviewed_workbook_and_canonical_unique_gene_ids':139,'checks':acceptance,
             'comparison_scope':'shards against completed produced metadata/marker tables; no fits or raw RNA aggregation',
             'browser_interaction_inspection_performed':False}
    path=BASE/'reviews/offline_dashboard_data_acceptance.json'
    path.write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps({'status':'PASS','receipt':str(path)}))

if __name__=='__main__':main()
