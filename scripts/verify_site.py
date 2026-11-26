"""Validate published transport identities without executing scientific code."""
from pathlib import Path
import gzip,hashlib,json,re,sys

ROOT=Path(__file__).resolve().parent.parent
SITE=ROOT/'site'
PREFIX='https://raw.githubusercontent.com/Ali-Mujahid/satb2_snRNA/atlas-v2.1.0/data_downloads/'
GENE_PREFIX='https://raw.githubusercontent.com/Ali-Mujahid/satb2_snRNA/atlas-v3.1.0/gene_expression/'
# Canonical v3.0 registry SHA, independently matched against the previous Git HEAD.
# These 3,140 descriptors and their scientific payloads must remain unchanged.
PREVIOUS_SHARD_REGISTRY_SHA256='2361052aad5957cf44ed5364fdf996c6295dd79207fa9ddbc59fda04df1923e4'
BRANCHES={'original_full','paired_original','paired_cellbender'}
MAX_SAFE_INTEGER=2**53-1

def sha(data):
    return hashlib.sha256(data).hexdigest()

def verify_gene_expression(manifest):
    """Check all outside-Pages per-nucleus counts and their immutable source axes."""
    catalog_path=SITE/'data/GENE_EXPRESSION_CATALOG.json'
    raw_catalog=catalog_path.read_bytes()
    catalog=json.loads(raw_catalog)
    js=(SITE/'data/gene_expression_catalog.js').read_text().strip()
    assignment='window.SATB2_GENE_EXPRESSION_CATALOG='
    assert js.startswith(assignment) and js.endswith(';'),'Gene catalog JS assignment'
    assert json.loads(js[len(assignment):-1])==catalog,'Gene JSON/JS catalogs differ'
    assert catalog['schema']=='satb2-nucleus-expression-catalog-1'
    assert catalog['release_tag']=='atlas-v3.1.0'
    assert catalog['scientific_model_fitting_performed'] is False
    assert catalog['genes']==manifest['genes'],'Gene registry order or symbols changed'
    gene_symbols={g['gene_id']:g['gene_symbol'] for g in catalog['genes']}
    assert len(gene_symbols)==len(catalog['genes'])==889
    assert set(catalog['branches'])==BRANCHES
    assert len(catalog['shards'])==2670
    norm=catalog['normalization']
    assert norm['target_sum']==10000 and norm['transform']=='log1p'
    assert norm['denominator']=='complete source library'
    receipt=json.loads((SITE/'data/GENE_EXPRESSION_PROVENANCE.json').read_text())
    assert receipt==json.loads((ROOT/'gene_expression/PROVENANCE.json').read_text())
    assert receipt['status']=='COMPLETE' and receipt['catalog_sha256']==sha(raw_catalog)
    assert receipt['scientific_model_fitting_performed'] is False
    assert receipt['source_count_objects_modified'] is False
    assert receipt['coordinate_membership_annotation_QC_DESeq2_modified'] is False
    assert receipt['descriptive_per_nucleus_count_export_performed'] is True
    assert receipt['n_shards']==2670 and receipt['n_genes']==889
    assert receipt['release_tag']==catalog['release_tag']
    branch_receipts={item['branch']:item for item in receipt['branches']}
    assert set(branch_receipts)==BRANCHES
    checked=set(); files=set(); total_bytes=0; coordinate_ids_checked=0; sparse_nonzero_counts=0

    def payload(ident):
        nonlocal total_bytes
        assert ident not in checked,('Repeated gene-expression shard',ident)
        checked.add(ident)
        desc=catalog['shards'][ident]
        name=ident+'.json.gz'
        assert re.fullmatch(r'[-A-Za-z0-9_.]+\.json\.gz',name),name
        assert desc['path']==GENE_PREFIX+name,('Gene source must use exact pinned release path',ident)
        assert desc['encoding']=='gzip-bytes'
        file=ROOT/'gene_expression'/name
        assert file.resolve().parent==(ROOT/'gene_expression').resolve()
        assert file.resolve().is_relative_to(ROOT.resolve()) and not file.resolve().is_relative_to(SITE.resolve())
        packed=file.read_bytes();raw=gzip.decompress(packed)
        assert sha(packed)==desc['sha256']==receipt['output_sha256'][name],ident
        assert sha(raw)==desc['payload_sha256'],ident
        assert len(packed)==desc['gzip_bytes'] and len(raw)==desc['uncompressed_bytes'],ident
        files.add(name);total_bytes+=len(packed)
        return json.loads(raw)

    for branch in sorted(BRANCHES):
        item=catalog['branches'][branch]
        assert set(item['genes'])==set(gene_symbols),branch
        assert item['axis_shard']=='axis__'+branch
        assert re.fullmatch(r'[0-9a-f]{64}',item['source_count_object_sha256'])
        assert item['source_count_object_release_asset']
        audit=branch_receipts[branch]
        assert audit['source_sha256']==item['source_count_object_sha256']
        assert audit['source_release_asset']==item['source_count_object_release_asset']
        assert audit['n_source_nuclei']==item['n_source_nuclei'] and audit['exported_genes']==889
        assert audit['all_full_libraries_equal_native_verified_libraries'] is True
        assert audit['source_889_gene_sums_and_detection_equal_decoded_vectors'] is True
        assert audit['independent_backed_native_cell_gene_agreement']=='EXACT'
        axis_desc=catalog['shards'][item['axis_shard']]
        axis=payload(item['axis_shard'])
        assert axis['schema']=='satb2-nucleus-expression-axis-1' and axis['branch']==branch
        nuclei=axis['nucleus_ids'];libraries=axis['library_counts'];n=item['n_source_nuclei']
        assert len(nuclei)==len(libraries)==n and n>0,branch
        assert all(isinstance(x,str) and x for x in nuclei),branch
        assert len(set(nuclei))==n,branch
        assert all(type(x) is int and 0<x<=MAX_SAFE_INTEGER for x in libraries),branch
        axis_ids=set(nuclei)
        for dataset in manifest['datasets']:
            if dataset['branch']!=branch:continue
            desc=manifest['shards'][dataset['coordinates_shard']]
            coords=json.loads(gzip.decompress((SITE/desc['path']).read_bytes()))
            coordinate_ids=[x['nucleus_id'] for x in coords['nuclei']]
            assert len(coordinate_ids)==len(set(coordinate_ids)),dataset['id']
            assert set(coordinate_ids)<=axis_ids,('Existing coordinates absent from source axis',dataset['id'])
            coordinate_ids_checked+=len(coordinate_ids)
        for gene_id in sorted(gene_symbols):
            ident=item['genes'][gene_id]
            assert ident=='gene__'+branch+'__'+gene_id,ident
            vector=payload(ident)
            assert vector['schema']=='satb2-nucleus-gene-counts-1'
            assert vector['branch']==branch and vector['gene_id']==gene_id
            assert vector['gene_symbol']==gene_symbols[gene_id]
            assert vector['n_source_nuclei']==n
            assert vector['axis_payload_sha256']==axis_desc['payload_sha256']
            indices=vector['indices'];counts=vector['counts']
            assert len(indices)==len(counts),ident
            previous=-1
            for i,count in zip(indices,counts):
                assert type(i) is int and previous<i<n,ident
                assert type(count) is int and 0<count<=MAX_SAFE_INTEGER and count<=libraries[i],ident
                previous=i
            sparse_nonzero_counts+=len(counts)
    assert checked==set(catalog['shards'])
    assert files=={p.name for p in (ROOT/'gene_expression').glob('*.json.gz')}
    assert files==set(receipt['output_sha256'])
    assert total_bytes==receipt['total_compressed_bytes']
    return {'n_gene_expression_shards':len(checked),'n_expression_branches':len(BRANCHES),
        'gene_expression_bytes_outside_pages':total_bytes,'coordinate_ids_covered':coordinate_ids_checked,
        'sparse_nonzero_counts_checked':sparse_nonzero_counts,'gene_json_js_catalogs_match':True}

def main():
    names=['manifest.json','DE_CATALOG.json','MOUSE_CATALOG.json']
    maps=[json.loads((SITE/'data'/name).read_text()) for name in names]
    previous_registry={name:registry['shards'] for name,registry in zip(names,maps)}
    canonical=json.dumps(previous_registry,sort_keys=True,separators=(',',':'),ensure_ascii=False).encode()
    assert sha(canonical)==PREVIOUS_SHARD_REGISTRY_SHA256,'The previous 3,140 scientific shard descriptors changed'
    count=0
    for registry in maps:
        for key,d in registry['shards'].items():
            path=d['path']
            if path.startswith(PREFIX):
                name=path[len(PREFIX):]
                assert re.fullmatch(r'[A-Za-z0-9_.-]+\.json\.gz',name),name
                p=ROOT/'data_downloads'/name
            else:
                assert path.startswith('data/shards/') and '..' not in path and '://' not in path,path
                p=SITE/path
            packed=p.read_bytes();assert hashlib.sha256(packed).hexdigest()==d['sha256'],key
            raw=gzip.decompress(packed);assert hashlib.sha256(raw).hexdigest()==d['payload_sha256'],key
            count+=1
    assert count==3140,count
    total=sum(p.stat().st_size for p in SITE.rglob('*') if p.is_file())
    assert total<950_000_000,total
    html=(SITE/'index.html').read_text();ids=re.findall(r'\bid="([^"]+)"',html)
    assert len(ids)==len(set(ids))
    for h in SITE.glob('*.html'):
        for ref in re.findall(r'(?:src|href)="([^"]+)"',h.read_text()):
            if ref.startswith(('https://','http://','#','data:','mailto:')):continue
            path=ref.split('#')[0].split('?')[0]
            if path:assert (h.parent/path).is_file(),(h.name,ref)
    assert len(maps[0]['datasets'])==24
    assert len(maps[0]['genes'])==889
    preset=next(p for p in maps[0]['gene_sets'] if p['id']=='screenshot_20261015_75')
    source=json.loads((ROOT/'analysis/dashboard/v3_marker_extension/manifest_patch.json').read_text())
    assert preset['gene_ids']==source['gene_sets']['Screenshot75_20261015']
    assert len(set(preset['gene_ids']))==75
    assert [g for b in preset['gene_blocks'] for g in b['gene_ids']]==preset['gene_ids']
    assert maps[0]['marker_extension']['new_scientific_model_fitting_performed'] is False
    assert len(maps[2]['genes'])==32285 and len(maps[2]['fits'])==1060
    expression_checks=verify_gene_expression(maps[0])
    print(json.dumps({'status':'PASS','shards':count,'previous_3140_shards_preserved':True,'site_bytes':total,'n_datasets':24,'n_marker_genes':889,'n_screenshot_genes':75,'n_mouse_genes':32285,'n_mouse_endpoints':1060,'scientific_model_fitting_performed':False,**expression_checks}))

if __name__=='__main__':main()
