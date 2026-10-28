#!/usr/bin/env python3
"""Synthetic third-component transport/gate QA; no study files or models."""
import argparse,base64,copy,gzip,hashlib,importlib.util,json,subprocess,sys,tempfile,zipfile
from pathlib import Path
from unittest.mock import patch

ROOT=Path(__file__).resolve().parent.parent
def module(name,path):
    spec=importlib.util.spec_from_file_location(name,path);m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m);return m
p=module('publication',ROOT/'scripts/build_regulatory_publication.py')
v=module('verification',ROOT/'scripts/verify_regulatory_atlas.py')
def reject(call,phrase):
    try:call()
    except ValueError as e:assert phrase in str(e),(phrase,str(e))
    else:raise AssertionError('Expected rejection: '+phrase)
def write(path,value):
    path.parent.mkdir(parents=True,exist_ok=True);path.write_bytes(p.encoded(value) if not isinstance(value,bytes) else value)
def rec(path):return p.receipt(path)
def native_files(folder):return {str(x.relative_to(folder)):p.receipt(x,folder) for x in p.all_files(folder)}

with tempfile.TemporaryDirectory(prefix='SYNTHETIC_PYSCENIC_TRANSPORT_ONLY_') as tmp:
    original_index_bytes=(ROOT/'site/regulatory/index.html').read_bytes()
    t=Path(tmp).resolve();endpoints=[{'id':ep,'family':a,'source_endpoint_id':b,'n_nuclei':c,'membership_sha256':d} for ep,(a,b,c,d) in v.EXPECTED.items()]
    mice=[{'id':m,'genotype':g} for m,g in v.MICE.items()]
    base={'schema':p.SCHEMA,'analysis_id':'SYNTHETIC_HD_TRANSPORT','source_fingerprint':'a'*64,'status':'partial','endpoints':endpoints,'samples':mice,
          'methods':[],'features':[],'values':[],'effects':[],'hubs':[],'nodes':[],'edges':[],'embeddings':[],'coverage':[],'downloads':[],'provenance':{'synthetic':True,'TF_activity':'NOT_RUN_HISTORICAL'}}
    hd=copy.deepcopy(base)
    hd['methods']=[{'id':'hdwgcna_'+lin+'_'+branch,'family':'hdWGCNA','status':'complete','endpoint_ids':[e['id'] for e in endpoints if e['family']==lin]} for lin in ('CA1','IN') for branch in ('pooled','consensus')]
    tf={**copy.deepcopy(base),'analysis_id':'SYNTHETIC_TF_TRANSPORT','source_fingerprint':'b'*64,'methods':[{'id':'collectri_ulm','family':'curated_tf_activity','status':'complete','endpoint_ids':list(v.EXPECTED)},{'id':'pyscenic_motif','family':'motif_regulon','status':'not_available'}],
        'features':[{'id':'collectri::SYNTHETIC','method_id':'collectri_ulm','label':'SYNTHETIC'}],'provenance':{'synthetic':True,'motif_supported_regulons':'BLOCKED_HISTORICAL'}}
    tf['values']=[{'method_id':'collectri_ulm','endpoint_id':e['id'],'feature_id':'collectri::SYNTHETIC','sample_id':m['id'],'value':float(i),'eligible':True} for e in endpoints for i,m in enumerate(mice)]
    ps={**copy.deepcopy(base),'analysis_id':'SYNTHETIC_PYSCENIC_TRANSPORT','source_fingerprint':'c'*64,'status':'complete','methods':[{'id':'pyscenic_motif_'+lin,'family':'motif_regulon','status':'complete','endpoint_ids':[e['id'] for e in endpoints if e['family']==lin]} for lin in ('CA1','IN')],
        'features':[{'id':'pyscenic_motif_'+lin+'::'+lin+'::Satb2(+)','method_id':'pyscenic_motif_'+lin,'native_name':'Satb2(+)','source_id':lin+'::Satb2(+)','n_targets':10} for lin in ('CA1','IN')],
        'native_supplements':{'Satb2_no_self_target':{'SYNTHETIC':'SEPARATE_UNCHANGED'}},'Satb2_status':{'SYNTHETIC':'EXACT_NATIVE_STATUS'},'provenance':{'synthetic':True,'motif_supported_regulons':'NATIVE_POSITIVE_PYSCENIC_COMPLETE'}}
    ps['values']=[{'method_id':'pyscenic_motif_'+e['family'],'endpoint_id':e['id'],'feature_id':'pyscenic_motif_'+e['family']+'::'+e['family']+'::Satb2(+)','sample_id':m['id'],'value':.1+i*.01,'eligible':True,'coverage_n':10} for e in endpoints for i,m in enumerate(mice)]
    before=copy.deepcopy((hd,tf,ps));combined=p.merge_payloads(hd,tf,ps)
    assert (hd,tf,ps)==before and len(combined['methods'])==7 and not any(m['id']=='pyscenic_motif' for m in combined['methods'])
    assert combined['provenance']['components']['curated_tf_activity']==tf['provenance'] and combined['provenance']['components']['motif_regulon']==ps['provenance']
    assert combined['values']==hd['values']+tf['values']+ps['values'] and combined['Satb2_status']==ps['Satb2_status']
    assert combined['component_supplements']['motif_regulon']==ps['native_supplements']
    for change,phrase in ((lambda x:x['endpoints'][0].update(membership_sha256='0'*64),'memberships'),(lambda x:x['samples'][0].update(id='UNKNOWN'),'roster'),(lambda x:x['features'][0].update(id='collectri::SYNTHETIC'),'overlapping'),(lambda x:x['features'][0].update(native_name='Satb2(-)'),'positive native'),(lambda x:x['values'][0].update(value=1.01),'[0,1]'),(lambda x:x['values'].pop(),'all twelve')):
        bad=copy.deepcopy(ps);change(bad);reject(lambda:p.merge_payloads(hd,tf,bad),phrase)
    empty=copy.deepcopy(ps);empty['features']=[];empty['values']=[]
    for method in empty['methods']:method['status']='not_available'
    empty_merged=p.merge_payloads(hd,tf,empty)
    assert empty_merged['status']=='partial' and not any(row['method_id'].startswith('pyscenic') for row in empty_merged['values'])
    print('PASS Actual motif methods supersede only empty top-level placeholder; native rows, Satb2 status and historical provenance remain immutable; incompatible identities/AUC rejected')

    hd_run=t/'HD';hr=hd_run/'results';hs=t/'HD_SOURCE';(hs/'scripts').mkdir(parents=True)
    for name in ('package_results.py','run_hdwgcna.R'):(hs/'scripts'/name).write_text('# SYNTHETIC, never executed\n')
    for lin in ('CA1','IN'):
        for filename in ('counts.mtx.gz','genes.tsv','cells.tsv'):write(hr/'prepared/lineages'/lin/filename,b'SYNTHETIC_SELECTED_INPUT\n'+(lin+filename).encode())
    write(hr/'prepared/COMPLETE.json',{'SYNTHETIC':'NO_SCIENTIFIC_EXECUTION'})
    setup=t/'SETUP';setup.mkdir();wrapper=setup/'run_R.sh';wrapper.write_text('# SYNTHETIC, never executed\n')
    samples=t/'SYNTHETIC_SAMPLES.tsv';samples.write_text('SYNTHETIC_INPUT_ONLY\n')
    ref=t/'SYNTHETIC_SHARED_INPUT.tsv';ref.write_text('SYNTHETIC_INPUT_ONLY\n')
    shared={name:{ep:rec(ref) for ep in p.TF_ENDPOINTS} for name in ('source_pseudobulks','existing_deseq2')}
    hc=t/'HD_CONFIG.json';write(hc,{'run_root':str(hd_run),'R_wrapper':rec(wrapper),'inputs':{'samples':rec(samples)},**shared})
    write(hr/'FINAL_VALIDATION.json',{'status':'SYNTHETIC_MOCKED','analysis_id':hd['analysis_id'],'source_fingerprint':hd['source_fingerprint']})
    write(hr/'atlas/data.json',p.encoded(hd))
    png=base64.b64decode('iVBORw0KGgoAAAANSUhEUgAAAAEAAAABCAQAAAC1HAwCAAAAC0lEQVR42mP8/x8AAusB9Wl6bAAAAABJRU5ErkJggg==')
    write(hr/'hdwgcna/CA1/pooled/figures/SYNTHETIC_HD.png',png)
    hv={'final':p.load(hr/'FINAL_VALIDATION.json'),'config':p.load(hc),'payload':hd,'manifest':{'analysis_id':hd['analysis_id'],'source_fingerprint':hd['source_fingerprint'],'branches':{lin+'_pooled':{'validated_complete':True} for lin in ('CA1','IN')}},'files':native_files(hr),'bindings':{'R_wrapper':rec(wrapper)}}
    tr=t/'TF/results';ts=t/'TF_SOURCE';(ts/'scripts').mkdir(parents=True)
    for name in ('run_tf_activity.py','analyze_tf_mouse.R','package_tf_results.py'):(ts/'scripts'/name).write_text('# SYNTHETIC, never executed\n')
    tc=t/'TF_CONFIG.json';write(tc,shared);write(tr/'FINAL_VALIDATION.json',{'status':'SYNTHETIC_MOCKED'})
    write(tr/'atlas/data.json',p.encoded(tf))
    write(tr/'tf_mouse/figures/SYNTHETIC_TF.png',png)
    tv={'final':{'analysis_id':tf['analysis_id'],'source_fingerprint':tf['source_fingerprint']},'config':shared,'payload':tf,'manifest':{'analysis_id':tf['analysis_id'],'source_fingerprint':tf['source_fingerprint'],'completed_primary':'curated_mouse_CollecTRI_ULM_and_independent_mouse_associations','blocked_optional':'pySCENIC_motif_cache_missing'},'files':native_files(tr),'bindings':{'shared':rec(ref)}}
    pr=t/'PS/results';src=t/'PS_SOURCE';(src/'scripts').mkdir(parents=True)
    for name in p.PYSCENIC_SOURCE:(src/name).write_text('# SYNTHETIC, never executed\n')
    pc=t/'PS_CONFIG.json';rankings=[];resources=[]
    for i,name in enumerate(('tf_names','motif_annotations','rankings0','rankings1')):
        path=t/('SYNTHETIC_'+name+('.feather' if name.startswith('rankings') else '.txt'));path.write_text('SYNTHETIC_REFERENCE_BYTES_'+name+'\n')
        row=rec(path);resources.append({**row,'requested_url':'https://resources.aertslab.org/cistarget/SYNTHETIC_TEST_ONLY/'+path.name})
        if name.startswith('rankings'):rankings.append(row)
        elif name=='tf_names':tf_names=row
        else:annotations=row
    ready=t/'READY_MANIFEST.json';write(ready,{'SYNTHETIC':True,'resources':resources})
    cfg={'schema':'satb2-pyscenic/v1','analysis_id':ps['analysis_id'],'run_root':str(pr.parent),'samples':rec(samples),'prepared_complete':rec(hr/'prepared/COMPLETE.json'),
         'lineages':{lin:{key:rec(hr/'prepared/lineages'/lin/filename) for key,filename in (('counts','counts.mtx.gz'),('genes','genes.tsv'),('cells','cells.tsv'))} for lin in ('CA1','IN')},
         'references':{'tf_names':tf_names,'motif_annotations':annotations,'rankings':rankings},'additional_hash_inputs':{'frozen_membership_audit':rec(ref),'mouse_coverage':rec(ref),'reference_ready':rec(ready)}}
    write(pc,cfg);source_hashes={name:p.digest(src/name) for name in p.PYSCENIC_SOURCE}
    for relative,schema,source in (('pyscenic_native/COMPLETE.json','satb2-pyscenic-native-complete/v1',p.PYSCENIC_SOURCE[0]),('regulon_mouse/COMPLETE.json','satb2-regulon-mouse-complete/v1',p.PYSCENIC_SOURCE[1])):
        write(pr/relative,{'schema':schema,'status':'COMPLETE','config':rec(pc),'source':rec(src/source)})
    axis={}
    for lin,n in (('CA1',20456),('IN',2687)):
        axis[lin]={'n_nuclei':n,'exact_source_nucleus_axis_verified':True,'no_score_recomputation_or_averaging':True,'auc_table':lin+'/aucell/nucleus_AUC.tsv.gz','metadata_table':lin+'/aucell/nucleus_axis_metadata.tsv.gz'}
        for key in ('auc_table','metadata_table'):write(pr/'pyscenic_native'/axis[lin][key],gzip.compress(b'SYNTHETIC_TEST_ONLY\n'))
    write(pr/'regulon_mouse/figures/SYNTHETIC_MOTIF.png',png)
    for folder,checks in (('pyscenic_native',p.PYSCENIC_NATIVE_CHECKS),('regulon_mouse',p.PYSCENIC_MOUSE_CHECKS)):
        path=pr/folder/'COMPLETE.json';marker=p.load(path);marker['checks']={k:True for k in checks}
        marker['output_files' if folder=='pyscenic_native' else 'outputs']={str(x.relative_to(path.parent)):rec(x) for x in p.all_files(path.parent) if x!=path}
        write(path,marker)
    ps['provenance']['packager_sha256']=p.digest(src/p.PYSCENIC_SOURCE[2]);raw=p.encoded(ps);atlas=pr/'atlas';write(atlas/'data.json',raw);write(atlas/'data.js',b'/* SYNTHETIC ONLY: decoder explicitly mocked in filesystem test */\n')
    manifest={'schema':'satb2-regulatory-transport/v1','analysis_id':ps['analysis_id'],'source_fingerprint':ps['source_fingerprint'],'status':'complete','completed_primary':'positive_motif_supported_regulons_native_nucleus_AUCell_and_independent_mouse_associations','transport':{name:p.receipt(atlas/name,atlas) for name in ('data.js','data.json')}}
    write(atlas/'manifest.json',manifest)
    flags=('TF_activity_run','motif_regulon_run','positive_motif_regulons_only','twelve_mice_per_endpoint_regulon','unique_mouse_keys','raw_mouse_values_finite','frozen_membership_hashes_verified','native_score_receipts_verified','exact_nucleus_axes_verified','Satb2_self_target_sensitivity_separate','all_configured_inputs_SHA256_before_after_verified','original_analysis_files_unchanged','no_scientific_fits_in_packager')
    final={'schema':'satb2-pyscenic-regulatory-final-validation/v1','status':'COMPLETE','analysis_id':ps['analysis_id'],'source_fingerprint':ps['source_fingerprint'],'configuration_sha256':p.digest(pc),'executed_source_sha256':source_hashes,'minimum_measured_targets':10,'activity_UMAP_packaged':False,'Satb2_official_TF_list_membership':False,'native_nucleus_axis_audit':axis,**{k:True for k in flags},'native_file_receipts':[r for k,r in native_files(pr).items() if not k.startswith('atlas/')],'package_file_receipts':[p.receipt(x,pr) for x in p.all_files(atlas)]}
    write(pr/'FINAL_VALIDATION.json',final)
    with patch.object(p,'decode_data',return_value=(ps,raw)):
        pv=p.validate_pyscenic_complete(pr,src,pc);p.validate_pyscenic_shared_inputs(hv,pv)
        write(pr/'.snakemake/log/SYNTHETIC_CONTROLLER.log',b'SYNTHETIC MUTABLE CONTROLLER METADATA\n')
        assert p.validate_pyscenic_complete(pr,src,pc)['files']==pv['files']
        extra=pr/'pyscenic_native/SYNTHETIC_UNRECEIPTED.tsv';write(extra,b'UNKNOWN SCIENTIFIC FILE\n');reject(lambda:p.validate_pyscenic_complete(pr,src,pc),'not exhaustive');extra.unlink()
        for flag in ('positive_motif_regulons_only','exact_nucleus_axes_verified','Satb2_self_target_sensitivity_separate'):
            bad=copy.deepcopy(final);bad[flag]=False;write(pr/'FINAL_VALIDATION.json',bad);reject(lambda:p.validate_pyscenic_complete(pr,src,pc),flag)
        bad=copy.deepcopy(final);bad['minimum_measured_targets']=5;write(pr/'FINAL_VALIDATION.json',bad);reject(lambda:p.validate_pyscenic_complete(pr,src,pc),'target eligibility')
        write(pr/'FINAL_VALIDATION.json',final)
        original=Path(rankings[0]['path']).read_bytes();Path(rankings[0]['path']).write_bytes(b'CHANGED REFERENCE\n');reject(lambda:p.validate_pyscenic_complete(pr,src,pc),'bytes / SHA256 differ');Path(rankings[0]['path']).write_bytes(original)
        assert p.validate_pyscenic_complete(pr,src,pc)['final']==final
    changed=copy.deepcopy(pv);changed['config']['lineages']['IN']['counts']['sha256']='0'*64;reject(lambda:p.validate_pyscenic_shared_inputs(hv,changed),'selected raw lineage')
    print('PASS Real motif completion/receipt/source/input/axis gate rejects incomplete support, min5, unreceipted scientific files and changed reference/raw-count bytes; controller metadata alone is ignored')

    out=t/'SYNTHETIC_THREE_COMPONENT_STAGE';args=argparse.Namespace(results=hr,source=hs,config=hc,tf_results=tr,tf_source=ts,tf_config=tc,pyscenic_results=pr,pyscenic_source=src,pyscenic_config=pc,outdir=out,site_dir=t/'DEPLOYED_SITE/regulatory',setup_dir=setup,tag='atlas-v3.8.0',stage_only=True)
    original_bytes={str(path):path.read_bytes() for path in (hr/'FINAL_VALIDATION.json',tr/'FINAL_VALIDATION.json',pr/'FINAL_VALIDATION.json',pc,src/p.PYSCENIC_SOURCE[0])}
    with patch.object(p,'validate_complete',return_value=hv),patch.object(p,'validate_tf_complete',return_value=tv),patch.object(p,'validate_pyscenic_complete',return_value=pv):report=p.build(args)
    assert report['motif_regulon_run'] and report['large_public_motif_references_external_hash_bound'] and report['native_figure_groups']==3 and p.verify_outdir(out)['status']=='PASS'
    merged_manifest=p.load(out/'site/regulatory/manifest.json')
    # Production decode rejects our synthetic identity; explicit synthetic helper
    # bypass only for examining generated bytes. CLI has no such bypass option.
    script=(out/'site/regulatory/data.js').read_text();payload=json.loads(gzip.decompress(base64.b64decode(json.loads(script.split('=',1)[1].strip().rstrip(';')))))
    assert payload['values']==tf['values']+ps['values'] and payload['Satb2_status']==ps['Satb2_status']
    assert v.validate_payload(payload,merged_manifest)['motif_component'] is True
    contents=p.load(out/'ARCHIVE_CONTENTS.json');locations=contents['pySCENIC_bound_input_locations']
    external=[r for r in locations if r['location']=='external_official_reference'];aliases=[r for r in locations if r['location']=='shared_byte_identical_archive_member']
    assert len(external)==3 and all(not r['included_in_archive'] for r in external) and len(aliases)>=6
    archived={r['archive_path'] for r in contents['files']}
    assert all(r['archive_path'] in archived for r in aliases) and any('reference_ready' in name for name in archived)
    assert not any(r['sha256'] in {x['sha256'] for x in external} for r in contents['files'])
    assert 'resource-blocked' not in (out/'site/regulatory/downloads.html').read_text() and 'native positive pySCENIC' in (out/'site/regulatory/gallery.html').read_text()
    for page in (out/'site/regulatory/index.html',out/'release_assets/SATB2_regulatory_interactive_offline.html'):
        text=page.read_text()
        assert 'Gene programs, TF activity &amp; motif regulons' in text
        assert 'analysis has completed' in text and 'no eligible positive regulons' in text
        assert 'curated CollecTRI analysis and its separate self-target sensitivity' in text and 'Read the motif/AUCell models separately for their actual genotype effects.' in text
        assert 'Earlier resource-blocked reports remain preserved as historical records' in text
        assert 'Motif-supported pySCENIC inference remains unavailable' not in text and '<b>SATB2 activity is not FDR-significant</b>' not in text
        assert 'blob/atlas-v3.7.0/analysis/reviews/20261026_regulatory_findings.md' not in text
        assert 'Show curated SATB2 activity in all five populations' in text
    base_index=(ROOT/'site/regulatory/index.html').read_text()
    assert (ROOT/'site/regulatory/index.html').read_bytes()==original_index_bytes,'The deployed/source index must remain unchanged during stage-only rendering'
    assert p.completed_motif_index(p.completed_motif_index(base_index))==p.completed_motif_index(base_index),'Completed motif wording must be stable if the next build starts from a prior completed template'
    extracted=t/'EXTRACTED'
    with zipfile.ZipFile(report['complete_native_zip']['path']) as z:z.extractall(extracted)
    helper=next(extracted.rglob('VALIDATE_EXTRACTED.py'));proc=subprocess.run([sys.executable,'-B',str(helper)],capture_output=True,text=True,check=True);assert proc.stdout.startswith('PASS ')
    readme=(helper.parent/'README.txt').read_text();assert p.digest(ready) in readme
    for row in external:
        assert row['sha256'] in readme and str(row['bytes']) in readme and row['download_provenance']['requested_url'] in readme
    assert (helper.parent/'results/atlas/data.json').read_bytes()==p.encoded(hd)
    assert (helper.parent/'tf_component/results/atlas/data.json').read_bytes()==p.encoded(tf)
    assert any(m['id']=='pyscenic_motif' and m['status']=='not_available' for m in json.loads((helper.parent/'tf_component/results/atlas/data.json').read_bytes())['methods'])
    assert all(Path(name).read_bytes()==content for name,content in original_bytes.items()) and not args.site_dir.exists()
    reject(lambda:p.decode_data(out/'site/regulatory/data.js'),'Synthetic')
    print('PASS Three-component native archive/figure gallery/offline transport: exact immutable native bytes, external official references, deduplicated selected raw inputs and extraction checks')

    bad=copy.deepcopy(payload);bad['methods'].append({'id':'pyscenic_motif','family':'motif_regulon','status':'not_available'});reject(lambda:v.validate_payload(bad,merged_manifest),'method IDs')
    bad=copy.deepcopy(payload);bad['values'][-1]['coverage_n']=5;reject(lambda:v.validate_payload(bad,merged_manifest),'minimum-target')
    bad=copy.deepcopy(merged_manifest);bad['component_validation']['motif_regulon']['positive_motif_regulons_only']=False;reject(lambda:v.validate_payload(payload,bad),'validation flags')
    bad=copy.deepcopy(payload);bad['features'][-1]['native_name']='Satb2(-)';reject(lambda:v.validate_payload(bad,merged_manifest),'positive native')
    bad=copy.deepcopy(payload);bad['values']=[r for r in bad['values'] if not(r['method_id']=='pyscenic_motif_CA1' and r['endpoint_id']=='CA1_DL')];reject(lambda:v.validate_payload(bad,merged_manifest),'entire matching')
    unmodelable=copy.deepcopy(payload)
    for m in unmodelable['methods']:
        if m['id']=='pyscenic_motif_IN':m['status']='not_available'
    for row in unmodelable['values']:
        if row['method_id']=='pyscenic_motif_IN':row.update(value=None,eligible=False,missing_reason='SYNTHETIC_LOW_TARGET_COUNT',coverage_n=2)
    assert v.validate_payload(unmodelable,merged_manifest)['motif_component']
    print('PASS Production verifier distinguishes completed native motif component from historical placeholder, validates positive lineage definitions/min10, and retains explicit unavailable lineage nulls')

    for overrides,phrase in (({'pyscenic_source':None},'together'),({'tf_results':None,'tf_source':None,'tf_config':None},'third pySCENIC')):
        bad=argparse.Namespace(**{**vars(args),**overrides,'outdir':t/'SYNTHETIC_REJECTED'})
        with patch.object(p,'validate_complete',return_value=hv),patch.object(p,'validate_tf_complete',return_value=tv):reject(lambda:p.build(bad),phrase)
    assert p.MAX_SITE==950_000_000 and p.MAX_NEW_PAGES==38_000_000
    old_tag=argparse.Namespace(**{**vars(args),'tag':'atlas-v3.7.0','outdir':t/'SYNTHETIC_OLD_TAG_REJECTED'})
    with patch.object(p,'validate_complete',return_value=hv),patch.object(p,'validate_tf_complete',return_value=tv),patch.object(p,'validate_pyscenic_complete',return_value=pv):reject(lambda:p.build(old_tag),'preserve v3.7')
    print('PASS All-or-none third-component CLI contract and prerequisite completed TF component; total950MB cap unchanged')

print(json.dumps({'status':'PASS','checks':5,'synthetic_only':True,'native_scientific_validation_mocked_only_for_transport_build':True,'real_outputs_opened':False,'scientific_models_run':False,'published':False},indent=2))
