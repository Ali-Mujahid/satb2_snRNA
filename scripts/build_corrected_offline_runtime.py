#!/usr/bin/env python3
"""Create a UI-only derivative of the preserved classic-script offline runtime."""
from pathlib import Path
import argparse, copy, datetime, hashlib, json, os, shutil, struct, subprocess, zipfile

def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        while block:=f.read(8*1024*1024): h.update(block)
    return h.hexdigest()

def compressed_sha(fp, info):
    fp.seek(info.header_offset)
    header=fp.read(30)
    fields=struct.unpack('<IHHHHHIIIHH',header)
    assert fields[0]==0x04034B50
    fp.seek(fields[9]+fields[10],1)
    h=hashlib.sha256();remaining=info.compress_size
    while remaining:
        block=fp.read(min(remaining,8*1024*1024));assert block
        h.update(block);remaining-=len(block)
    return h.hexdigest()

def replace_once(text,old,new):
    assert text.count(old)==1, 'Expected unique preserved source fragment'
    return text.replace(old,new)

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--repo',type=Path,default=Path(__file__).resolve().parents[1]);parser.add_argument('--node',type=Path,required=True);parser.add_argument('--replace-owned-derivative',action='store_true');args=parser.parse_args()
    repo=args.repo.resolve();archives=repo/'release_staging/archives';source=archives/'SATB2_atlas_v2_runtime.zip';target=archives/'SATB2_atlas_v2_1_runtime.zip';work=repo/'release_staging/offline_ui_v2_1_build'
    if target.exists() or work.exists():
        assert args.replace_owned_derivative, 'Refusing existing derivative or work directory'
        prior=json.loads((repo/'docs/offline_ui_update.json').read_text());assert prior['derivative_archive']==target.name and sha(target)==prior['derivative_archive_sha256'], 'Existing target is not the verified owned derivative'
        backup=archives/'superseded_ui_patches'/datetime.datetime.now(datetime.timezone.utc).strftime('%Y%m%dT%H%M%S%fZ');backup.mkdir(parents=True)
        shutil.move(str(target),str(backup/target.name));shutil.copy2(repo/'docs/offline_ui_update.json',backup/'offline_ui_update.json')
        if work.exists():shutil.move(str(work),str(backup/'build_work'))
    assets=json.loads((archives/'ARCHIVE_ASSETS.json').read_text());expected=next(r for r in assets['assets'] if r['asset_name']==source.name)
    source_sha=sha(source);assert source_sha==expected['sha256']
    work.mkdir(parents=True);folder=work/'offline_dashboard_v2';folder.mkdir()
    with zipfile.ZipFile(source) as original:
        assert len(original.infolist())==3294
        original_app=original.read('offline_dashboard_v2/app.js').decode();original_index=original.read('offline_dashboard_v2/index.html').decode()
        app=original_app
        old="  function fmt(value,digits=4) { return number(value)?(value===0?'0':Math.abs(value)>=1000?value.toLocaleString(undefined,{maximumFractionDigits:2}):value.toPrecision(digits).replace(/\\.?0+(?=e|$)/,'')):'NA'; }"
        new="  function fmt(value,digits=4) {if(!number(value))return 'NA';if(value===0)return '0';if(Math.abs(value)>=1000)return value.toLocaleString(undefined,{maximumFractionDigits:2});const parts=value.toPrecision(digits).split('e');if(parts[0].includes('.'))parts[0]=parts[0].replace(/0+$/,'').replace(/\\.$/,'');return parts.join('e');}"
        app=replace_once(app,old,new)
        app=replace_once(app,"if($('overlay-gene'))$('overlay-gene').value=state.overlayGene||manifest.genes[0]?.gene_id;","if($('overlay-gene'))filterOverlayOptions();")
        app=replace_once(app,"function filterOverlayOptions() {const search=$('overlay-search').value.trim().toLowerCase(),genes=manifest.genes.filter","function filterOverlayOptions() {if(!manifest.genes.some(g=>g.gene_id===state.overlayGene))state.overlayGene=manifest.genes[0]?.gene_id||'';const search=$('overlay-search').value.trim().toLowerCase(),genes=manifest.genes.filter")
        app=replace_once(app,"gene=(preferred[state.population]||[]).map(symbol=>manifest.genes.find(g=>g.gene_symbol===symbol)).find(Boolean)||manifest.genes[0];state.overlayGene=gene.gene_id;","gene=manifest.genes.find(g=>g.gene_id===state.overlayGene)||(preferred[state.population]||[]).map(symbol=>manifest.genes.find(g=>g.gene_symbol===symbol)).find(Boolean)||manifest.genes[0];if(!gene){toast('No saved overlay genes are available.');return;}state.overlayGene=gene.gene_id;")
        app=replace_once(app,"state.annotation=b.dataset.category==='CA1_programs'?'layer_program_label':'reviewed_label';syncControls();await changeScope();","state.annotation=b.dataset.category==='CA1_programs'?'layer_program_label':'reviewed_label';syncControls();persist();await changeScope();")
        index=replace_once(original_index,'Show gene expression example','Show selected gene on UMAP')
        (folder/'app.js').write_text(app);(folder/'index.html').write_text(index)
        manifest_path='offline_dashboard_v2/BUNDLE_MANIFEST.json';manifest=json.loads(original.read(manifest_path));before=copy.deepcopy(manifest)
        changed=[]
        for name,old_bytes in [('app.js',original_app.encode()),('index.html',original_index.encode())]:
            row=next(r for r in manifest['files'] if r['path']==name);assert row['sha256']==hashlib.sha256(old_bytes).hexdigest();new_bytes=(folder/name).read_bytes();changed.append({'path':name,'source_sha256':row['sha256'],'source_bytes':row['bytes'],'derivative_sha256':hashlib.sha256(new_bytes).hexdigest(),'derivative_bytes':len(new_bytes)});row.update(bytes=len(new_bytes),sha256=hashlib.sha256(new_bytes).hexdigest())
        assert before['bytes']==sum(r['bytes'] for r in before['files'])
        manifest['bytes']=sum(r['bytes'] for r in manifest['files']);manifest['generated_at']=datetime.datetime.now(datetime.timezone.utc).isoformat()
        manifest['ui_revision']={'version':'2.1','source_archive':source.name,'source_archive_sha256':source_sha,'derivative_archive':target.name,'scope':'UI-only number formatting, selected-overlay retention and category-shortcut persistence; scientific payloads and classic local transport unchanged','changed_ui_files':changed,'companion_saved_results_archive':'SATB2_atlas_v2_saved_results.zip','companion_folder':'offline_dashboard_v2','new_scientific_fitting':False,'raw_nucleus_aggregation':False}
        assert manifest['file_count']==len(manifest['files'])==before['file_count']
        (folder/'BUNDLE_MANIFEST.json').write_text(json.dumps(manifest,indent=2)+'\n')
        for name in ['v2_features.js','vendor/fflate.min.js','vendor/jspdf.umd.min.js']:
            q=folder/name;q.parent.mkdir(parents=True,exist_ok=True);q.write_bytes(original.read('offline_dashboard_v2/'+name))
        # The loader block must remain byte-for-byte the original classic loader.
        loader=lambda t:t[t.index('  async function loadShard(id) {'):t.index('  function descriptor(tool)')]
        assert loader(app)==loader(original_app)
    print('UI patches, full-folder manifest and unchanged classic loader verified',flush=True)
    qa=work/'qa';qa.mkdir();env={**os.environ,'SATB2_SITE_DIR':str(folder)}
    for name in ['number_format.cjs','overlay_selection.cjs','category_persistence.cjs']:
        shutil.copy2(repo/'tests'/name,qa/name)
        with (qa/(name+'.log')).open('w') as log:subprocess.run([str(args.node),str(qa/name)],env=env,check=True,stdout=log,stderr=subprocess.STDOUT)
    (work/'original_app.js').write_text(original_app)
    js="const fs=require('fs'),vm=require('vm'),assert=require('assert');function keys(p){const c={SATB2_NO_AUTO_INIT:true};c.window=c;c.globalThis=c;vm.createContext(c);vm.runInContext(fs.readFileSync(p,'utf8'),c);return{api:Object.keys(c.SATB2_APP).sort(),helpers:Object.keys(c.SATB2_APP.helpers).sort()};}assert.deepStrictEqual(keys(process.argv[1]),keys(process.argv[2]));console.log('PASS original API/helper contracts');"
    subprocess.run([str(args.node),'-e',js,str(work/'original_app.js'),str(folder/'app.js')],check=True)
    print('Number-format, overlay-selection and original API checks pass',flush=True)
    shutil.copyfile(source,target)
    subprocess.run(['/usr/bin/zip','-q',str(target),'offline_dashboard_v2/app.js','offline_dashboard_v2/index.html','offline_dashboard_v2/BUNDLE_MANIFEST.json'],cwd=work,check=True)
    modified={'offline_dashboard_v2/app.js','offline_dashboard_v2/index.html',manifest_path};unchanged=0
    with zipfile.ZipFile(source) as a,zipfile.ZipFile(target) as b,source.open('rb') as af,target.open('rb') as bf:
        assert len(a.infolist())==len(b.infolist())==3294;assert set(a.namelist())==set(b.namelist())
        for info in a.infolist():
            other=b.getinfo(info.filename)
            if info.filename in modified:continue
            assert (info.CRC,info.file_size,info.compress_size,info.compress_type)==(other.CRC,other.file_size,other.compress_size,other.compress_type)
            assert compressed_sha(af,info)==compressed_sha(bf,other);unchanged+=1
        assert b.testzip() is None
        for name in modified:assert b.read(name)==(work/name).read_bytes()
    assert sha(source)==source_sha
    target_sha=sha(target);target_bytes=target.stat().st_size;assert target_bytes<2*1024**3
    public={'schema':'satb2-offline-ui-update-1','status':'PREPARED_AND_VERIFIED','release_tag':'atlas-v2.1.0','source_archive':source.name,'source_archive_sha256':source_sha,'derivative_archive':target.name,'derivative_archive_sha256':target_sha,'derivative_archive_bytes':target_bytes,'zip_entries':3294,'unchanged_compressed_entries_verified':unchanged,'changed_ui_files':changed,'bundle_manifest_sha256':sha(folder/'BUNDLE_MANIFEST.json'),'bundle_manifest_scope':'Combined offline_dashboard_v2 folder after extracting runtime and unchanged saved-results companion; manifest excludes itself as in original package','companion_saved_results_archive':'SATB2_atlas_v2_saved_results.zip','companion_saved_results_sha256':next(r['sha256'] for r in assets['assets'] if r['asset_name']=='SATB2_atlas_v2_saved_results.zip'),'transport':'UNCHANGED_CLASSIC_LOCAL_SCRIPT_GZIP_BASE64','changes':['Correct scientific notation by trimming fractional mantissa zeros without editing exponent digits','Show selected gene on UMAP; retain the chosen gene when the search filters other options','Persist category-shortcut selection before loading so reload restores the chosen population/programme'],'validation':{'number_format_examples':20,'exponent_preservation_cases':8,'overlay_selection_checks':6,'category_shortcut_checks':10,'original_api_helper_keys':'PASS','classic_loader_source_identical':True,'all_unchanged_entry_compressed_bytes':'PASS','zip_CRC_all_entries':'PASS','original_archive_preserved':'PASS'},'scientific_fitting_performed':False,'raw_nucleus_count_aggregation_performed':False,'scientific_payloads_modified':False,'original_or_previous_v2_files_modified':False}
    docs=repo/'docs';docs.mkdir(exist_ok=True);(docs/'offline_ui_update.json').write_text(json.dumps(public,indent=2)+'\n')
    private={'status':'PREPARED_FOR_UPLOAD','assets':[{'asset_name':target.name,'local_path':str(target),'bytes':target_bytes,'sha256':target_sha,'kind':'corrected_offline_runtime'}],'public_provenance':str(docs/'offline_ui_update.json'),'public_provenance_sha256':sha(docs/'offline_ui_update.json')}
    (archives/'OFFLINE_UI_UPDATE_UPLOAD_MANIFEST.json').write_text(json.dumps(private,indent=2)+'\n');print(json.dumps(public,indent=2),flush=True)

if __name__=='__main__':main()
