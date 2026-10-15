/* Actual produced-data/state checks in a JavaScript VM; no browser navigation. */
'use strict';
const fs=require('fs'),path=require('path'),vm=require('vm'),assert=require('assert'),zlib=require('zlib'),crypto=require('crypto');
const base=path.dirname(__dirname), bundle=path.join(base,'reports/offline_dashboard_v1');
const source=path.join(base,'offline_dashboard_source');
const checks=[];
function check(name,fn){fn();checks.push(name);console.log(JSON.stringify({check:name,status:'PASS'}));}
function readJson(name){return JSON.parse(fs.readFileSync(path.join(bundle,name),'utf8'));}
function hash(bytes){return crypto.createHash('sha256').update(bytes).digest('hex');}
function readPayload(id,desc){
  assert(desc&&desc.path,'Missing descriptor: '+id);
  const bytes=fs.readFileSync(path.join(bundle,desc.path));
  assert.equal(hash(bytes),desc.sha256,'Shard integrity '+id);
  const text=bytes.toString('utf8');
  const m=text.match(/(?:\bdata\s*:|"data"\s*:)\s*("[A-Za-z0-9+/=]+")/);
  assert(m,'No gzip data: '+id);
  return {encoding:'gzip-base64',data:JSON.parse(m[1])};
}
const context={console,Uint8Array,Uint16Array,Uint32Array,Int32Array,Float32Array,Float64Array,ArrayBuffer,TextEncoder,TextDecoder,atob,btoa,Blob,setTimeout,clearTimeout,performance,
  localStorage:{getItem(){return null;},setItem(){}},SATB2_NO_AUTO_INIT:true};
context.window=context;context.self=context;context.globalThis=context;vm.createContext(context);
const appBytes=fs.readFileSync(path.join(source,'app.js'));
vm.runInContext(fs.readFileSync(path.join(source,'vendor/fflate.min.js'),'utf8'),context,{filename:'fflate.min.js'});
vm.runInContext(appBytes.toString('utf8'),context,{filename:'app.js'});
const app=context.SATB2_APP;
const plain=value=>JSON.parse(JSON.stringify(value));
function key(tool){return tool+'|'+app.state.population+'|original_frozen';}
function clusterKey(tool,id){return key(tool)+'|'+id;}
async function scope(pop,branch='original_full',tab='de'){
  app.state.population=pop;app.state.branch=branch;app.state.tab=tab;
  app.state.deEndpoint='PRIMARY';app.state.exclusions={};app.state.labels={};
  await app.changeScope();
  assert(app.active.scanpy&&app.active.seurat,'Scope failed '+pop+'/'+branch);
}

(async()=>{
  const manifest=readJson('data/manifest.json'),catalog=readJson('DE_CATALOG.json');
  assert.equal(catalog.status,'COMPLETE','Final acceptance requires the completed catalog');
  assert.equal(catalog.native_comparison_completed,true);
  await app.init({manifest,catalog,headless:true,ignoreSavedState:true,shardLoader:readPayload});
  check('Default coordinates load for both tools',()=>{
    assert.equal(app.active.scanpy.descriptor.population,manifest.defaults.population);
    assert.equal(app.active.seurat.coords.nuclei.length,app.active.scanpy.coords.nuclei.length);
  });
  for(const pop of ['allcells','Excitatory','Inhibitory','CA1_ProS']){
    await scope(pop);
    for(const tool of ['scanpy','seurat']){
      const fit=app.resolveFit(tool);
      check(pop+'/'+tool+' baseline total binds exact source fit',()=>{
        assert.equal(fit.status,'EXACT_MATCH');assert.equal(fit.fit.analysis_kind,'total');
        assert.equal(fit.fit.members.length,app.active[tool].coords.nuclei.length);
        assert(app.active[tool].deRows.length>0);
      });
      const originalId=fit.fit.id;
      const clusters=[...new Set(app.active[tool].coords.nuclei.map(n=>String(n.cluster_id)))];
      app.state.exclusions[key(tool)]=[clusters[0]];
      await app.selectionChanged();
      const edited=app.resolveFit(tool);
      check(pop+'/'+tool+' edited total cannot borrow whole-population inference',()=>{
        assert(!edited.fit||edited.fit.id!==originalId);
        if(edited.status!=='EXACT_MATCH')assert.equal(app.active[tool].deRows,null);
        const request=app.request(tool);
        assert.equal(request.n_nuclei,app.targetMembers(tool).length);
        assert(request.no_fit_performed_in_dashboard);
        assert.equal(request.selected_nucleus_ids.length,request.n_nuclei);
      });
      app.state.deTargets[tool]='cluster:'+clusters[1];
      await app.selectionChanged();
      const individual=app.resolveFit(tool);
      check(pop+'/'+tool+' unchanged cluster remains source-bound after other exclusion',()=>{
        assert(['EXACT_MATCH','NOT_TESTED'].includes(individual.status));
        assert(individual.members.every(i=>String(app.active[tool].coords.nuclei[i].cluster_id)===clusters[1]));
      });
      app.state.exclusions[key(tool)]=clusters;
      app.state.deTargets[tool]='total';await app.selectionChanged();
      check(pop+'/'+tool+' exclude-all has no inference or stale rows',()=>{
        const x=app.resolveFit(tool);assert.equal(x.members.length,0);
        assert(['NO_NUCLEI','NOT_TESTED'].includes(x.status));
        assert.equal(app.active[tool].deRows,null);
      });
      app.state.exclusions[key(tool)]=[];app.state.deTargets[tool]='total';await app.selectionChanged();
    }
  }
  await scope('CA1_ProS');await app.ensureMarkers();
  const tool='scanpy',clusters=[...new Set(app.active[tool].coords.nuclei.map(n=>String(n.cluster_id)))];
  app.state.labels[clusterKey(tool,clusters[0])]='Displayed name';
  check('Cosmetic rename preserves default total fit',()=>assert.equal(app.resolveFit(tool).status,'EXACT_MATCH'));
  app.state.deTargets[tool]='alias:Displayed name';await app.selectionChanged();
  check('One-cluster alias resolves exact cluster membership',()=>{
    const x=app.resolveFit(tool);assert(['EXACT_MATCH','NOT_TESTED'].includes(x.status));
    assert(x.members.every(i=>String(app.active[tool].coords.nuclei[i].cluster_id)===clusters[0]));
  });
  app.state.labels[clusterKey(tool,clusters[1])]='Displayed name';await app.selectionChanged();
  check('Merged label uses exact union or no matching fit',()=>{
    const x=app.resolveFit(tool),want=app.active[tool].coords.nuclei.reduce((n,r)=>n+Number([clusters[0],clusters[1]].includes(String(r.cluster_id))),0);
    assert.equal(x.members.length,want);
    if(x.status==='EXACT_MATCH')assert.deepEqual(plain(x.members),plain(x.fit.members));
    else assert.equal(app.active[tool].deRows,null);
  });
  app.state.markerGrouping='alias';app.state.markerCondition='ALL';
  const gene=manifest.defaults.marker_ids[0];
  const sourceRows=app.active[tool].markers.rows.filter(r=>[clusters[0],clusters[1]].includes(String(r.cluster_id))&&r.gene_id===gene);
  const custom=app.dotRows(tool,[gene]).find(r=>r.display_label==='Displayed name');
  check('Custom marker union retains correct additive numerators and denominator',()=>{
    assert(custom);for(const field of ['n_nuclei','sum_counts','n_detected'])assert.equal(custom[field],sourceRows.reduce((s,r)=>s+r[field],0));
    assert(Math.abs(custom.fraction_detected-custom.n_detected/custom.n_nuclei)<1e-12);
  });
  app.state.markerMetric='scaled_expression_z';
  const neg=app.drawDotValues({scaled_expression_z:-1,fraction_detected:.8});
  check('Negative-z dots retain the scientific value and receive requested appearance',()=>{
    assert.equal(neg.value,-1);assert.equal(neg.area,20);assert.equal(neg.alpha,.18);
    app.state.negativeZ=false;assert.equal(app.drawDotValues({scaled_expression_z:-1,fraction_detected:.8}).area,80);app.state.negativeZ=true;
  });
  app.state.markerCondition='cKO';
  check('Unavailable cKO scaled values are missing rather than zero',()=>{
    const rr=app.dotRows(tool,[gene]);assert(rr.length);assert(rr.every(r=>r.scaled_expression_z===null));
    assert.equal(app.drawDotValues(rr[0]).value,null);
  });
  app.state.markerCondition='ALL';app.state.markerMetric='mean_log1p_normalized';app.state.sharedScale=true;
  check('Actual marker values share one numeric color range across both tools',()=>{
    const a=app.markerRange('scanpy',[gene]),b=app.markerRange('seurat',[gene]);
    assert.deepEqual(plain(a),plain(b));
    const vals=['scanpy','seurat'].flatMap(t=>app.dotRows(t,[gene]).map(r=>r.mean_log1p_normalized)).filter(Number.isFinite);
    assert.strictEqual(a[0],Math.min(...vals));assert.strictEqual(a[1],Math.max(...vals));
    const midpoint=(a[0]+a[1])/2;assert.strictEqual(app.helpers.continuousColor(midpoint,...a),app.helpers.continuousColor(midpoint,...b));
  });
  app.state.annotation='gene';app.state.overlayGene=gene;
  check('Actual gene-overlay values use the shared UMAP range',()=>{
    assert.deepEqual(plain(app.numericRange('scanpy')),plain(app.numericRange('seurat')));
    app.state.sharedScale=false;const left=app.numericRange('scanpy'),right=app.numericRange('seurat');assert(left&&right);
    app.state.sharedScale=true;
  });
  app.state.deShowNA=false;app.state.deFdr=.05;app.state.deLfc=0;app.state.deSearch='';
  const sample=[{gene:'Duplicated',gene_id:'id1',log2FoldChange:1,padj:.01},
    {gene:'Duplicated',gene_id:'id2',log2FoldChange:1,padj:null},{gene:'Other',gene_id:'id3',log2FoldChange:1,padj:.5}];
  check('DE filters preserve duplicate symbols and distinguish unavailable FDR',()=>{
    assert.equal(app.filteredDE(sample).length,1);app.state.deShowNA=true;
    const z=app.filteredDE(sample);assert.equal(z.length,2);assert.equal(z[1].padj,null);
  });
  const fit=app.active.seurat.fits.fits.find(f=>f.status==='TESTED_REVIEW_REQUIRED'&&f.endpoint==='PRIMARY'&&f.download_shard);
  const payload=await app.loadShard(fit.download_shard);
  check('Seven-column source download preserves exact bytes and symbols',()=>{
    assert.equal(hash(Buffer.from(payload.text,'utf8')),payload.sha256);
    assert.equal(payload.text.split('\n')[0],'gene\tbaseMean\tlog2FoldChange\tlfcSE\tstat\tpvalue\tpadj');
    assert.equal(payload.text.split('\n').filter(Boolean).slice(1).every(r=>r.split('\t').length===7),true);
  });
  app.state.deTargets.seurat='cluster:'+String(app.active.seurat.coords.nuclei[0].cluster_id);
  const saved=plain(app.exportState());await app.importState(saved);
  check('Imported settings restore exact selected DE target and custom edits',()=>{
    assert.equal(app.state.deTargets.seurat,saved.state.deTargets.seurat);
    assert.deepEqual(plain(app.state.labels),saved.state.labels);
  });
  const wrong={...saved,source_fingerprint:'DIFFERENT'};
  check('Foreign source presets are rejected',()=>assert.throws(()=>app.importState(wrong),/different source/));
  await app.resetState();
  check('Reset restores baseline without altering source data',()=>{
    assert.equal(Object.keys(app.state.exclusions).length,0);assert.equal(Object.keys(app.state.labels).length,0);
    assert.equal(app.active.scanpy.coords.nuclei.length,manifest.datasets.find(d=>d.id===app.active.scanpy.descriptor.id).n_nuclei);
  });
  for(const branch of ['original_full','paired_original','paired_cellbender']){
    for(const pop of ['allcells','Excitatory','Inhibitory','CA1_ProS']){
      await scope(pop,branch);
      check(branch+'/'+pop+' both tools bind exact completed totals and equal source statistics',()=>{
        for(const t of ['scanpy','seurat']){
          const f=app.resolveFit(t);assert.equal(f.status,'EXACT_MATCH');
          assert.equal(f.fit.branch,branch);assert.equal(f.fit.population,pop);
          assert.equal(f.fit.n_WT,6);assert.equal(f.fit.n_cKO,6);
          assert.equal(f.fit.coverage.length,12);
        }
        const right=new Map(app.active.seurat.deRows.map(r=>[r.gene_id,r]));
        assert.equal(app.active.scanpy.deRows.length,right.size);
        for(const a of app.active.scanpy.deRows){const b=right.get(a.gene_id);assert(b);assert.equal(a.gene,b.gene);
          for(const field of ['baseMean','log2FoldChange','lfcSE','stat','pvalue','padj'])assert.strictEqual(a[field],b[field]);}
      });
    }
  }
  assert.equal(hash(fs.readFileSync(path.join(source,'app.js'))),hash(appBytes),'Frontend changed during QA');
  const receipt={status:'PASS',validation_mode:'Pure JavaScript VM on actual produced data, without browser navigation',browser_interaction_inspection_performed:false,
    source_app_sha256:hash(appBytes),catalog_sha256:hash(fs.readFileSync(path.join(bundle,'DE_CATALOG.json'))),checks};
  fs.writeFileSync(path.join(base,'reviews/offline_dashboard_app_logic_acceptance.json'),JSON.stringify(receipt,null,2)+'\n');
  console.log(JSON.stringify({status:'PASS',checks:checks.length}));
})().catch(err=>{console.error(err);process.exitCode=1;});
