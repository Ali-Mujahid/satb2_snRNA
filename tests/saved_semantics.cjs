/* Actual produced-data/state checks in a JavaScript VM; no browser navigation. */
'use strict';
const fs=require('fs'),path=require('path'),vm=require('vm'),assert=require('assert'),zlib=require('zlib'),crypto=require('crypto');
const bundle=(process.env.SATB2_SITE_DIR||path.resolve(__dirname,'../site'));
const source=(process.env.SATB2_SITE_DIR||path.resolve(__dirname,'../site'));
const checks=[];
function check(name,fn){fn();checks.push(name);console.log(JSON.stringify({check:name,status:'PASS'}));}
function readJson(name){return JSON.parse(fs.readFileSync(path.join(bundle,name),'utf8'));}
function hash(bytes){return crypto.createHash('sha256').update(bytes).digest('hex');}
function readPayload(id,desc){
  assert(desc&&desc.path,'Missing descriptor: '+id);
  const remote='https://raw.githubusercontent.com/Ali-Mujahid/satb2_snRNA/atlas-v2.1.0/data_downloads/';
  const filename=desc.path.startsWith(remote)?path.join(path.dirname(bundle),'data_downloads',desc.path.slice(remote.length)):path.join(bundle,desc.path);
  const bytes=fs.readFileSync(filename);assert.equal(hash(bytes),desc.sha256,'Binary shard integrity '+id);
  return {encoding:'gzip-bytes',data:new Uint8Array(bytes),payload_sha256:desc.payload_sha256};
}

const context={console,Uint8Array,Uint16Array,Uint32Array,Int32Array,Float32Array,Float64Array,ArrayBuffer,TextEncoder,TextDecoder,atob,btoa,Blob,setTimeout,clearTimeout,performance,
  crypto:crypto.webcrypto,localStorage:{getItem(){return null;},setItem(){}},SATB2_NO_AUTO_INIT:true};
context.window=context;context.self=context;context.globalThis=context;vm.createContext(context);
const appBytes=fs.readFileSync(path.join(source,'app.js'));
vm.runInContext(fs.readFileSync(path.join(source,'vendor/fflate.min.js'),'utf8'),context,{filename:'fflate.min.js'});
vm.runInContext(fs.readFileSync(path.join(source,'v2_features.js'),'utf8'),context,{filename:'v2_features.js'});
vm.runInContext(appBytes.toString('utf8'),context,{filename:'app.js'});
for(const name of ['data/mouse_catalog.js','data/annotation_suggestions.js']){const p=path.join(bundle,name);if(fs.existsSync(p))vm.runInContext(fs.readFileSync(p,'utf8'),context,{filename:name});}
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
  // This fixture loads saved cluster summaries only; nucleus-vector behavior
  // and its measured-expression default are covered by plot_controls_v31.cjs.
  app.state.annotation='gene';app.state.overlayGene=gene;app.state.expressionResolution='cluster';
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

  await scope('CA1_ProS','original_full','markers');await app.ensureMarkers();
  const sourceSummaryBytes=JSON.stringify(app.active.scanpy.markers),id0=manifest.genes[0].gene_id;
  app.state.markers=manifest.genes.map(g=>g.gene_id);
  for(const count of [60,90,120,180]){app.state.genesPerPage=count;check('V2 '+count+' markers per page preserves every selected gene',()=>{assert.equal(app.pageGenes(),count);assert.equal(app.state.markers.length,manifest.genes.length);assert.equal(app.dotRows('scanpy',app.state.markers.slice(0,count)).length,app.active.scanpy.coords.clusters.filter(c=>app.active.scanpy.coords.nuclei.some(n=>String(n.cluster_id)===String(c.cluster_id))).length*count);});}
  app.state.markerMetric='scaled_expression_z';app.state.dotRangeMode='manual';app.state.dotMin=-1;app.state.dotMax=1;
  check('V2 manual z range changes display limits only',()=>{assert.deepEqual(plain(app.markerRange('scanpy',[id0])),[-1,1]);assert.equal(JSON.stringify(app.active.scanpy.markers),sourceSummaryBytes);});
  app.state.dotMin=2;app.state.dotMax=1;check('V2 invalid manual ranges remain explicitly invalid',()=>assert.equal(app.markerRange('scanpy',[id0]),null));
  app.state.dotRangeMode='source';app.state.dotArea=200;app.state.negativeArea=.1;app.state.negativeAlpha=.05;
  check('V2 configured negative glyph policy is explicit and numeric-preserving',()=>{const g=app.drawDotValues({scaled_expression_z:-1,fraction_detected:.8});assert.equal(g.value,-1);assert.equal(g.area,16);assert.equal(g.alpha,.05);});
  app.state.labels[clusterKey('scanpy',app.active.scanpy.coords.nuclei[0].cluster_id)]='Display label <script> harmless';
  check('V2 label export preserves source identities and marks display-only annotation',()=>{const q=app.exportLabels();assert(q.text.includes('Display label <script> harmless'));assert(q.text.includes('scientific_annotation_modified'));assert(q.text.includes('false'));});
  app.state.spacing=.5;app.state.umapLabelMode='both';app.state.genesPerPage=120;const preset=plain(app.exportState());await app.importState(preset);
  check('V2 presets restore60+ gene/spacing/label controls with schema2',()=>{assert.equal(preset.schema_version,'satb2-offline-2');assert.equal(app.state.genesPerPage,120);assert.equal(app.state.spacing,.5);assert.equal(app.state.umapLabelMode,'both');});
  await scope('Inhibitory','original_full','umap');const sourceCoordinates=JSON.stringify(app.active.scanpy.coords.nuclei),original=context.SATB2_V2.transformCoordinates(app.active.scanpy.coords.nuclei,'harmony',1),altered=context.SATB2_V2.transformCoordinates(app.active.scanpy.coords.nuclei,'harmony',.5);
  check('Actual IN spacing preserves original defaults and tags changed display coordinates',()=>{assert(original.every((p,i)=>p.x===app.active.scanpy.coords.nuclei[i].x_harmony&&!p.display_only));assert(altered.some((p,i)=>p.x!==original[i].x));assert(altered.every(p=>p.source_coordinates_modified===false));assert.equal(JSON.stringify(app.active.scanpy.coords.nuclei),sourceCoordinates);});
  app.state.umapCondition='WT';const targetBefore=app.targetMembers('scanpy').length;check('UMAP-only condition visibility does not modify inferential membership',()=>{assert(app.umapIndices('scanpy').every(i=>app.active.scanpy.coords.nuclei[i].condition==='WT'));assert.equal(app.targetMembers('scanpy').length,targetBefore);});app.state.umapCondition='ALL';
  const mouse=context.SATB2_MOUSE_CATALOG;if(!mouse||mouse.status!=='COMPLETE')throw Error('Closed mouse catalog required for final v2 integration QA');
  check('Full-gene boxplot catalog is source-bound and complete',()=>{assert.equal(mouse.genes.length,32285);assert.equal(Object.keys(mouse.fits).length,1060);assert.equal(mouse.de_catalog_sha256,hash(fs.readFileSync(path.join(bundle,'DE_CATALOG.json'))));});
  await scope('Inhibitory','original_full','boxplots');await app.ensureFits();await app.ensureMouseData();app.state.boxGenes=mouse.genes.filter(g=>['Cck','Gad1','Satb2'].includes(g.gene_symbol)).map(g=>g.gene_id);
  for(const tool of ['scanpy','seurat']){
    const exact=app.resolveFit(tool),meta=mouse.fits[exact.fit.id],source=app.active[tool].mouseData;
    app.state.boxTransform='normalized';const points=app.boxValues(tool),lookup=new Map(source.rows);
    check(tool+' exact mouse boxplot normalized values equal source integer/SF',()=>{assert.equal(points.length,app.state.boxGenes.length*12);for(const p of points){const gi=mouse.genes.findIndex(g=>g.gene_id===p.gene_id),i=meta.sample_order.indexOf(p.sample),count=lookup.get(gi)?.[i]??0;assert.equal(p.source_count,count);assert.equal(p.display_value,meta.size_factors[i]>0?count/meta.size_factors[i]:null);assert.equal(p.fit_id,exact.fit.id);assert.equal(p.biological_unit,'independent_mouse');assert.equal(p.statistical_test_performed,false);}});
    app.state.boxTransform='raw';check(tool+' raw boxplot values preserve all source zeros exactly',()=>{for(const p of app.boxValues(tool))assert.equal(p.display_value,p.source_count);});
    app.state.boxTransform='log1p_normalized';check(tool+' log boxplots transform finite source-normalized values only',()=>{for(const p of app.boxValues(tool))assert.equal(p.display_value,p.normalized_value===null?null:Math.log1p(p.normalized_value));});
    const zeroIndex=mouse.genes.findIndex((g,i)=>!lookup.has(i));assert(zeroIndex>=0);app.state.boxGenes=[mouse.genes[zeroIndex].gene_id];app.state.boxTransform='raw';check(tool+' explicit all-zero gene stays0, not missing',()=>assert(app.boxValues(tool).every(p=>p.source_count===0&&p.display_value===0)));
    app.state.boxGenes=mouse.genes.filter(g=>['Cck','Gad1','Satb2'].includes(g.gene_symbol)).map(g=>g.gene_id);
  }
  await scope('CA1_ProS','original_full','boxplots');await app.ensureFits();
  for(const tool of ['scanpy','seurat']){
    const rare=app.active[tool].fits.fits.find(f=>f.endpoint==='PRIMARY'&&f.status==='NOT_TESTED'&&f.analysis_kind==='cluster'&&f.members.length>0);assert(rare,'Need actual source unsupported rare group');app.state.deTargets[tool]='fit:'+rare.id;await app.selectionChanged();app.state.boxGenes=[mouse.genes.find(g=>g.gene_symbol==='Cck').gene_id];app.state.boxTransform='normalized';
    check(tool+' unsupported source SF stays unavailable, never borrowed',()=>{const points=app.boxValues(tool);assert.equal(points.length,12);assert(points.every(p=>p.normalized_value===null&&p.display_value===null));});app.state.boxTransform='raw';app.state.boxShowIneligible=true;
    check(tool+' unsupported population still exposes exact raw all12 counts without inference',()=>{const points=app.boxValues(tool);assert.equal(points.length,12);assert(points.every(p=>p.display_value===p.source_count));assert(points.some(p=>p.source_count>0)||app.active[tool].mouseData.omitted_rows==='EXPLICIT_ALL_ZERO');});
    app.state.deTargets[tool]='total';await app.selectionChanged();const before=app.resolveFit(tool).fit.id,clusters=[...new Set(app.active[tool].coords.nuclei.map(n=>String(n.cluster_id)))];app.state.exclusions[key(tool)]=[clusters[0]];
    check(tool+' edited membership immediately rejects stale mouse boxplot data',()=>{assert(!app.resolveFit(tool).fit||app.resolveFit(tool).fit.id!==before);assert.equal(app.boxValues(tool).length,0);});await app.selectionChanged();if(app.resolveFit(tool).status!=='EXACT_MATCH'&&app.resolveFit(tool).status!=='NOT_TESTED')assert.equal(app.active[tool].mouseData,null);
    app.state.exclusions[key(tool)]=[];app.state.deTargets[tool]='total';await app.selectionChanged();
  }
  assert.equal(hash(fs.readFileSync(path.join(source,'app.js'))),hash(appBytes),'Frontend changed during QA');

  const receipt={status:'PASS',validation_mode:'Pure JavaScript VM on actual produced data, without browser navigation',browser_interaction_inspection_performed:false,
    source_app_sha256:hash(appBytes),catalog_sha256:hash(fs.readFileSync(path.join(bundle,'DE_CATALOG.json'))),mouse_catalog_sha256:hash(fs.readFileSync(path.join(bundle,'data/MOUSE_CATALOG.json'))),source_helper_sha256:hash(fs.readFileSync(path.join(source,'v2_features.js'))),checks};
  fs.writeFileSync(path.join(__dirname,'hosted_binary_semantics_acceptance.json'),JSON.stringify(receipt,null,2)+'\n');
  console.log(JSON.stringify({status:'PASS',checks:checks.length}));
})().catch(err=>{console.error(err);process.exitCode=1;});
