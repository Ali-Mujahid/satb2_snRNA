#!/usr/bin/env node
/* Read-only actual-payload UI QA: no browser, native analysis, fit or publish. */
'use strict';
const fs=require('node:fs'),vm=require('node:vm'),zlib=require('node:zlib'),crypto=require('node:crypto'),assert=require('node:assert/strict'),{performance}=require('node:perf_hooks');
const [appPath,dataPath,outPath]=process.argv.slice(2);
assert(appPath&&dataPath&&outPath,'Usage: node smoke.cjs APP_JS STAGED_DATA_JS NEW_REPORT_JSON');
assert(!fs.existsSync(outPath),'Refusing to replace a QA report');
const timings=[],time=(label,fn)=>{const t=performance.now(),value=fn(),ms=performance.now()-t;timings.push({label,milliseconds:Math.round(ms*100)/100});assert(ms<5000,label+' exceeded bounded 5-second UI budget: '+ms);return value;};
const hash=x=>crypto.createHash('sha256').update(x).digest('hex'),appBytes=fs.readFileSync(appPath),dataBytes=fs.readFileSync(dataPath);
const appBox={window:{setTimeout:()=>1,clearTimeout(){}}};vm.createContext(appBox);vm.runInContext(appBytes.toString(),appBox,{timeout:5000});
const A=appBox.window.SATB2Regulatory,dataBox={window:{}};vm.runInNewContext(dataBytes.toString(),dataBox,{timeout:5000});
const original=time('decode_production_payload',()=>JSON.parse(zlib.gunzipSync(Buffer.from(dataBox.window.SATB2_REGULATORY_DATA_GZIP,'base64'))));
assert(!original.provenance?.synthetic&&!original.analysis_id.includes('SYNTHETIC'));
const d=time('validatePayload_production',()=>A.validatePayload(original)),before=hash(JSON.stringify(d));
const state=patch=>A.normalizeState({...A.DEFAULT,...patch},d),ready=d.methods.filter(m=>['complete','completed','partial'].includes(m.status));
assert.equal(ready.filter(m=>m.family==='hdWGCNA').length,4);assert.equal(ready.filter(m=>m.family==='curated_tf_activity').length,1);
const resultRows=[],views=['heatmap','dotplot','replicates','network','hubs','effects','coverage'];
for(const method of ready){
 const feature=d.features.find(f=>f.method_id===method.id&&(method.family!=='curated_tf_activity'||f.label==='Satb2')),base={methodId:method.id,featureId:feature.id};
 for(const view of views){
  const s=state({...base,view}),p=time(method.id+'_'+view,()=>A.buildView(d,s));
  assert(Array.isArray(p.rows));assert(!p.rows.some(r=>r.method_id&&r.method_id!==method.id));
  if(['heatmap','dotplot','replicates','network'].includes(view)){assert(p.svg,method.id+' '+view);assert(!/NaN|Infinity|undefined/.test(p.svg));}
  if(view==='network'){
   const g=A.graph(d,s);assert(g.nodes.length<=s.limit&&g.edges.length<=s.edgeLimit);assert(g.edges.every(e=>e.directed===(method.family==='curated_tf_activity')));
   assert.equal(p.svg.includes('class="reg-arrow"'),method.family==='curated_tf_activity'&&g.edges.length>0);
  }
  if(view==='effects'){
   const wt=A.buildView(d,state({...base,view,genotype:'WT'})),ko=A.buildView(d,state({...base,view,genotype:'cKO'}));
   assert.equal(JSON.stringify(wt.rows),JSON.stringify(ko.rows),'Display genotype must not change saved effects');
  }
  resultRows.push({method:method.id,view,rows:p.rows.length,svg_bytes:p.svg?.length||0});
 }
 const all=time(method.id+'_all_feature_heatmap',()=>A.buildView(d,state({methodId:method.id,view:'heatmap',scoreMode:'z'})));
 assert(all.svg);assert(!/NaN|Infinity|undefined/.test(all.svg));
 const graph=time(method.id+'_all_feature_network',()=>A.buildView(d,state({methodId:method.id,view:'network',limit:60,edgeLimit:200})));
 assert(graph.svg);
}
const focus=A.satb2Focus(d,state({search:'irrelevant',genotype:'WT',endpointFamily:'IN'}));
assert(focus&&focus.methodId==='collectri_ulm'&&focus.featureId==='collectri_ulm::Satb2');
assert.equal(focus.genotype,'all');assert.equal(focus.endpointFamily,'all');assert.equal(focus.endpointId,'');assert.equal(focus.search,'');assert.equal(focus.scoreMode,'raw');
const satb2=time('SATB2_one_click_all_five_endpoints',()=>A.buildView(d,focus));assert(satb2.svg);assert.equal(satb2.rows.length,60);assert.equal(new Set(satb2.rows.map(r=>r.endpoint_id)).size,5);assert(satb2.rows.every(r=>r.eligible&&Number.isFinite(r.value)));
const zAll=A.mouseRows(d,state({...focus,scoreMode:'z'})),zWt=A.mouseRows(d,state({...focus,scoreMode:'z',genotype:'WT'}));
assert.equal(zWt.length,30);for(const r of zWt){const all=zAll.find(x=>x.endpoint_id===r.endpoint_id&&x.sample_id===r.sample_id);assert.equal(r.display_value,all.display_value);assert.equal(r.display_reference_n,12);}
const satb2Net=time('SATB2_prior_network_60_nodes_arc',()=>A.buildView(d,state({...focus,view:'network',limit:60,networkLayout:'arc'})));assert(satb2Net.svg.includes('class="reg-arrow"'));
const motif=d.methods.find(m=>m.family==='motif_regulon');assert.equal(motif.status,'not_available');const blocked=A.buildView(d,state({methodId:motif.id,view:'heatmap'}));assert.equal(blocked.svg,null);assert.equal(blocked.rows.length,0);assert.match(blocked.note,/unavailable/);
for(const method of ready){const p=A.buildView(d,state({methodId:method.id,view:'umap'}));assert.equal(p.svg,null);assert.equal(p.rows.length,0);assert.match(p.note,/not substituted/);}
assert.equal(hash(JSON.stringify(d)),before,'Pure UI views changed source payload');assert.equal(hash(fs.readFileSync(dataPath)),hash(dataBytes),'Staged data bytes changed');
const report={schema:'satb2-regulatory-production-ui-smoke/v1',status:'PASS',created_at_utc:new Date().toISOString(),analysis_id:d.analysis_id,source_fingerprint:d.source_fingerprint,app:{path:appPath,sha256:hash(appBytes)},data:{path:dataPath,sha256:hash(dataBytes)},counts:Object.fromEntries(['values','features','edges','nodes','effects','coverage','hubs'].map(k=>[k,d[k].length])),timings,maximum_timed_call_ms:Math.max(...timings.map(x=>x.milliseconds)),rendered_views:resultRows,satb2_rows:satb2.rows.length,satb2_endpoints:5,mouse_z_reference_preserved_across_genotype_filters:true,saved_effects_unchanged_by_display_filter:true,HD_undirected_TF_prior_directed:true,blocked_motif_and_activity_umap_unavailable:true,source_data_immutable:true,browser_test_performed:false,scientific_models_run:false,published:false};
fs.writeFileSync(outPath,JSON.stringify(report,null,2)+'\n');console.log(JSON.stringify({status:report.status,calls:timings.length,max_call_ms:report.maximum_timed_call_ms,rendered_views:resultRows.length,report:outPath},null,2));
