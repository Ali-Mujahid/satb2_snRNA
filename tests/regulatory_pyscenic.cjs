#!/usr/bin/env node
/* Synthetic positive-motif UI/export compatibility; no browser or models. */
'use strict';
const assert=require('node:assert/strict'),fs=require('node:fs'),path=require('node:path'),vm=require('node:vm');
const base=path.join(__dirname,'../site/regulatory'),source=fs.readFileSync(path.join(base,'app.js'),'utf8'),html=fs.readFileSync(path.join(base,'index.html'),'utf8');
const mice=['WT','cKO'].flatMap(genotype=>Array.from({length:6},(_,i)=>({id:'SYNTH_'+genotype+(i+1),genotype})));
const endpoints=[['SYNTH_CA1_DL','CA1'],['SYNTH_CA1_SL','CA1'],['SYNTH_IN_C1','IN'],['SYNTH_IN_C6','IN'],['SYNTH_IN_C8','IN']].map(([id,family])=>({id,family,label:id,n_nuclei:12,membership_sha256:'SYNTHETIC_MEMBERSHIP_ONLY'}));
const F=lin=>'pyscenic_motif_'+lin+'::'+lin+'::Satb2(+)',low='pyscenic_motif_IN::IN::SYNTH_LOW(+)';
const fixture={schema:'satb2-regulatory-atlas/v1',analysis_id:'SYNTHETIC_POSITIVE_MOTIF_UI_ONLY',source_fingerprint:'SYNTHETIC_NOT_A_SCIENTIFIC_RESULT',status:'complete',endpoints,samples:mice,
 methods:['CA1','IN'].map(lin=>({id:'pyscenic_motif_'+lin,family:'motif_regulon',status:'complete',endpoint_ids:endpoints.filter(e=>e.family===lin).map(e=>e.id),label:lin+' native positive motif regulons',interpretation:'Positive native motif-pruned candidate sets; arrows do not establish binding or causality.'})),
 features:['CA1','IN'].map(lin=>({id:F(lin),method_id:'pyscenic_motif_'+lin,label:'Satb2(+)',native_name:'Satb2(+)',source_id:lin+'::Satb2(+)',kind:'regulon',n_targets:lin==='CA1'?12:10,units:'Native AUCell AUC'})),
 values:[],coverage:[],effects:[],hubs:[],nodes:[],edges:[],embeddings:[],downloads:[],
 native_supplements:{Satb2_no_self_target:{primary_full_regulons_unchanged:true,coverage:endpoints.map(e=>({lineage:e.family,endpoint_id:e.id,regulon_id:e.family+'::Satb2(+)_NO_SELF_TARGET',n_targets:e.family==='CA1'?11:9,eligible:e.family==='CA1',reason:e.family==='CA1'?'':'FEWER_THAN_10_REMAINING_TARGETS'})),effects:[{lineage:'CA1',regulon_id:'CA1::Satb2(+)_NO_SELF_TARGET',effect:.03,padj:.4}]}},
 Satb2_status:{official_TF_list_membership:true,manual_additions:false,native_inferred_sets:['CA1::Satb2(+)','IN::Satb2(+)']},provenance:{synthetic:true,native_positive_motif_regulons:true,scientific_fits_in_UI:false}};
fixture.features.push({id:low,method_id:'pyscenic_motif_IN',label:'SYNTH_LOW(+)',native_name:'SYNTH_LOW(+)',source_id:'IN::SYNTH_LOW(+)',kind:'regulon',n_targets:3});
for(const [ei,e] of endpoints.entries()){
 const method='pyscenic_motif_'+e.family,feature=F(e.family);
 for(const [i,m] of mice.entries())fixture.values.push({method_id:method,endpoint_id:e.id,feature_id:feature,sample_id:m.id,value:.1+ei*.1+i*.01,eligible:true,coverage_n:e.family==='CA1'?12:10,n_nuclei:1});
 fixture.coverage.push({method_id:method,endpoint_id:e.id,feature_id:feature,eligible:true,measured_targets:e.family==='CA1'?12:10,required_targets:10});
 fixture.effects.push({method_id:method,endpoint_id:e.id,feature_id:feature,tested:true,effect:.03,statistic:1.2,pvalue:.2,padj:.4});
 const tf=e.id+'::native_Satb2',target=e.id+'::native_target';
 fixture.nodes.push({id:tf,method_id:method,endpoint_id:e.id,feature_id:feature,label:'Satb2',kind:'tf'},{id:target,method_id:method,endpoint_id:e.id,feature_id:feature,label:e.family==='CA1'?'SYNTH_Gria1':'SYNTH_Gad1',kind:'target'});
 fixture.edges.push({method_id:method,endpoint_id:e.id,feature_id:feature,source:tf,target:tf,weight:1.25,sign:'positive',directed:true,self_target:true,edge_basis:'native_positive_motif_pruned_gene2weight_not_binding_strength'},
                    {method_id:method,endpoint_id:e.id,feature_id:feature,source:tf,target,weight:2.75,sign:'positive',directed:true,self_target:false,edge_basis:'native_positive_motif_pruned_gene2weight_not_binding_strength'});
 if(e.family==='IN'){
  fixture.coverage.push({method_id:method,endpoint_id:e.id,feature_id:low,eligible:false,measured_targets:3,required_targets:10,reason:'FEWER_THAN_10_MEASURED_TARGETS'});
  for(const m of mice)fixture.values.push({method_id:method,endpoint_id:e.id,feature_id:low,sample_id:m.id,value:null,eligible:false,coverage_n:3,missing_reason:'FEWER_THAN_10_MEASURED_TARGETS'});
 }
}
const box={window:{setTimeout:()=>1,clearTimeout(){}}};vm.createContext(box);vm.runInContext(source,box);const A=box.window.SATB2Regulatory,d=A.validatePayload(fixture),original=JSON.stringify(d),checks=[];
const state=patch=>A.normalizeState({...A.DEFAULT,methodId:'pyscenic_motif_CA1',...patch},d);
function check(name,fn){fn();checks.push(name);console.log('PASS '+name);}
check('Native positive regulons retain separate CA1/IN definitions and exact matching endpoint scope',()=>{
 assert.notEqual(F('CA1'),F('IN'));
 const ca=A.mouseRows(d,state({featureId:F('CA1')})),inh=A.mouseRows(d,state({methodId:'pyscenic_motif_IN',featureId:F('IN')}));assert.equal(ca.length,24);assert.equal(inh.length,36);
 assert(ca.every(r=>r.endpoint_id.startsWith('SYNTH_CA1'))&&inh.every(r=>r.endpoint_id.startsWith('SYNTH_IN')));
 const wrong=state({featureId:F('IN'),endpointId:'SYNTH_IN_C1'});assert.equal(wrong.featureId,'');assert.equal(wrong.endpointId,'');
 const changed=JSON.parse(JSON.stringify(fixture));changed.values[0].feature_id=F('IN');assert.throws(()=>A.validatePayload(changed),/different method/);
});
check('All12-mouse display-z reference remains fixed when genotype display is restricted',()=>{
 const all=A.mouseRows(d,state({featureId:F('CA1'),scoreMode:'z'})),wt=A.mouseRows(d,state({featureId:F('CA1'),scoreMode:'z',genotype:'WT'}));assert.equal(wt.length,12);
 for(const r of wt){const x=all.find(a=>a.endpoint_id===r.endpoint_id&&a.sample_id===r.sample_id);assert.equal(r.display_reference_n,12);assert.equal(r.display_value,x.display_value);assert.equal(r.display_reference_mean,x.display_reference_mean);}
 const effects=A.buildView(d,state({featureId:F('CA1'),view:'effects',genotype:'WT'}));assert.equal(JSON.stringify(effects.rows),JSON.stringify(A.buildView(d,state({featureId:F('CA1'),view:'effects',genotype:'cKO'})).rows));
});
check('Positive motif edges retain native weights, candidate direction and explicit selfgene flags',()=>{
 const s=state({featureId:F('CA1'),endpointId:'SYNTH_CA1_DL',view:'network'}),g=A.graph(d,s),view=A.buildView(d,s);assert.equal(g.edges.length,2);assert(g.edges.every(e=>e.directed&&e.sign==='positive'));assert(g.edges.some(e=>e.self_target&&e.source===e.target));assert.deepEqual(Array.from(g.edges.map(e=>e.weight)).sort(),[1.25,2.75]);assert.match(view.svg,/reg-arrow/);assert.match(view.svg,/not_binding_strength/);
 const styled=A.buildView(d,state({...s,edgeWidth:12,arrowSize:20,curvature:.7,arcWidth:20,edgeColorMode:'custom',edgeColor:'#102030'}));assert.equal(JSON.stringify(styled.rows),JSON.stringify(view.rows));assert.match(styled.svg,/#102030/);assert.match(A.csv(view.rows),/self_target/);
});
const empty=JSON.parse(JSON.stringify(fixture));for(const m of empty.methods)if(m.id==='pyscenic_motif_IN'){m.status='not_available';m.interpretation='Native analysis completed; no eligible positive motif-supported regulons. No activity score is available.';m.resources={native_lineage_status:{status:'EMPTY_NO_SUPPORTED_REGULONS'}};}
for(const field of ['features','values','coverage','effects','hubs','nodes','edges'])empty[field]=empty[field].filter(r=>r.method_id!=='pyscenic_motif_IN');
check('Completed-empty lineage and low-target sets remain unavailable without zero or RNA substitution',()=>{
 const dd=A.validatePayload(empty),s=A.normalizeState({...A.DEFAULT,methodId:'pyscenic_motif_IN',view:'heatmap'},dd),view=A.buildView(dd,s);assert.equal(view.svg,null);assert.equal(view.rows.length,0);assert.match(view.note,/Native analysis completed; no eligible/);assert.doesNotMatch(view.note,/resource-blocked|awaiting/i);
 const lowRows=A.mouseRows(d,state({methodId:'pyscenic_motif_IN',featureId:low}));assert.equal(lowRows.length,36);assert(lowRows.every(r=>r.value===null&&r.display_value===null&&r.eligible===false));
 const umap=A.buildView(d,state({methodId:'pyscenic_motif_IN',featureId:F('IN'),view:'umap'}));assert.equal(umap.svg,null);assert.match(umap.note,/not substituted/);
});
check('Satb2 no-self IDs and sensitivity remain separate from primary motif definitions and curated shortcut',()=>{
 assert.equal(JSON.stringify(d.native_supplements),JSON.stringify(fixture.native_supplements));assert.equal(d.native_supplements.Satb2_no_self_target.coverage.length,5);assert(d.native_supplements.Satb2_no_self_target.coverage.every(r=>r.regulon_id.endsWith('_NO_SELF_TARGET')));assert(d.features.every(f=>!f.id.endsWith('_NO_SELF_TARGET')));assert.equal(A.satb2Focus(d,state({})),null,'Curated Satb2 shortcut must not substitute motif activity');
 const absent=JSON.parse(JSON.stringify(fixture));for(const name of ['features','values','coverage','effects','nodes','edges'])absent[name]=absent[name].filter(r=>r.method_id!=='pyscenic_motif_CA1');absent.methods.find(m=>m.id==='pyscenic_motif_CA1').status='not_available';absent.Satb2_status.native_inferred_sets=['IN::Satb2(+)'];
 for(const r of absent.native_supplements.Satb2_no_self_target.coverage)if(r.lineage==='CA1'){r.eligible=false;r.reason='NO_NATIVE_SATB2_REGULON';}
 const missing=A.validatePayload(absent);assert.equal(missing.Satb2_status.official_TF_list_membership,true);assert(!missing.features.some(f=>f.id===F('CA1')),'Official TF-list membership must not create a motif regulon');assert.equal(A.satb2Focus(missing,state({})),null);
 const curated=JSON.parse(JSON.stringify(fixture));curated.methods.push({id:'collectri_ulm',family:'curated_tf_activity',status:'complete',endpoint_ids:endpoints.map(e=>e.id)});curated.features.push({id:'collectri_ulm::Satb2',method_id:'collectri_ulm',label:'Satb2'});
 for(const e of endpoints)for(const m of mice)curated.values.push({method_id:'collectri_ulm',endpoint_id:e.id,feature_id:'collectri_ulm::Satb2',sample_id:m.id,value:1.2,eligible:true});
 const cc=A.validatePayload(curated),focus=A.satb2Focus(cc,state({methodId:'pyscenic_motif_IN',featureId:F('IN')}));assert.equal(focus.methodId,'collectri_ulm');assert.equal(focus.featureId,'collectri_ulm::Satb2');assert.equal(A.buildView(cc,focus).rows.length,60);
 curated.methods.find(m=>m.id==='collectri_ulm').status='not_available';assert.equal(A.satb2Focus(A.validatePayload(curated),state({})),null);assert.equal(JSON.stringify(d),original);
});

// The same DOM/image/PDF wrapper interface as the existing regulatory tests.
const elements=new Map(),downloads=[],pdfs=[],blobs=new Map();let urlCount=0,decodeWait=null;
class El{constructor(id,tag='div',attrs=''){this.id=id;this.tagName=tag.toUpperCase();this.type=(attrs.match(/type="([^"]+)"/)||[])[1]||'';this.value='';this.checked=false;this.hidden=false;this.disabled=false;this.dataset={};this.listeners={};this.classList={toggle(){}};this._html='';}set innerHTML(v){this._html=v;}get innerHTML(){return this._html;}addEventListener(k,fn){this.listeners[k]=fn;}append(el){if(this.tagName==='HEAD')el.onerror();}click(){if(this.tagName==='A')downloads.push({name:this.download,blob:blobs.get(this.href)});else if(this.onclick)this.onclick();}remove(){}}
for(const m of html.matchAll(/<([\w-]+)\b([^>]*\bid="([^"]+)"[^>]*)>/g))elements.set(m[3],new El(m[3],m[1],m[2]));const tabs=[...html.matchAll(/data-view="([^"]+)"/g)].map(m=>{const e=new El('tab');e.dataset.view=m[1];return e;});
const document={readyState:'loading',addEventListener(){},getElementById:id=>elements.get(id),querySelectorAll:()=>tabs,head:new El('head','head'),body:new El('body','body'),createElement:tag=>{const el=new El('',tag);if(tag==='canvas'){el.getContext=()=>({fillRect(){},drawImage(){}});el.toDataURL=()=> 'data:image/png;base64,AA==';el.toBlob=fn=>fn(new Blob(['SYNTHETIC PNG'],{type:'image/png'}));}return el;}};
class ImageMock{async decode(){if(decodeWait)await decodeWait;}}
class Parser{parseFromString(svg){return {documentElement:{getAttribute:k=>(svg.match(new RegExp(k+'="([^"]+)"'))||[])[1]}};}}
class Pdf{constructor(o){this.options=o;pdfs.push(this);}setProperties(p){this.properties=p;}addImage(){}save(name){this.name=name;}}
const store=new Map(),window={document,SATB2_REGULATORY_DATA:fixture,localStorage:{getItem:k=>store.get(k)||null,setItem:(k,v)=>store.set(k,v)},setTimeout:()=>1,clearTimeout(){},jspdf:{jsPDF:Pdf}},ui={window,Blob,Image:ImageMock,DOMParser:Parser,URL:{createObjectURL:b=>{const url='blob:synthetic/'+ ++urlCount;blobs.set(url,b);return url;},revokeObjectURL(){}},TextDecoder,Uint8Array};
vm.runInNewContext(source,ui);const U=window.SATB2Regulatory;
(async()=>{
 await U.start();assert.equal(U.state.methodId,'pyscenic_motif_CA1');assert.equal(elements.get('focus-satb2').disabled,true);assert(U.plot.svg);
 U.change('featureId',F('CA1'));U.change('endpointId','SYNTH_CA1_DL');U.change('view','network');U.render();const originalPlot=U.plot;assert.equal(originalPlot.rows.length,2);
 await U.exportPlot('svg');assert.equal(await downloads.at(-1).blob.text(),originalPlot.svg);await U.exportPlot('csv');assert.match(await downloads.at(-1).blob.text(),/self_target/);await U.exportPlot('json');const record=JSON.parse(await downloads.at(-1).blob.text());assert.equal(record.method.family,'motif_regulon');assert(record.rows.some(r=>r.self_target));await U.exportPlot('png');assert.equal(downloads.at(-1).blob.type,'image/png');
 let release;decodeWait=new Promise(r=>{release=r;});const pending=U.exportPlot('pdf');U.change('methodId','pyscenic_motif_IN');U.change('featureId',F('IN'));U.change('view','replicates');U.render();release();decodeWait=null;await pending;
 assert.equal(pdfs.length,1);assert(pdfs[0].name.includes('pyscenic_motif_CA1')&&pdfs[0].name.includes('network'));assert.equal(JSON.parse(pdfs[0].properties.subject).method,'pyscenic_motif_CA1');assert.equal(JSON.stringify(U.data),JSON.stringify(d));checks.push('Generic PDF/SVG/PNG/CSV/JSON exports preserve motif lineage, native selfgene flags, weights and captured source during navigation');
 U.setData(empty,{ignoreSaved:true});U.change('methodId','pyscenic_motif_IN');U.change('view','heatmap');U.render();assert.equal(U.plot.svg,null);assert.equal(elements.get('plot-pdf').disabled,true);assert.equal(elements.get('plot-svg').disabled,true);assert.equal(elements.get('plot-png').disabled,true);await U.exportPlot('pdf');assert.equal(pdfs.length,1);assert.match(elements.get('view-note').textContent,/no eligible/);checks.push('Unavailable empty motif branch clears stale plots and disables PDF/SVG/PNG export without substituting another method');
 console.log(JSON.stringify({status:'PASS',checks:checks.length,names:checks,synthetic_fixture_only:true,scientific_analysis_performed:false,browser_launched:false},null,2));
})().catch(e=>{console.error(e);process.exitCode=1;});
