#!/usr/bin/env node
/* Synthetic source-preservation and UI/export tests. No browser or analysis is run. */
'use strict';
const assert=require('node:assert/strict'),fs=require('node:fs'),path=require('node:path'),vm=require('node:vm');
const fixture=require('./fixtures/regulatory_synthetic.cjs'),base=path.join(__dirname,'../site/regulatory'),source=fs.readFileSync(path.join(base,'app.js'),'utf8'),html=fs.readFileSync(path.join(base,'index.html'),'utf8');
const sandbox={window:{setTimeout:()=>1,clearTimeout(){}}};vm.createContext(sandbox);vm.runInContext(source,sandbox);const A=sandbox.window.SATB2Regulatory,d=A.validatePayload(fixture),before=JSON.stringify(d),s=patch=>A.normalizeState({...A.DEFAULT,methodId:'SYNTH_HD_CA1',...patch},d),checks=[];
function check(name,fn){fn();checks.push(name);console.log('PASS '+name);}
check('No synthetic results are bundled in site/regulatory; real saved data remains optional',()=>{for(const p of fs.readdirSync(base))assert(!/synthetic|fixture/i.test(p));if(fs.existsSync(path.join(base,'data.js'))){const c={window:{}};vm.runInNewContext(fs.readFileSync(path.join(base,'data.js'),'utf8'),c);const payload=c.window.SATB2_REGULATORY_DATA||JSON.parse(require('node:zlib').gunzipSync(Buffer.from(c.window.SATB2_REGULATORY_DATA_GZIP,'base64')));assert.notEqual(payload.analysis_id,fixture.analysis_id);assert.notEqual(payload.provenance?.synthetic,true);}assert.match(html,/Awaiting verified analysis outputs/);assert.match(source,/SATB2_REGULATORY_DATA_GZIP/);assert.match(source,/script\.src='data\.js'/);});
check('Exact references, independent mouse uniqueness and eligibility are validated',()=>{
 const invalid=patch=>{const x=JSON.parse(JSON.stringify(fixture));patch(x);assert.throws(()=>A.validatePayload(x));};invalid(x=>x.values.push({...x.values[0]}));invalid(x=>x.values[0].sample_id='UNKNOWN');invalid(x=>x.values[0].value=NaN);invalid(x=>delete x.values[0].eligible);invalid(x=>x.values[0].feature_id='SYNTH_MOD_IN');invalid(x=>x.edges[0].directed=true);invalid(x=>x.edges[0].target='SYNTH_TF_NODE');invalid(x=>x.embeddings[0].activity.SYNTH_TF_A.pop());invalid(x=>x.embeddings[0].coordinates_kind='artificial_schematic');
});
check('Mouse z reference includes all eligible mice across genotypes and excludes ineligible outliers',()=>{
 const rows=A.mouseRows(d,s({featureId:'SYNTH_MOD_A',scoreMode:'z'}));assert.equal(rows.length,4);assert.equal(rows[0].display_reference_n,3);assert.equal(rows[0].display_reference_mean,2);assert.equal(rows[0].display_reference_sample_sd,1);assert.equal(rows[0].display_value,-1);assert.equal(rows[2].display_value,1);assert.equal(rows[3].display_value,null);
 const wt=A.mouseRows(d,s({featureId:'SYNTH_MOD_A',scoreMode:'z',genotype:'WT'}));assert.equal(wt.length,2);assert.equal(wt[0].display_reference_mean,2);assert.equal(wt[0].display_value,-1);
 const include=A.mouseRows(d,s({featureId:'SYNTH_MOD_A',scoreMode:'z',showIneligible:true}));assert.equal(include[3].display_value,98);assert.equal(include[3].display_reference_mean,2);
 const constant=A.mouseRows(d,s({featureId:'SYNTH_MOD_CONSTANT',scoreMode:'z'}));assert(constant.every(r=>r.display_value===null));
});
check('Genotype means exclude ineligible points and missing rows never become zero',()=>{
 const st=s({featureId:'SYNTH_MOD_A',heatColumns:'genotype',showIneligible:true}),m=A.matrix(d,st,A.mouseRows(d,st));assert.equal(m.cells.find(c=>c.col==='cKO').value,3);assert.equal(m.cells.find(c=>c.col==='cKO').n,1);
 const m2=A.matrix(d,s({featureId:'SYNTH_MOD_CONSTANT'}),A.mouseRows(d,s({featureId:'SYNTH_MOD_CONSTANT'})));assert.equal(m2.cells.find(c=>c.col==='SYNTH_K2').value,null);
});
check('All plot types use native rows; undirected hdWGCNA never gains TF arrows',()=>{
 for(const view of ['heatmap','dotplot','replicates','network']){const p=A.buildView(d,s({view}));assert(p.svg,view);assert.doesNotMatch(p.svg,/NaN|Infinity|undefined/);}
 const module=A.buildView(d,s({view:'network',arrowSize:24}));assert.doesNotMatch(module.svg,/class="reg-arrow"/);assert.match(module.svg,/weighted coexpression/);
 const tf=A.buildView(d,s({methodId:'SYNTH_TF',view:'network',featureId:'SYNTH_TF_A'}));assert.match(tf.svg,/class="reg-arrow"/);assert.match(tf.svg,/curated prior/);
 const limited=A.buildView(d,s({view:'network',limit:2,edgeLimit:1}));assert.equal(limited.rows.length,2,'CSV rows retain every filtered native edge.');assert.equal(A.graph(d,s({limit:2,edgeLimit:1})).edges.length,1);
 const targetSearch=A.graph(d,s({search:'SYNTHETIC GENE B'}));assert.equal(targetSearch.edges.length,1);assert.equal(targetSearch.nodes.length,2,'A matching edge retains its connected partner even when that partner does not match the search.');
});
check('Unavailable method or activity feature is never replaced by expression or another score',()=>{
 assert.equal(A.buildView(d,s({methodId:'SYNTH_MOTIF_PENDING',view:'heatmap'})).svg,null);
 assert.equal(A.buildView(d,s({view:'umap'})).svg,null);
 const actual=A.buildView(d,s({methodId:'SYNTH_TF',featureId:'SYNTH_TF_A',view:'umap'}));assert(actual.svg);assert.equal(actual.rows.length,3);assert.equal(actual.rows[2].display_value,null);assert.match(actual.svg,/SYNTH_N3/);
 const missing=A.buildView(d,s({methodId:'SYNTH_TF',featureId:'SYNTH_TF_UNSCORED',view:'umap'}));assert.equal(missing.svg,null);assert.match(missing.note,/not substituted/);
});
check('Shared pooled-lineage hubs and networks are retained once without endpoint refitting',()=>{
 const x=JSON.parse(JSON.stringify(fixture));x.endpoints.push({...x.endpoints[0],id:'SYNTH_CA1_SL',label:'SYNTHETIC CA1 SL'});x.methods[0].endpoint_ids=['SYNTH_CA1_DL','SYNTH_CA1_SL'];delete x.features.find(f=>f.id==='SYNTH_MOD_A').endpoint_id;
 for(const name of ['hubs','nodes','edges'])for(const r of x[name])if(r.method_id==='SYNTH_HD_CA1')r.endpoint_id=null;
 const shared=A.validatePayload(x),state=ep=>A.normalizeState({...A.DEFAULT,methodId:'SYNTH_HD_CA1',endpointId:ep,view:'network'},shared),all=A.graph(shared,state('')),dl=A.graph(shared,state('SYNTH_CA1_DL')),sl=A.graph(shared,state('SYNTH_CA1_SL'));
 assert.equal(all.edges.length,2);assert.equal(all.nodes.length,3);assert.equal(JSON.stringify(dl.edges),JSON.stringify(sl.edges));assert.equal(JSON.stringify(all.edges),JSON.stringify(dl.edges));assert.match(A.buildView(shared,state('SYNTH_CA1_SL')).note,/endpoint selection does not refit edges/);assert.equal(A.buildView(shared,{...state('SYNTH_CA1_SL'),view:'hubs'}).rows.length,2);
 const missing=JSON.parse(JSON.stringify(x));delete missing.methods[0].endpoint_ids;assert.throws(()=>A.validatePayload(missing),/absent method, endpoint/);
});
check('Saved mouse-level comparisons preserve the original cohort and source values',()=>{
 const wt=A.buildView(d,s({view:'effects',genotype:'WT'})),ko=A.buildView(d,s({view:'effects',genotype:'cKO'}));assert.equal(JSON.stringify(wt.rows),JSON.stringify(ko.rows));assert.equal(wt.rows[0].padj,.4);assert.match(A.csv(A.mouseRows(d,s({featureId:'SYNTH_MOD_CONSTANT'}))),/NA/);assert.equal(JSON.stringify(d),before);
});
check('Display controls are present and invalid numeric colors are rejected',()=>{
 for(const id of Object.keys(A.DEFAULT).filter(k=>k!=='view'))assert.match(html,new RegExp('id="'+id+'"'),id);assert.throws(()=>A.buildView(d,s({view:'heatmap',scaleMode:'manual',colorMin:2,colorMax:1})),/maximum/);const normalized=s({lowColor:'javascript:bad',limit:9000});assert.equal(normalized.lowColor,A.DEFAULT.lowColor);assert.equal(normalized.limit,120);
});
check('SATB2 focus requires a completed curated method and selects its saved mouse scores',()=>{
 assert.equal(A.satb2Focus(d,s({})),null);
 const x=JSON.parse(JSON.stringify(fixture));x.features.find(f=>f.id==='SYNTH_TF_A').label='Satb2';const available=A.validatePayload(x),original=JSON.stringify(available),focus=A.satb2Focus(available,s({search:'unrelated',endpointFamily:'CA1',genotype:'WT',scaleMode:'manual'}));
 assert.equal(focus.methodId,'SYNTH_TF');assert.equal(focus.featureId,'SYNTH_TF_A');assert.equal(focus.endpointId,'');assert.equal(focus.endpointFamily,'all');assert.equal(focus.genotype,'all');assert.equal(focus.view,'replicates');assert.equal(focus.scoreMode,'raw');assert.equal(focus.scaleMode,'auto');assert.equal(focus.search,'');assert.equal(JSON.stringify(available),original);
});

// DOM, image and PDF mocks exercise current-view export snapshots without a browser.
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
 await U.start();assert.equal(U.state.methodId,'SYNTH_HD_CA1');assert(U.plot.svg);for(const id of Object.keys(U.DEFAULT).filter(k=>k!=='view'))assert(elements.get(id).listeners.input||elements.get(id).listeners.change,id+' is unbound');checks.push('Every saved control is bound, and loaded methods retain distinct families');
 U.change('view','network');U.render();const original=U.plot;await U.exportPlot('svg');assert.equal(await downloads.at(-1).blob.text(),original.svg);await U.exportPlot('csv');assert.equal((await downloads.at(-1).blob.text()).trimEnd().split('\n').length,3);await U.exportPlot('json');assert.equal(JSON.parse(await downloads.at(-1).blob.text()).analysis_id,fixture.analysis_id);await U.exportPlot('png');assert.equal(downloads.at(-1).blob.type,'image/png');
 let release;decodeWait=new Promise(r=>{release=r;});const pending=U.exportPlot('pdf');U.change('methodId','SYNTH_TF');U.change('view','umap');U.render();release();decodeWait=null;await pending;assert.equal(pdfs.length,1);assert(pdfs[0].name.includes('SYNTH_HD_CA1')&&pdfs[0].name.includes('network'));assert.equal(JSON.parse(pdfs[0].properties.subject).method,'SYNTH_HD_CA1');checks.push('PDF, SVG, PNG, CSV and JSON exports preserve the captured source/view during navigation');
 const saved=U.settingsRecord();assert.throws(()=>U.importSettings({...saved,analysis_id:'OTHER_ANALYSIS'}),/different source/);U.importSettings({...saved,state:{...saved.state,pointSize:999,nodeColor:'javascript:bad'}});assert.equal(U.state.pointSize,8);assert.equal(U.state.nodeColor,U.DEFAULT.nodeColor);checks.push('Settings import enforces source identity and bounded appearance controls');
 U.change('view','heatmap');U.change('scaleMode','manual');U.change('colorMin',5);U.change('colorMax',1);U.render();assert.equal(U.plot,null);assert.equal(elements.get('plot-pdf').disabled,true);assert.match(elements.get('plot-status').textContent,/maximum/);checks.push('An invalid color range clears stale plots and disables old exports');
 const satb2=JSON.parse(JSON.stringify(fixture));satb2.features.find(f=>f.id==='SYNTH_TF_A').label='Satb2';U.setData(satb2,{ignoreSaved:true});assert.equal(elements.get('focus-satb2').disabled,false);elements.get('focus-satb2').click();assert.equal(U.state.methodId,'SYNTH_TF');assert.equal(U.state.featureId,'SYNTH_TF_A');assert.equal(U.state.view,'replicates');assert(U.plot.svg);assert.equal(JSON.stringify(U.data.values),JSON.stringify(satb2.values));checks.push('One-click SATB2 view is bound and preserves the native values');
 console.log(JSON.stringify({status:'PASS',checks:checks.length,names:checks,synthetic_fixture_only:true,scientific_analysis_performed:false,browser_launched:false},null,2));
})().catch(e=>{console.error(e);process.exitCode=1;});
