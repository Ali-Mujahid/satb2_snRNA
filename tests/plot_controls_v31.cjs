'use strict';
// VM integration: actual source Cck vectors and an unchanged native-coordinate
// subset, plus adversarial tiny source-axis fixtures. No browser/network/fits.
const fs=require('fs'),path=require('path'),vm=require('vm'),assert=require('assert'),zlib=require('zlib'),crypto=require('crypto');
const root=path.resolve(__dirname,'..'),site=process.env.SATB2_SITE_DIR||path.join(root,'site'),checks=[];
const readJSON=file=>JSON.parse(fs.readFileSync(file,'utf8')),plain=x=>JSON.parse(JSON.stringify(x)),hash=x=>crypto.createHash('sha256').update(x).digest('hex');
const source=fs.readFileSync(path.join(site,'app.js'),'utf8'),nativeManifest=readJSON(path.join(site,'data/manifest.json')),geneCatalog=readJSON(path.join(site,'data/GENE_EXPRESSION_CATALOG.json'));
const cck=geneCatalog.genes.find(g=>g.gene_symbol==='Cck').gene_id,branches=Object.keys(geneCatalog.branches),payloads=new Map();
function rawPayload(id){if(!payloads.has(id)){const d=geneCatalog.shards[id],bytes=fs.readFileSync(path.join(root,'gene_expression',path.basename(d.path)));assert.equal(hash(bytes),d.sha256);const raw=zlib.gunzipSync(bytes);assert.equal(hash(raw),d.payload_sha256);payloads.set(id,JSON.parse(raw));}return payloads.get(id);}
const axes=Object.fromEntries(branches.map(b=>[b,rawPayload(geneCatalog.branches[b].axis_shard)]));
const vectors=Object.fromEntries(branches.map(b=>[b,rawPayload(geneCatalog.branches[b].genes[cck])]));
const axisMaps=Object.fromEntries(branches.map(b=>[b,new Map(axes[b].nucleus_ids.map((id,i)=>[id,i]))]));
const counts=Object.fromEntries(branches.map(b=>[b,new Map(vectors[b].indices.map((i,j)=>[i,vectors[b].counts[j]]))]));
function expected(id,b,raw=false){const i=axisMaps[b].get(id);if(i===undefined)return null;const n=counts[b].get(i)||0;return raw?n:Math.log1p(10000*n/axes[b].library_counts[i]);}
const d=nativeManifest.datasets.find(d=>d.tool==='scanpy'&&d.branch==='original_full'&&d.population==='CA1_ProS'),coords=JSON.parse(zlib.gunzipSync(fs.readFileSync(path.join(site,nativeManifest.shards[d.coordinates_shard].path))));
const groups=new Map();for(const n of coords.nuclei){if(!branches.every(b=>axisMaps[b].has(n.nucleus_id)))continue;const key=String(n.cluster_id);if(!groups.has(key))groups.set(key,[]);groups.get(key).push(n);}
const group=[...groups.values()].find(g=>g.some(n=>expected(n.nucleus_id,'original_full',true)===0)&&new Set(g.map(n=>expected(n.nucleus_id,'original_full',true))).size>=4);
assert(group,'Need a real source cluster with zero and varying positive Cck counts');
const selected=[];for(const n of group){const ncount=expected(n.nucleus_id,'original_full',true);if(!selected.some(x=>expected(x.nucleus_id,'original_full',true)===ncount))selected.push(n);if(selected.length===7)break;}
const nativeSubset={nuclei:selected.map(n=>({...n})),clusters:coords.clusters.filter(c=>String(c.cluster_id)===String(group[0].cluster_id))};
function context(cat=geneCatalog){const c={console,URL,Uint8Array,Float64Array,ArrayBuffer,TextEncoder,TextDecoder,atob,Blob,Response,DecompressionStream,crypto:crypto.webcrypto,setTimeout,clearTimeout,SATB2_NO_AUTO_INIT:true,localStorage:{getItem:()=>null,setItem(){}},SATB2_GENE_EXPRESSION_CATALOG:cat};c.window=c;c.globalThis=c;vm.createContext(c);for(const f of ['vendor/fflate.min.js','v2_features.js','display_controls.js','volcano_controls.js'])if(fs.existsSync(path.join(site,f)))vm.runInContext(fs.readFileSync(path.join(site,f),'utf8'),c);vm.runInContext(source.replace('const api={','const api={testCategory:category,testCategoryColor:categoricalColor,testOverlayValue:overlayValue,'),c);return c;}
function fixtureManifest(nuclei=nativeSubset,cat=geneCatalog){return {schema_version:'satb2-offline-1',study:{source_fingerprint:'v31-actual-count-fixture'},defaults:{population:'CA1_ProS',branch:'original_full',marker_ids:[]},genes:cat.genes||[{gene_id:'g',gene_symbol:'Cck'}],shards:Object.fromEntries(branches.flatMap(b=>['scanpy','seurat'].map(t=>[t+'__'+b+'__coords',{}]))),datasets:branches.flatMap(branch=>['scanpy','seurat'].map(tool=>({tool,branch,population:'CA1_ProS',cluster_scope:'original_frozen',coordinates_shard:tool+'__'+branch+'__coords'})))};}
function check(name,fn){fn();checks.push(name);console.log('PASS '+name);}
async function rejects(name,job,pattern){await assert.rejects(job,pattern);checks.push(name);console.log('PASS '+name);}
(async()=>{
 let deferredId=null,releaseDeferred=null,deferredStarted=null;
 const c=context(),app=c.SATB2_APP,manifest=fixtureManifest(),fitCatalog={fits:branches.flatMap(branch=>['scanpy','seurat'].map(tool=>({id:tool+'_'+branch,status:'TESTED_REVIEW_REQUIRED',endpoint:'PRIMARY',branch,tool,population:'CA1_ProS',cluster_scope:'original_frozen',analysis_kind:'total',analysis_label:'Fixed source fixture',members:selected.map((_,i)=>i)})))};
 const loader=async(id,desc)=>{if(id.endsWith('__coords'))return {encoding:'json',data:plain(nativeSubset)};if(id===deferredId)await new Promise(resolve=>{releaseDeferred=resolve;deferredStarted();});return {encoding:'json',data:rawPayload(id)};};
 await app.init({manifest,catalog:fitCatalog,headless:true,ignoreSavedState:true,shardLoader:loader});
 check('New atlas defaults to individual-nucleus expression when catalog is present',()=>assert.equal(app.state.expressionResolution,'nucleus'));
 check('A true expression default does not eagerly load all gene vectors',()=>assert.equal(app.nucleusGeneReady(),false));
 app.state.annotation='gene';app.state.overlayGene=cck;await app.ensureNucleusGene();
 const frozen=JSON.stringify(app.active),members=JSON.stringify(app.targetMembers('scanpy')),fitId=app.resolveFit('scanpy').fit.id;
 check('Actual Cck vector uses full-source log1p libraries for every displayed nucleus',()=>{for(const t of ['scanpy','seurat'])for(const n of app.active[t].coords.nuclei)assert.equal(app.nucleusGeneValue(n),expected(n.nucleus_id,'original_full'));});
 check('The same stable nucleus ID gives the same expression in both tools',()=>{for(let i=0;i<selected.length;i++)assert.equal(app.nucleusGeneValue(app.active.scanpy.coords.nuclei[i]),app.nucleusGeneValue(app.active.seurat.coords.nuclei[i]));});
 check('Nuclei from one native cluster retain measured within-cluster expression variation',()=>{assert.equal(new Set(selected.map(n=>n.cluster_id)).size,1);assert(new Set(selected.map(n=>app.nucleusGeneValue(n))).size>=4);});
 check('Measured zeros remain exact zero and an unknown source nucleus is unavailable',()=>{const zero=selected.find(n=>expected(n.nucleus_id,'original_full',true)===0);assert.equal(app.nucleusGeneValue(zero),0);assert.equal(app.nucleusGeneValue({nucleus_id:'NOT_A_SOURCE_NUCLEUS'}),null);});
 app.state.overlayMetric='mean_source_counts';
 check('Raw mode returns exact native integer counts instead of normalized or cluster means',()=>{for(const n of selected)assert.equal(app.nucleusGeneValue(n),expected(n.nucleus_id,'original_full',true));});
 app.state.overlayMetric='mean_log1p_normalized';app.state.branch='paired_cellbender';
 check('Changing the source branch immediately makes stale old values unavailable',()=>{assert.equal(app.nucleusGeneReady(),false);for(const n of selected)assert.equal(app.nucleusGeneValue(n),null);});
 await app.changeScope();
 check('CellBender overlay uses its own exact counts and complete corrected library',()=>{for(const n of selected)assert.equal(app.nucleusGeneValue(n),expected(n.nucleus_id,'paired_cellbender'));assert(selected.some(n=>expected(n.nucleus_id,'paired_cellbender')!==expected(n.nucleus_id,'original_full')));});
 app.state.branch='paired_original';app.state.annotation='reviewed_label';await app.changeScope();app.state.annotation='gene';
 deferredId=geneCatalog.branches.paired_original.genes[cck];const pendingStarted=new Promise(resolve=>deferredStarted=resolve),oldJob=app.ensureNucleusGene();await pendingStarted;
 app.state.branch='original_full';await app.changeScope();releaseDeferred();await oldJob;deferredId=null;
 check('A late previous-branch response never replaces the current branch expression',()=>{assert.equal(app.nucleusGeneReady(),true);for(const n of selected)assert.equal(app.nucleusGeneValue(n),expected(n.nucleus_id,'original_full'));});
 app.state.labelStyle='short';app.state.labelAliases={'Superficial-associated expression program enriched':'My SL'};app.state.categoryColors={reviewed_label:{'Superficial-associated expression program enriched':'#112233'}};app.state.categoricalPalette='bright';
 const sourceCategory=selected[0].reviewed_label;
 check('Automatic compact labels expose SL-like and DL-like display names',()=>{assert.equal(app.displayLabel('Superficial-associated expression program enriched'),'SL-like');assert.equal(app.displayLabel('Deep-associated expression program enriched'),'DL-like');});
 app.state.annotation='reviewed_label';
 check('Category colors remain exact original-label keyed rather than shortened-label keyed',()=>assert.equal(app.testCategoryColor('Superficial-associated expression program enriched'),'#112233'));
 app.state.annotation='gene';app.state.umapCondition='WT';app.state.umapProgram='ALL';app.state.umapLabelMode='both';app.state.markerRowOrder='descending';app.state.markerOrder='focus';
 check('Palette, compact labels, gene overlay and display-only filters preserve coordinates and DE members',()=>{assert.equal(JSON.stringify(app.active),frozen);assert.equal(JSON.stringify(app.targetMembers('scanpy')),members);assert.equal(app.resolveFit('scanpy').fit.id,fitId);assert.equal(app.active.scanpy.coords.nuclei[0].reviewed_label,sourceCategory);});
 const saved=plain(app.exportState());app.state.labelStyle='full';app.state.expressionResolution='cluster';app.state.volcanoLabelCount=0;await app.importState(saved);
 check('Settings import round-trips individual gene, colors, compact labels and new plot controls',()=>{assert.equal(app.state.overlayGene,cck);assert.equal(app.state.expressionResolution,'nucleus');assert.equal(app.state.labelStyle,'short');assert.equal(app.state.categoryColors.reviewed_label['Superficial-associated expression program enriched'],'#112233');assert.equal(app.state.volcanoLabelCount,20);assert.equal(app.state.umapCondition,'WT');assert.equal(app.state.markerRowOrder,'descending');});
 check('Settings reload keeps per-nucleus source values and inferential membership intact',()=>{assert.equal(JSON.stringify(app.targetMembers('scanpy')),members);for(const n of selected)assert.equal(app.nucleusGeneValue(n),expected(n.nucleus_id,'original_full'));assert.equal(app.resolveFit('scanpy').fit.id,fitId);});
 const normal={expressionResolution:'nucleus',overlayMetric:'fraction_detected',volcanoLabelCount:999,labelStyle:'invalid'};app.normalizePlotState(normal);
 check('Unsupported nucleus summary statistics normalize to measured per-nucleus expression',()=>{assert.equal(normal.overlayMetric,'mean_log1p_normalized');assert.equal(normal.volcanoLabelCount,100);assert.equal(normal.labelStyle,'compact');});
 await app.resetState();
 check('Reset preserves the measured-expression default and restores plot-control defaults',()=>{assert.equal(app.state.expressionResolution,'nucleus');assert.equal(app.state.labelStyle,'compact');assert.equal(app.state.volcanoLabelCount,20);assert.equal(JSON.stringify(app.targetMembers('scanpy')),members);});
 const legacy=context({}),legacyApp=legacy.SATB2_APP;await legacyApp.init({manifest:fixtureManifest(nativeSubset,{genes:geneCatalog.genes}),catalog:fitCatalog,headless:true,ignoreSavedState:true,shardLoader:async()=>({encoding:'json',data:plain(nativeSubset)})});
 check('A legacy atlas without source vectors explicitly falls back to cluster summaries',()=>assert.equal(legacyApp.state.expressionResolution,'cluster'));
 const advertised=context({}),advertisedApp=advertised.SATB2_APP,advertisedManifest=fixtureManifest(nativeSubset,{genes:geneCatalog.genes});advertisedManifest.measurement_contract={source_nucleus_expression_extracted:true};
 await advertisedApp.init({manifest:advertisedManifest,catalog:fitCatalog,headless:true,ignoreSavedState:true,shardLoader:async()=>({encoding:'json',data:plain(nativeSubset)})});
 check('An atlas advertising measured nucleus expression never silently falls back when its catalog is absent',()=>assert.equal(advertisedApp.state.expressionResolution,'nucleus'));
 await rejects('An advertised measured-expression request fails closed when its catalog is unavailable',advertisedApp.ensureNucleusGene(),/not packaged/);
 legacyApp.state.expressionResolution='nucleus';legacyApp.normalizePlotState(legacyApp.state);
 check('An explicit individual-nucleus selection survives normalization when its catalog is absent',()=>assert.equal(legacyApp.state.expressionResolution,'nucleus'));
 await rejects('An explicit measured-expression request in a legacy atlas fails instead of substituting cluster means',legacyApp.ensureNucleusGene(),/not packaged/);
 const tinyBase={schema:'satb2-nucleus-expression-catalog-1',genes:[{gene_id:'g',gene_symbol:'Cck'}],branches:{original_full:{axis_shard:'axis',genes:{g:'vector'}}},shards:{axis:{payload_sha256:'a'.repeat(64)},vector:{}}},tinyCoords={nuclei:[{nucleus_id:'n1',cluster_id:'0',condition:'WT',x_harmony:0,y_harmony:1,x_unintegrated:0,y_unintegrated:1}],clusters:[{cluster_id:'0',reviewed_label:'Fixture'}]};
 async function tinyCase(axis,vector,extra=[]){const tc=context(tinyBase),ta=tc.SATB2_APP,tm=fixtureManifest(tinyCoords,tinyBase);await ta.init({manifest:tm,catalog:{fits:[]},headless:true,ignoreSavedState:true,shardLoader:async id=>({encoding:'json',data:id.endsWith('__coords')?{...tinyCoords,nuclei:[...tinyCoords.nuclei,...extra]}:id==='axis'?axis:vector})});ta.state.overlayGene='g';return ta;}
 const validAxis={branch:'original_full',nucleus_ids:['n1','n2'],library_counts:[10,20]},validVector={branch:'original_full',gene_id:'g',n_source_nuclei:2,axis_payload_sha256:'a'.repeat(64),indices:[0],counts:[2]};
 let tiny=await tinyCase({...validAxis,library_counts:[10]},validVector);await rejects('Mismatched nucleus/library axis lengths are rejected',tiny.ensureNucleusGene(),/axis mismatch/);
 tiny=await tinyCase({...validAxis,nucleus_ids:['n1','n1']},validVector);await rejects('Duplicate source nucleus IDs are rejected',tiny.ensureNucleusGene(),/Invalid full-library/);
 tiny=await tinyCase({...validAxis,library_counts:[10,0]},validVector);await rejects('A zero full-source library is rejected rather than divided by',tiny.ensureNucleusGene(),/Invalid full-library/);
 tiny=await tinyCase(validAxis,{...validVector,axis_payload_sha256:'b'.repeat(64)});await rejects('A sparse gene vector cannot bind to another source-axis hash',tiny.ensureNucleusGene(),/does not match/);
 for(const [suffix,change] of [['duplicate indices',{indices:[0,0],counts:[1,2]}],['fractional count',{counts:[1.5]}],['unsafe count',{counts:[Number.MAX_SAFE_INTEGER+1]}],['out-of-range index',{indices:[2]}]]){tiny=await tinyCase(validAxis,{...validVector,...change});await rejects('Invalid sparse input rejected: '+suffix,tiny.ensureNucleusGene(),/Invalid sparse source/);}
 tiny=await tinyCase(validAxis,validVector,[{...tinyCoords.nuclei[0],nucleus_id:'absent'}]);await rejects('A displayed nucleus missing from its measured source axis fails closed',tiny.ensureNucleusGene(),/displayed nucleus is absent/);
 const pdfContext=context(),pdfApp=pdfContext.SATB2_APP;let pdfStartedResolve,pdfRelease,constructedPDF=0,savedPDF=0;
 const pdfStarted=new Promise(resolve=>pdfStartedResolve=resolve);
 await pdfApp.init({manifest:fixtureManifest(),catalog:fitCatalog,headless:true,ignoreSavedState:true,shardLoader:async id=>{if(id.endsWith('__coords'))return {encoding:'json',data:plain(nativeSubset)};if(id===geneCatalog.branches.original_full.genes[cck]){pdfStartedResolve();await new Promise(resolve=>pdfRelease=resolve);}return {encoding:'json',data:rawPayload(id)};}});
 pdfApp.state.annotation='gene';pdfApp.state.overlayGene=cck;
 pdfContext.jspdf={jsPDF:class {constructor(){constructedPDF++;this.internal={pageSize:{getWidth:()=>297,getHeight:()=>210}};}save(){savedPDF++;}}};
 // A changed selection must abort before a canvas read or PDF construction. This
 // intentionally has no canvas fixture: an old/stale plot must not be exported.
 pdfContext.document={getElementById(){throw Error('Stale PDF path reached a canvas before its pending source expression resolved');}};
 const pendingPDF=pdfApp.exportPDF('umap','scanpy');
 let timer;await Promise.race([pdfStarted,new Promise((_,reject)=>timer=setTimeout(()=>reject(Error('PDF export never requested its required source expression')),2000))]);clearTimeout(timer);
 check('PDF waits for the requested measured-expression vector before constructing a document',()=>{assert.equal(constructedPDF,0);assert.equal(savedPDF,0);});
 pdfApp.state.palette='magma';pdfRelease();const canceledPDF=await pendingPDF;
 check('A changed plot selection cancels a pending PDF before drawing or saving stale content',()=>{assert.equal(canceledPDF,undefined);assert.equal(constructedPDF,0);assert.equal(savedPDF,0);});
 const report={status:'PASS',scope:'VM_INTEGRATION_ACTUAL_CCK_VECTORS_AND_NATIVE_COORDINATE_SUBSET_PLUS_SYNTHETIC_EDGE_CASES_NO_BROWSER_NO_NEW_FITS',checks:checks.length,check_names:checks,source_vector_catalog_sha256:hash(fs.readFileSync(path.join(site,'data/GENE_EXPRESSION_CATALOG.json'))),app_sha256:hash(source),selected_actual_nucleus_ids:selected.map(n=>n.nucleus_id),source_native_cluster:String(selected[0].cluster_id)};
 fs.writeFileSync(path.join(__dirname,'plot_controls_v31_acceptance.json'),JSON.stringify(report,null,2)+'\n');console.log(JSON.stringify({status:'PASS',checks:checks.length}));
})().catch(e=>{console.error(e);process.exitCode=1;});
