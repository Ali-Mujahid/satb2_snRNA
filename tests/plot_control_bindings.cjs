'use strict';
// Actual control handlers in a tiny DOM; no browser, network, fitting or study data.
const fs=require('fs'),path=require('path'),vm=require('vm'),assert=require('assert'),crypto=require('crypto');
const site=process.env.SATB2_SITE_DIR||path.resolve(__dirname,'../site'),checks=[];
class Element {constructor(tag='input'){this.tagName=tag;this.value='';this.checked=false;this.children=[];this.style={};this.dataset={};this.textContent='';}append(...x){this.children.push(...x);}replaceChildren(...x){this.children=[...x];}addEventListener(name,fn){this[name]=fn;}get options(){return this.children;}}
const controls={},document={activeElement:null,getElementById:id=>controls[id]||(controls[id]=new Element()),createElement:tag=>new Element(tag),querySelectorAll:()=>[]};
for(const match of fs.readFileSync(path.join(site,'index.html'),'utf8').matchAll(/\bid="([^"]+)"/g))document.getElementById(match[1]);
const context={console,document,setTimeout,clearTimeout,localStorage:{getItem:()=>null,setItem(){}},SATB2_NO_AUTO_INIT:true,SATB2_GENE_EXPRESSION_CATALOG:{branches:{original_full:{}},genes:[]}};context.window=context;context.globalThis=context;vm.createContext(context);
for(const file of ['v2_features.js','display_controls.js','volcano_controls.js'])vm.runInContext(fs.readFileSync(path.join(site,file),'utf8'),context);
const source=fs.readFileSync(path.join(site,'app.js'),'utf8');vm.runInContext(source.replace('const api={','const api={testBindPlot:bindPlotControls,testCategoryColor:categoricalColor,'),context);
const app=context.SATB2_APP,plain=x=>JSON.parse(JSON.stringify(x)),label='Superficial-associated expression program enriched',genes=[{gene_id:'g1',gene_symbol:'Cck'},{gene_id:'g2',gene_symbol:'Calb1'},{gene_id:'g3',gene_symbol:'Apoe'}];
const coords={nuclei:['0','1'].map((cluster_id,i)=>({nucleus_id:'n'+i,sample:'mouse'+i,condition:i?'cKO':'WT',cluster_id,reviewed_label:label,x_harmony:i,y_harmony:i,x_unintegrated:i,y_unintegrated:i})),clusters:['0','1'].map(cluster_id=>({cluster_id,reviewed_label:label}))};
const manifest={schema_version:'satb2-offline-1',study:{source_fingerprint:'bindings-fixture'},defaults:{population:'CA1_ProS',branch:'original_full',marker_ids:[]},genes,shards:{scanpycoords:{},seuratcoords:{}},datasets:['scanpy','seurat'].map(tool=>({tool,branch:'original_full',population:'CA1_ProS',cluster_scope:'original_frozen',coordinates_shard:tool+'coords'}))};
const catalog={fits:['scanpy','seurat'].map(tool=>({id:tool+'fit',tool,branch:'original_full',population:'CA1_ProS',cluster_scope:'original_frozen',status:'TESTED_REVIEW_REQUIRED',endpoint:'PRIMARY',analysis_kind:'total',members:[0,1]}))};
function check(name,fn){fn();checks.push(name);console.log('PASS '+name);}
async function input(id,value){controls[id].value=value;await controls[id].oninput();}
(async()=>{
 await app.init({manifest,catalog,headless:true,ignoreSavedState:true,shardLoader:async()=>({encoding:'json',data:plain(coords)})});app.testBindPlot();
 const frozen=JSON.stringify(app.active),members=JSON.stringify(['scanpy','seurat'].map(t=>app.targetMembers(t))),fits=JSON.stringify(['scanpy','seurat'].map(t=>app.resolveFit(t).fit.id));
 check('Every new plot-binding control exists exactly once in actual HTML',()=>{const html=fs.readFileSync(path.join(site,'index.html'),'utf8');for(const id of Object.keys(controls)){const matches=[...html.matchAll(new RegExp('id="'+id+'"','g'))];assert.equal(matches.length,1,id);}});
 app.state.markers=['g1','g2','g3'];app.state.markerOrder='selection';app.state.markerPreset='fixture';app.state.markerPage=2;controls['marker-reverse'].onclick();
 check('Reverse button reverses displayed genes and switches to editable selection order',()=>{assert.deepEqual(plain(app.state.markers),['g3','g2','g1']);assert.equal(app.state.markerOrder,'selection');assert.equal(app.state.markerPreset,'');assert.equal(app.state.markerPage,0);});
 app.state.markers=['g1','g3','g2'];controls['marker-alphabetical'].onclick();check('A–Z button sorts actual symbols without adding or losing genes',()=>assert.deepEqual(plain(app.state.markers),['g3','g2','g1']));
 controls['label-target'].value=label;controls['label-alias'].value='  Sup.custom  ';controls['apply-label-alias'].onclick();
 check('Applying display alias trims text and activates custom display only',()=>{assert.equal(app.state.labelAliases[label],'Sup.custom');assert.equal(app.state.labelStyle,'custom');assert.equal(app.displayLabel('C0 · '+label),'C0 · Sup.custom');assert.deepEqual(plain(app.state.labels),{});});
 controls['label-target'].onchange();check('Selecting alias target loads existing display alias into editor',()=>assert.equal(controls['label-alias'].value,'Sup.custom'));
 controls['reset-label-alias'].onclick();check('Reset display alias removes its exact source-keyed entry',()=>assert(!Object.hasOwn(app.state.labelAliases,label)));
 controls['label-alias'].value=' ';controls['apply-label-alias'].onclick();check('Blank alias is ignored without source-label changes',()=>assert(!Object.hasOwn(app.state.labelAliases,label)));
 controls['category-color-target'].value=label;controls['category-color'].value='#112233';controls['apply-category-color'].onclick();
 check('Category Apply uses current field and original unshortened category key',()=>{assert.equal(app.state.categoryColors.reviewed_label[label],'#112233');assert.equal(app.testCategoryColor(label),'#112233');});
 controls['category-color'].value='url(invalid)';controls['apply-category-color'].onclick();check('Invalid color cannot replace a valid category override',()=>assert.equal(app.state.categoryColors.reviewed_label[label],'#112233'));
 app.state.categoryColors.reviewed_label.Other='#445566';controls['reset-category-color'].onclick();check('Category Reset removes only the selected override',()=>{assert(!Object.hasOwn(app.state.categoryColors.reviewed_label,label));assert.equal(app.state.categoryColors.reviewed_label.Other,'#445566');});
 await input('volcano-label-count','20');await input('volcano-up','#aa0000');await input('volcano-down','#0000aa');await input('volcano-ns','#777777');
 check('Volcano numeric/count/color controls update actual saved settings',()=>{assert.equal(app.state.volcanoLabelCount,20);assert.equal(app.state.volcanoUp,'#aa0000');assert.equal(app.state.volcanoDown,'#0000aa');assert.equal(app.state.volcanoNS,'#777777');});
 await input('volcano-label-count','999');check('Volcano label maximum is clamped to100',()=>assert.equal(app.state.volcanoLabelCount,100));
 controls['volcano-labels'].checked=false;controls['volcano-labels'].onchange();controls['volcano-thresholds'].checked=false;controls['volcano-thresholds'].onchange();
 check('Volcano label/threshold checkbox handlers toggle actual settings',()=>{assert.equal(app.state.volcanoLabels,false);assert.equal(app.state.volcanoThresholds,false);});
 controls['auto-gene'].checked=false;controls['auto-gene'].onchange();check('Automatic gene checkbox uses actual handler',()=>assert.equal(app.state.autoGene,false));
 await input('expression-resolution','cluster');check('Expression resolution handler activates cluster mode',()=>assert.equal(app.state.expressionResolution,'cluster'));
 app.state.overlayMetric='fraction_detected';await input('expression-resolution','nucleus');check('Nucleus mode normalizes unsupported summary metric without inventing individual detection fractions',()=>{assert.equal(app.state.expressionResolution,'nucleus');assert.equal(app.state.overlayMetric,'mean_log1p_normalized');});
 await input('umap-range-mode','manual');await input('umap-min','0.5');await input('umap-max','4');await input('gene-zero-color','#eeeeee');await input('categorical-palette','bright');await input('marker-label-mode','id');
 check('Manual expression ranges, zero color, categorical palette and row-label mode are wired',()=>{assert.equal(app.state.umapRangeMode,'manual');assert.equal(app.state.umapMin,.5);assert.equal(app.state.umapMax,4);assert.equal(app.state.geneZeroColor,'#eeeeee');assert.equal(app.state.categoricalPalette,'bright');assert.equal(app.state.markerLabelMode,'id');});
 check('Every display-handler operation preserves source annotations, coordinates, exact DE members and fit IDs',()=>{assert.equal(JSON.stringify(app.active),frozen);assert.equal(JSON.stringify(['scanpy','seurat'].map(t=>app.targetMembers(t))),members);assert.equal(JSON.stringify(['scanpy','seurat'].map(t=>app.resolveFit(t).fit.id)),fits);assert.deepEqual(plain(app.state.labels),{});});
 const report={status:'PASS',scope:'ACTUAL_PRIVATE_HANDLER_VM_BINDINGS_NO_BROWSER_NO_SCIENTIFIC_CHANGES',checks:checks.length,check_names:checks,app_sha256:crypto.createHash('sha256').update(source).digest('hex')};fs.writeFileSync(path.join(__dirname,'plot_control_bindings_acceptance.json'),JSON.stringify(report,null,2)+'\n');console.log(JSON.stringify({status:'PASS',checks:checks.length}));
})().catch(e=>{console.error(e);process.exitCode=1;});
