'use strict';
// Actual application export path, shared canvas renderer and a recording PDF.
// Tiny fixtures only: no scientific model, browser, network or study-data writes.
const fs=require('node:fs'),path=require('node:path'),vm=require('node:vm'),assert=require('node:assert/strict');
const site=process.env.SATB2_SITE_DIR||path.resolve(__dirname,'../site');
const checks=[],plain=value=>JSON.parse(JSON.stringify(value));
function check(name,fn){fn();checks.push(name);console.log('PASS '+name);}
class CanvasContext {
  constructor(){this.calls=[];this.font='10px system-ui';this.stack=[];this.fillStyle='';this.strokeStyle='';this.globalAlpha=1;}
  save(){this.stack.push({font:this.font,fillStyle:this.fillStyle,strokeStyle:this.strokeStyle,globalAlpha:this.globalAlpha});}
  restore(){Object.assign(this,this.stack.pop());}
  measureText(text){return {width:String(text).length*parseFloat(this.font)*.55};}
}
for(const method of ['setTransform','clearRect','fillRect','strokeRect','beginPath','moveTo','lineTo','stroke','arc','fill','translate','rotate','rect','clip','setLineDash','fillText'])CanvasContext.prototype[method]=function(...args){this.calls.push({method,args,fillStyle:this.fillStyle,strokeStyle:this.strokeStyle});};
class Element {
  constructor(tag='input'){this.tagName=tag;this.value='';this.checked=false;this.children=[];this.style={};this.dataset={};this.textContent='';this.classList={add(){},toggle(){}};this.clientWidth=700;}
  append(...values){this.children.push(...values);}
  replaceChildren(...values){this.children=[...values];}
  getContext(){return this.context||(this.context=new CanvasContext());}
  toDataURL(){return 'data:image/png;base64,recording';}
  addEventListener(){}
  get options(){return this.children;}
}
const types={allcells:'All nuclei',Excitatory:'Excitatory',Inhibitory:'Inhibitory',CA1_ProS:'CA1 / ProS'};
const tools=['scanpy','seurat'],branches=['original_full','paired_cellbender'];
const coords={nuclei:Array.from({length:4},(_,i)=>({nucleus_id:'mouse_'+i,cluster_id:String(i%2),sample:'mouse'+i,condition:i<2?'WT':'cKO',x_harmony:i,y_harmony:i,x_unintegrated:i,y_unintegrated:i})),clusters:[{cluster_id:'0',reviewed_label:'Pvalb Gaba'},{cluster_id:'1',reviewed_label:'Sst Gaba'}]};
function rowsFor(id){const seurat=id.startsWith('seurat');return [{gene_id:'g1',gene:seurat?'SeuratHigh':'ScanpyHigh',padj:1e-80,log2FoldChange:seurat?2:-2},{gene_id:'g2',gene:'LowGene',padj:.01,log2FoldChange:1.5},{gene_id:'g3',gene:'GrayGene',padj:.8,log2FoldChange:-.1},{gene_id:'g4',gene:'NAGene',padj:null,log2FoldChange:2}];}
function fixture(population='Inhibitory'){
  const descriptors=branches.flatMap(branch=>tools.map(tool=>({tool,branch,population,cluster_scope:'original_frozen',coordinates_shard:tool+'_'+branch+'_coords'})));
  const fits=descriptors.flatMap(d=>['total','C0'].map(scope=>({...d,id:d.tool+'_'+d.branch+'_'+scope,status:'TESTED_REVIEW_REQUIRED',endpoint:'PRIMARY',analysis_kind:scope==='total'?'total':'cluster',analysis_label:scope==='total'?'Total nuclei':'C0 Pvalb Gaba',members:scope==='total'?[0,1,2,3]:[0,2],rows_shard:d.tool+'_'+d.branch+'_'+scope+'_rows',n_nuclei:scope==='total'?4:2,n_WT:6,n_cKO:6,residual_df:10,membership_sha256:'fixture_membership_'+scope})));
  const shards=Object.fromEntries([...descriptors.map(d=>d.coordinates_shard),...fits.map(f=>f.rows_shard)].map(id=>[id,{}]));
  return {manifest:{schema_version:'satb2-offline-1',study:{source_fingerprint:'clean-pdf-fixture'},defaults:{population,branch:'original_full',marker_ids:[]},genes:[{gene_id:'g1',gene_symbol:'Gene'}],datasets:descriptors,shards},catalog:{fits}};
}
async function context(population='Inhibitory',gate=null){
  const controls={},created=[],documents=[],renders=[];
  const document={activeElement:null,getElementById:id=>controls[id]||(controls[id]=new Element()),createElement:tag=>{const e=new Element(tag);created.push(e);return e;},querySelectorAll:()=>[]};
  class PDF {
    constructor(options){this.options=options;this.texts=[];this.images=[];this.pages=1;this.properties=null;this.filename=null;this.internal={pageSize:{getWidth:()=>options.format[0],getHeight:()=>options.format[1]}};documents.push(this);}
    setFontSize(){}text(value){this.texts.push(String(value));}addImage(...args){this.images.push(args);}addPage(){this.pages++;}splitTextToSize(value){return [String(value)];}setFillColor(){}rect(){}setProperties(value){this.properties=value;}save(filename){this.filename=filename;}
  }
  const c={console,document,setTimeout,clearTimeout,localStorage:{getItem:()=>null,setItem(){}},SATB2_NO_AUTO_INIT:true,jspdf:{jsPDF:PDF}};c.window=c;c.globalThis=c;vm.createContext(c);
  for(const file of ['v2_features.js','display_controls.js','volcano_controls.js','app.js'])vm.runInContext(fs.readFileSync(path.join(site,file),'utf8'),c);
  const render=c.SATB2_VOLCANO.draw;c.SATB2_VOLCANO.draw=(plot,rows,settings,options)=>{const start=plot.c.calls.length,result=render(plot,rows,settings,options);renders.push({plot,rows:plain(rows),options:plain(options||{}),result,calls:plot.c.calls.slice(start)});return result;};
  const f=fixture(population),app=c.SATB2_APP;
  await app.init({...f,headless:true,ignoreSavedState:true,shardLoader:async id=>{if(id.endsWith('_coords'))return {encoding:'json',data:plain(coords)};if(gate&&id===gate.id){gate.start();await gate.promise;}return {encoding:'json',data:{rows:rowsFor(id)}};}});
  return {app,c,documents,renders,controls,created};
}
function frozenScience(app){return JSON.stringify(tools.map(tool=>({coords:app.active[tool].coords,fit:app.resolveFit(tool).fit,members:app.targetMembers(tool),rows:app.boundDERows(tool)})));}
(async()=>{
  const {app,documents,renders}=await context();await app.ensureFits();await app.refreshDERows();
  check('Clean PDF is enabled by default',()=>assert.equal(app.state.volcanoCleanPDF,true));
  const sourceBefore=frozenScience(app);
  for(const branch of branches){app.state.branch=branch;await app.changeScope();await app.ensureFits();await app.refreshDERows();const branchBefore=frozenScience(app);
    for(const tool of tools)for(const view of ['full','low','both']){
      app.state.volcanoView=view;app.state.volcanoCleanPDF=true;const start=renders.length,pdf=await app.exportPDF('volcano',tool),draws=renders.slice(start),screen=draws[0],clean=draws.at(-1),fit=app.resolveFit(tool).fit;
      check(`${branch} ${tool} ${view}: clean PDF shows only type title outside plot`,()=>{assert.deepEqual(pdf.texts,['Inhibitory']);assert.equal(pdf.images.length,1);assert.equal(pdf.pages,1);assert.ok(pdf.filename.includes('_'+tool+'_Inhibitory_Total_nuclei_'+branch+'.pdf'));assert.equal(pdf.properties.title,'Inhibitory');const meta=JSON.parse(pdf.properties.subject);assert.equal(meta.source_fit,fit.id);assert.equal(meta.tool,tool);assert.equal(meta.branch,branch);assert.equal(meta.view,view);assert.equal(meta.FDR,.05);assert.equal(meta.abs_log2FC,0);assert.equal(meta.source_fingerprint,'clean-pdf-fixture');});
      check(`${branch} ${tool} ${view}: clean shared canvas retains exact points, axes and names without notes or leaders`,()=>{assert.equal(screen.result.summaryText,true);assert.equal(clean.result.summaryText,false);assert.equal(clean.result.labelLeaderLines,false);assert.deepEqual(clean.rows,rowsFor(fit.rows_shard));assert.equal(clean.result.validRows,3);assert.equal(clean.result.counts.nonsignificant,1);assert.equal(clean.result.lowPanel?.visibleRows,view==='full'?undefined:2);const texts=clean.calls.filter(x=>x.method==='fillText').map(x=>String(x.args[0]));assert.ok(texts.includes('log2 fold change (cKO / WT)'));assert.ok(texts.includes('−log10 adjusted P'));assert.ok(texts.includes(view==='low'?'LowGene':tool==='seurat'?'SeuratHigh':'ScanpyHigh'));assert.ok(!texts.some(t=>/shown|gene labels|Full range|Low-y zoom|Higher in cKO|Lower in cKO|Below thresholds|omitted|above view/.test(t)));assert.ok(!clean.calls.some(x=>x.method==='stroke'&&x.strokeStyle==='#687e8d'));assert.equal(clean.calls.filter(x=>x.method==='arc').length,screen.calls.filter(x=>x.method==='arc').length);});
      app.state.volcanoCleanPDF=false;const detailed=await app.exportPDF('volcano',tool),last=renders.at(-1);
      check(`${branch} ${tool} ${view}: optional detailed PDF retains provenance`,()=>{assert.ok(detailed.texts[0].includes(tool+' | Inhibitory | '+branch));assert.ok(detailed.texts.some(t=>t.includes('Saved fit '+fit.id)));assert.ok(detailed.texts.some(t=>t.includes('independent mice')));assert.equal(last.result.summaryText,true);assert.equal(last.result.labelLeaderLines,false);});
    }
    check(branch+': all PDF presentation operations preserve both exact memberships, source fits and rows',()=>assert.equal(frozenScience(app),branchBefore));
  }
  app.state.branch='original_full';await app.changeScope();await app.ensureFits();await app.refreshDERows();
  check('Returning to original branch restores exactly the source science snapshot',()=>assert.equal(frozenScience(app),sourceBefore));
  for(const tool of tools){app.state.deTargets[tool]='cluster:0';app.state.volcanoCleanPDF=true;const pdf=await app.exportPDF('volcano',tool);check(tool+': clean cluster PDF identifies the exact cluster fit in filename and hidden provenance',()=>{assert.ok(pdf.filename.includes('_C0_Pvalb_Gaba_'));assert.equal(JSON.parse(pdf.properties.subject).source_fit,tool+'_original_full_C0');assert.deepEqual(pdf.texts,['Inhibitory']);assert.deepEqual(plain(app.targetMembers(tool)),[0,2]);});}
  const savedCount=()=>documents.filter(doc=>doc.filename).length,oldSaved=savedCount();app.state.exclusions={'scanpy|Inhibitory|original_frozen':['0']};
  check('Excluded cluster makes old bound DE rows unavailable immediately',()=>assert.equal(app.boundDERows('scanpy'),null));
  const unmatched=await app.exportPDF('volcano','scanpy');check('Unmatched selection never downloads a clean volcano with stale source rows',()=>{assert.equal(unmatched,undefined);assert.equal(savedCount(),oldSaved);});
  for(const [population,title]of Object.entries(types)){const x=await context(population),pdf=await x.app.exportPDF('volcano','scanpy');check(population+': correct nucleus-type-only PDF title',()=>assert.deepEqual(pdf.texts,[title]));}
  let release,startedResolve;const started=new Promise(resolve=>startedResolve=resolve),gate={id:'scanpy_original_full_total_rows',start:startedResolve,promise:new Promise(resolve=>release=resolve)};
  const waiting=await context('Inhibitory',gate),pending=waiting.app.exportPDF('volcano','scanpy');await started;
  check('PDF construction waits until the selected exact DE source rows finish loading',()=>assert.equal(waiting.documents.length,0));
  waiting.app.state.volcanoLowYMax=2;release();const canceled=await pending;
  check('Changing a plot control during loading cancels export before any PDF construction',()=>{assert.equal(canceled,undefined);assert.equal(waiting.documents.length,0);});
  console.log(JSON.stringify({status:'PASS',checks:checks.length,scope:'ACTUAL_APP_EXPORT_AND_SHARED_RENDERER_VM_ONLY_NO_NEW_SCIENCE'}));
})().catch(error=>{console.error(error);process.exitCode=1;});
