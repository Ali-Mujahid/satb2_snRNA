'use strict';
// Regression for the actual-browser Calb1 selection defect. Pure VM; no network/browser.
const fs=require('fs'),path=require('path'),vm=require('vm'),assert=require('assert'),crypto=require('crypto');
const site=process.env.SATB2_SITE_DIR||path.resolve(__dirname,'../site'),checks=[];
class Element{
 constructor(tag='input'){this.tagName=tag;this.children=[];this.listeners={};this._value='';this.textContent='';this.style={};}
 set value(x){this._value=this.tagName==='select'&&!this.children.some(c=>c.value===String(x))?'':String(x);}
 get value(){return this._value;}
 append(...x){this.children.push(...x);}
 replaceChildren(...x){this.children=[...x];if(this.tagName==='select'&&!this.children.some(c=>c.value===this._value))this._value=this.children[0]?.value||'';}
 addEventListener(k,f){this.listeners[k]=f;}
 get options(){return this.children;}
}
const controls={};const document={getElementById:id=>controls[id]||(controls[id]=new Element(id==='overlay-gene'?'select':'input')),createElement:tag=>new Element(tag),querySelectorAll:()=>[]};
const c={console,document,setTimeout,clearTimeout,localStorage:{getItem:()=>null,setItem(){}},SATB2_NO_AUTO_INIT:true};c.window=c;c.globalThis=c;vm.createContext(c);
vm.runInContext(fs.readFileSync(path.join(site,'v2_features.js'),'utf8'),c);
const original=fs.readFileSync(path.join(site,'app.js'),'utf8');
// Test-only private hooks inserted into the VM string, never written to production source.
vm.runInContext(original.replace('const api={','const api={testBindV2:bindV2Controls,testSync:syncControls,testFilter:filterOverlayOptions,'),c);
const app=c.SATB2_APP,coords={nuclei:[{nucleus_id:'sample_A',sample:'sample',condition:'WT',cluster_id:'0',x_harmony:1,y_harmony:2,x_unintegrated:1,y_unintegrated:2}],clusters:[{cluster_id:'0',reviewed_label:'Fixture'}]};
const genes=[{gene_id:'calb1',gene_symbol:'Calb1'},{gene_id:'etv1',gene_symbol:'Etv1'},{gene_id:'cck',gene_symbol:'Cck'}];
const manifest={schema_version:'satb2-offline-1',defaults:{population:'CA1_ProS',branch:'original_full',marker_ids:[]},study:{source_fingerprint:'fixture'},genes,shards:{scanpycoords:{},seuratcoords:{},scanpymarkers:{},seuratmarkers:{}},datasets:['scanpy','seurat'].map(t=>({tool:t,branch:'original_full',population:'CA1_ProS',cluster_scope:'original_frozen',coordinates_shard:t+'coords',markers_shard:t+'markers'}))};
function check(name,fn){fn();checks.push(name);console.log('PASS '+name);}
(async()=>{
 await app.init({manifest,catalog:{datasets:[],shards:{}},headless:true,ignoreSavedState:true,shardLoader:async id=>({encoding:'json',data:id.endsWith('coords')?coords:{rows:[]}})});assert(app.active.scanpy?.coords.nuclei.length===1&&app.active.seurat?.coords.nuclei.length===1,'Fixture coordinates must load before UI regression');app.testBindV2();
 document.getElementById('overlay-search').value='Calb1';app.state.overlayGene='calb1';app.testFilter();
 check('Filtered dropdown selects Calb1',()=>{assert.equal(controls['overlay-gene'].value,'calb1');assert.equal(controls['overlay-gene'].options.length,1);});
 const members=JSON.stringify(app.targetMembers('scanpy'));
 await controls['show-gene-overlay'].onclick();
 check('Show-selected button retains Calb1 rather than preferred Etv1',()=>{assert.equal(app.state.overlayGene,'calb1');assert.equal(controls['overlay-gene'].value,'calb1');assert.equal(app.state.annotation,'gene');assert.equal(app.state.tab,'umap');});
 app.state.overlayGene='etv1';app.testSync();
 check('Sync keeps current gene visible outside search matches',()=>{assert.equal(controls['overlay-gene'].value,'etv1');assert(controls['overlay-gene'].options.some(o=>o.value==='etv1'));assert.equal(controls['overlay-search'].value,'Calb1');});
 app.state.overlayGene='invalid';await controls['show-gene-overlay'].onclick();
 check('Invalid selection falls back to an available preferred gene without blank select',()=>{assert.equal(app.state.overlayGene,'etv1');assert.equal(controls['overlay-gene'].value,'etv1');});
 controls['overlay-search'].value='NoMatches';app.testFilter();
 check('No search matches still retain the selected current overlay',()=>{assert.equal(controls['overlay-gene'].value,'etv1');assert.equal(controls['overlay-gene'].options.length,1);});
 check('Gene display changes preserve inference membership',()=>assert.equal(JSON.stringify(app.targetMembers('scanpy')),members));
 const report={status:'PASS',scope:'PURE_VM_NATIVE_SELECT_SEMANTICS_NO_BROWSER',checks:checks.length,check_names:checks,app_sha256:crypto.createHash('sha256').update(original).digest('hex')};fs.writeFileSync(path.join(__dirname,'overlay_selection_acceptance.json'),JSON.stringify(report,null,2)+'\n');
})().catch(e=>{console.error(e);process.exitCode=1;});
