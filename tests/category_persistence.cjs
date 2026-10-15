'use strict';
// Five category shortcuts must persist before async loading and survive reload.
const fs=require('fs'),path=require('path'),vm=require('vm'),assert=require('assert'),crypto=require('crypto');
const site=process.env.SATB2_SITE_DIR||path.resolve(__dirname,'../site'),checks=[],trace=[];let stored=null;
class Element{constructor(tag='input'){this.tagName=tag;this.children=[];this.value='';this.textContent='';this.dataset={};this.style={};}append(...x){this.children.push(...x);}replaceChildren(...x){this.children=x;}addEventListener(){}get options(){return this.children;}}
const ids={},categories=['allcells','Excitatory','CA1_ProS','CA1_programs','Inhibitory'].map(id=>{const b=new Element('button');b.dataset.category=id;return b;});
const c={console,setTimeout,clearTimeout,SATB2_NO_AUTO_INIT:true,document:{getElementById:id=>ids[id]||(ids[id]=new Element()),createElement:t=>new Element(t),querySelectorAll:q=>q==='[data-category]'?categories:[]},localStorage:{getItem:()=>stored,setItem:(key,value)=>{stored=value;trace.push('persist:'+JSON.parse(value).state.population);}}};c.window=c;c.globalThis=c;vm.createContext(c);
vm.runInContext(fs.readFileSync(path.join(site,'v2_features.js'),'utf8'),c);const source=fs.readFileSync(path.join(site,'app.js'),'utf8');vm.runInContext(source.replace('const api={','const api={testBindV2:bindV2Controls,'),c);const app=c.SATB2_APP;
const pops=['allcells','Excitatory','CA1_ProS','Inhibitory'],datasets=pops.flatMap(pop=>['scanpy','seurat'].map(tool=>({tool,population:pop,branch:'original_full',cluster_scope:'original_frozen',coordinates_shard:tool+'_'+pop}))),shards=Object.fromEntries(datasets.map(d=>[d.coordinates_shard,{}]));
const manifest={schema_version:'satb2-offline-1',study:{source_fingerprint:'category_fixture'},defaults:{population:'CA1_ProS',branch:'original_full',marker_ids:[]},genes:[{gene_id:'g1',gene_symbol:'Calb1'}],datasets,shards};
const loader=async id=>{trace.push('load:'+id);return {encoding:'json',data:{nuclei:[{nucleus_id:'mouse_B',cluster_id:'0',x_harmony:1,y_harmony:2,x_unintegrated:1,y_unintegrated:2}],clusters:[{cluster_id:'0',reviewed_label:'Fixture'}]}};};
function check(name,fn){fn();checks.push(name);console.log('PASS '+name);}
(async()=>{
 await app.init({manifest,catalog:{datasets:[],shards:{}},headless:true,ignoreSavedState:true,shardLoader:loader});app.testBindV2();
 for(const b of categories){trace.length=0;const expected=b.dataset.category==='CA1_programs'?'CA1_ProS':b.dataset.category,promise=b.onclick();
  check(b.dataset.category+' persists before loading',()=>{const saved=JSON.parse(stored);assert.equal(saved.state.population,expected);assert.equal(saved.state.annotation,b.dataset.category==='CA1_programs'?'layer_program_label':'reviewed_label');assert.equal(trace[0],'persist:'+expected);});
  await promise;await app.init({manifest,catalog:{datasets:[],shards:{}},headless:true,shardLoader:loader});
  check(b.dataset.category+' reload restores chosen scope',()=>{assert.equal(app.state.population,expected);assert.equal(app.state.annotation,b.dataset.category==='CA1_programs'?'layer_program_label':'reviewed_label');assert.equal(app.active.scanpy.descriptor.population,expected);assert.equal(app.active.seurat.descriptor.population,expected);});
 }
 const report={status:'PASS',scope:'PURE_VM_CATEGORY_STATE_ASYNC_ORDER_NO_BROWSER',checks:checks.length,check_names:checks,app_sha256:crypto.createHash('sha256').update(source).digest('hex')};fs.writeFileSync(path.join(__dirname,'category_persistence_acceptance.json'),JSON.stringify(report,null,2)+'\n');
})().catch(e=>{console.error(e);process.exitCode=1;});
