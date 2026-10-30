process.chdir(require('path').join(__dirname, '../site/communication/ca1_partner'));
/* Dependency-free synthetic DOM smoke test. This is not a browser or PDF-file test. */
const assert = require('assert'), fs = require('fs'), vm = require('vm');
let downloads = 0, prints = 0, windows = 0;
class Element {
  constructor(tag='div') {this.tag=tag;this.children=[];this.listeners={};this.attrs={};this.value='';this.checked=false;this.style={};this.textContent='';}
  setAttribute(k,v) {this.attrs[k]=String(v)} removeAttribute(k){delete this.attrs[k]}
  appendChild(x){this.children.push(x);return x} replaceChildren(...x){this.children=x}
  addEventListener(k,f){this.listeners[k]=f} remove(){} click(){if(this.tag==='a')downloads++;if(this.onclick)this.onclick();if(this.listeners.click)this.listeners.click({target:this});}
  fire(name){if(this.listeners[name])this.listeners[name]({target:this})}
  cloneNode(){const e=new Element(this.tag);e.attrs={...this.attrs};return e}
}
const html=fs.readFileSync('index.html','utf8'), els={};
for(const m of html.matchAll(/id="([^"]+)"/g)) els[m[1]]=new Element();
Object.assign(els.threshold,{value:'5'}); Object.assign(els['mouse-metric'],{value:'detection_fraction'});Object.assign(els.ca1color,{value:'#245e9e'});Object.assign(els.incolor,{value:'#b64f48'});Object.assign(els.dotsize,{value:'5'});Object.assign(els.majority,{type:'checkbox'});Object.assign(els.replicated,{type:'checkbox'});
const body=new Element('body'), doc={getElementById:k=>{assert.ok(els[k],`unknown ID ${k}`);return els[k]},createElement:t=>new Element(t),createElementNS:(_,t)=>new Element(t),body};
const makePair=(id,direction,p5,p10,majority5,majority10)=>({record_id:id,source:'CA1',target:'IN',ligand_complex:id==='a'?'Slit1':'B',receptor_complex:id==='a'?'Robo1':'C',direction,ca1_population:'CA1_Deep',in_population:'IN_C6',ca1_significant_genes:'Slit1',in_significant_genes:'',in_status:'CA1 change only detected',passes_5pct:p5,passes_10pct:p10,majority_mouse_support_5pct:majority5,majority_mouse_support_10pct:majority10,wt_supported_mice:'4',wt_supported_mice_10pct:'2',wt_eligible_mice:'6',cko_supported_mice:'3',cko_supported_mice_10pct:'1',cko_eligible_mice:'6'});
const data={metadata:{state:'reviewed_complete',provenance:'Synthetic only'},pairs:[makePair('a','CA1 to IN','True','False','True','False'),makePair('b','IN to CA1','True','True','False','True')],subunits:[{record_id:'a',endpoint:'CA1',population:'CA1_Deep',role:'ligand',gene:'Slit1',log2FoldChange:'-1',padj:'1.19752181884482e-29',de_status:'significant',WT_detection:'0.3',cKO_detection:'0.2',WT_mean:'0.1',cKO_mean:'0.05'},{record_id:'a',endpoint:'IN',population:'IN_C6',role:'receptor',gene:'<img src=x>',log2FoldChange:'0.2',padj:'',de_status:'unavailable',WT_detection:'0.2',cKO_detection:'0.1',WT_mean:'0.04',cKO_mean:'0.02'},{record_id:'a',endpoint:'IN',population:'IN_C6',role:'receptor',gene:'Zero',log2FoldChange:'0.2',padj:'0',de_status:'significant',WT_detection:'0.2',cKO_detection:'0.1',WT_mean:'0.04',cKO_mean:'0.02'}],mouse_expression:[{record_id:'a',endpoint:'CA1',gene:'Slit1',condition:'WT',sample:'w1',n_nuclei:'100',detection_fraction:'0.3',mean_expression:'0.1'}]};
const win={SATB2_CA1_IN_DATA:data,open:()=>{windows++;return {document:{title:'',body:new Element('body'),createElement:t=>new Element(t),importNode:x=>x,close(){}},focus(){},print(){prints++}}}};
const context={window:win,document:doc,globalThis:win,Map,Set,Number,String,Object,Math,Array,Blob,URL:{createObjectURL:()=>{return 'blob:synthetic'},revokeObjectURL(){}},XMLSerializer:class {serializeToString(){return '<svg/>'}},setTimeout:f=>f(),alert:msg=>{throw Error('Unexpected alert: '+msg)}};
vm.createContext(context);vm.runInContext(fs.readFileSync('model.js','utf8'),context);vm.runInContext(fs.readFileSync('app.js','utf8'),context);
assert.match(els['pair-summary'].textContent,/2 of 2/);assert.equal(els['pairs-body'].children.length,2);
els.threshold.value='10';els.threshold.fire('change');assert.match(els['pair-summary'].textContent,/1 of 2/);assert.equal(els['pairs-body'].children.length,1);
els.majority.checked=true;els.majority.fire('change');assert.match(els['pair-summary'].textContent,/1 of 2/);
els.reset.click();assert.match(els['pair-summary'].textContent,/2 of 2/);
els['pairs-body'].children[0].children[0].children[0].click();assert.match(els['selected-title'].textContent,/Slit1/);assert.equal(els['subunit-body'].children.length,3);assert.equal(els['mouse-body'].children.length,1);assert.equal(els['subunit-body'].children[0].children[5].textContent,'1.198e-29');assert.equal(els['subunit-body'].children[1].children[5].textContent,'—');assert.equal(els['subunit-body'].children[2].children[5].textContent,'0');
assert.equal(els['subunit-body'].children[1].children[3].textContent,'<img src=x>'); // textContent, no HTML execution
const findNodes=(node,tag)=>[...(node.tag===tag?[node]:[]),...node.children.flatMap(child=>findNodes(child,tag))];assert.ok(findNodes(els['effect-plot'].children[0],'title').some(t=>t.textContent.includes('FDR 1.198e-29')));assert.ok(findNodes(els['effect-plot'].children[0],'title').some(t=>t.textContent.includes('FDR 0;')));
assert.equal(els['mouse-plot'].children[0].tag,'svg');els['mouse-metric'].value='mean_expression';els['mouse-metric'].fire('change');assert.equal(els['mouse-plot'].children[0].tag,'svg');
els['download-pairs'].click();els['download-subunits'].click();els['download-mice'].click();els['svg-main'].click();els['svg-mice'].click();assert.equal(downloads,5);
els['pdf-main'].click();els['pdf-mice'].click();assert.equal(windows,2);assert.equal(prints,2);
els.search.value='unmatched';els.search.fire('input');assert.match(els['pair-summary'].textContent,/0 of 2/);assert.equal(els['pairs-body'].children.length,0);
// Exercise the actual reviewed Lrrc4c record through the DOM event handlers.
const realWindow={};vm.runInNewContext(fs.readFileSync('data.js','utf8'),{window:realWindow});const real=realWindow.SATB2_CA1_IN_DATA, rid='0061c003e4432d6b';
win.SATB2_CA1_IN_DATA={metadata:real.metadata,pairs:real.pairs.filter(r=>r.record_id===rid),subunits:real.subunits.filter(r=>r.record_id===rid),mouse_expression:real.mouse_expression.filter(r=>r.record_id===rid)};
els.threshold.value='5';els.majority.checked=false;els.replicated.checked=false;els.search.value='';vm.runInContext(fs.readFileSync('app.js','utf8'),context);
assert.equal(els['pairs-body'].children.length,1);els['pairs-body'].children[0].children[0].children[0].click();const actualRow=els['subunit-body'].children.find(tr=>tr.children[3].textContent==='Lrrc4c');assert.ok(actualRow);assert.equal(actualRow.children[5].textContent,'1.198e-29');assert.ok(findNodes(els['effect-plot'].children[0],'title').some(t=>t.textContent.includes('Lrrc4c')&&t.textContent.includes('FDR 1.198e-29')));
console.log('synthetic and actual-record DOM behavior smoke tests passed (not a real browser or saved PDF test)');
