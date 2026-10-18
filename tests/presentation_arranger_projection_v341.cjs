'use strict';
// Execute the actual v3.4 arranger handlers and canvas/PDF path with frozen
// fixtures. This is not a browser, screenshot, statistical fit or data export.
const fs = require('fs'), path = require('path'), vm = require('vm'), assert = require('assert'), crypto = require('crypto');
const site = process.env.SATB2_SITE_DIR || path.resolve(__dirname, '../site'), html = fs.readFileSync(path.join(site, 'index.html'), 'utf8'), source = fs.readFileSync(path.join(site, 'app.js'), 'utf8');
const checks = [], plain = x => JSON.parse(JSON.stringify(x)), near = (a, b, eps = 1e-8) => assert(Math.abs(a - b) < eps, `${a} != ${b}`);
class Paint {
  constructor() {this.texts = []; this.rects = []; this.strokes = []; this.path = []; this.globalAlpha = 1;}
  setTransform() {} clearRect() {this.texts = []; this.rects = []; this.strokes = [];}
  fillRect(x,y,w,h) {this.rects.push({x,y,w,h,color:this.fillStyle});} fillText(text,x,y) {this.texts.push({text:String(text),x,y});} strokeText() {} save() {} restore() {} translate() {} rotate() {} setLineDash() {}
  beginPath() {this.path = [];} moveTo(x,y) {this.path.push({x,y});} lineTo(x,y) {this.path.push({x,y});} stroke() {this.strokes.push(this.path.slice());} arc() {} fill() {} measureText(text) {return {width:String(text).length*5};}
}
class Element {
  constructor(tag = 'input', id = '') {this.tagName=tag;this.id=id;this.value='';this.checked=false;this.disabled=false;this.children=[];this.style={};this.dataset={};this.attributes={};this._text='';this.clientWidth=700;this.paint=new Paint();}
  append(...nodes) {this.children.push(...nodes);} replaceChildren(...nodes) {this.children=[...nodes];this._text='';}
  addEventListener(name,fn) {this[name]=fn;} setAttribute(name,value) {this.attributes[name]=String(value);} getAttribute(name) {return this.attributes[name]??null;}
  get textContent() {return this._text + this.children.map(c=>c.textContent||'').join('');} set textContent(value) {this._text=String(value);this.children=[];}
  get options() {return this.children;} getContext() {return this.paint;} getBoundingClientRect() {return {left:0,top:0,width:this.clientWidth,height:430};} setPointerCapture() {} toDataURL() {return 'data:image/png;base64,fixture';}
}
const controls={}, presets=[], nudges=[], resets=[];
for (const match of html.matchAll(/<([a-z]+)\b([^>]*)>/gi)) {
  const attrs=match[2], id=attrs.match(/\bid="([^"]+)"/)?.[1], preset=attrs.match(/\bdata-layout-preset="([^"]+)"/)?.[1], nudge=attrs.match(/\bdata-layout-nudge="([^"]+)"/)?.[1], reset=attrs.match(/\bdata-reset-zoom="([^"]+)"/)?.[1];
  if (!id && !preset && !nudge && !reset) continue;
  const element=new Element(match[1],id||''); if(id) controls[id]=element;
  const value=attrs.match(/\bvalue="([^"]*)"/)?.[1];if(value!==undefined)element.value=value;
  if(reset){element.dataset.resetZoom=reset;resets.push(element);}if(preset){element.dataset.layoutPreset=preset;presets.push(element);}if(nudge){element.dataset.layoutNudge=nudge;nudges.push(element);}
}
const document={activeElement:null,getElementById:id=>controls[id]||(controls[id]=new Element('input',id)),createElement:tag=>new Element(tag),querySelectorAll:selector=>selector==='[data-layout-preset]'?presets:selector==='[data-layout-nudge]'?nudges:selector==='[data-reset-zoom]'?resets:[]};
class PDF {constructor(options){this.options=options;this.texts=[];this.images=[];this.pages=1;this.internal={pageSize:{getWidth:()=>options.format[0],getHeight:()=>options.format[1]}};}setFontSize(){}text(value){this.texts.push(String(value));}addImage(...args){this.images.push(args);}addPage(){this.pages++;}splitTextToSize(value){return [String(value)];}setFillColor(){}rect(){}save(filename){this.filename=filename;}}
const context={console,document,setTimeout,clearTimeout,localStorage:{getItem:()=>null,setItem(){}},SATB2_NO_AUTO_INIT:true,jspdf:{jsPDF:PDF},SATB2_TEST_DOWNLOADS:[]};context.window=context;context.globalThis=context;vm.createContext(context);
for(const file of ['v2_features.js','display_controls.js','volcano_controls.js','presentation_layouts.js','marker_styling.js'])vm.runInContext(fs.readFileSync(path.join(site,file),'utf8'),context);
const instrumented=source.replace('const api={','const api={testRenderUMAP:renderUMAP,testBindCanvas:bindCanvas,testBindControls:bindControls,testSetHeadless:value=>headless=value,testSyncView:syncViewControls,')
 .replace("function download(name,text,mime='text/plain;charset=utf-8') {","function download(name,text,mime='text/plain;charset=utf-8') {global.SATB2_TEST_DOWNLOADS.push({name,text,mime});");
assert.notEqual(instrumented,source);vm.runInContext(instrumented,context);
const app=context.SATB2_APP,zlib=require('zlib'),nativeManifest=JSON.parse(fs.readFileSync(path.join(site,'data/manifest.json'))),nativeCatalog=JSON.parse(fs.readFileSync(path.join(site,'data/DE_CATALOG.json'))),tools=['scanpy','seurat'];
const testManifest={...nativeManifest,defaults:{...nativeManifest.defaults,population:'Inhibitory',branch:'original_full',annotation:'cluster_id',tab:'umap'}};
const loaded={};for(const tool of tools){const ds=nativeManifest.datasets.find(d=>d.tool===tool&&d.population==='Inhibitory'&&d.branch==='original_full');loaded[ds.coordinates_shard]=JSON.parse(zlib.gunzipSync(fs.readFileSync(path.join(site,nativeManifest.shards[ds.coordinates_shard].path))));}
function check(name,fn){fn();checks.push(name);console.log('PASS '+name);}function preset(mode){presets.find(b=>b.dataset.layoutPreset===mode).onclick();}function nudge(direction){nudges.find(b=>b.dataset.layoutNudge===direction).onclick();}
function render(tool='seurat'){app.testRenderUMAP(tool);return new Map(app.active[tool].screenPoints.map(p=>[p.index,{...p}]));}
function hashes(){return tools.map(tool=>crypto.createHash('sha256').update(JSON.stringify(app.active[tool].coords)).digest('hex'));}
function compareDelta(before,after,expectedX,expectedY,source){let c4=0,other=0;for(const[index,p]of before){const q=after.get(index);if(!q)continue;if(String(source[index].cluster_id)==='4'){near(q.x-p.x,expectedX);near(q.y-p.y,expectedY);c4++;}else{near(q.x,p.x);near(q.y,p.y);other++;}}assert(c4>500,'Full C4 body remains measurable');assert(other>5000,'Unselected nuclei retain their screen positions');return{C4_checked:c4,other_nuclei_checked:other};}
const metrics=[];
(async()=>{
 await app.init({manifest:testManifest,catalog:nativeCatalog,headless:true,ignoreSavedState:true,shardLoader:async id=>{assert(loaded[id],'Only frozen native coordinate shards are needed');return{encoding:'json',data:plain(loaded[id])};}});app.testBindControls();app.bindPresentationControls();controls['layout-tool'].value='seurat';controls['layout-cluster'].value='4';controls['layout-nudge-step'].value='5';const frozen=hashes(),sourceNuclei=app.active.seurat.coords.nuclei;
 check('Regression uses the original native 8,911 inhibitory nuclei in each tool with source Seurat C4',()=>{assert.equal(sourceNuclei.length,8911);assert.equal(app.active.scanpy.coords.nuclei.length,8911);assert.equal(sourceNuclei.filter(n=>String(n.cluster_id)==='4').length,593);assert.equal(resets.length,2);});
 for(const mode of ['original','packed','rounded']){
  preset(mode);controls['layout-cluster'].value='4';const layout=app.presentationCoordinates('seurat'),b=layout.bounds,unit=Math.min((700-56)/(b.xmax-b.xmin),(430-56)/(b.ymax-b.ymin)),before=render();const otherBefore=render('scanpy');nudge('up');const after=render(),expected=-.05*layout.source_span*unit;
  check(mode+': actual C4 glyphs move upward exactly 5% of source span without auto-fit cancellation or X drift',()=>{const count=compareDelta(before,after,0,expected,sourceNuclei);metrics.push({mode,action:'up',expected_pixel_delta_y:expected,actual_first_C4_delta_y:after.get([...before.keys()].find(i=>String(sourceNuclei[i].cluster_id)==='4')).y-before.get([...before.keys()].find(i=>String(sourceNuclei[i].cluster_id)==='4')).y,...count});});
  check(mode+': moving Seurat C4 preserves all Scanpy glyph screen positions',()=>{const otherAfter=render('scanpy');assert.equal(otherAfter.size,otherBefore.size);for(const[index,p]of otherBefore){const q=otherAfter.get(index);near(p.x,q.x);near(p.y,q.y);}});
  const beforeRight=render();nudge('right');const afterRight=render();check(mode+': actual right nudge follows screen X while retaining every other cluster position',()=>{compareDelta(beforeRight,afterRight,.05*layout.source_span*unit,0,sourceNuclei);});
 }
 preset('original');controls['layout-cluster'].value='4';const layout=app.presentationCoordinates('seurat'),unit=Math.min((700-56)/(layout.bounds.xmax-layout.bounds.xmin),(430-56)/(layout.bounds.ymax-layout.bounds.ymin)),numericBefore=render();controls['layout-dx'].value='0';controls['layout-dy'].value='5';controls['layout-scale'].value='100';controls['apply-layout-cluster'].onclick();const numericAfter=render();
 check('Visible numeric Apply uses the same fixed projection as arrow movement for boundary C4',()=>{compareDelta(numericBefore,numericAfter,0,-.05*layout.source_span*unit,sourceNuclei);});
 preset('original');controls['layout-cluster'].value='4';const canvas=controls['umap-seurat'];canvas.wheel({deltaY:-Math.log(.8)/.001,preventDefault(){}});canvas.pointerdown({clientX:0,clientY:0,pointerId:1});canvas.pointermove({clientX:7,clientY:-8,pointerId:1});canvas.pointerup({clientX:7,clientY:-8,pointerId:1});const pannedBefore=render();nudge('up');const pannedAfter=render();
 check('Manual movement preserves the current pan and zoom rather than resetting either viewport control',()=>{compareDelta(pannedBefore,pannedAfter,0,-.05*layout.source_span*unit*.8,sourceNuclei);});
 const resetButton=resets.find(b=>b.dataset.resetZoom==='seurat');resetButton.onclick();const refit=render(),current=app.presentationCoordinates('seurat'),cb=current.bounds,currentUnit=Math.min((700-56)/(cb.xmax-cb.xmin),(430-56)/(cb.ymax-cb.ymin));
 check('Actual Reset zoom handler releases the frozen frame and fits the current adjusted display',()=>{for(const[index,p]of refit){const q=current.points[index];near(p.x,350+(q.x-(cb.xmin+cb.xmax)/2)*currentUnit);near(p.y,215-(q.y-(cb.ymin+cb.ymax)/2)*currentUnit);}});
 // A second move freezes the refitted display; rotation changes the frame key.
 nudge('up');controls['layout-rotation'].value='35';controls['layout-rotation'].oninput();const rotated=app.presentationCoordinates('seurat'),rb=rotated.bounds,rotatedUnit=Math.min((700-56)/(rb.xmax-rb.xmin),(430-56)/(rb.ymax-rb.ymin)),rotatedScreen=render();
 check('Changing rotation invalidates the old projection frame and fits the rotated layout',()=>{for(const[index,p]of rotatedScreen){const q=rotated.points[index];near(p.x,350+(q.x-(rb.xmin+rb.xmax)/2)*rotatedUnit);near(p.y,215-(q.y-(rb.ymin+rb.ymax)/2)*rotatedUnit);}});
 preset('original');const reset=app.presentationCoordinates('seurat');check('Restore original recovers exact native geometry after all screen-level movement tests',()=>{assert.equal(reset.display_only,false);for(let i=0;i<sourceNuclei.length;i++){assert.equal(reset.points[i].x,sourceNuclei[i].x_harmony);assert.equal(reset.points[i].y,sourceNuclei[i].y_harmony);}assert.deepEqual(hashes(),frozen);});
 const report={status:'PASS',scope:'ACTUAL_NATIVE_BOUNDARY_C4_CANVAS_PROJECTION_REGRESSION_NO_SOURCE_DATA_EDITS',checks:checks.length,check_names:checks,native_source_hashes:frozen,metrics,app_sha256:crypto.createHash('sha256').update(source).digest('hex')};fs.writeFileSync(path.join(__dirname,'presentation_arranger_projection_v341_acceptance.json'),JSON.stringify(report,null,2)+'\n');console.log(JSON.stringify(report));
})().catch(error=>{console.error(error);process.exitCode=1;});
