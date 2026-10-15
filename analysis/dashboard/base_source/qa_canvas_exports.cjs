/* Simulated DOM and canvas-command artifact rendering; never launches a browser. */
'use strict';
const fs=require('fs'),path=require('path'),vm=require('vm'),assert=require('assert'),crypto=require('crypto');
const deps='/Users/mujahid357/.cache/codex-runtimes/codex-primary-runtime/dependencies/node/node_modules';
const sharp=require(path.join(deps,'sharp'));
const base=path.dirname(__dirname),src=path.join(base,'offline_dashboard_source'),bundle=path.join(base,'reports/offline_dashboard_v1');
const out=path.join(base,'reviews/offline_canvas_export_verification');fs.mkdirSync(out,{recursive:true});
const html=fs.readFileSync(path.join(src,'index.html'),'utf8'),appBytes=fs.readFileSync(path.join(src,'app.js'));
const elements=new Map(),all=[],pngCache=new Map(),saved=[];
const esc=x=>String(x).replace(/[&<>"']/g,c=>({'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;',"'":'&apos;'}[c]));
function color(x){const m=String(x).match(/^hsl\((\d+) (\d+)% (\d+)%\)$/);if(!m)return x;const h=+m[1]/360,s=+m[2]/100,l=+m[3]/100;const f=n=>{const k=(n+h*12)%12,a=s*Math.min(l,1-l);return Math.round(255*(l-a*Math.max(-1,Math.min(k-3,9-k,1))));};return `rgb(${f(0)},${f(8)},${f(4)})`;}
class Drawing {
 constructor(canvas){this.canvas=canvas;this.parts=[];this.path=[];this.stack=[];this.transform='';this.globalAlpha=1;this.fillStyle='#000';this.strokeStyle='#000';this.font='12px sans-serif';this.textAlign='left';this.lineWidth=1;}
 setTransform(a,b,c,d,e,f){this.transform=`matrix(${a} ${b} ${c} ${d} ${e} ${f})`;}
 clearRect(){this.parts=[];}
 style(stroke=false){return `transform="${this.transform}" opacity="${this.globalAlpha}" ${stroke?'fill="none" stroke="'+color(this.strokeStyle)+'" stroke-width="'+this.lineWidth+'"':'fill="'+color(this.fillStyle)+'"'}`;}
 fillRect(x,y,w,h){this.parts.push(`<rect x="${x}" y="${y}" width="${w}" height="${h}" ${this.style()}/>`);}
 strokeRect(x,y,w,h){this.parts.push(`<rect x="${x}" y="${y}" width="${w}" height="${h}" ${this.style(true)}/>`);}
 text(text,x,y,stroke){const px=Number(this.font.match(/([0-9.]+)px/)?.[1]||12),anchor=this.textAlign==='right'?'end':this.textAlign==='center'?'middle':'start';this.parts.push(`<text x="${x}" y="${y}" font-size="${px}" font-weight="${this.font.includes('bold')?'bold':'normal'}" font-family="Arial,sans-serif" text-anchor="${anchor}" ${this.style(stroke)}>${esc(text)}</text>`);}
 fillText(text,x,y){this.text(text,x,y,false);}strokeText(text,x,y){this.text(text,x,y,true);}
 beginPath(){this.path=[];}
 arc(x,y,r){this.path.push(`<circle cx="${x}" cy="${y}" r="${r}" STYLE/>`);}
 moveTo(x,y){this.path.push(`M${x} ${y}`);}
 lineTo(x,y){this.path.push(`L${x} ${y}`);}
 finish(stroke){if(this.path[0]?.startsWith('<'))for(const p of this.path)this.parts.push(p.replace('STYLE',this.style(stroke)));else this.parts.push(`<path d="${this.path.join(' ')}" ${this.style(stroke)}/>`);}
 fill(){this.finish(false);}stroke(){this.finish(true);}closePath(){this.path.push('Z');}
 save(){this.stack.push({transform:this.transform,globalAlpha:this.globalAlpha,fillStyle:this.fillStyle,strokeStyle:this.strokeStyle,font:this.font,textAlign:this.textAlign});}
 restore(){Object.assign(this,this.stack.pop()||{});}
 translate(x,y){this.transform+=` translate(${x} ${y})`;}
 rotate(r){this.transform+=` rotate(${r*180/Math.PI})`;}
 measureText(t){return {width:String(t).length*(parseFloat(this.font)||12)*.55};}
}
class Element {
 constructor(tag,attrs={}){this.tagName=tag.toUpperCase();this.attrs=attrs;this.children=[];this.style={};this.dataset={};this.listeners={};this.hidden='hidden'in attrs;this.checked='checked'in attrs;this.disabled=false;this.value=attrs.value||'';this._text='';this.clientWidth=700;this.width=700;this.height=430;this.files=[];this.className=attrs.class||'';this.classList={toggle(){},add(){},remove(){}};for(const[k,v]of Object.entries(attrs))if(k.startsWith('data-'))this.dataset[k.slice(5).replace(/-([a-z])/g,(_,x)=>x.toUpperCase())]=v;if(attrs.id){this.id=attrs.id;elements.set(attrs.id,this);}all.push(this);}
 append(...items){this.children.push(...items);for(const x of items)if(x&&typeof x==='object')x.parentElement=this;}
 appendChild(x){this.append(x);return x;}
 replaceChildren(...items){this.children=[];this.append(...items);}
 remove(){}setAttribute(k,v){this.attrs[k]=v;}
 get options(){const result=[];const walk=n=>{for(const c of n.children||[]){if(c.tagName==='OPTION')result.push(c);walk(c);}};walk(this);return result;}
 set textContent(v){this._text=String(v);this.children=[];}get textContent(){return this._text+this.children.map(x=>x.textContent||'').join('');}
 addEventListener(k,fn){(this.listeners[k]??=[]).push(fn);}
 async dispatch(k,extra={}){const e={target:this,preventDefault(){},...extra};if(this['on'+k])await this['on'+k](e);for(const fn of this.listeners[k]||[])await fn(e);}
 click(){return this.dispatch('click');}
 getBoundingClientRect(){return {left:0,top:0,width:this.clientWidth,height:this.height};}
 setPointerCapture(){}getContext(){return this.drawing??=new Drawing(this);}
 svg(){return `<svg xmlns="http://www.w3.org/2000/svg" width="${this.width}" height="${this.height}" viewBox="0 0 ${this.width} ${this.height}">${this.getContext().parts.join('')}</svg>`;}
 toDataURL(){const key=crypto.createHash('sha256').update(this.svg()).digest('hex');assert(pngCache.has(key),'Canvas rendering must be prepared before PDF test: '+this.id);return 'data:image/png;base64,'+pngCache.get(key).toString('base64');}
}
for(const m of html.matchAll(/<([a-z][\w-]*)\b([^>]*)>/gi)){const a={};for(const k of m[2].matchAll(/([\w-]+)(?:="([^"]*)")?/g))a[k[1]]=k[2]??'';if(a.id||Object.keys(a).some(k=>k.startsWith('data-')))new Element(m[1],a);}
for(const[id,val]of Object.entries({'marker-search':'','de-search':'','marker-paste':'','marker-metric':'scaled_expression_z','marker-condition':'ALL','marker-grouping':'cluster'}))elements.get(id).value=val;
const document={readyState:'complete',getElementById:id=>elements.get(id)||null,createElement:tag=>new Element(tag),querySelectorAll:selector=>{const m=selector.match(/^\[([^\]]+)\]$/);return m?all.filter(x=>m[1]in x.attrs):[];},head:new Element('head'),body:new Element('body'),addEventListener(){},createEvent(){return {initMouseEvent(){}};}};
const context={console,document,Uint8Array,Uint16Array,Uint32Array,Int32Array,Float32Array,Float64Array,ArrayBuffer,TextEncoder,TextDecoder,atob,btoa,Blob,URL,setTimeout,clearTimeout,performance,
 localStorage:{getItem(){return null;},setItem(){}},navigator:{userAgent:'Pure simulated DOM artifact verification'},devicePixelRatio:1,addEventListener(){},SATB2_NO_AUTO_INIT:true};
context.window=context;context.self=context;context.globalThis=context;vm.createContext(context);
for(const name of ['vendor/fflate.min.js','vendor/jspdf.umd.min.js'])vm.runInContext(fs.readFileSync(path.join(src,name),'utf8'),context,{filename:name});
vm.runInContext(appBytes.toString('utf8'),context,{filename:'app.js'});
const NativePDF=context.jspdf.jsPDF;context.jspdf.jsPDF=function(options){const doc=new NativePDF(options);doc.save=name=>{const b=Buffer.from(doc.output('arraybuffer'));fs.writeFileSync(path.join(out,name),b);saved.push({name,bytes:b.length,pages:doc.getNumberOfPages()});};return doc;};
const app=context.SATB2_APP;
function loader(id,desc){const text=fs.readFileSync(path.join(bundle,desc.path),'utf8'),m=text.match(/(?:\bdata\s*:|"data"\s*:)\s*("[A-Za-z0-9+/=]+")/);return {encoding:'gzip-base64',data:JSON.parse(m[1])};}
async function prepareCanvas(id,label){const canvas=elements.get(id),svg=canvas.svg(),key=crypto.createHash('sha256').update(svg).digest('hex');const data=await sharp(Buffer.from(svg)).png().toBuffer();pngCache.set(key,data);fs.writeFileSync(path.join(out,label+'.png'),data);return {id,width:canvas.width,height:canvas.height,draw_commands:canvas.getContext().parts.length};}
(async()=>{
 const manifest=JSON.parse(fs.readFileSync(path.join(bundle,'data/manifest.json'))),catalog=JSON.parse(fs.readFileSync(path.join(bundle,'DE_CATALOG.json')));
 await app.init({manifest,catalog,ignoreSavedState:true,shardLoader:loader});assert(app.active.scanpy&&app.active.seurat);
 const canvases=[];
 for(const tool of ['scanpy','seurat']){canvases.push(await prepareCanvas('umap-'+tool,'actual_CA1_UMAP_'+tool));await app.exportPDF('umap',tool);}
 app.state.tab='markers';await app.ensureMarkers();app.render();
 for(const tool of ['scanpy','seurat']){canvases.push(await prepareCanvas('dots-'+tool,'actual_CA1_marker_dots_'+tool));await app.exportPDF('dots',tool);}
 app.state.tab='de';await app.ensureFits();await app.refreshDERows();app.render();
 for(const tool of ['scanpy','seurat']){canvases.push(await prepareCanvas('volcano-'+tool,'actual_CA1_volcano_'+tool));await app.exportPDF('volcano',tool);}
 canvases.push(await prepareCanvas('de-compare','actual_CA1_effect_compare'));await app.exportPDF('compare','both');
 app.state.tab='clusters';app.render();assert(elements.get('clusters-scanpy').children.length>0);
 const none=all.find(x=>x.dataset.clusters==='scanpy:none');await none.click();
 assert.equal(app.retainedIndices('scanpy').length,0);
 app.state.tab='umap';app.render();await prepareCanvas('umap-scanpy','empty_selection_message');
 assert(elements.get('umap-scanpy').getContext().parts.some(x=>x.includes('No nuclei selected')));
 const restored=all.find(x=>x.dataset.clusters==='scanpy:all');await restored.click();assert(app.retainedIndices('scanpy').length>0);
 assert.equal(crypto.createHash('sha256').update(fs.readFileSync(path.join(src,'app.js'))).digest('hex'),crypto.createHash('sha256').update(appBytes).digest('hex'),'Frontend changed during QA');
 const receipt={status:'PASS',validation_mode:'Simulated DOM plus canvas-command-to-PNG rendering and actual bundled jsPDF export',browser_rendering_performed:false,
  source_app_sha256:crypto.createHash('sha256').update(appBytes).digest('hex'),canvases,pdfs:saved,limitations:['This checks the application drawing commands and export artifacts, not real browser/CSS rendering.']};
 fs.writeFileSync(path.join(out,'COMPLETE.json'),JSON.stringify(receipt,null,2)+'\n');console.log(JSON.stringify(receipt));
})().catch(e=>{console.error(e);process.exitCode=1;});
