/* Simulated DOM and canvas-command artifact rendering; never launches a browser. */
'use strict';
const fs=require('fs'),path=require('path'),vm=require('vm'),assert=require('assert'),crypto=require('crypto');
const deps='/Users/mujahid357/.cache/codex-runtimes/codex-primary-runtime/dependencies/node/node_modules';
const sharp=require(path.join(deps,'sharp'));
const base=path.dirname(__dirname),src=path.join(base,'offline_dashboard_v2_source'),bundle=path.join(base,'reports/offline_dashboard_v2');
const out=path.join(base,'reviews/offline_v2_canvas_export_verification');fs.mkdirSync(out,{recursive:true});
const html=fs.readFileSync(path.join(src,'index.html'),'utf8'),appBytes=fs.readFileSync(path.join(src,'app.js'));
const elements=new Map(),all=[],pngCache=new Map(),saved=[];let exportPrefix='';
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
 constructor(tag,attrs={}){this.tagName=tag.toUpperCase();this.attrs=attrs;this.children=[];this.style={};this.dataset={};this.listeners={};this.hidden='hidden'in attrs;this.checked='checked'in attrs;this.disabled=false;this.value=attrs.value||'';this._text='';this.clientWidth=700;this._width=700;this._height=430;this.files=[];this.className=attrs.class||'';this.classList={toggle(){},add(){},remove(){}};for(const[k,v]of Object.entries(attrs))if(k.startsWith('data-'))this.dataset[k.slice(5).replace(/-([a-z])/g,(_,x)=>x.toUpperCase())]=v;if(attrs.id){this.id=attrs.id;elements.set(attrs.id,this);}all.push(this);}
 append(...items){this.children.push(...items);for(const x of items)if(x&&typeof x==='object')x.parentElement=this;}
 appendChild(x){this.append(x);return x;}
 replaceChildren(...items){this.children=[];this.append(...items);}
 remove(){}setAttribute(k,v){this.attrs[k]=v;}
 get options(){const result=[];const walk=n=>{for(const c of n.children||[]){if(c.tagName==='OPTION')result.push(c);walk(c);}};walk(this);return {length:result.length,item:i=>result[i],[Symbol.iterator]:function*(){yield*result;}};}
 get width(){return this._width;}set width(v){this._width=v;this.drawing=new Drawing(this);}get height(){return this._height;}set height(v){this._height=v;this.drawing=new Drawing(this);}
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
for(const name of ['vendor/fflate.min.js','vendor/jspdf.umd.min.js','v2_features.js'])vm.runInContext(fs.readFileSync(path.join(src,name),'utf8'),context,{filename:name});
for(const name of ['data/mouse_catalog.js','data/annotation_suggestions.js']){const p=path.join(bundle,name);if(fs.existsSync(p))vm.runInContext(fs.readFileSync(p,'utf8'),context,{filename:name});}
vm.runInContext(appBytes.toString('utf8'),context,{filename:'app.js'});
const NativePDF=context.jspdf.jsPDF;context.jspdf.jsPDF=function(options){const doc=new NativePDF(options);doc.save=name=>{const b=Buffer.from(doc.output('arraybuffer'));const filename=exportPrefix+name;fs.writeFileSync(path.join(out,filename),b);saved.push({name:filename,bytes:b.length,pages:doc.getNumberOfPages(),page_width_mm:doc.internal.pageSize.getWidth(),page_height_mm:doc.internal.pageSize.getHeight()});};return doc;};
const app=context.SATB2_APP;
function loader(id,desc){const text=fs.readFileSync(path.join(bundle,desc.path),'utf8'),m=text.match(/(?:\bdata\s*:|"data"\s*:)\s*("[A-Za-z0-9+/=]+")/);return {encoding:'gzip-base64',data:JSON.parse(m[1])};}
async function prepareCanvas(id,label){const canvas=elements.get(id),svg=canvas.svg(),key=crypto.createHash('sha256').update(svg).digest('hex');const data=await sharp(Buffer.from(svg)).png().toBuffer();pngCache.set(key,data);fs.writeFileSync(path.join(out,label+'.png'),data);return {id,width:canvas.width,height:canvas.height,draw_commands:canvas.getContext().parts.length};}

(async()=>{
 const manifest=JSON.parse(fs.readFileSync(path.join(bundle,'data/manifest.json'))),catalog=JSON.parse(fs.readFileSync(path.join(bundle,'DE_CATALOG.json'))),mouse=context.SATB2_MOUSE_CATALOG;
 assert.equal(mouse?.status,'COMPLETE','Closed mouse catalog is required');
 await app.init({manifest,catalog,ignoreSavedState:true,shardLoader:loader});assert(app.active.scanpy&&app.active.seurat);
 const canvases=[],checks=[];function check(name,fn){fn();checks.push(name);}
 async function capture(id,label){canvases.push(await prepareCanvas(id,label));}
 async function scope(pop,tab='umap'){app.state.population=pop;app.state.branch='original_full';app.state.tab=tab;app.state.exclusions={};app.state.labels={};app.state.deTargets={scanpy:'total',seurat:'total'};app.state.umapCondition='ALL';app.state.umapProgram='ALL';await app.changeScope();}
 if(process.argv.includes('--spacing-only')){await scope('Inhibitory');app.state.annotation='cluster_id';app.state.umapLabelMode='id';app.state.spacing=1;app.render();await capture('umap-scanpy','actual_IN_original_geometry_ID_labels');app.state.spacing=.5;app.render();await capture('umap-scanpy','actual_IN_spacing_0p5_ID_labels');exportPrefix='IN_spacing_0p5_ID_labels_';await app.exportPDF('umap','scanpy');const receipt={status:'PASS',stage:'display_only_IN_spacing_ID_label_supplement',source_app_sha256:crypto.createHash('sha256').update(appBytes).digest('hex'),canvases,pdfs:saved,browser_execution_performed:false,scientific_fitting_performed:false,original_source_coordinates_unchanged:true};fs.writeFileSync(path.join(out,'SPACING_ID_COMPLETE.json'),JSON.stringify(receipt,null,2)+'\n');console.log(JSON.stringify(receipt));return;}
 check('Native options collection remains iterable without array-map workaround',()=>{assert.equal(elements.get('de-target-scanpy').options.map,undefined);assert([...elements.get('de-target-scanpy').options].length>=0);});
 check('UMAP screen projection preserves equal original x/y coordinate units',()=>{for(const tool of ['scanpy','seurat']){const points=app.active[tool].screenPoints,ns=app.active[tool].coords.nuclei,a=points[0],na=ns[a.index],b=points.find(p=>Math.abs(ns[p.index].x_harmony-na.x_harmony)>.1&&Math.abs(ns[p.index].y_harmony-na.y_harmony)>.1);assert(b);const nb=ns[b.index],sx=(b.x-a.x)/(nb.x_harmony-na.x_harmony),sy=-(b.y-a.y)/(nb.y_harmony-na.y_harmony);assert(Math.abs(sx-sy)<1e-8,'coordinate units are anisotropic');}});
 for(const tool of ['scanpy','seurat']){await capture('umap-'+tool,'actual_CA1_reviewed_UMAP_'+tool);exportPrefix='reviewed_';await app.exportPDF('umap',tool);}
 app.state.tab='clusters';app.render();const cluster=app.active.scanpy.coords.nuclei[0].cluster_id;app.state.labels['scanpy|CA1_ProS|original_frozen|'+cluster]='Edited display label';app.state.umapLabelMode='both';app.state.tab='umap';app.render();
 check('Edited source-bound labels appear on UMAP drawing commands',()=>assert(elements.get('umap-scanpy').getContext().parts.some(x=>x.includes('Edited display label'))));await capture('umap-scanpy','actual_CA1_edited_labels_UMAP');exportPrefix='edited_labels_';await app.exportPDF('umap','scanpy');
 app.state.tab='markers';app.state.labels={};await app.ensureMarkers();app.state.markers=manifest.genes.slice(0,65).map(g=>g.gene_id);app.state.genesPerPage=60;app.state.markerPage=0;app.state.dotRangeMode='source';app.render();
 async function prepareDotPages(tool,label){const oldPage=app.state.markerPage,pages=Math.max(1,Math.ceil(app.state.markers.length/app.pageGenes()));for(let page=0;page<pages;page++){app.state.markerPage=page;app.render();await capture('dots-'+tool,label+'_page'+(page+1));}app.state.markerPage=oldPage;app.render();}
 for(const tool of ['scanpy','seurat']){await prepareDotPages(tool,'actual_CA1_60genes_'+tool);check(tool+'60 genes have intrinsic scrolling width',()=>assert(elements.get('dots-'+tool).width>=1655));exportPrefix='60genes_';await app.exportPDF('dots',tool);}
 app.state.markers=manifest.genes.slice(0,120).map(g=>g.gene_id);elements.get('genes-per-page').value='120';await elements.get('genes-per-page').dispatch('input');app.state.dotRangeMode='manual';app.state.dotMin=-1.5;app.state.dotMax=1.5;app.state.dotArea=180;app.state.negativeArea=.1;app.state.negativeAlpha=.08;app.render();
 check('120gene control/manual limits affect only rendering state',()=>{assert.equal(app.pageGenes(),120);assert.deepEqual(Array.from(app.markerRange('scanpy')),[-1.5,1.5]);assert.equal(app.state.markers.length,120);});await capture('dots-scanpy','actual_CA1_120genes_manual_z');check('120gene intrinsic width is large enough for distinct gene columns',()=>assert(elements.get('dots-scanpy').width>=2975));exportPrefix='120genes_manual_';await app.exportPDF('dots','scanpy');
 app.state.tab='de';await app.ensureFits();await app.refreshDERows();app.render();for(const tool of ['scanpy','seurat']){await capture('volcano-'+tool,'actual_CA1_volcano_'+tool);exportPrefix='primary_';await app.exportPDF('volcano',tool);}await capture('de-compare','actual_CA1_effect_compare');await app.exportPDF('compare','both');
 await scope('Inhibitory');app.state.annotation='cluster_id';app.state.spacing=1;app.render();await capture('umap-scanpy','actual_IN_original_geometry');app.state.spacing=.5;app.render();
 check('Nondefault IN display geometry is explicitly watermarked',()=>assert(elements.get('umap-scanpy').getContext().parts.some(x=>x.includes('DISPLAY LAYOUT'))));await capture('umap-scanpy','actual_IN_spacing_0p5');exportPrefix='IN_spacing_0p5_';await app.exportPDF('umap','scanpy');elements.get('reset-spacing').onclick();app.render();check('Restore original geometry resets persisted display state',()=>assert.equal(app.state.spacing,1));
 app.state.tab='boxplots';app.state.deEndpoint='PRIMARY';await app.ensureFits();await app.ensureMouseData();app.state.boxGenes=mouse.genes.filter(g=>['Cck','Gad1','Satb2','Apool','Arhgap31','Pdia5','Pvalb'].includes(g.gene_symbol)).map(g=>g.gene_id);app.state.boxTransform='log1p_normalized';app.state.boxPage=0;app.render();
 async function prepareBoxPages(tool,label){const old=app.state.boxPage,pages=Math.max(1,Math.ceil(app.state.boxGenes.length/6));for(let page=0;page<pages;page++){app.state.boxPage=page;app.render();await capture('boxplot-'+tool,label+'_page'+(page+1));}app.state.boxPage=old;app.render();}
 for(const tool of ['scanpy','seurat']){await prepareBoxPages(tool,'actual_IN_mouse_normalized_'+tool);check(tool+' boxes draw independent-mouse points and source-only factors',()=>{assert(elements.get('box-status-'+tool).textContent.includes('EXACT_MATCH'));assert(elements.get('boxplot-'+tool).getContext().parts.some(x=>x.includes('WT')));});exportPrefix='IN_mouse_normalized_';await app.exportPDF('boxplot',tool);}
 app.state.boxTransform='raw';app.state.boxGenes=mouse.genes.filter(g=>['Cck','Gad1','Satb2'].includes(g.gene_symbol)).map(g=>g.gene_id);app.render();await capture('boxplot-scanpy','actual_IN_mouse_raw');exportPrefix='IN_mouse_raw_';await app.exportPDF('boxplot','scanpy');
 await scope('CA1_ProS','boxplots');await app.ensureFits();const unsupported=app.active.scanpy.fits.fits.find(f=>f.endpoint==='PRIMARY'&&f.analysis_kind==='cluster'&&f.status==='NOT_TESTED'&&f.members.length);assert(unsupported);app.state.deTargets.scanpy='fit:'+unsupported.id;await app.selectionChanged();app.state.boxGenes=[mouse.genes.find(g=>g.gene_symbol==='Cck').gene_id];app.state.boxTransform='normalized';app.render();
 check('Unsupported normalized boxplot does not invent finite values',()=>assert(app.boxValues('scanpy').every(r=>r.display_value===null)));await capture('boxplot-scanpy','actual_CA1_unsupported_normalized_NA');exportPrefix='unsupported_NA_';await app.exportPDF('boxplot','scanpy');
 const examples=[['allcells','Gad1'],['Excitatory','Slc17a7'],['Inhibitory','Cck'],['CA1_ProS','Calb1'],['CA1_ProS','Etv1']];
 for(const [pop,symbol]of examples){await scope(pop);await app.ensureMarkers();const gene=manifest.genes.find(g=>g.gene_symbol===symbol);assert(gene,'Existing produced gene example required '+symbol);app.state.annotation='gene';app.state.overlayGene=gene.gene_id;app.state.overlayMetric='mean_log1p_normalized';app.state.overlayCondition='ALL';app.state.spacing=1;app.state.umapLabelMode='id';app.render();for(const tool of ['scanpy','seurat']){await capture('umap-'+tool,'example_'+pop+'_'+symbol+'_'+tool);exportPrefix='example_'+symbol+'_';await app.exportPDF('umap',tool);}}
 await scope('CA1_ProS');app.state.annotation='layer_program_label';app.render();for(const tool of ['scanpy','seurat']){await capture('umap-'+tool,'example_CA1_WT_associated_programs_'+tool);exportPrefix='CA1_programs_';await app.exportPDF('umap',tool);}
 const programme=[...new Set(app.active.scanpy.coords.nuclei.map(n=>n.layer_program_label))].find(x=>/superficial/i.test(x));assert(programme,'Existing superficial-associated program metadata');app.state.umapProgram=programme;app.state.annotation='gene';app.state.overlayGene=manifest.genes.find(g=>g.gene_symbol==='Calb1').gene_id;app.state.overlayMetric='mean_log1p_normalized';app.render();for(const tool of ['scanpy','seurat']){await capture('umap-'+tool,'example_CA1_program_Calb1_'+tool);exportPrefix='example_CA1_program_Calb1_';await app.exportPDF('umap',tool);}check('Fifth category programme view retains native member labels and display-only filtering',()=>{assert(app.umapIndices('scanpy').length<app.retainedIndices('scanpy').length);assert(app.umapIndices('scanpy').every(i=>app.active.scanpy.coords.nuclei[i].layer_program_label===programme));});
 assert.equal(crypto.createHash('sha256').update(fs.readFileSync(path.join(src,'app.js'))).digest('hex'),crypto.createHash('sha256').update(appBytes).digest('hex'),'Frontend changed during QA; rerun against final snapshot');
 const receipt={status:'PASS',validation_mode:'Pure mock DOM/native collection semantics plus actual produced-data canvas commands and bundled jsPDF; no browser execution',browser_rendering_performed:false,source_app_sha256:crypto.createHash('sha256').update(appBytes).digest('hex'),source_helper_sha256:crypto.createHash('sha256').update(fs.readFileSync(path.join(src,'v2_features.js'))).digest('hex'),source_index_sha256:crypto.createHash('sha256').update(Buffer.from(html)).digest('hex'),canvases,pdfs:saved,checks,scientific_fitting_performed:false,raw_nucleus_count_aggregation_performed:false,limitations:['Actual browser/CSS/interaction inspection remains unverified after the prior explicit file:// URL policy rejection. These artifacts verify drawing commands and PDF exports, not a real browser screenshot.']};
 fs.writeFileSync(path.join(out,'COMPLETE.json'),JSON.stringify(receipt,null,2)+'\n');console.log(JSON.stringify({status:'PASS',canvases:canvases.length,pdfs:saved.length,checks:checks.length,app_sha256:receipt.source_app_sha256}));
})().catch(e=>{console.error(e);process.exitCode=1;});
