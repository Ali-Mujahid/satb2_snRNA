'use strict';
/* Pure VM/data and DOM-mock acceptance checks. This does not launch a browser. */
const assert=require('node:assert/strict'),fs=require('node:fs'),path=require('node:path'),vm=require('node:vm'),zlib=require('node:zlib');
const base=path.join(__dirname,'../site/communication'),html=fs.readFileSync(path.join(base,'index.html'),'utf8'),elements=new Map(),downloads=[],pdfs=[],errors=[];
let nextUrl=0,decodeBarrier=null,decodeEntered=0;const blobs=new Map(),drawnImages=[];
class Element{
 constructor(id,tag='div',attrs=''){this.id=id;this.tagName=tag.toUpperCase();this.type=(attrs.match(/\btype="([^"]*)"/)||[])[1]||'';this.multiple=/\bmultiple\b/.test(attrs);this.value=(attrs.match(/\bvalue="([^"]*)"/)||[])[1]||'';this.dataset={};this.hidden=false;this.disabled=false;this.options=[];this.listeners={};this.textContent='';this._html='';this.attributes=[];this.classList={toggle(){},add(){},remove(){}};}
 set innerHTML(value){this._html=value;for(const m of value.matchAll(/<([\w-]+)\b([^>]*\bid="([^"]+)"[^>]*)>/g))if(!elements.has(m[3]))elements.set(m[3],new Element(m[3],m[1],m[2]));if(this.tagName==='SELECT'){this.options=[...value.matchAll(/<option\b([^>]*)>([\s\S]*?)<\/option>/g)].map(m=>({value:(m[1].match(/value="([^"]*)"/)||[])[1]||'',selected:/\bselected\b/.test(m[1]),textContent:m[2]}));if(!this.value&&this.options.length)this.value=this.options[0].value;}}
 get innerHTML(){return this._html;}
 get selectedOptions(){return this.options.filter(o=>o.selected);}
 addEventListener(type,fn){this.listeners[type]=fn;}
 setAttribute(k,v){this[k]=String(v);}
 getAttribute(k){return this[k]??null;}
 querySelectorAll(){return [];}
 scrollIntoView(){}
 append(){}
 remove(){}
 click(){if(this.tagName==='A')downloads.push({name:this.download,url:this.href,blob:blobs.get(this.href)});else if(this.onclick)this.onclick();}
}
for(const m of html.matchAll(/<([\w-]+)\b([^>]*\bid="([^"]+)"[^>]*)>/g))elements.set(m[3],new Element(m[3],m[1],m[2]));
for(const m of html.matchAll(/<select\b[^>]*id="([^"]+)"[^>]*>([\s\S]*?)<\/select>/g))elements.get(m[1]).innerHTML=m[2];
const tabs=[...html.matchAll(/<button\b[^>]*data-view="([^"]+)"[^>]*>/g)].map(m=>{const e=new Element('tab-'+m[1],'button');e.dataset.view=m[1];return e;});
const document={getElementById:id=>elements.get(id)||null,querySelectorAll:q=>q==='[data-view]'?tabs:[],querySelector:()=>null,body:new Element('body'),createElement:tag=>{
 const el=new Element('',tag);if(tag==='canvas'){el.getContext=()=>({fillRect(){},drawImage(img){drawnImages.push(blobs.get(img.src));},fillStyle:''});el.toDataURL=()=> 'data:image/png;base64,AA==';el.toBlob=fn=>fn(new Blob(['PNG'],{type:'image/png'}));}return el;
}};
class ImageMock{async decode(){assert.ok(blobs.has(this.src));if(decodeBarrier){decodeEntered++;await decodeBarrier.promise;}}}
class DOMParserMock{parseFromString(s){const attrs={};for(const m of s.slice(0,s.indexOf('>')).matchAll(/(width|height)="([^"]+)"/g))attrs[m[1]]=m[2];return {documentElement:{getAttribute:k=>attrs[k]}};}}
class PDFMock{constructor(opts){this.opts=opts;pdfs.push(this);}setProperties(p){this.properties=p;}addImage(...args){this.image=args;}save(name){this.name=name;}}
const window={jspdf:{jsPDF:PDFMock}},store=new Map(),sandbox={window,document,Blob,Image:ImageMock,DOMParser:DOMParserMock,URL:{createObjectURL:b=>{const u='blob:mock/'+ ++nextUrl;blobs.set(u,b);return u;},revokeObjectURL(){}},localStorage:{getItem:k=>store.get(k)||null,setItem:(k,v)=>store.set(k,v)},setTimeout:()=>1,clearTimeout(){},console:{error:e=>errors.push(e.message)},TextDecoder,Uint8Array};
for(const file of ['model.js','charts.js','networks.js','data/data.js'])vm.runInNewContext(fs.readFileSync(path.join(base,file),'utf8'),sandbox,{filename:file});
const packed=JSON.parse(zlib.gunzipSync(Buffer.from(window.CCC_DATA_GZIP,'base64'))),data=window.CCCModel.decode(packed);
for(const [key,d]of Object.entries(data.datasets)){const[scope,condition]=key.split('|');d.rows.forEach(r=>{r.scope=scope;r.condition=condition;});}
const app=fs.readFileSync(path.join(base,'app.js'),'utf8');
assert.match(app,/\nstart\(\);\n\}\)\(\);\s*$/,'Only the startup call is replaced for pure VM testing.');
const instrument=app.replace(/\nstart\(\);\n\}\)\(\);\s*$/,"\nwindow.__audit={DEFAULT,NUMBERS,valid,nodes,selected,getPairs,commonContext,bind,sync,change,render,exportPlot,scopeSVG,renderTable,renderDetail,populationOptions,setup(d,s){data=d;state=valid(s);tableCondition=state.condition;filteredCache={};pairPage=0;tablePage=0;},get state(){return state;},get displayed(){return displayed;},get tableCondition(){return tableCondition;}};\n})();");
vm.runInNewContext(instrument,sandbox,{filename:'app.js'});const A=window.__audit,M=window.CCCModel,P=window.CCCPlots,N=window.CCCNetworks;
const checks=[];function check(name,fn){fn();checks.push(name);console.log('PASS '+name);}
const original=JSON.stringify(data),baseState={...A.DEFAULT,sources:[],targets:[],nodeColors:{}};
check('Every visible data/appearance control has a binding; all bound IDs exist',()=>{
 A.setup(data,baseState);A.populationOptions();A.bind();
 for(const key of Object.keys(A.DEFAULT).filter(k=>!['version','nodeColors','view'].includes(k))){const el=elements.get(key);assert.ok(el,'Missing control '+key);assert.ok(Object.keys(el.listeners).length,'Unbound control '+key);}
 for(const m of app.matchAll(/\$\('([^']+)'\)/g)){if(['close-detail','download-component'].includes(m[1]))continue;assert.ok(elements.has(m[1]),'app references missing HTML ID '+m[1]);}
 assert.equal(tabs.length,8);assert.ok(tabs.every(b=>b.listeners.click));
});
check('All six fits reproduce saved row totals and both-direction unions',()=>{
 const expected={'IN_clusters|pooled':20471,'IN_clusters|WT':19467,'IN_clusters|cKO':21680,'IN_total|pooled':1108,'IN_total|WT':1040,'IN_total|cKO':1170};
 for(const[key,d]of Object.entries(data.datasets)){const s={...baseState,maxRank:1,tier:'all_candidates'};assert.equal(M.filter(d.rows,s).length,expected[key]);const forward=M.filter(d.rows,{...s,direction:'CA1_to_IN'}),reverse=M.filter(d.rows,{...s,direction:'IN_to_CA1'});assert.equal(forward.length+reverse.length,expected[key]);assert.ok(forward.every(r=>r.source.startsWith('CA1_')&&r.target.startsWith('IN_')));assert.ok(reverse.every(r=>r.source.startsWith('IN_')&&r.target.startsWith('CA1_')));}
 const strong=M.filter(data.datasets['IN_clusters|pooled'].rows,{...baseState,tier:'replicated_specific'});assert.equal(strong.length,324);assert.equal(strong.filter(r=>r.source.startsWith('CA1_')).length,112);assert.equal(strong.filter(r=>r.source.startsWith('IN_')).length,212);assert.equal(strong.filter(r=>r.both_sides_deg_supported).length,2);assert.equal(M.filter(data.datasets['IN_total|pooled'].rows,{...baseState,tier:'replicated_specific'}).length,0);
});
check('One-sided population filters retain every eligible opposite endpoint and accurate heatmaps',()=>{
 A.setup(data,{...baseState,sources:['CA1_Deep'],view:'heatmap'});const ns=A.nodes();assert.equal(ns.length,21);assert.ok(ns.some(n=>n.id==='IN_C20'));assert.ok(!ns.some(n=>n.id==='CA1_Superficial'));
 const rows=A.selected('pooled'),ctx=A.commonContext({pooled:rows});const cells=P.helpers.aggregateHeatmap(rows,{...ctx,eligiblePairs:A.getPairs('pooled')});assert.equal(cells.filter(c=>c.tested).length,20);assert.equal(cells.reduce((n,c)=>n+(c.count||0),0),rows.length);
 A.setup(data,{...baseState,targets:['CA1_Deep']});assert.equal(A.nodes().length,21);assert.ok(A.selected('pooled').every(r=>r.target==='CA1_Deep'));
 A.setup(data,{...baseState,sources:['IN_C1'],targets:['IN_C2']});assert.equal(A.getPairs('pooled').length,0);assert.equal(A.nodes().length,2);
});
check('WT/cKO comparisons share LR identities, population columns, node order and numeric expression scales',()=>{
 A.setup(data,{...baseState,compare:true,view:'dotplot',dotMetric:'lr_gmeans'});const sets={WT:A.selected('WT'),cKO:A.selected('cKO')},ctx=A.commonContext(sets);
 const selections=Object.keys(sets).map(c=>P.helpers.selectDotplot(sets[c],{...ctx,condition:c}));assert.equal(JSON.stringify(selections[0].lrs.map(x=>x.key)),JSON.stringify(selections[1].lrs.map(x=>x.key)));assert.equal(JSON.stringify(selections[0].pairs.map(x=>x.key)),JSON.stringify(selections[1].pairs.map(x=>x.key)));
 const lrs=new Set(ctx.lrOrder),pairs=new Set(ctx.pairOrder),gmeans=Object.values(sets).flat().filter(r=>lrs.has(P.helpers.lrKey(r))&&pairs.has(P.helpers.pairKey(r))).map(r=>r.lr_gmeans);assert.equal(ctx.colorDomain[1],Math.max(...gmeans));
 A.setup(data,{...baseState,compare:true,view:'geneCircle'});const geneSets={WT:A.selected('WT'),cKO:A.selected('cKO')},gctx=A.commonContext(geneSets),selected=Object.values(geneSets).flat().filter(r=>gctx.geneKeys.includes(N.geneKey(r)));assert.equal(gctx.edgeMax,Math.max(...selected.map(r=>r.lr_gmeans)),'Gene circle widths use expression units, not candidate count units.');
});
check('Imported settings are validated; custom color controls select their active modes',()=>{
 A.setup(data,{...baseState,scope:'bad',fontSize:99,pairLimit:-10,nodeColors:{IN_C1:'#123456',fake:'#ffffff'},lowColor:'javascript:x',compare:'yes'});assert.equal(A.state.scope,'IN_clusters');assert.equal(A.state.fontSize,20);assert.equal(A.state.pairLimit,1);assert.equal(A.state.compare,false);assert.equal(A.state.nodeColors.fake,undefined);
 A.change('lowColor','#001122');assert.equal(A.state.palette,'custom');A.change('edgeColor','#112233');assert.equal(A.state.edgeColorMode,'custom');A.change('arcColor','#223344');assert.equal(A.state.arcColorMode,'custom');A.change('nodeColor','#334455');assert.equal(A.state.nodeColorMode,'custom');
 A.change('condition','WT');assert.equal(A.tableCondition,'WT');A.change('compare',true);assert.equal(A.tableCondition,'WT');
});
check('All charts render across all six fits with exact current filters and immutable native evidence',()=>{
 for(const scope of ['IN_clusters','IN_total'])for(const condition of ['pooled','WT','cKO'])for(const view of ['dotplot','heatmap','circle','circos','geneCircle']){
  A.setup(data,{...baseState,scope,condition,view,topN:5,pairLimit:5});A.render();assert.equal(A.displayed.length,1,scope+' '+condition+' '+view);assert.equal(A.displayed[0].rows.length,M.filter(data.datasets[scope+'|'+condition].rows,baseState).length);assert.doesNotMatch(A.displayed[0].svg,/NaN|Infinity|undefined/);
 }
 assert.equal(JSON.stringify(data),original,'No filters, styles or renders may change scientific rows.');
});
check('Invalid numeric limits clear prior chart/export state and changing context dismisses stale details',()=>{
 A.setup(data,{...baseState});A.render();assert.equal(A.displayed.length,1);elements.get('selected-detail').hidden=false;A.change('scope','IN_total');A.render();assert.equal(elements.get('selected-detail').hidden,true);
 A.setup(data,{...baseState,scaleMode:'manual',colorMin:2,colorMax:1});A.render();assert.equal(A.displayed.length,0);assert.match(elements.get('render-status').textContent,/maximum must be greater/);
});
check('SVG style scoping isolates text/classes in adjacent plots and rewrites local URL IDs',()=>{
 const style={textContent:'text{fill:red}.title,.small{font-size:12px}'},child={id:'clip',attributes:[]},pathEl={attributes:[{name:'clip-path',value:'url(#clip)'}],setAttribute(k,v){this[k]=v;}},svg={setAttribute(k,v){this[k]=v;},querySelectorAll(q){return q==='style'?[style]:q==='[id]'?[child]:q==='*'?[pathEl]:[];}};A.scopeSVG(svg,'plot-x');assert.equal(svg.id,'plot-x');assert.equal(style.textContent,'#plot-x text{fill:red}#plot-x .title,#plot-x .small{font-size:12px}');assert.equal(child.id,'plot-x-clip');assert.equal(pathEl['clip-path'],'url(#plot-x-clip)');
});
(async()=>{
 A.setup(data,{...baseState,view:'circle',edgeColorMode:'custom',edgeColor:'#123456',arrowSize:18,edgeWidth:11,opacity:.4});A.render();const plot=A.displayed[0];assert.ok(plot.svg.includes('#123456'));
 await A.exportPlot(0,'svg');assert.equal(await downloads.at(-1).blob.text(),plot.svg);assert.ok(downloads.at(-1).name.endsWith('.svg'));
 await A.exportPlot(0,'json');const spec=JSON.parse(await downloads.at(-1).blob.text());assert.equal(spec.settings.edgeWidth,11);assert.equal(spec.settings.arrowSize,18);assert.equal(spec.condition,'pooled');assert.equal(spec.records.length,1589);
 await A.exportPlot(0,'png');assert.equal(downloads.at(-1).blob.type,'image/png');
 await A.exportPlot(0,'pdf');assert.equal(pdfs.length,1);assert.equal(pdfs[0].opts.format[0],1280);assert.equal(pdfs[0].opts.format[1],1080);assert.equal(pdfs[0].image[1],'PNG');assert.equal(pdfs[0].image[4],1280);assert.equal(pdfs[0].image[5],1080);assert.ok(pdfs[0].name.endsWith('.pdf'));
 checks.push('Current styled SVG, JSON, PNG and PDF export code paths complete with finite dimensions (DOM/canvas/PDF mocks)');console.log('PASS '+checks.at(-1));
 A.setup(data,{...baseState,scope:'IN_clusters',condition:'WT',direction:'CA1_to_IN',view:'circle',edgeColorMode:'custom',edgeColor:'#789abc'});A.render();const initial=A.displayed[0],oldPdfCount=pdfs.length;
 let releaseDecode;decodeBarrier={promise:new Promise(resolve=>{releaseDecode=resolve;})};const pending=A.exportPlot(0,'pdf');assert.equal(decodeEntered,1,'PDF export must yield during image decoding.');
 A.change('scope','IN_total');A.change('condition','cKO');A.change('direction','IN_to_CA1');A.change('view','heatmap');A.render();assert.equal(A.state.scope,'IN_total');assert.equal(A.state.condition,'cKO');assert.equal(A.state.view,'heatmap');assert.equal(A.displayed[0].context.scope,'IN_total');
 decodeBarrier=null;releaseDecode();await pending;
 assert.equal(pdfs.length,oldPdfCount+1);const resumed=pdfs.at(-1);assert.equal(resumed.name,'SATB2_LR_IN_clusters_WT_CA1_to_IN_circle.pdf','An export begun before navigation retains its original scope, condition, direction and view.');assert.equal(resumed.properties.title,'SATB2 Directed population circle · WT','PDF metadata must describe the original chart, not the current tab.');assert.equal(resumed.opts.format[0],1280);assert.equal(resumed.opts.format[1],1080);assert.equal(await drawnImages.at(-1).text(),initial.svg,'The rasterized SVG must retain the exact initial plot data and style.');
 checks.push('An in-flight PDF export retains its original plot, title and filename after changing scope, condition, direction and view during image decoding');console.log('PASS '+checks.at(-1));
 A.setup(data,{...baseState,scope:'IN_clusters',condition:'WT',direction:'CA1_to_IN',view:'circle',edgeColorMode:'custom',edgeColor:'#abcdef',nodeColors:{IN_C1:'#102030'}});A.render();const lastRendered=A.displayed[0];
 A.change('scope','IN_total');A.change('condition','cKO');A.change('direction','IN_to_CA1');A.change('view','heatmap');A.change('edgeColor','#998877');
 assert.equal(A.displayed[0],lastRendered,'The debounced redraw has not yet replaced the visible old plot.');assert.equal(A.state.scope,'IN_total');assert.equal(lastRendered.settings.scope,'IN_clusters');
 await A.exportPlot(0,'json');const beforeRedraw=JSON.parse(await downloads.at(-1).blob.text());assert.equal(beforeRedraw.settings.scope,'IN_clusters');assert.equal(beforeRedraw.settings.condition,'WT');assert.equal(beforeRedraw.settings.view,'circle');assert.equal(beforeRedraw.settings.edgeColor,'#abcdef');assert.equal(beforeRedraw.settings.nodeColors.IN_C1,'#102030');assert.equal(beforeRedraw.scope,'IN_clusters');assert.equal(beforeRedraw.condition,'WT');assert.ok(beforeRedraw.records.every(r=>r.scope==='IN_clusters'&&r.condition==='WT'));
 await A.exportPlot(0,'pdf');assert.equal(pdfs.at(-1).name,'SATB2_LR_IN_clusters_WT_CA1_to_IN_circle.pdf');assert.equal(pdfs.at(-1).properties.title,'SATB2 Directed population circle · WT');assert.equal(await drawnImages.at(-1).text(),lastRendered.svg);
 checks.push('Export during a pending debounced redraw uses the visible plot’s deep settings snapshot, exact native records and original PDF provenance');console.log('PASS '+checks.at(-1));
 assert.equal(JSON.stringify(data),original);assert.equal(errors.filter(e=>!/maximum must be greater/.test(e)).length,0,errors.join('\n'));
 console.log(JSON.stringify({status:'PASS',checks:checks.length,mode:'PURE_VM_DATA_DOM_MOCKS_NO_BROWSER',check_names:checks},null,2));
})().catch(e=>{console.error(e);process.exitCode=1;});
