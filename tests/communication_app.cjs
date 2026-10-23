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
for(const file of ['model.js','charts.js','networks.js','data/deg_only_data.js'])vm.runInNewContext(fs.readFileSync(path.join(base,file),'utf8'),sandbox,{filename:file});
const packed=JSON.parse(zlib.gunzipSync(Buffer.from(window.CCC_DATA_GZIP,'base64'))),data=window.CCCModel.decode(packed);
for(const [key,d]of Object.entries(data.datasets)){const[scope,condition]=key.split('|');d.rows.forEach(r=>{r.scope=scope;r.condition=condition;});}
const app=fs.readFileSync(path.join(base,'app.js'),'utf8');
assert.match(app,/\nstart\(\);\n\}\)\(\);\s*$/,'Only the startup call is replaced for pure VM testing.');
const instrument=app.replace(/\nstart\(\);\n\}\)\(\);\s*$/,"\nwindow.__audit={DEFAULT,NUMBERS,valid,nodes,selected,getPairs,commonContext,bind,sync,change,render,exportPlot,scopeSVG,renderTable,renderDetail,populationOptions,tableDownload,setup(d,s){data=d;state=valid(s);tableCondition=state.condition;filteredCache={};pairPage=0;tablePage=0;},get state(){return state;},get displayed(){return displayed;},get tableCondition(){return tableCondition;}};\n})();");
vm.runInNewContext(instrument,sandbox,{filename:'app.js'});const A=window.__audit,M=window.CCCModel,P=window.CCCPlots,N=window.CCCNetworks;
const checks=[];function check(name,fn){fn();checks.push(name);console.log('PASS '+name);}
const original=JSON.stringify(data),baseState={...A.DEFAULT,sources:[],targets:[],nodeColors:{}};
const allowed=new Set(['CA1_Deep','CA1_Superficial','IN_C1','IN_C6','IN_C8']);
const totals={pooled:16,WT:15,cKO:16};
check('Every visible data/appearance control has a binding; all bound IDs exist',()=>{
 A.setup(data,baseState);A.populationOptions();A.bind();
 for(const key of Object.keys(A.DEFAULT).filter(k=>!['version','nodeColors','view','scope','bothDEG'].includes(k))){const el=elements.get(key);assert.ok(el,'Missing control '+key);assert.ok(Object.keys(el.listeners).length,'Unbound control '+key);}
 for(const m of app.matchAll(/\$\('([^']+)'\)/g)){if(['close-detail','download-component'].includes(m[1]))continue;assert.ok(elements.has(m[1]),'app references missing HTML ID '+m[1]);}
 assert.equal(tabs.length,8);assert.ok(tabs.every(b=>b.listeners.click));
});
check('Only the three saved strict fits are loaded; defaults include all native ranks and both DEG endpoints',()=>{
 assert.equal(Object.keys(data.datasets).length,3);assert.equal(data.nodes.length,5);assert(data.nodes.every(n=>allowed.has(n.id)));
 assert.equal(A.DEFAULT.version,2);assert.equal(A.DEFAULT.maxRank,1);assert.equal(A.DEFAULT.bothDEG,true);assert.equal(A.DEFAULT.scope,'IN_clusters');assert.equal(A.DEFAULT.tier,'DEG_only');
 for(const condition of ['pooled','WT','cKO']){const rows=data.datasets['IN_clusters|'+condition].rows;assert.equal(rows.length,totals[condition]);assert.equal(M.filter(rows,baseState).length,totals[condition]);const f=M.filter(rows,{...baseState,direction:'CA1_to_IN'}),r=M.filter(rows,{...baseState,direction:'IN_to_CA1'});assert.equal(f.length+r.length,totals[condition]);assert(rows.every(x=>allowed.has(x.source)&&allowed.has(x.target)&&x.ligand_all_subunits_deg&&x.receptor_all_subunits_deg&&x.both_sides_deg_supported));}
 const strong=M.filter(data.datasets['IN_clusters|pooled'].rows,{...baseState,tier:'replicated_specific'});assert.equal(strong.length,2);assert(strong.every(r=>r.source==='IN_C6'&&r.ligand_complex==='Cntn3'&&r.receptor_complex==='Ptprg'));
});
check('Population filters retain C8 empty context and accurate strict five-population heatmaps',()=>{
 A.setup(data,{...baseState,sources:['CA1_Deep'],view:'heatmap'});const ns=A.nodes();assert.equal(ns.length,4);assert.ok(ns.some(n=>n.id==='IN_C8'));assert.ok(!ns.some(n=>n.id==='CA1_Superficial'));assert.ok(ns.every(n=>allowed.has(n.id)));
 const rows=A.selected('pooled'),ctx=A.commonContext({pooled:rows});const cells=P.helpers.aggregateHeatmap(rows,{...ctx,eligiblePairs:A.getPairs('pooled')});assert.equal(cells.filter(c=>c.tested).length,3);assert.equal(cells.reduce((n,c)=>n+(c.count||0),0),rows.length);
 const c8=cells.find(c=>c.source==='CA1_Deep'&&c.target==='IN_C8');assert(c8&&c8.tested&&c8.count===0);
 A.setup(data,{...baseState,targets:['CA1_Deep']});assert.equal(A.nodes().length,4);assert(A.selected('pooled').every(r=>r.target==='CA1_Deep'));
 A.setup(data,{...baseState,sources:['IN_C1'],targets:['IN_C6']});assert.equal(A.getPairs('pooled').length,0);assert.equal(A.nodes().length,2);assert.equal(A.selected('pooled').length,0);
 A.setup(data,{...baseState,sources:['IN_C8']});assert.equal(A.nodes().length,3);assert.equal(A.getPairs('pooled').length,2);assert.equal(A.selected('pooled').length,0);
});
check('Endpoint DEG direction selectors enforce four signed combinations without changing native scores',()=>{
 const expectations={pooled:{'up→up':3,'up→down':2,'down→up':7,'down→down':4},WT:{'up→up':2,'up→down':2,'down→up':7,'down→down':4},cKO:{'up→up':3,'up→down':2,'down→up':7,'down→down':4}};
 for(const condition of ['pooled','WT','cKO'])for(const [combo,n]of Object.entries(expectations[condition])){const [l,r]=combo.split('→');A.setup(data,{...baseState,condition,directionPair:combo});const rows=A.selected(condition);assert.equal(rows.length,n,condition+' '+combo);assert(rows.every(x=>x.ligand_deg_direction===l&&x.receptor_deg_direction===r));A.setup(data,{...baseState,condition,ligandDirection:l,receptorDirection:r});assert.equal(A.selected(condition).length,n);}
 A.setup(data,{...baseState,directionPair:'down→up',ligandDirection:'up'});assert.equal(A.selected('pooled').length,0,'Conflicting filters combine with AND.');
 A.setup(data,{...baseState,directionPair:'mixed'});assert.equal(A.selected('pooled').length,0);
});
check('WT/cKO comparisons share LR identities, population columns, node order and numeric expression scales',()=>{
 A.setup(data,{...baseState,compare:true,view:'dotplot',dotMetric:'lr_gmeans'});const sets={WT:A.selected('WT'),cKO:A.selected('cKO')},ctx=A.commonContext(sets);
 const selections=Object.keys(sets).map(c=>P.helpers.selectDotplot(sets[c],{...ctx,condition:c}));assert.equal(JSON.stringify(selections[0].lrs.map(x=>x.key)),JSON.stringify(selections[1].lrs.map(x=>x.key)));assert.equal(JSON.stringify(selections[0].pairs.map(x=>x.key)),JSON.stringify(selections[1].pairs.map(x=>x.key)));assert.equal(ctx.nodes.length,5);
 const lrs=new Set(ctx.lrOrder),pairs=new Set(ctx.pairOrder),gmeans=Object.values(sets).flat().filter(r=>lrs.has(P.helpers.lrKey(r))&&pairs.has(P.helpers.pairKey(r))).map(r=>r.lr_gmeans);assert.equal(ctx.colorDomain[1],Math.max(...gmeans));
 A.setup(data,{...baseState,compare:true,view:'geneCircle'});const geneSets={WT:A.selected('WT'),cKO:A.selected('cKO')},gctx=A.commonContext(geneSets),selected=Object.values(geneSets).flat().filter(r=>gctx.geneKeys.includes(N.geneKey(r)));assert.equal(gctx.edgeMax,Math.max(...selected.map(r=>r.lr_gmeans)),'Gene circle widths use expression units, not candidate count units.');
});
check('Legacy policy settings cannot loosen scope, cohorts, all-subunit support or DEG-only evidence',()=>{
 A.setup(data,{...baseState,version:1,scope:'IN_total',tier:'all_candidates',bothDEG:false,sources:['IN_C2'],targets:['IN_Total'],fontSize:99,pairLimit:-10,nodeColors:{IN_C1:'#123456',IN_C20:'#aaaaaa',fake:'#ffffff'},lowColor:'javascript:x',compare:'yes'});assert.equal(A.state.version,2);assert.equal(A.state.scope,'IN_clusters');assert.equal(A.state.bothDEG,true);assert.equal(A.state.tier,'DEG_only');assert.equal(A.state.sources.length,0);assert.equal(A.state.targets.length,0);assert.equal(A.state.fontSize,20);assert.equal(A.state.pairLimit,1);assert.equal(A.state.compare,false);assert.equal(A.state.nodeColors.fake,undefined);assert.equal(A.state.nodeColors.IN_C20,undefined);assert.equal(A.selected('pooled').length,16);
 A.change('bothDEG',false);assert.equal(A.state.bothDEG,true);A.change('scope','IN_total');assert.equal(A.state.scope,'IN_clusters');A.change('tier','all_candidates');assert.equal(A.state.tier,'DEG_only');
 assert.match(app,/satb2-communication-deg-only-v2/);assert.doesNotMatch(app,/getItem\(['"]satb2-communication-v1/);assert.match(app,/satb2-communication-settings\/v2/);
 A.change('lowColor','#001122');assert.equal(A.state.palette,'custom');A.change('edgeColor','#112233');assert.equal(A.state.edgeColorMode,'custom');A.change('arcColor','#223344');assert.equal(A.state.arcColorMode,'custom');A.change('nodeColor','#334455');assert.equal(A.state.nodeColorMode,'custom');
 A.change('condition','WT');assert.equal(A.tableCondition,'WT');A.change('compare',true);assert.equal(A.tableCondition,'WT');
});
check('Reset, lead and context-specific shortcuts retain only the corrected strict cohort',()=>{
 A.setup(data,{...baseState,condition:'cKO',directionPair:'up→down',sources:['IN_C1'],compare:true});A.populationOptions();A.bind();
 elements.get('show-strong').onclick();assert.equal(A.state.tier,'replicated_specific');assert.equal(A.state.scope,'IN_clusters');assert.equal(A.state.bothDEG,true);assert.equal(A.state.directionPair,'any');assert.equal(A.selected('pooled').length,2);assert.equal(A.state.compare,false);
 elements.get('show-lead').onclick();assert.equal(A.state.tier,'DEG_only');assert.equal(A.state.view,'table');assert.equal(A.state.search,'Cntn3');assert.equal(A.selected('pooled').length,2);assert(A.selected('pooled').every(r=>r.source==='IN_C6'&&r.ligand_deg_direction==='down'&&r.receptor_deg_direction==='up'));
 elements.get('reset-settings').onclick();assert.equal(A.state.condition,'pooled');assert.equal(A.state.tier,'DEG_only');assert.equal(A.state.maxRank,1);assert.equal(A.state.search,'');assert.equal(A.state.sources.length,0);assert.equal(A.state.targets.length,0);assert.equal(A.selected('pooled').length,16);assert.equal(A.state.showDEGDirections,true);
 assert.equal(store.has('satb2-communication-v1'),false);const saved=JSON.parse(store.get('satb2-communication-deg-only-v2'));assert.equal(saved.version,2);assert.equal(saved.bothDEG,true);
});
check('All five charts render each strict fit with exact signed filters and immutable scientific evidence',()=>{
 for(const condition of ['pooled','WT','cKO'])for(const view of ['dotplot','heatmap','circle','circos','geneCircle']){
  A.setup(data,{...baseState,condition,view,topN:5,pairLimit:5});A.render();assert.equal(A.displayed.length,1,condition+' '+view);assert.equal(A.displayed[0].rows.length,totals[condition]);assert.doesNotMatch(A.displayed[0].svg,/NaN|Infinity|undefined/);
  A.setup(data,{...baseState,condition,view,directionPair:'down→up'});A.render();assert.equal(A.displayed.length,1);assert(A.displayed[0].rows.every(r=>r.ligand_deg_direction==='down'&&r.receptor_deg_direction==='up'));assert.doesNotMatch(A.displayed[0].svg,/NaN|Infinity|undefined/);
 }
 assert.equal(JSON.stringify(data),original,'Filters, styles or renders may not change saved scientific rows.');
});
check('Up/down annotations survive SVG export geometry and can be recolored or hidden independently of native numeric fills',()=>{
 for(const view of ['dotplot','heatmap','circle','circos','geneCircle']){
  A.setup(data,{...baseState,view,showDEGDirections:true});A.render();assert.match(A.displayed[0].svg,/data-deg-summary="true"/,view+' needs a signed RNA endpoint summary');
  A.setup(data,{...baseState,view,showDEGDirections:false});A.render();assert.doesNotMatch(A.displayed[0].svg,/data-deg-summary="true"/,view+' show/hide changes appearance only');
 }
 // Hiding the summary may reclaim header space and move the plot vertically.
 // Compare the scientific encodings independently of that layout offset.
 const geneCircles=svg=>[...svg.matchAll(/<circle\b[^>]*data-lr="[^>]*>/g)].map(m=>m[0].replace(/\s(?:cx|cy)="[^"]*"/g,''));
 A.setup(data,{...baseState,view:'dotplot',showDEGDirections:true,degUpColor:'#aa0011',degDownColor:'#0033cc'});A.render();const signed=A.displayed[0].svg;
 assert.match(signed,/ccc-deg-dot-rim/);assert.match(signed,/data-role="ligand"/);assert.match(signed,/data-role="receptor"/);assert(signed.includes('#aa0011')&&signed.includes('#0033cc'));
 A.setup(data,{...baseState,view:'dotplot',showDEGDirections:false,degUpColor:'#aa0011',degDownColor:'#0033cc'});A.render();const plain=A.displayed[0].svg;
 assert.doesNotMatch(plain,/ccc-deg-dot-rim/);assert.equal(JSON.stringify(geneCircles(signed)),JSON.stringify(geneCircles(plain)),'Signed rims must not alter native metric fill, detected-fraction area or plotted rows.');
 A.setup(data,{...baseState,view:'geneCircle',search:'Cntn3',showDEGDirections:true,degUpColor:'#aa0011',degDownColor:'#0033cc'});A.render();const lead=A.displayed[0].svg;
 assert.match(lead,/data-deg-direction="down"/);assert.match(lead,/data-deg-direction="up"/);assert.match(lead,/Cntn3 ↓/);assert.match(lead,/Ptprg ↑/);assert(lead.includes('#aa0011')&&lead.includes('#0033cc'));
 assert.equal(JSON.stringify(data),original,'Up/down display options cannot alter saved DESeq2 effects or LIANA scores.');
});
check('Invalid numeric limits clear stale exports and changing filters dismisses old component details',()=>{
 A.setup(data,{...baseState});A.render();assert.equal(A.displayed.length,1);elements.get('selected-detail').hidden=false;A.change('condition','WT');A.render();assert.equal(elements.get('selected-detail').hidden,true);
 A.setup(data,{...baseState,scaleMode:'manual',colorMin:2,colorMax:1});A.render();assert.equal(A.displayed.length,0);assert.match(elements.get('render-status').textContent,/maximum must be greater/);
});
check('SVG style scoping isolates adjacent plots and rewrites local URL IDs',()=>{
 const style={textContent:'text{fill:red}.title,.small{font-size:12px}'},child={id:'clip',attributes:[]},pathEl={attributes:[{name:'clip-path',value:'url(#clip)'}],setAttribute(k,v){this[k]=v;}},svg={setAttribute(k,v){this[k]=v;},querySelectorAll(q){return q==='style'?[style]:q==='[id]'?[child]:q==='*'?[pathEl]:[];}};A.scopeSVG(svg,'plot-x');assert.equal(svg.id,'plot-x');assert.equal(style.textContent,'#plot-x text{fill:red}#plot-x .title,#plot-x .small{font-size:12px}');assert.equal(child.id,'plot-x-clip');assert.equal(pathEl['clip-path'],'url(#plot-x-clip)');
});
(async()=>{
 A.setup(data,{...baseState,view:'circle',edgeColorMode:'custom',edgeColor:'#123456',arrowSize:18,edgeWidth:11,opacity:.4});A.populationOptions();A.bind();A.render();const plot=A.displayed[0];assert.ok(plot.svg.includes('#123456'));
 await A.exportPlot(0,'svg');assert.equal(await downloads.at(-1).blob.text(),plot.svg);assert.ok(downloads.at(-1).name.endsWith('.svg'));
 await A.exportPlot(0,'json');const spec=JSON.parse(await downloads.at(-1).blob.text());assert.equal(spec.settings.edgeWidth,11);assert.equal(spec.settings.arrowSize,18);assert.equal(spec.condition,'pooled');assert.equal(spec.records.length,16);assert(spec.records.every(r=>r.ligand_all_subunits_deg&&r.receptor_all_subunits_deg));
 await A.exportPlot(0,'png');assert.equal(downloads.at(-1).blob.type,'image/png');
 await A.exportPlot(0,'pdf');assert.equal(pdfs.length,1);assert.equal(pdfs[0].opts.format[0],1280);assert.equal(pdfs[0].opts.format[1],1080);assert.equal(pdfs[0].image[1],'PNG');assert.equal(pdfs[0].image[4],1280);assert.equal(pdfs[0].image[5],1080);assert.ok(pdfs[0].name.endsWith('.pdf'));
 checks.push('Current styled SVG, signed JSON, PNG and PDF exports complete with finite dimensions (DOM/canvas/PDF mocks)');console.log('PASS '+checks.at(-1));
 A.tableDownload(false,'pooled');let text=await downloads.at(-1).blob.text();assert.equal(text.trimEnd().split('\n').length,17);assert.equal(text.split('\n')[0],Array.from(M.SEVEN).join('\t'));
 A.tableDownload(true,'pooled');text=await downloads.at(-1).blob.text();assert.match(text.split('\n')[0],/ligand_deg_direction/);assert.match(text.split('\n')[0],/receptor_deg_direction/);assert.match(text.split('\n')[0],/ligand_deseq2_log2fc/);assert.match(text.split('\n')[0],/receptor_deseq2_log2fc/);
 checks.push('Seven-column exports retain the native format; extended exports include matching signed endpoint effects');console.log('PASS '+checks.at(-1));
 A.setup(data,{...baseState,view:'table',search:'Cntn3'});A.render();const lead=A.selected('pooled').find(r=>r.target==='CA1_Deep');A.renderDetail(lead);const detail=elements.get('selected-detail').innerHTML;
 assert.match(detail,/Ligand source[\s\S]*Down in cKO/);assert.match(detail,/Receptor target[\s\S]*Up in cKO/);assert.match(detail,/log2\(cKO \/ WT\)/);assert.match(detail,/group-versus-rest lr_logfc/);
 elements.get('download-component').onclick();const component=JSON.parse(await downloads.at(-1).blob.text());assert.equal(component.record.source,'IN_C6');assert.equal(component.record.target,'CA1_Deep');assert(component.ligand_components[0].stable_id_results[0].log2FoldChange<0);assert(component.receptor_components[0].stable_id_results[0].log2FoldChange>0);assert.equal(component.record.ligand_deseq2_log2fc,component.ligand_components[0].stable_id_results[0].log2FoldChange);assert.equal(component.record.receptor_deseq2_log2fc,component.receptor_components[0].stable_id_results[0].log2FoldChange);
 checks.push('Component details and exact JSON identify down-regulated C6 ligand and up-regulated CA1 receptor independently of native LR specificity');console.log('PASS '+checks.at(-1));
 A.setup(data,{...baseState,condition:'WT',direction:'CA1_to_IN',view:'circle',edgeColorMode:'custom',edgeColor:'#789abc',directionPair:'down→down'});A.render();const initial=A.displayed[0],oldPdfCount=pdfs.length;
 let releaseDecode;decodeBarrier={promise:new Promise(resolve=>{releaseDecode=resolve;})};const pending=A.exportPlot(0,'pdf');assert.equal(decodeEntered,1,'PDF export must yield during image decoding.');
 A.change('scope','IN_total');A.change('condition','cKO');A.change('direction','IN_to_CA1');A.change('directionPair','down→up');A.change('view','heatmap');A.render();assert.equal(A.state.scope,'IN_clusters');assert.equal(A.state.condition,'cKO');assert.equal(A.state.view,'heatmap');
 decodeBarrier=null;releaseDecode();await pending;
 assert.equal(pdfs.length,oldPdfCount+1);const resumed=pdfs.at(-1);assert(resumed.name.includes('WT')&&resumed.name.includes('CA1_to_IN')&&resumed.name.includes('circle'));assert.equal(resumed.properties.title,'SATB2 Directed population circle · WT');assert.equal(await drawnImages.at(-1).text(),initial.svg,'The rasterized SVG retains the initial signed filter and style.');
 checks.push('In-flight PDF retains its original dataset, signed filter, title and appearance during navigation');console.log('PASS '+checks.at(-1));
 A.setup(data,{...baseState,condition:'WT',direction:'CA1_to_IN',view:'circle',edgeColorMode:'custom',edgeColor:'#abcdef',nodeColors:{IN_C1:'#102030'},directionPair:'down→down'});A.render();const lastRendered=A.displayed[0];
 A.change('condition','cKO');A.change('direction','IN_to_CA1');A.change('directionPair','down→up');A.change('view','heatmap');A.change('edgeColor','#998877');
 assert.equal(A.displayed[0],lastRendered,'The debounced redraw has not yet replaced the visible old plot.');
 await A.exportPlot(0,'json');const beforeRedraw=JSON.parse(await downloads.at(-1).blob.text());assert.equal(beforeRedraw.settings.condition,'WT');assert.equal(beforeRedraw.settings.view,'circle');assert.equal(beforeRedraw.settings.edgeColor,'#abcdef');assert.equal(beforeRedraw.settings.nodeColors.IN_C1,'#102030');assert.equal(beforeRedraw.settings.directionPair,'down→down');assert.equal(beforeRedraw.condition,'WT');assert(beforeRedraw.records.every(r=>r.condition==='WT'&&r.ligand_deg_direction==='down'&&r.receptor_deg_direction==='down'));
 await A.exportPlot(0,'pdf');assert.equal(pdfs.at(-1).properties.title,'SATB2 Directed population circle · WT');assert.equal(await drawnImages.at(-1).text(),lastRendered.svg);
 checks.push('Export before a pending redraw uses the visible chart’s exact saved DEG subset and signed settings');console.log('PASS '+checks.at(-1));
 A.setup(data,{...baseState,condition:'cKO',view:'table',directionPair:'down→up'});A.populationOptions();A.bind();const before=JSON.stringify(A.state);elements.get('settings-file').files=[{text:async()=>JSON.stringify({schema:'satb2-communication-settings/v1',settings:{scope:'IN_total',tier:'all_candidates',bothDEG:false}})}];await elements.get('settings-file').onchange();assert.equal(JSON.stringify(A.state),before,'Legacy settings imports must be rejected, not silently restored.');assert(elements.get('render-status').textContent.length>0);
 const validSettings={...baseState,condition:'WT',directionPair:'up→down',bothDEG:false,scope:'IN_total',tier:'all_candidates'};elements.get('settings-file').files=[{text:async()=>JSON.stringify({schema:'satb2-communication-settings/v2',settings:validSettings})}];await elements.get('settings-file').onchange();assert.equal(A.state.scope,'IN_clusters');assert.equal(A.state.bothDEG,true);assert.equal(A.state.tier,'DEG_only');assert.equal(A.state.directionPair,'up→down');assert.equal(A.selected('WT').length,2);
 checks.push('Old v1 settings imports are rejected; v2 imports cannot bypass the strict DEG-only invariant');console.log('PASS '+checks.at(-1));
 assert.equal(JSON.stringify(data),original);assert.equal(errors.filter(e=>!/maximum must be greater/.test(e)).length,0,errors.join('\n'));
 console.log(JSON.stringify({status:'PASS',checks:checks.length,mode:'PURE_VM_DATA_DOM_MOCKS_NO_BROWSER',check_names:checks},null,2));
})().catch(e=>{console.error(e);process.exitCode=1;});
