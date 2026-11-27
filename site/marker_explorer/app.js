/* Shared, read-only Seurat marker editor. Guides and borders are display settings. */
(() => {
 'use strict';
 const D=window.SATB2_MARKER_EXPLORER_DATA,root=document.querySelector('[data-marker-explorer],[data-ca1-top-markers]');
 if(!D||!root)return;
 const site=new URL('../',document.currentScript.src),NS='http://www.w3.org/2000/svg';
 const esc=s=>String(s).replace(/[&<>"']/g,c=>({'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;',"'":'&#39;'}[c]));
 const $=s=>root.querySelector(s),all=s=>[...root.querySelectorAll(s)];
 const catalog=new Map(D.catalog.map(([id,symbol])=>[id,{id,symbol}])),names=new Map(D.catalog.flatMap(([id,symbol])=>[[id.toLowerCase(),id],[symbol.toLowerCase(),id]]));
 const preset=id=>D.presets.find(p=>p.id===id)||D.presets[0];
 const initial=()=>({condition:'ALL',mode:'dot',metric:'z',order:'peak',low:'#2166ac',mid:'#f7f7f7',high:'#b2182b',min:-.7,max:2,radius:8,font:15,opacity:1,cellWidth:34,groups:true,rowLabel:D.default_row_label,plotBorder:false,dotBorder:false,cellBorder:false,hLines:false,vLines:false,borderColor:'#203c50',dotBorderColor:'#203c50',guideColor:'#a5b5c0',borderWidth:1,dotBorderWidth:.7,guideWidth:.6,guideOpacity:.6,page:0,perPage:100,genes:[...preset(D.default_preset).genes],customGenes:[],clusters:[...D.default_clusters]});
 let state=initial(),current=null,ready=false,payloads=[],byCondition={},references=new Map(),peaks=new Map();
 root.classList.add('mxe');
 root.innerHTML=`<p class="mx-intro">All <strong>${D.screenshot_markers.length} screenshot markers</strong> are available, plus ${D.catalog.length.toLocaleString()} genes in the existing saved marker catalog. Add or remove genes and cluster rows below. ${D.id==='ca1_markers'?'C7 starts hidden and can be restored; C8 is CA1-prosubiculum.':'All 20 original inhibitory clusters start visible; subtype names are working labels.'}</p>
 <div class="mx-controls">
 <label>Plot type<select data-setting="mode"><option value="dot">Dotplot</option><option value="heat">Heatmap · pheatmap style</option></select></label>
 <label>Expression displayed<select data-setting="condition"><option value="ALL">Pooled · all 12 mice</option><option value="WT">WT · 6 mice</option><option value="cKO">cKO · 6 mice</option></select></label>
 <label>Colour measurement<select data-setting="metric"><option value="z">Shared-reference expression z</option><option value="mean">Mean log1p CP10k</option><option value="detected">Detected fraction</option></select></label>
 <label>Gene-column order<select data-setting="order"><option value="peak">Peak expression · top to bottom</option><option value="selection">Screenshot / manual order</option><option value="alpha">Gene names A–Z</option></select></label>
 <label>Colour minimum<input data-setting="min" type="number" step=".1" value="-.7"></label><label>Colour maximum<input data-setting="max" type="number" step=".1" value="2"></label>
 <label>Low / midpoint / high colours<span class="mx-colours"><input data-setting="low" type="color" value="#2166ac" aria-label="Low colour"><input data-setting="mid" type="color" value="#f7f7f7" aria-label="Midpoint colour"><input data-setting="high" type="color" value="#b2182b" aria-label="High colour"></span></label>
 <label>Maximum dot radius<input data-setting="radius" type="range" min="3" max="13" step=".5" value="8"></label><label>Label size<input data-setting="font" type="range" min="12" max="23" value="15"></label><label>Mark opacity<input data-setting="opacity" type="range" min=".2" max="1" step=".05" value="1"></label><label>Column spacing<input data-setting="cellWidth" type="range" min="28" max="70" value="34"></label>
 <label>Cluster row labels<select data-setting="rowLabel"><option value="id">Cluster IDs</option><option value="both">IDs + short labels</option><option value="label">Short labels</option></select></label>
 </div>
 <details class="mx-box" open><summary>Plot boundary, dot outlines and connector lines</summary><div class="mx-toggles">
 <label><input data-setting="plotBorder" type="checkbox">Boundary around plot</label><label><input data-setting="dotBorder" type="checkbox">Boundary around each dot</label><label><input data-setting="cellBorder" type="checkbox">Heatmap cell boundaries</label><label><input data-setting="hLines" type="checkbox">Horizontal connector lines</label><label><input data-setting="vLines" type="checkbox">Vertical connector lines</label><label><input data-setting="groups" type="checkbox" checked>Peak-expression blocks</label></div>
 <div class="mx-style-controls"><label>Plot / cell border colour<input data-setting="borderColor" type="color" value="#203c50"></label><label>Plot / cell border width<input data-setting="borderWidth" type="number" min=".2" max="4" step=".2" value="1"></label><label>Dot outline colour<input data-setting="dotBorderColor" type="color" value="#203c50"></label><label>Dot outline width<input data-setting="dotBorderWidth" type="number" min=".2" max="3" step=".1" value=".7"></label><label>Connector colour<input data-setting="guideColor" type="color" value="#a5b5c0"></label><label>Connector width<input data-setting="guideWidth" type="number" min=".2" max="3" step=".1" value=".6"></label><label>Connector opacity<input data-setting="guideOpacity" type="range" min=".1" max="1" step=".05" value=".6"></label></div>
 <p class="mx-note">Connector lines are visual guides through dot centres. In heatmap view they form row/column grid lines. They do not represent biological connections.</p></details>
 <div class="mx-toolbar">${D.presets.map(p=>`<button data-preset="${p.id}">${esc(p.label)}</button>`).join('')}<button data-action="auto-range">Auto colour limits</button><button data-action="reset">Reset plot</button></div>
 <p class="mx-status" aria-live="polite">Loading verified saved Seurat summaries…</p><div class="mx-scroll" aria-label="Scrollable Seurat marker plot"></div>
 <div class="mx-toolbar"><button data-action="prev">← Previous genes</button><span data-page-note></span><button data-action="next">Next genes →</button><label>Genes per page<select data-setting="perPage"><option>40</option><option>60</option><option selected>100</option><option>180</option></select></label></div>
 <div class="mx-toolbar"><button data-export="svg">SVG · this page</button><button data-export="pdf">PDF · all gene pages</button><button data-export="png">PNG · high resolution</button><button data-export="values">All selected values TSV</button><button data-export="genes">Gene list TSV</button><button data-export="state">Plot settings JSON</button><button data-action="load-state">Load settings</button><input data-state-file type="file" accept="application/json,.json" hidden></div><p class="mx-message" aria-live="polite"></p>
 <details class="mx-box" open><summary>Add or remove cluster rows</summary><div class="mx-toolbar"><button data-action="all-clusters">Add all clusters</button><button data-action="no-clusters">Clear clusters</button><button data-action="default-clusters">Restore default rows</button></div><div class="mx-clusters"></div></details>
 <details class="mx-box"><summary>Add or remove marker genes</summary><p>Uncheck any gene to remove it; check it again to restore it. Screenshot markers retain their original positions in manual order.</p><div class="mx-toolbar"><button data-preset="screenshot">Select all screenshot markers</button><button data-action="no-genes">Clear genes</button><label>Filter these marker choices<input data-gene-filter type="search" placeholder="Gene symbol"></label></div><div class="mx-genes"></div>
 <label>Find more genes in the existing saved catalog<input data-catalog-filter type="search" placeholder="Search any of the ${D.catalog.length} available genes"></label><select data-catalog-choice size="6" aria-label="Available saved marker genes"></select><button data-action="add-gene">Add selected gene</button>
 <label>Paste gene names or stable IDs<textarea data-gene-paste rows="3" placeholder="e.g. Satb2 Ano3 Calb1"></textarea></label><div class="mx-toolbar"><button data-action="add-paste">Add pasted genes</button><button data-action="use-paste">Use pasted genes only</button></div></details>
 <p class="mx-note">Dot area is the measured nonzero fraction. Heatmap tiles show the selected colour measurement using the same genes, row order and colour scale. Shared-reference z uses the saved pooled gene mean and SD; filtering rows never rescales it. Gene selection and WT/cKO views are descriptive displays, not new differential-expression tests. Original counts and classifications are preserved.</p>
 <p class="mx-note"><a data-provenance href="provenance.json" download>Source and screenshot marker mapping</a></p>`;
 $('[data-provenance]').href=new URL(D.provenance_path,site).href;
 const colorKeys=['low','mid','high','borderColor','dotBorderColor','guideColor'];
 const booleans=['groups','plotBorder','dotBorder','cellBorder','hLines','vLines'];
 const ranges={min:[-100,100],max:[-100,100],radius:[3,13],font:[12,23],opacity:[.2,1],cellWidth:[28,70],borderWidth:[.2,4],dotBorderWidth:[.2,3],guideWidth:[.2,3],guideOpacity:[.1,1],perPage:[40,180]};
 const strings={condition:['ALL','WT','cKO'],mode:['dot','heat'],metric:['z','mean','detected'],order:['peak','selection','alpha'],rowLabel:['id','both','label']};
 const hex=c=>[1,3,5].map(i=>parseInt(c.slice(i,i+2),16));
 function colour(v){const t=Math.max(0,Math.min(1,(v-state.min)/(state.max-state.min))),k=t<=.5?0:1,f=t<=.5?t*2:(t-.5)*2,p=[state.low,state.mid,state.high],a=hex(p[k]),b=hex(p[k+1]);return '#'+a.map((x,i)=>Math.round(x+(b[i]-x)*f).toString(16).padStart(2,'0')).join('');}
 const displayZ=r=>{const ref=references.get(r.gene_id);return ref&&ref.scale_reference_sd>0?(r.mean_log1p_normalized-ref.scale_reference_mean)/ref.scale_reference_sd:null;};
 const value=r=>!r?null:state.metric==='mean'?r.mean_log1p_normalized:state.metric==='detected'?r.fraction_detected:displayZ(r);
 const measureLabel=()=>state.metric==='mean'?'Mean log1p CP10k':state.metric==='detected'?'Detected fraction':'Shared-reference expression z';
 function markerOrder(){const a=state.genes.filter(g=>catalog.has(g));if(state.order==='alpha')a.sort((x,y)=>catalog.get(x).symbol.localeCompare(catalog.get(y).symbol));if(state.order==='peak')a.sort((x,y)=>Number(peaks.get(x)||999)-Number(peaks.get(y)||999)||state.genes.indexOf(x)-state.genes.indexOf(y));return a;}
 const chosenClusters=()=>D.clusters.filter(c=>state.clusters.includes(c.id));
 const rowLabel=c=>state.rowLabel==='id'?'C'+c.id:state.rowLabel==='label'?c.label:'C'+c.id+' · '+c.label;
 function choices(){
  $('.mx-clusters').innerHTML=D.clusters.map(c=>`<label><input type="checkbox" data-cluster="${c.id}" ${state.clusters.includes(c.id)?'checked':''}>C${c.id} · ${esc(c.label)} <span class="mx-note">(${c.n.toLocaleString()})</span></label>`).join('');
  const pool=[...new Set([...D.presets.flatMap(p=>p.genes),...state.customGenes])],filter=$('[data-gene-filter]').value.trim().toLowerCase();
  $('.mx-genes').innerHTML=pool.filter(id=>!filter||catalog.get(id).symbol.toLowerCase().includes(filter)).map(id=>`<label><input type="checkbox" data-gene="${id}" ${state.genes.includes(id)?'checked':''}>${esc(catalog.get(id).symbol)}</label>`).join('');
  const q=$('[data-catalog-filter]').value.trim().toLowerCase();
  $('[data-catalog-choice]').innerHTML=D.catalog.filter(([id,s])=>!q||s.toLowerCase().includes(q)||id.toLowerCase().includes(q)).slice(0,q?150:50).map(([id,s])=>`<option value="${id}">${esc(s)} · ${id}</option>`).join('');
 }
 function makeChart(page=state.page){
  const ordered=markerOrder(),pages=Math.max(1,Math.ceil(ordered.length/state.perPage)),pageIndex=Math.max(0,Math.min(pages-1,page)),genes=ordered.slice(pageIndex*state.perPage,(pageIndex+1)*state.perPage),clusters=chosenClusters(),lookup=byCondition[state.condition];
  const ctx=document.createElement('canvas').getContext('2d');ctx.font=`600 ${state.font}px Arial`;
  const left=Math.max(94,...clusters.map(c=>ctx.measureText(rowLabel(c)).width+46)),dx=Math.max(state.cellWidth,state.font*1.3),dy=Math.max(30,state.font*1.65,state.radius*2+6),top=108,gridX=left-dx/2,gridY=top-dy/2,W=Math.max(730,left+genes.length*dx+30),bottom=top+clusters.length*dy,H=bottom+Math.max(125,state.font*8)+90,ops=[];
  const text=(x,y,t,size=state.font,fill='#183448',weight=400,align='start',angle=0)=>ops.push({type:'text',x,y,t,size,fill,weight,align,angle});
  const rect=(x,y,w,h,fill,stroke=null,width=0,opacity=1,cell=null,title='')=>ops.push({type:'rect',x,y,w,h,fill,stroke,width,opacity,cell,title});
  const line=(x,y,x2,y2,stroke=state.guideColor,width=state.guideWidth,opacity=state.guideOpacity)=>ops.push({type:'line',x,y,x2,y2,stroke,width,opacity});
  const circle=(x,y,r,fill,opacity=1,stroke=null,width=0,cell=null,title='')=>ops.push({type:'circle',x,y,r,fill,opacity,stroke,width,cell,title});
  rect(0,0,W,H,'#ffffff');text(24,31,`Seurat · ${D.title} · ${state.mode==='dot'?'marker dotplot':'marker heatmap'}`,23,'#183448',700);
  const first=D.catalog[0][0],n=clusters.reduce((sum,c)=>sum+(lookup.get(c.id+'|'+first)?.n_nuclei||0),0);
  text(24,57,`${state.condition==='ALL'?'Pooled 12 mice':state.condition+' · 6 mice'} · ${n.toLocaleString()} nuclei · ${ordered.length} selected genes · ${clusters.length} clusters${pages>1?' · page '+(pageIndex+1)+'/'+pages:''}`,14,'#526c7c');
  if(state.groups&&state.order==='peak'){
   for(let j=0;j<genes.length;){const group=peaks.get(genes[j]);let k=j+1;while(k<genes.length&&peaks.get(genes[k])===group)k++;const x=left+(j-.5)*dx,w=(k-j)*dx;rect(x,77,w-1,12,'#e9f1f5');text(x+w/2,88,group?'C'+group:'—',12,'#526c7c',600,'middle');j=k;}
  }
  clusters.forEach((c,i)=>{const y=top+i*dy;if(state.mode==='dot'&&i%2===0)rect(12,y-dy/2,W-28,dy,'#f5f8fa');text(left-27,y+state.font*.34,rowLabel(c),state.font,'#183448',600,'end');rect(14,y-9,5,18,c.color);});
  if(state.mode==='dot'&&genes.length&&clusters.length){if(state.hLines)clusters.forEach((c,i)=>line(left,top+i*dy,left+(genes.length-1)*dx,top+i*dy));if(state.vLines)genes.forEach((g,j)=>line(left+j*dx,top,left+j*dx,top+(clusters.length-1)*dy));}
  clusters.forEach((c,i)=>genes.forEach((g,j)=>{
   const r=lookup.get(c.id+'|'+g),v=value(r),x=left+j*dx,y=top+i*dy,cell={cluster:c.id,gene:g},title=`C${c.id} · ${catalog.get(g).symbol}\n${r?`Detected ${(r.fraction_detected*100).toFixed(1)}% (${r.n_detected}/${r.n_nuclei})\nMean log1p CP10k ${r.mean_log1p_normalized.toFixed(3)}\n${measureLabel()}: ${Number.isFinite(v)?v.toFixed(3):'unavailable'}`:'Source value unavailable'}`;
   if(r&&r.n_nuclei>0&&Number.isFinite(v)){
    if(state.mode==='heat')rect(x-dx/2,y-dy/2,dx,dy,colour(v),state.cellBorder?state.borderColor:null,state.cellBorder?state.borderWidth:0,state.opacity,cell,title);
    else circle(x,y,state.radius*Math.sqrt(Math.max(0,Math.min(1,r.fraction_detected))),colour(v),state.opacity,state.dotBorder?state.dotBorderColor:null,state.dotBorder?state.dotBorderWidth:0,cell,title);
   }else{if(state.mode==='heat')rect(x-dx/2,y-dy/2,dx,dy,'#e2e7eb',null,0,1,cell,title);line(x-3,y-3,x+3,y+3,'#687b88',1,1);line(x-3,y+3,x+3,y-3,'#687b88',1,1);}
  }));
  if(state.mode==='heat'&&genes.length&&clusters.length){if(state.hLines)for(let i=0;i<=clusters.length;i++)line(gridX,gridY+i*dy,gridX+genes.length*dx,gridY+i*dy);if(state.vLines)for(let j=0;j<=genes.length;j++)line(gridX+j*dx,gridY,gridX+j*dx,gridY+clusters.length*dy);}
  if(state.plotBorder&&genes.length&&clusters.length)rect(gridX,gridY,genes.length*dx,clusters.length*dy,'none',state.borderColor,state.borderWidth);
  genes.forEach((g,j)=>text(left+j*dx+4,bottom+6,catalog.get(g).symbol,state.font,'#183448',400,'end',-55));
  const ly=H-67;if(state.mode==='dot'){text(24,ly-2,'Detected nuclei · dot area',13,'#526c7c');[.1,.5,1].forEach((f,i)=>{circle(234+i*80,ly-6,state.radius*Math.sqrt(f),'#69869a',state.opacity,state.dotBorder?state.dotBorderColor:null,state.dotBorder?state.dotBorderWidth:0);text(250+i*80,ly-2,(f*100)+'%',12,'#526c7c');});}
  const gx=W-230;for(let i=0;i<180;i++)rect(gx+i,ly-13,1.1,8,colour(state.min+i/179*(state.max-state.min)));text(gx,ly-22,measureLabel(),13,'#526c7c');text(gx,ly+13,String(state.min),12,'#526c7c');text(gx+180,ly+13,String(state.max),12,'#526c7c',400,'end');
  text(24,H-20,`Saved Seurat summaries · ${state.order==='peak'?'pooled peak-expression column order':'custom column order'} · borders and guide lines are display settings`,12,'#526c7c');
  return {W,H,ops,genes,allGenes:ordered,clusters,n,pages,page:pageIndex,dx,dy,lookup};
 }
 function svg(c){
  const cell=o=>o.cell?` data-cell="true" data-cluster="${o.cell.cluster}" data-gene="${o.cell.gene}"`:'';
  const stroke=o=>o.stroke?` stroke="${o.stroke}" stroke-width="${o.width}"`:'';
  return `<svg xmlns="${NS}" width="${c.W}" height="${c.H}" viewBox="0 0 ${c.W} ${c.H}" role="img" aria-label="Seurat ${esc(D.title)} ${state.mode==='dot'?'dotplot':'heatmap'}"><title>${esc(D.title)} marker expression</title><desc>Saved Seurat means and detected fractions. No clustering, classification or differential model is refitted.</desc><style>text{font-family:Arial,Helvetica,sans-serif}</style>${c.ops.map(o=>o.type==='text'?`<text x="${o.x}" y="${o.y}" font-size="${o.size}" fill="${o.fill}" font-weight="${o.weight}" text-anchor="${o.align}"${o.angle?` transform="rotate(${o.angle} ${o.x} ${o.y})"`:''}>${esc(o.t)}</text>`:o.type==='line'?`<line x1="${o.x}" y1="${o.y}" x2="${o.x2}" y2="${o.y2}" stroke="${o.stroke}" stroke-width="${o.width}" opacity="${o.opacity}"/>`:o.type==='circle'?`<circle cx="${o.x}" cy="${o.y}" r="${o.r}" fill="${o.fill}" opacity="${o.opacity}"${stroke(o)}${cell(o)}><title>${esc(o.title)}</title></circle>`:`<rect x="${o.x}" y="${o.y}" width="${o.w}" height="${o.h}" fill="${o.fill}" opacity="${o.opacity}"${stroke(o)}${cell(o)}>${o.title?'<title>'+esc(o.title)+'</title>':''}</rect>`).join('')}</svg>`;
 }
 function render(){
  const valid=ready&&state.min<state.max;all('[data-export]').forEach(b=>b.disabled=!valid);
  if(!ready)return;if(!valid){$('.mx-status').textContent='Colour minimum must be below the maximum. The last valid plot is retained.';return;}
  current=makeChart();state.page=current.page;$('.mx-scroll').innerHTML=svg(current);$('.mx-scroll').dataset.genes=String(current.allGenes.length);$('.mx-scroll').dataset.clusters=String(current.clusters.length);
  $('.mx-status').textContent=`${current.allGenes.length} selected genes · ${current.clusters.length} / ${D.clusters.length} clusters · ${state.mode==='dot'?'dotplot':'heatmap'} · ${state.order==='peak'?'peak expression from top to bottom':'selected gene order'}`;
  $('[data-page-note]').textContent=`Genes ${current.allGenes.length?state.page*state.perPage+1:0}–${Math.min((state.page+1)*state.perPage,current.allGenes.length)} of ${current.allGenes.length}`;
  $('[data-action="prev"]').disabled=state.page===0;$('[data-action="next"]').disabled=state.page>=current.pages-1;
 }
 function sync(){all('[data-setting]').forEach(el=>{const v=state[el.dataset.setting];if(el.type==='checkbox')el.checked=v;else el.value=v;});choices();render();}
 async function digest(bytes){return [...new Uint8Array(await crypto.subtle.digest('SHA-256',bytes))].map(x=>x.toString(16).padStart(2,'0')).join('');}
 async function source(desc){
  const url=desc.path.startsWith('https://')?desc.path:new URL(desc.path,site).href,response=await fetch(url);if(!response.ok)throw Error('Source HTTP '+response.status);
  const packed=await response.arrayBuffer();if(await digest(packed)!==desc.sha256)throw Error('Source checksum mismatch');
  let raw;if(typeof DecompressionStream==='function')raw=await new Response(new Blob([packed]).stream().pipeThrough(new DecompressionStream('gzip'))).arrayBuffer();else{if(!window.fflate)await new Promise((ok,bad)=>{const s=document.createElement('script');s.src=new URL('vendor/fflate.min.js',site).href;s.onload=ok;s.onerror=bad;document.head.append(s);});raw=window.fflate.gunzipSync(new Uint8Array(packed));}
  if(await digest(raw)!==desc.payload_sha256)throw Error('Decoded source checksum mismatch');return JSON.parse(new TextDecoder().decode(raw));
 }
 async function load(){
  payloads=await Promise.all(D.sources.map(source));
  for(const [condition,key] of [['ALL','rows'],['WT','WT_rows'],['cKO','cKO_rows']]){const rows=payloads.flatMap(p=>p[key]);byCondition[condition]=new Map();for(const r of rows){const k=r.cluster_id+'|'+r.gene_id;if(byCondition[condition].has(k))throw Error('Duplicate source cell');byCondition[condition].set(k,r);}}
  const pooled=payloads.flatMap(p=>p.rows),peakValues=new Map();for(const r of pooled){references.set(r.gene_id,r);const last=peakValues.get(r.gene_id);if(!last||r.mean_log1p_normalized>last.mean_log1p_normalized||(r.mean_log1p_normalized===last.mean_log1p_normalized&&Number(r.cluster_id)<Number(last.cluster_id)))peakValues.set(r.gene_id,r);}
  for(const [id,r] of peakValues)peaks.set(id,r.cluster_id);
  ready=true;sync();
 }
 const safeFile=()=>`Seurat_${D.id}_${state.mode}_${state.condition}`;
 function save(name,data,type){const url=URL.createObjectURL(new Blob([data],{type})),a=document.createElement('a');a.href=url;a.download=name;a.click();setTimeout(()=>URL.revokeObjectURL(url),2000);$('.mx-message').textContent='Downloaded '+name;}
 function tsv(rows){if(!rows.length)return '';const keys=[...new Set(rows.flatMap(r=>Object.keys(r)))],q=v=>v==null?'NA':String(v).replace(/[\t\r\n]/g,' ');return [keys.join('\t'),...rows.map(r=>keys.map(k=>q(r[k])).join('\t'))].join('\n')+'\n';}
 function values(){return chosenClusters().flatMap(c=>markerOrder().map(g=>{const r=byCondition[state.condition].get(c.id+'|'+g);return {condition:state.condition,display_cluster_label:c.label,...r,display_shared_reference_z:r?displayZ(r):null,display_measurement:state.metric,display_value:value(r),pooled_peak_cluster:peaks.get(g),plot_mode:state.mode};}));}
 function pdf(){
  if(!window.jspdf?.jsPDF)throw Error('PDF library is unavailable; please reload or use SVG.');let p;
  const ascii=s=>String(s).replace(/–|—/g,'-').replace(/†/g,'*').replace(/‡/g,'**');
  for(let page=0;page<current.pages;page++){
   const c=makeChart(page),format=[c.W,c.H],orientation=c.W>=c.H?'landscape':'portrait';if(!p)p=new window.jspdf.jsPDF({orientation,unit:'pt',format,compress:true});else p.addPage(format,orientation);
   for(const o of c.ops){if(o.opacity!=null&&o.opacity!==1)p.setGState(new p.GState({opacity:o.opacity}));
    if(o.type==='line'){p.setDrawColor(o.stroke);p.setLineWidth(o.width);p.line(o.x,o.y,o.x2,o.y2);}
    else if(o.type==='rect'||o.type==='circle'){if(o.fill!=='none')p.setFillColor(o.fill);if(o.stroke){p.setDrawColor(o.stroke);p.setLineWidth(o.width);}const style=o.fill==='none'?'S':o.stroke?'FD':'F';if(o.type==='rect')p.rect(o.x,o.y,o.w,o.h,style);else if(o.r>0)p.circle(o.x,o.y,o.r,style);}
    else{p.setFont('helvetica',o.weight>=600?'bold':'normal');p.setFontSize(o.size);p.setTextColor(o.fill);const t=ascii(o.t);if(o.angle){const width=p.getTextWidth(t),offset=o.align==='end'?width:o.align==='middle'?width/2:0,rad=o.angle*Math.PI/180;p.text(t,o.x-offset*Math.cos(rad),o.y-offset*Math.sin(rad),{align:'left',angle:-o.angle});}else p.text(t,o.x,o.y,{align:o.align==='middle'?'center':o.align==='end'?'right':'left'});}
    if(o.opacity!=null&&o.opacity!==1)p.setGState(new p.GState({opacity:1}));
   }
  }
  p.setProperties({title:'Seurat '+D.title+' marker '+state.mode,subject:'Saved expression summaries; guide lines and borders are display only'});p.save(safeFile()+'.pdf');$('.mx-message').textContent='Downloaded vector PDF · '+current.pages+' gene page(s).';
 }
 async function png(){const url=URL.createObjectURL(new Blob([svg(current)],{type:'image/svg+xml'}));try{const img=new Image();img.src=url;await img.decode();const scale=Math.min(3,16384/Math.max(current.W,current.H),Math.sqrt(60e6/(current.W*current.H))),c=document.createElement('canvas');c.width=Math.ceil(current.W*scale);c.height=Math.ceil(current.H*scale);const ctx=c.getContext('2d');ctx.scale(scale,scale);ctx.drawImage(img,0,0);save(safeFile()+'.png',await new Promise(ok=>c.toBlob(ok,'image/png')),'image/png');}finally{URL.revokeObjectURL(url);}}
 async function exportFile(kind){if(!ready||!(state.min<state.max))return;if(kind==='pdf')return pdf();if(kind==='png')return png();if(kind==='svg')return save(safeFile()+'.svg',svg(current),'image/svg+xml');if(kind==='values')return save(safeFile()+'_values.tsv',tsv(values()),'text/tab-separated-values');if(kind==='genes')return save(safeFile()+'_genes.tsv',tsv(markerOrder().map((id,i)=>({display_order:i+1,gene_id:id,gene_symbol:catalog.get(id).symbol,pooled_peak_cluster:peaks.get(id),screenshot_position:D.screenshot_markers.indexOf(catalog.get(id).symbol)+1||null}))), 'text/tab-separated-values');if(kind==='state')save(safeFile()+'_settings.json',JSON.stringify({schema:D.schema,dataset:D.id,display_only:true,...state},null,2),'application/json');}
 function paste(replace){const words=$('[data-gene-paste]').value.split(/[\s,;]+/).filter(Boolean),valid=[],missing=[];for(const word of words){const id=names.get(word.toLowerCase());if(id)valid.push(id);else missing.push(word);}if(valid.length){state.customGenes=[...new Set([...state.customGenes,...valid])];state.genes=[...new Set([...(replace?[]:state.genes),...valid])];state.page=0;sync();}$('.mx-message').textContent=(valid.length?valid.length+' names matched. ':'')+(missing.length?'Unavailable names: '+missing.join(', '):'');}
 async function action(a){
  if(a==='load-state')return $('[data-state-file]').click();
  if(a==='reset'){state=initial();return sync();}
  if(a==='all-clusters')state.clusters=D.clusters.map(c=>c.id);
  if(a==='no-clusters')state.clusters=[];
  if(a==='default-clusters')state.clusters=[...D.default_clusters];
  if(a==='no-genes'){state.genes=[];state.page=0;}
  if(a==='prev')state.page--;
  if(a==='next')state.page++;
  if(a==='add-gene'){const id=$('[data-catalog-choice]').value;if(id){state.customGenes=[...new Set([...state.customGenes,id])];state.genes=[...new Set([...state.genes,id])];state.page=0;}}
  if(a==='add-paste'||a==='use-paste')return paste(a==='use-paste');
  if(a==='auto-range'&&ready){const a=values().map(r=>r.display_value).filter(Number.isFinite);if(a.length){state.min=Math.min(...a);state.max=Math.max(...a);if(state.min===state.max)state.max=state.min+1;}}
  sync();
 }
 root.addEventListener('input',e=>{
  const k=e.target.dataset.setting;if(k){state[k]=e.target.type==='checkbox'?e.target.checked:strings[k]||colorKeys.includes(k)?e.target.value:Number(e.target.value);if(k==='metric'){[state.min,state.max]=state.metric==='z'?[-.7,2]:state.metric==='mean'?[0,3]:[0,1];sync();return;}if(k==='perPage')state.page=0;render();}
  if(e.target.hasAttribute('data-gene-filter')||e.target.hasAttribute('data-catalog-filter'))choices();
 });
 root.addEventListener('change',e=>{if(e.target.dataset.cluster){state.clusters=all('[data-cluster]:checked').map(x=>x.dataset.cluster);render();}if(e.target.dataset.gene){const id=e.target.dataset.gene;state.genes=e.target.checked?[...new Set([...state.genes,id])]:state.genes.filter(g=>g!==id);state.page=0;render();}});
 root.addEventListener('click',e=>{const p=e.target.closest('[data-preset]'),a=e.target.closest('[data-action]'),ex=e.target.closest('[data-export]');if(p){state.genes=[...preset(p.dataset.preset).genes];state.page=0;sync();}if(a)action(a.dataset.action).catch(err=>$('.mx-message').textContent=err.message);if(ex)exportFile(ex.dataset.export).catch(err=>$('.mx-message').textContent=err.message);});
 $('[data-state-file]').addEventListener('change',async e=>{try{const s=JSON.parse(await e.target.files[0].text());const legacy=D.id==='ca1_markers'&&s.schema==='satb2-ca1-top2-seurat-v2';if(!legacy&&(s.schema!==D.schema||s.dataset!==D.id))throw Error('Settings belong to another marker editor.');const n=initial();for(const k of colorKeys)if(/^#[\da-f]{6}$/i.test(s[k]))n[k]=s[k];for(const [k,allowed] of Object.entries(strings))if(allowed.includes(s[k]))n[k]=s[k];for(const [k,[lo,hi]] of Object.entries(ranges))if(Number.isFinite(s[k]))n[k]=Math.max(lo,Math.min(hi,s[k]));for(const k of booleans)if(typeof s[k]==='boolean')n[k]=s[k];for(const k of ['genes','customGenes'])if(Array.isArray(s[k]))n[k]=[...new Set(s[k].filter(g=>catalog.has(g)))];if(Array.isArray(s.clusters))n.clusters=D.clusters.map(c=>c.id).filter(id=>s.clusters.includes(id));if(!(n.min<n.max))throw Error('Invalid colour limits.');state=n;sync();}catch(err){$('.mx-message').textContent=err.message;}});
 const api={get state(){return JSON.parse(JSON.stringify(state));},get chart(){return current;},get ready(){return ready;},values,svg,render,makeChart};window.SATB2_MARKER_EXPLORER_VIEW=api;if(D.id==='ca1_markers')window.SATB2_CA1_TOP_MARKERS_VIEW=api;
 sync();load().catch(err=>{$('.mx-status').textContent='Could not load verified marker data: '+err.message;});
})();
