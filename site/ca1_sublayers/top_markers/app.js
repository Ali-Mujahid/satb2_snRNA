/* Read-only views of saved Seurat summaries. No fits, classifications or source values change. */
(() => {
  'use strict';
  const D=window.SATB2_CA1_TOP_MARKERS,root=document.querySelector('[data-ca1-top-markers]');
  if(!D||!root)return;
  const base=new URL('.',document.currentScript.src),NS='http://www.w3.org/2000/svg';
  const esc=s=>String(s).replace(/[&<>"']/g,c=>({'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;',"'":'&#39;'}[c]));
  const initial=()=>({condition:'ALL',low:'#2166ac',mid:'#f7f7f7',high:'#b2182b',min:-0.7,max:2,radius:8,font:15,opacity:1,supportedOnly:false,groups:true,genes:D.markers.map(m=>m.gene_id),clusters:D.clusters.map(c=>c.id)});
  let state=initial(),current;
  root.classList.add('ctm');
  root.innerHTML=`<p>Two candidates for each original Seurat CA1/ProS cluster, using the saved marker ranking. Rows are all 11 clusters; each two-gene block is labelled by its selection cluster. <strong>C8 is CA1-prosubiculum.</strong></p>
  <p class="ctm-warning"><strong>C7:</strong> no gene passed the original representative-marker criteria. Its two displayed candidates (†) are exploratory, with weak specificity; they are not established C7 or CA1 markers. C11 remains an identity-conflicting group.</p>
  <div class="ctm-controls">
    <label>Expression displayed<select data-setting="condition"><option value="ALL">Pooled · all 12 mice</option><option value="WT">WT · 6 mice</option><option value="cKO">cKO · 6 mice</option></select></label>
    <label>Colour minimum<input data-setting="min" type="number" step="0.1" value="-0.7"></label>
    <label>Colour maximum<input data-setting="max" type="number" step="0.1" value="2"></label>
    <label>Low / midpoint / high colours<span class="ctm-palette"><input data-setting="low" type="color" value="#2166ac" aria-label="Low expression colour"><input data-setting="mid" type="color" value="#f7f7f7" aria-label="Midpoint expression colour"><input data-setting="high" type="color" value="#b2182b" aria-label="High expression colour"></span></label>
    <label>Maximum dot radius<input data-setting="radius" type="range" min="3" max="13" step="0.5" value="8"></label>
    <label>Gene / cluster label size<input data-setting="font" type="range" min="12" max="23" value="15"></label>
    <label>Dot opacity<input data-setting="opacity" type="range" min="0.2" max="1" step="0.05" value="1"></label>
  </div>
  <div class="ctm-toolbar"><label class="ctm-check"><input data-setting="supportedOnly" type="checkbox">Supported representatives only (hide C7 candidates)</label><label class="ctm-check"><input data-setting="groups" type="checkbox" checked>Show selection-cluster blocks</label><button data-action="reset">Reset plot</button></div>
  <p class="ctm-status" aria-live="polite"></p><div class="ctm-scroll" aria-label="Scrollable Seurat CA1 top-two marker dotplot"></div>
  <div class="ctm-toolbar"><button data-action="svg">SVG</button><button data-action="pdf">PDF · vector</button><button data-action="png">PNG · 3×</button><button data-action="values">Displayed data TSV</button><button data-action="genes">Gene list TSV</button><button data-action="state">Plot settings JSON</button><button data-action="load-state">Load settings</button><input data-state-file type="file" accept="application/json,.json" hidden></div><p class="ctm-action-message" aria-live="polite"></p>
  <p class="ctm-note">Dot <strong>area</strong> is the fraction of nuclei with nonzero source counts. Colour is the saved gene-wise expression z-score, clipped to the chosen limits. Low-expression dots keep their measured detection area; no negative-z shrinking is applied. Gene selection is fixed from the pooled saved marker panel, not genome-wide discovery or a WT–cKO DEG test. Hiding rows does not recalculate z-scores.</p>
  <details class="ctm-selection"><summary>Choose genes to keep or remove</summary><div class="ctm-toolbar"><button data-action="all-genes">All 22 candidates</button><button data-action="no-genes">Clear genes</button></div><div class="ctm-selections" data-gene-choices></div></details>
  <details class="ctm-selection"><summary>Choose cluster rows</summary><div class="ctm-selections" data-cluster-choices></div></details>
  <details class="ctm-selection"><summary>Marker ranking and support evidence</summary><div class="ctm-tablewrap"><table><thead><tr><th>Selection cluster</th><th>Rank</th><th>Gene</th><th>Native support</th><th>Detection difference vs rest</th><th>Mean log1p difference</th><th>Mice supporting expression</th></tr></thead><tbody>${D.markers.map(m=>`<tr><td>C${m.source_cluster}</td><td>${m.source_rank}</td><td>${esc(m.gene_symbol)}</td><td>${m.supported?'Supported candidate':'Exploratory †'}</td><td>${(m.detection_difference_vs_rest*100).toFixed(1)} percentage points</td><td>${m.mean_log1p_difference_vs_rest.toFixed(3)}</td><td>${m.mice_expression_supported} / ${m.mice_eligible}</td></tr>`).join('')}</tbody></table></div></details>
  <p class="ctm-note ctm-downloads"><a data-source="marker_selection.tsv" download>Full marker selection TSV</a><a data-source="all_expression_values.tsv" download>All conditions: expression TSV</a><a data-source="provenance.json" download>Source and method JSON</a><a data-source="figures/Seurat_CA1_top2_ALL.svg" download>Default SVG</a><a data-source="figures/Seurat_CA1_top2_ALL.pdf" download>Default PDF</a></p>`;
  const $=s=>root.querySelector(s),all=s=>[...root.querySelectorAll(s)];
  all('[data-source]').forEach(a=>a.href=new URL(a.dataset.source,base).href);
  const geneChoices=$('[data-gene-choices]'),clusterChoices=$('[data-cluster-choices]');
  function choices(){
    geneChoices.innerHTML=D.markers.map(m=>`<label class="ctm-check"><input type="checkbox" data-gene="${m.gene_id}" ${state.genes.includes(m.gene_id)?'checked':''}>C${m.source_cluster} · ${esc(m.gene_symbol)}${m.supported?'':' †'}</label>`).join('');
    clusterChoices.innerHTML=D.clusters.map(c=>`<label class="ctm-check"><input type="checkbox" data-cluster="${c.id}" ${state.clusters.includes(c.id)?'checked':''}>C${c.id} · ${esc(c.program)}</label>`).join('');
  }
  const hex=c=>[1,3,5].map(i=>parseInt(c.slice(i,i+2),16));
  function colour(v){const t=Math.max(0,Math.min(1,(v-state.min)/(state.max-state.min))),k=t<=.5?0:1,f=t<=.5?t*2:(t-.5)*2,p=[state.low,state.mid,state.high],a=hex(p[k]),b=hex(p[k+1]);return '#'+a.map((x,i)=>Math.round(x+(b[i]-x)*f).toString(16).padStart(2,'0')).join('');}
  function chart(){
    const markers=D.markers.filter(m=>state.genes.includes(m.gene_id)&&(!state.supportedOnly||m.supported)),clusters=D.clusters.filter(c=>state.clusters.includes(c.id));
    const lookup=new Map(D.rows[state.condition].map(r=>[r.cluster_id+'|'+r.gene_id,r]));
    const step=Math.max(39,state.font*2.3),dy=Math.max(30,state.font*1.65,state.radius*2+6),left=Math.max(92,45+state.font*3),top=108,font=state.font,W=Math.max(710,left+markers.length*step+32),bottom=top+clusters.length*dy,H=bottom+Math.max(122,font*8)+86,ops=[];
    const text=(x,y,t,size=font,fill='#183448',weight=400,align='start',angle=0)=>ops.push({type:'text',x,y,t,size,fill,weight,align,angle});
    const rect=(x,y,w,h,fill)=>ops.push({type:'rect',x,y,w,h,fill});
    const line=(x,y,x2,y2,stroke='#dce5ea',width=1)=>ops.push({type:'line',x,y,x2,y2,stroke,width});
    const circle=(x,y,r,fill,opacity=1,title='')=>ops.push({type:'circle',x,y,r,fill,opacity,title});
    rect(0,0,W,H,'#ffffff');text(24,31,'Seurat · CA1 / ProS · top-two marker candidates',23,'#183448',700);
    const n=clusters.reduce((s,c)=>s+(lookup.get(c.id+'|'+D.markers[0].gene_id)?.n_nuclei||0),0);
    text(24,57,`${state.condition==='ALL'?'Pooled 12 mice':state.condition+' · 6 mice'} · ${n.toLocaleString()} nuclei · ${markers.length} genes · ${clusters.length} cluster rows`,14,'#526c7c');
    const groupMap=new Map();markers.forEach((m,j)=>{if(!groupMap.has(m.source_cluster))groupMap.set(m.source_cluster,[]);groupMap.get(m.source_cluster).push(j);});
    if(state.groups)for(const [c,indices] of groupMap){const x=left+(indices[0]-.5)*step,w=indices.length*step;rect(x,77,w-1,12,c==='7'?'#fff1d4':'#e9f1f5');text(x+w/2,88,'C'+c+(c==='7'?' †':''),12,'#526c7c',600,'middle');if(indices[0]>0)line(x,77,x,bottom-12);}
    clusters.forEach((c,i)=>{const y=top+i*dy;if(i%2===0)rect(13,y-dy/2,W-33,dy,'#f5f8fa');text(left-27,y+font*.34,'C'+c.id+(c.id==='11'?' ‡':''),font,'#183448',600,'end');
      const rowColor=c.program==='Superficial-associated'?'#d53c57':c.program==='Deep-associated'?'#3789c2':c.id==='8'?'#9065b6':c.id==='11'?'#a5b3bc':'#f4a35f';rect(14,y-9,5,18,rowColor);
      markers.forEach((m,j)=>{const r=lookup.get(c.id+'|'+m.gene_id),x=left+j*step;if(r&&Number.isFinite(r.fraction_detected)&&r.n_nuclei>0&&Number.isFinite(r.scaled_expression_z)){
        const rad=state.radius*Math.sqrt(Math.max(0,Math.min(1,r.fraction_detected)));circle(x,y,rad,colour(r.scaled_expression_z),state.opacity,`C${c.id} · ${m.gene_symbol}\nDetected ${(r.fraction_detected*100).toFixed(1)}% (${r.n_detected}/${r.n_nuclei})\nMean log1p CP10k: ${r.mean_log1p_normalized.toFixed(3)}\nSaved z: ${r.scaled_expression_z.toFixed(3)}`);
      }else{line(x-3,y-3,x+3,y+3,'#6e7b82');line(x-3,y+3,x+3,y-3,'#6e7b82');}});
    });
    markers.forEach((m,j)=>text(left+j*step+4,bottom+6,m.gene_symbol+(m.supported?'':' †'),font,m.supported?'#183448':'#81591a',400,'end',-55));
    const ly=H-68;text(24,ly-2,'Detected nuclei · dot area',13,'#526c7c');[.1,.5,1].forEach((f,i)=>{circle(235+i*80,ly-6,state.radius*Math.sqrt(f),'#69869a',state.opacity);text(250+i*80,ly-2,(f*100)+'%',12,'#526c7c');});
    const gx=W-225;for(let i=0;i<180;i++)rect(gx+i,ly-13,1.1,8,colour(state.min+i/179*(state.max-state.min)));text(gx,ly-22,'Saved expression z-score',13,'#526c7c');text(gx,ly+13,state.min.toFixed(1),12,'#526c7c');text(gx+180,ly+13,state.max.toFixed(1),12,'#526c7c',400,'end');
    text(24,H-21,'† C7 exploratory only.  ‡ C11 identity-conflicting.  Selection fixed from saved pooled markers.',12,'#526c7c');
    return {W,H,ops,markers,clusters,lookup,n};
  }
  function toSVG(c){return `<svg xmlns="${NS}" width="${c.W}" height="${c.H}" viewBox="0 0 ${c.W} ${c.H}" role="img" aria-label="Seurat CA1 top-two marker dotplot"><title>Seurat CA1 marker candidates, top two per cluster</title><desc>Dot area: detected fraction. Colour: saved gene-wise expression z-score. C7 candidates are exploratory. C11 identity remains conflicting.</desc><style>text{font-family:Arial,Helvetica,sans-serif}</style>${c.ops.map(o=>o.type==='rect'?`<rect x="${o.x}" y="${o.y}" width="${o.w}" height="${o.h}" fill="${o.fill}"/>`:o.type==='line'?`<line x1="${o.x}" y1="${o.y}" x2="${o.x2}" y2="${o.y2}" stroke="${o.stroke}" stroke-width="${o.width}"/>`:o.type==='circle'?`<circle cx="${o.x}" cy="${o.y}" r="${o.r}" fill="${o.fill}" opacity="${o.opacity}"><title>${esc(o.title)}</title></circle>`:`<text x="${o.x}" y="${o.y}" font-size="${o.size}" fill="${o.fill}" font-weight="${o.weight}" text-anchor="${o.align}"${o.angle?` transform="rotate(${o.angle} ${o.x} ${o.y})"`:''}>${esc(o.t)}</text>`).join('')}</svg>`;}
  function render(){
    for(const b of all('[data-action]'))if(['svg','pdf','png','values','genes','state'].includes(b.dataset.action))b.disabled=!(state.min<state.max);
    if(!(state.min<state.max)){$('.ctm-status').textContent='Colour minimum must be below the maximum. The last valid plot remains visible.';return;}
    current=chart();$('.ctm-scroll').innerHTML=toSVG(current);$('.ctm-status').textContent=`${current.markers.length} selected genes · ${current.clusters.length} clusters · ${state.supportedOnly?'supported representatives only':'C7 exploratory candidates flagged †'} · C8: CA1-prosubiculum`;
    $('.ctm-scroll').dataset.geneCount=String(current.markers.length);$('.ctm-scroll').dataset.clusterCount=String(current.clusters.length);
  }
  function sync(){all('[data-setting]').forEach(el=>{const v=state[el.dataset.setting];if(el.type==='checkbox')el.checked=v;else el.value=v;});choices();render();}
  function save(name,data,type){const url=URL.createObjectURL(new Blob([data],{type})),a=document.createElement('a');a.href=url;a.download=name;a.click();setTimeout(()=>URL.revokeObjectURL(url),2000);$('.ctm-action-message').textContent='Downloaded '+name;}
  function tsv(rows){if(!rows.length)return '';const keys=Object.keys(rows[0]),q=v=>v==null?'':String(v).replace(/[\t\r\n]/g,' ');return [keys.join('\t'),...rows.map(r=>keys.map(k=>q(r[k])).join('\t'))].join('\n')+'\n';}
  function values(){return current.clusters.flatMap(c=>current.markers.map(m=>({condition:state.condition,cluster:'C'+c.id,program:c.program,selection_cluster:'C'+m.source_cluster,selection_rank:m.source_rank,supported_representative:m.supported,...current.lookup.get(c.id+'|'+m.gene_id)})));}
  function pdf(){
    if(!window.jspdf?.jsPDF)throw Error('PDF library did not load. Please reload or use SVG.');
    const c=current,p=new window.jspdf.jsPDF({orientation:c.W>=c.H?'landscape':'portrait',unit:'pt',format:[c.W,c.H],compress:true});
    const pdfText=s=>String(s).replace(/†/g,'*').replace(/‡/g,'**');
    for(const o of c.ops){if(o.type==='rect'){p.setFillColor(o.fill);p.rect(o.x,o.y,o.w,o.h,'F');}else if(o.type==='line'){p.setDrawColor(o.stroke);p.setLineWidth(o.width);p.line(o.x,o.y,o.x2,o.y2);}else if(o.type==='circle'){if(o.opacity!==1)p.setGState(new p.GState({opacity:o.opacity}));p.setFillColor(o.fill);if(o.r>0)p.circle(o.x,o.y,o.r,'F');if(o.opacity!==1)p.setGState(new p.GState({opacity:1}));}else{p.setFont('helvetica',o.weight>=600?'bold':'normal');p.setFontSize(o.size);p.setTextColor(o.fill);if(o.angle){const width=p.getTextWidth(pdfText(o.t)),offset=o.align==='end'?width:o.align==='middle'?width/2:0,rad=o.angle*Math.PI/180;p.text(pdfText(o.t),o.x-offset*Math.cos(rad),o.y-offset*Math.sin(rad),{align:'left',angle:-o.angle});}else p.text(pdfText(o.t),o.x,o.y,{align:o.align==='middle'?'center':o.align==='end'?'right':'left'});}}
    p.setProperties({title:'Seurat CA1 top-two marker candidates',subject:'Read-only saved Seurat marker summaries; C7 exploratory; C11 identity-conflicting'});p.save('Seurat_CA1_top2_'+state.condition+'.pdf');$('.ctm-action-message').textContent='Downloaded vector PDF. * C7 exploratory; ** C11 identity-conflicting.';
  }
  async function png(){const svg=toSVG(current),url=URL.createObjectURL(new Blob([svg],{type:'image/svg+xml'}));try{const image=new Image();image.src=url;await image.decode();const canvas=document.createElement('canvas');canvas.width=current.W*3;canvas.height=current.H*3;const ctx=canvas.getContext('2d');ctx.scale(3,3);ctx.drawImage(image,0,0);const blob=await new Promise(resolve=>canvas.toBlob(resolve,'image/png'));save('Seurat_CA1_top2_'+state.condition+'.png',blob,'image/png');}finally{URL.revokeObjectURL(url);}}
  async function action(a){switch(a){case 'reset':state=initial();sync();break;case 'all-genes':state.genes=D.markers.map(m=>m.gene_id);sync();break;case 'no-genes':state.genes=[];sync();break;case 'svg':save('Seurat_CA1_top2_'+state.condition+'.svg',toSVG(current),'image/svg+xml');break;case 'pdf':pdf();break;case 'png':await png();break;case 'values':save('Seurat_CA1_top2_'+state.condition+'_values.tsv',tsv(values()),'text/tab-separated-values');break;case 'genes':save('Seurat_CA1_top2_selected_genes.tsv',tsv(current.markers),'text/tab-separated-values');break;case 'state':save('Seurat_CA1_top2_plot_settings.json',JSON.stringify({schema:D.schema,display_only:true,...state},null,2),'application/json');break;case 'load-state':$('[data-state-file]').click();}}
  root.addEventListener('input',e=>{const k=e.target.dataset.setting;if(k){state[k]=e.target.type==='checkbox'?e.target.checked:['condition','low','mid','high'].includes(k)?e.target.value:Number(e.target.value);render();}});
  root.addEventListener('change',e=>{if(e.target.dataset.gene){state.genes=all('[data-gene]:checked').map(x=>x.dataset.gene);render();}if(e.target.dataset.cluster){state.clusters=all('[data-cluster]:checked').map(x=>x.dataset.cluster);render();}});
  root.addEventListener('click',e=>{const a=e.target.closest('[data-action]');if(a)action(a.dataset.action).catch(err=>$('.ctm-action-message').textContent=err.message);});
  $('[data-state-file]').addEventListener('change',async e=>{try{const s=JSON.parse(await e.target.files[0].text());if(s.schema!==D.schema)throw Error('Settings are not for this dotplot.');const n=initial();for(const k of ['low','mid','high'])if(/^#[\da-f]{6}$/i.test(s[k]))n[k]=s[k];if(['ALL','WT','cKO'].includes(s.condition))n.condition=s.condition;for(const [k,min,max] of [['min',-100,100],['max',-100,100],['radius',3,13],['font',12,23],['opacity',.2,1]])if(Number.isFinite(s[k]))n[k]=Math.max(min,Math.min(max,s[k]));if(n.min>=n.max)throw Error('Colour range is invalid.');for(const k of ['supportedOnly','groups'])if(typeof s[k]==='boolean')n[k]=s[k];for(const k of ['genes','clusters'])if(Array.isArray(s[k]))n[k]=n[k].filter(v=>s[k].includes(v));state=n;sync();}catch(err){$('.ctm-action-message').textContent=err.message;}});
  window.SATB2_CA1_TOP_MARKERS_VIEW={get state(){return JSON.parse(JSON.stringify(state));},get chart(){return current;},values,render,toSVG};
  sync();
})();
