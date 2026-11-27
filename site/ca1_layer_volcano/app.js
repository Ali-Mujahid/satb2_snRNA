/* Exact saved Seurat fits. All controls below affect presentation/selection only. */
(function(root){'use strict';
const defaults={fit:'superficial',fdr:.05,fc:0,metric:'padj',xmin:'',xmax:'',ymin:'',ymax:'',xtick:'',ytick:'','fc-lines':true,'guide-fc':1,'fdr-line':true,boundary:true,'line-color':'#667782',up:'#ad3030',down:'#6c91e5',ns:'#c1c7ca','dot-size':2.4,opacity:.75,'label-size':15,width:900,height:650,labels:true,'top-n':10,custom:''};
const finite=x=>typeof x==='number'&&Number.isFinite(x),clamp=(x,a,b)=>Math.max(a,Math.min(b,x)),fallback=(x,v)=>x===''||x==null||!Number.isFinite(Number(x))?v:Number(x);
const esc=x=>String(x??'').replace(/[&<>"']/g,c=>({'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;',"'":'&apos;'}[c]));
const nice=span=>{const p=10**Math.floor(Math.log10(span/5));return (span/(5*p)<=1?1:span/(5*p)<=2?2:span/(5*p)<=5?5:10)*p;};
function ticks(a,b,step){if(!(step>0)||(b-a)/step>80)throw Error('Tick intervals must be positive and give at most 80 marks per axis.');const out=[];for(let v=Math.ceil((a-1e-9)/step)*step;v<=b+1e-9&&out.length<82;v+=step)out.push(Math.abs(v)<1e-9?0:v);return out;}
const fmt=v=>v===null||v===undefined?'NA':typeof v==='number'?(Math.abs(v)>0&&Math.abs(v)<.0001?v.toExponential(3):Number(v.toPrecision(5))):v;
const significant=(r,s)=>finite(r.padj)&&r.padj<s.fdr&&finite(r.log2FoldChange)&&Math.abs(r.log2FoldChange)>=s.fc&&r.log2FoldChange!==0;
function records(fit){return fit.rows.map(r=>Object.fromEntries(fit.columns.map((c,i)=>[c,r[i]])));}
function build(fit,s){
 s={...defaults,...s};
 if(!['padj','pvalue'].includes(s.metric))throw Error('Choose adjusted P or raw P for the Y-axis.');
 for(const [key,min,max] of [['fdr',.000001,1],['fc',0,30],['guide-fc',0,30],['dot-size',.5,7],['opacity',.1,1],['label-size',12,28],['top-n',0,50],['width',650,1600],['height',500,1200]])if(!finite(s[key])||s[key]<min||s[key]>max)throw Error(`${key} must be between ${min} and ${max}.`);
 const rows=records(fit),W=clamp(Math.round(s.width),650,1600),H=clamp(Math.round(s.height),500,1200),left=95,right=W-35,top=105,bottom=H-120,ops=[];
 const add=(type,o)=>ops.push({type,...o}),text=(t,x,y,size=14,o={})=>add('text',{t,x,y,size,fill:'#203c4c',align:'start',...o}),line=(x,y,x2,y2,o={})=>add('line',{x,y,x2,y2,stroke:s['line-color'],width:1,...o});
 const valid=rows.filter(r=>finite(r.log2FoldChange)&&finite(r[s.metric])&&r[s.metric]>=0&&r[s.metric]<=1).map(r=>({...r,x:r.log2FoldChange,y:-Math.log10(Math.max(r[s.metric],1e-300)),sig:significant(r,s)}));
 const maxX=valid.reduce((a,r)=>Math.max(a,Math.abs(r.x)),1),maxY=valid.reduce((a,r)=>Math.max(a,r.y),1);
 const xmin=fallback(s.xmin,-Math.ceil(maxX)),xmax=fallback(s.xmax,Math.ceil(maxX)),ymin=fallback(s.ymin,0),ymax=fallback(s.ymax,Math.ceil(maxY*1.08));
 if(!(xmin<xmax&&ymin>=0&&ymin<ymax))throw Error('Use X minimum < X maximum and 0 ≤ Y minimum < Y maximum.');
 if(!(s.fdr>0&&s.fdr<=1&&s.fc>=0))throw Error('Adjusted P threshold must be >0 and ≤1; fold-change cutoff must be ≥0.');
 const dx=fallback(s.xtick,nice(xmax-xmin)),dy=fallback(s.ytick,nice(ymax-ymin)),xt=ticks(xmin,xmax,dx),yt=ticks(ymin,ymax,dy),xp=x=>left+(x-xmin)/(xmax-xmin)*(right-left),yp=y=>bottom-(y-ymin)/(ymax-ymin)*(bottom-top);
 add('rect',{x:0,y:0,w:W,h:H,fill:'#ffffff'});
 text(fit.metadata.label+' · cKO vs WT',24,30,22,{weight:600});text('Seurat mouse pseudobulk · 6 WT + 6 cKO · C8 excluded',24,53,13,{fill:'#526c79'});
 const degs=rows.filter(r=>significant(r,s)),up=degs.filter(r=>r.log2FoldChange>0).length,down=degs.length-up;
 [['up',`${up} higher in cKO`],['down',`${down} lower in cKO`],['ns','Below cutoffs / unavailable FDR']].forEach(([k,name],i)=>{const x=left+i*(right-left)/3;add('circle',{x,y:79,r:4,fill:s[k]});text(name,x+10,83,12);});
 xt.forEach(v=>{line(xp(v),bottom,xp(v),bottom+6);text(fmt(v),xp(v),bottom+25,13,{align:'middle'});});
 yt.forEach(v=>{line(left-6,yp(v),left,yp(v));text(fmt(v),left-12,yp(v)+4,13,{align:'end'});});
 if(s['fc-lines'])for(const v of new Set([-s['guide-fc'],s['guide-fc']]))if(v>=xmin&&v<=xmax)line(xp(v),top,xp(v),bottom,{dash:'2 4',width:1.2,role:'fc-guide',value:v});
 if(s['fdr-line']&&s.metric==='padj'){const v=-Math.log10(s.fdr);if(v>=ymin&&v<=ymax)line(left,yp(v),right,yp(v),{dash:'5 4',role:'fdr-guide',value:v});}
 const visible=valid.filter(r=>r.x>=xmin&&r.x<=xmax&&r.y>=ymin&&r.y<=ymax).sort((a,b)=>Number(a.sig)-Number(b.sig));
 for(const r of visible)add('circle',{x:xp(r.x),y:yp(r.y),r:s['dot-size'],fill:r.sig?(r.x>0?s.up:s.down):s.ns,opacity:s.opacity,gene_id:r.gene_id,gene:r.gene_symbol,status:r.sig?(r.x>0?'up':'down'):'below',tip:`${r.gene_symbol||r.gene_id} · log2FC ${r.x} · adjusted P ${r.padj??'NA'} · raw P ${r.pvalue??'NA'}`});
 line(left,top,left,bottom,{width:1.3});line(left,bottom,right,bottom,{width:1.3});
 if(s.boundary){line(left,top,right,top,{role:'plot-boundary',width:1.3});line(right,top,right,bottom,{role:'plot-boundary',width:1.3});}
 text('log2 fold change (cKO / WT)',(left+right)/2,H-66,16,{align:'middle'});text(s.metric==='padj'?'−log10 adjusted P':'−log10 raw P',25,(top+bottom)/2,16,{align:'middle',angle:-90});
 const ranked=visible.filter(r=>r.sig).sort((a,b)=>a.padj-b.padj||Math.abs(b.x)-Math.abs(a.x)||a.gene_id.localeCompare(b.gene_id)),wanted=new Map((s.labels?ranked.slice(0,Math.floor(s['top-n'])):[]).map(r=>[r.gene_id,r]));
 const custom=new Set(s.custom.split(/[\s,;]+/).filter(Boolean).map(g=>g.toLowerCase()));for(const r of visible)if(custom.has(r.gene_id.toLowerCase())||custom.has((r.gene_symbol||'').toLowerCase()))wanted.set(r.gene_id,r);
 const boxes=[],labelRows=[],size=s['label-size'];let omitted=0;
 // Search increasingly distant positions; avoid text collisions and preserve point identity with leaders.
 for(const r of wanted.values()){
  const name=r.gene_symbol||r.gene_id,width=name.length*size*.62+6,px=xp(r.x),py=yp(r.y);let found=null;
  for(const offset of [0,1,-1,2,-2,3,-3,4,-4,5,-5,6,-6,7,-7,8,-8,9,-9,10,-10]){
   for(const side of [r.x>=0?1:-1,r.x>=0?-1:1]){
    const x=clamp(side>0?px+9:px-width-9,left+4,right-width-4),y=clamp(py-9+offset*(size+5),top+size+4,bottom-8),box={x,y:y-size,w:width,h:size+5};
    if(!boxes.some(b=>box.x<b.x+b.w+4&&box.x+box.w+4>b.x&&box.y<b.y+b.h+2&&box.y+box.h+2>b.y)){found={x,y,box};break;}
   }if(found)break;
  }
  if(!found){omitted++;continue;}boxes.push(found.box);line(px,py,clamp(px,found.x,found.x+width),found.y-size/2,{stroke:'#78848b',width:.6,role:'label-leader'});text(name,found.x,found.y,size,{weight:600,halo:true,role:'gene-label',gene_id:r.gene_id});labelRows.push(r);
 }
 text(`Adjusted P < ${s.fdr}; |log2FC| ≥ ${s.fc}. Dotted FC guides: ±${s['guide-fc']} log2FC.`,24,H-38,12,{fill:'#526c79'});
 text(`${visible.length.toLocaleString()} / ${rows.length.toLocaleString()} model genes shown · ${labelRows.length} gene labels · ${valid.length-visible.length} outside axes · ${rows.length-valid.length} missing coordinates`,24,H-18,12,{fill:'#526c79'});
 return {W,H,ops,rows,visible,degs,up,down,labelRows,requested:wanted.size,omitted,missing:rows.length-valid.length,outside:valid.length-visible.length,axes:{xmin,xmax,ymin,ymax,dx,dy},zeros:valid.filter(r=>r[s.metric]===0).length};
}
function svg(c){return `<svg xmlns="http://www.w3.org/2000/svg" width="${c.W}" height="${c.H}" viewBox="0 0 ${c.W} ${c.H}" role="img" aria-label="CA1 layer volcano plot">`+c.ops.map(o=>{
 const role=o.role?` data-role="${o.role}"`:'',attrs=o.gene_id?` data-gene-id="${esc(o.gene_id)}"`:'';
 if(o.type==='rect')return `<rect x="${o.x}" y="${o.y}" width="${o.w}" height="${o.h}" fill="${o.fill}"/>`;
 if(o.type==='line')return `<line x1="${o.x}" y1="${o.y}" x2="${o.x2}" y2="${o.y2}" stroke="${o.stroke}" stroke-width="${o.width}"${o.dash?` stroke-dasharray="${o.dash}"`:''}${role}/>`;
 if(o.type==='circle')return `<circle cx="${o.x}" cy="${o.y}" r="${o.r}" fill="${o.fill}" opacity="${o.opacity??1}"${attrs}${o.status?` data-status="${o.status}"`:''}>${o.tip?`<title>${esc(o.tip)}</title>`:''}</circle>`;
 return `<text x="${o.x}" y="${o.y}" font-family="Arial,Helvetica,sans-serif" font-size="${o.size}" font-weight="${o.weight||400}" text-anchor="${o.align}" fill="${o.fill}"${o.angle?` transform="rotate(${o.angle} ${o.x} ${o.y})"`:''}${o.halo?' stroke="white" stroke-width="3" paint-order="stroke" stroke-linejoin="round"':''}${role}${attrs}>${esc(o.t)}</text>`;
 }).join('')+'</svg>';}
const API={defaults,build,svg,records,significant};if(typeof module!=='undefined')module.exports=API;
if(typeof document==='undefined')return;
const $=id=>document.getElementById(id),M=root.CA1_LAYER_VOLCANO_MANIFEST,cache={},bools=Object.keys(defaults).filter(k=>typeof defaults[k]==='boolean'),strings=['fit','metric','custom','line-color','up','down','ns'],axes=['xmin','xmax','ymin','ymax','xtick','ytick'];let current=null,selected=null,sequence=0;
function state(){return Object.fromEntries(Object.keys(defaults).map(k=>[k,bools.includes(k)?$(k).checked:strings.includes(k)?$(k).value:axes.includes(k)&&$(k).value===''?'':fallback($(k).value,defaults[k])]));}
function setState(s){Object.entries(s).forEach(([k,v])=>{if(!(k in defaults))return;if(bools.includes(k))$(k).checked=v;else $(k).value=v;});}
function ready(v){document.querySelectorAll('[data-export]').forEach(b=>b.disabled=!v);}
function table(rows,cols){return '<table><thead><tr>'+cols.map(c=>`<th>${esc(c)}</th>`).join('')+'</tr></thead><tbody>'+rows.map(r=>'<tr>'+cols.map(c=>`<td>${esc(fmt(r[c]))}</td>`).join('')+'</tr>').join('')+'</tbody></table>';}
function render(){if(!selected)return;try{
 const s=state(),fit=cache[selected];current=build(fit,s);$('error').textContent='';$('figure').innerHTML=svg(current);$('layer-title').textContent=fit.metadata.label;
 $('summary').textContent=`Saved FDR < 0.05: ${fit.metadata.DEGs.toLocaleString()} DEGs (${fit.metadata.up} up / ${fit.metadata.down} down). Current display cutoffs: ${current.degs.length.toLocaleString()} DEGs (${current.up} up / ${current.down} down). ${fit.metadata.n_nuclei.toLocaleString()} nuclei; 6 WT + 6 cKO mice.`;
 $('fc-ratio').textContent=`|log2FC| ${s.fc} corresponds to ≥${fmt(2**s.fc)}-fold in either direction.`;
 $('status').textContent=`${current.visible.length.toLocaleString()} model genes visible; ${current.outside} outside axes; ${current.missing} unavailable coordinates; ${current.zeros} zero ${s.metric} values displayed at 10⁻³⁰⁰. ${current.labelRows.length}/${current.requested} requested labels shown; ${current.omitted} omitted to avoid overlap. ${s.metric==='pvalue'?'Dot significance and ranking still use adjusted P; the FDR guide is hidden on the raw-P axis.':''}`;
 $('top-table').innerHTML=table(current.labelRows,['gene_symbol','gene_id','log2FoldChange','padj']);$('mouse-table').innerHTML=table(fit.mice,['sample','condition','n_nuclei','total_UMIs']);ready(true);
 }catch(e){current=null;ready(false);$('figure').textContent='';$('error').textContent=e.message;}}
async function digest(buf){return Array.from(new Uint8Array(await crypto.subtle.digest('SHA-256',buf))).map(x=>x.toString(16).padStart(2,'0')).join('');}
async function fetchFit(key){if(cache[key])return;const d=M.fits[key],response=await fetch(d.url);if(!response.ok)throw Error('Could not load saved results ('+response.status+').');const packed=await response.arrayBuffer();if(await digest(packed)!==d.sha256)throw Error('Saved-fit checksum mismatch.');let raw;
 if(typeof DecompressionStream==='function')raw=await new Response(new Blob([packed]).stream().pipeThrough(new DecompressionStream('gzip'))).arrayBuffer();else{if(!root.fflate)await new Promise((ok,bad)=>{const s=document.createElement('script');s.src='../vendor/fflate.min.js';s.onload=ok;s.onerror=bad;document.head.append(s);});raw=root.fflate.gunzipSync(new Uint8Array(packed));}
 if(await digest(raw)!==d.payload_sha256)throw Error('Decoded-fit checksum mismatch.');const fit=JSON.parse(new TextDecoder().decode(raw));if(fit.metadata.id!==key||fit.rows.length!==d.n_rows)throw Error('Unexpected fit identity.');cache[key]=fit;
}
async function select(){const key=$('fit').value;if(!M.fits[key])return;const token=++sequence;selected=null;current=null;ready(false);$('figure').textContent='';$('layer-title').textContent='Loading '+M.fits[key].label+'…';$('status').textContent='Verifying saved results…';try{await fetchFit(key);if(token!==sequence)return;selected=key;const d=M.fits[key];$('source').textContent=`${d.label}; source SHA256 ${d.source_sha256}. ${d.n_rows} modeled genes; cluster membership ${d.clusters.map(k=>'C'+k).join(' + ')}. Existing ${d.fit_date} fit. No model refitted for this page.`;history.replaceState(null,'','#'+key);render();}catch(e){if(token===sequence){$('error').textContent=e.message;$('status').textContent='Results unavailable; missing values are not replaced with zeros.';}}}
const filename=()=>`Seurat_CA1_${selected==='superficial'?'Superficial_C2_C3':'Deep_C1_C5'}_volcano`;
function save(name,data,type){const u=URL.createObjectURL(new Blob([data],{type})),a=document.createElement('a');a.href=u;a.download=name;a.click();setTimeout(()=>URL.revokeObjectURL(u),3000);$('message').textContent='Downloaded '+name;}
function tsv(rows){const cols=cache[selected].columns;return [cols.join('\t'),...rows.map(r=>cols.map(k=>r[k]==null?'NA':String(r[k]).replace(/[\t\r\n]/g,' ')).join('\t'))].join('\n')+'\n';}
function pdf(){const p=new root.jspdf.jsPDF({unit:'pt',format:[current.W,current.H],orientation:current.W>=current.H?'landscape':'portrait',compress:true});let opacity=1;const ascii=s=>String(s).replace(/[−–—]/g,'-').replace(/≥/g,'>=').replace(/±/g,'+/-');
 for(const o of current.ops){const op=o.opacity??1;if(op!==opacity){p.setGState(new p.GState({opacity:op}));opacity=op;}
  if(o.type==='line'){p.setDrawColor(o.stroke);p.setLineWidth(o.width);p.setLineDashPattern(o.dash?o.dash.split(' ').map(Number):[],0);p.line(o.x,o.y,o.x2,o.y2);}
  else if(o.type==='circle'||o.type==='rect'){p.setFillColor(o.fill);if(o.type==='circle')p.circle(o.x,o.y,o.r,'F');else p.rect(o.x,o.y,o.w,o.h,'F');}
  else{p.setFont('helvetica',o.weight>=600?'bold':'normal');p.setFontSize(o.size);p.setTextColor(o.fill);const t=ascii(o.t);if(o.angle){const width=p.getTextWidth(t),offset=o.align==='middle'?width/2:0,rad=o.angle*Math.PI/180;p.text(t,o.x-offset*Math.cos(rad),o.y-offset*Math.sin(rad),{angle:-o.angle});}else{if(o.halo){p.setFillColor('#ffffff');p.rect(o.x-1,o.y-o.size*.82,p.getTextWidth(t)+2,o.size,'F');}p.text(t,o.x,o.y,{align:o.align==='middle'?'center':o.align==='end'?'right':'left'});}}
 }p.setProperties({title:cache[selected].metadata.label+' cKO vs WT',subject:'Existing Seurat mouse-level DESeq2; display controls only'});p.save(filename()+'.pdf');$('message').textContent='Downloaded vector PDF';}
async function exportFile(type){if(!current)return;const name=filename();if(type==='pdf')return pdf();if(type==='svg')return save(name+'.svg',svg(current),'image/svg+xml');if(type==='degs'||type==='all'||type==='labels')return save(name+'_'+type+'.tsv',tsv(type==='degs'?current.degs:type==='labels'?current.labelRows:current.rows),'text/tab-separated-values');if(type==='state')return save(name+'_settings.json',JSON.stringify({schema:M.schema,display_only:true,source_sha256:M.fits[selected].source_sha256,state:state()},null,2),'application/json');
 if(type==='png'){const u=URL.createObjectURL(new Blob([svg(current)],{type:'image/svg+xml'}));try{const im=new Image();im.src=u;await im.decode();const c=document.createElement('canvas');c.width=current.W*3;c.height=current.H*3;c.getContext('2d').drawImage(im,0,0,c.width,c.height);save(name+'.png',await new Promise(ok=>c.toBlob(ok,'image/png')),'image/png');}finally{URL.revokeObjectURL(u);}}
}
let timer;Object.keys(defaults).filter(k=>k!=='fit').forEach(k=>$(k).addEventListener('input',()=>{clearTimeout(timer);timer=setTimeout(render,120);}));$('fit').addEventListener('change',select);window.addEventListener('hashchange',()=>{const key=location.hash.slice(1);if(M.fits[key]&&$('fit').value!==key){$('fit').value=key;select();}});
$('auto-axes').onclick=()=>{axes.forEach(k=>$(k).value='');render();};$('reset').onclick=()=>{setState({...defaults,fit:$('fit').value});render();};document.querySelectorAll('[data-export]').forEach(b=>b.onclick=()=>exportFile(b.dataset.export).catch(e=>$('error').textContent=e.message));$('load-state').onclick=()=>$('state-file').click();
$('state-file').onchange=async e=>{try{const d=JSON.parse(await e.target.files[0].text());if(d.schema!==M.schema||!M.fits[d.state?.fit])throw Error('Not a CA1 revised-layer plot settings file.');if(d.source_sha256!==M.fits[d.state.fit].source_sha256)throw Error('Settings source does not match the saved fit.');for(const [k,v] of Object.entries(d.state)){if(!(k in defaults))continue;if(bools.includes(k)&&typeof v!=='boolean')throw Error('Invalid checkbox setting.');if(['up','down','ns','line-color'].includes(k)&&!/^#[a-f0-9]{6}$/i.test(v))throw Error('Invalid colour setting.');if(!strings.includes(k)&&!bools.includes(k)&&v!==''&&!finite(v))throw Error('Invalid numeric setting.');}setState({...defaults,...d.state});await select();$('message').textContent='Loaded plot settings';}catch(err){$('error').textContent=err.message;}};
root.CA1_LAYER_VOLCANO={...API,get current(){return current;},get state(){return state();},get ready(){return !!selected;},cache,render};const start=location.hash.slice(1);if(M.fits[start])$('fit').value=start;select();
})(typeof window==='undefined'?globalThis:window);
