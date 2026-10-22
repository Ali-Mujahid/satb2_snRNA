(function(){
'use strict';
const $=id=>document.getElementById(id),M=window.CCCModel,P=window.CCCPlots,N=window.CCCNetworks;
const esc=s=>String(s==null?'':s).replace(/[&<>"']/g,c=>({'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;',"'":'&#39;'}[c]));
const DEFAULT={version:1,scope:'IN_clusters',condition:'pooled',compare:false,direction:'both',tier:'DEG_supported',maxRank:.05,bothDEG:false,search:'',sources:[],targets:[],view:'dotplot',dotMetric:'magnitude_rank',heatMetric:'count',palette:'viridis',lowColor:'#440154',midColor:'#21918c',highColor:'#fde725',reverseScale:false,scaleMode:'auto',colorMin:0,colorMax:.05,opacity:.8,fontSize:12,topN:24,pairLimit:20,order:'rank',dotRadius:9,labelMode:'compact',edgeWidth:5,arrowSize:10,curvature:.35,arcWidth:17,arcGap:3,ribbonScale:1,nodeScale:1,edgeColorMode:'source',edgeColor:'#397eb5',arcColorMode:'node',arcColor:'#9b6cb0',geneEdges:12,nodeColorMode:'palette',nodeColor:'#438a9c',nodeColors:{}};
const NUMBERS={maxRank:[0,1],colorMin:[-1e9,1e9],colorMax:[-1e9,1e9],opacity:[.1,1],fontSize:[8,20],topN:[1,120],pairLimit:[1,80],dotRadius:[2,16],edgeWidth:[.5,16],arrowSize:[0,24],curvature:[0,.9],arcWidth:[2,40],arcGap:[0,10],ribbonScale:[.1,1],nodeScale:[.3,3],geneEdges:[1,30]};
const BOOLEANS=['compare','bothDEG','reverseScale'],COLOR_KEYS=['lowColor','midColor','highColor','edgeColor','arcColor','nodeColor'];
const ENUMS={scope:['IN_clusters','IN_total'],condition:['pooled','WT','cKO'],direction:['both','CA1_to_IN','IN_to_CA1'],tier:['DEG_supported','replicated_specific','all_candidates'],view:['dotplot','heatmap','circos','circle','geneCircle','table','findings','downloads'],dotMetric:['magnitude_rank','lr_gmeans','lr_logfc','cellphone_pvals'],heatMetric:['count','mean_gmeans'],palette:['viridis','blue','rdBu','custom'],scaleMode:['auto','manual'],order:['rank','alphabetical'],labelMode:['compact','annotation','id'],edgeColorMode:['source','target','metric','custom'],arcColorMode:['node','custom'],nodeColorMode:['palette','custom']};
const PALETTE=['#b86740','#438db5','#6a40b7','#33a77e','#bb7bd4','#4f4bbc','#7eb336','#cf485b','#9ba6ae','#b243a9','#16aa94','#5443b1','#af88bd','#b95c3b','#4c74c0','#be4c87','#6ea291','#8a84b3','#a39f7f','#b68592','#60919b','#77873f'];
const STORE='satb2-communication-v1';let data,state,tableCondition='pooled',tablePage=0,pairPage=0,renderToken=0,timer,displayed=[],filteredCache={},currentTable=[];
const format=v=>v===null||v===undefined?'NA':typeof v==='boolean'?(v?'yes':'no'):!Number.isFinite(Number(v))?String(v):Number(v)===0?'0':Math.abs(Number(v))<.001?Number(v).toExponential(3):Number.isInteger(Number(v))?Number(v).toLocaleString('en-US'):Number(v).toPrecision(4);
function valid(s){
 const out={...DEFAULT,sources:[],targets:[],nodeColors:{}};
 if(!s||typeof s!=='object')return out;
 for(const [k,range] of Object.entries(NUMBERS)){if(M.finite(s[k]))out[k]=Math.max(range[0],Math.min(range[1],Number(s[k])));}
 ['topN','pairLimit','geneEdges','fontSize'].forEach(k=>out[k]=Math.round(out[k]));
 for(const k of BOOLEANS)if(typeof s[k]==='boolean')out[k]=s[k];
 for(const [k,choices] of Object.entries(ENUMS))if(choices.includes(s[k]))out[k]=s[k];
 for(const k of COLOR_KEYS)if(/^#[0-9a-f]{6}$/i.test(s[k]||''))out[k]=s[k];
 if(typeof s.search==='string')out.search=s.search.slice(0,200);
 for(const k of ['sources','targets'])if(Array.isArray(s[k]))out[k]=s[k].filter(x=>typeof x==='string'&&data.nodes.some(n=>n.id===x));
 if(s.nodeColors&&typeof s.nodeColors==='object')for(const [k,v] of Object.entries(s.nodeColors))if(data.nodes.some(n=>n.id===k)&&/^#[0-9a-f]{6}$/i.test(v))out.nodeColors[k]=v;
 return out;
}
function save(){try{localStorage.setItem(STORE,JSON.stringify(state));}catch(_){}}
function baseNodes(){return [...data.nodes].filter(n=>state.scope==='IN_total'?!/^IN_C\d+$/.test(n.id):n.id!=='IN_Total').sort((a,b)=>nodeOrder(a.id)-nodeOrder(b.id));}
function nodeOrder(id){return id==='CA1_Deep'?0:id==='CA1_Superficial'?1:id==='IN_Total'?2:2+Number(id.replace('IN_C',''));}
function label(n){const compact=n.id==='CA1_Deep'?'CA1 DL':n.id==='CA1_Superficial'?'CA1 SL':n.id==='IN_Total'?'IN total':n.id.replace('IN_','');return state.labelMode==='id'?n.id:state.labelMode==='annotation'?compact+' · '+n.annotation+(n.needs_review?' *':''):compact+(n.needs_review?' *':'');}
function nodes(){const ids=new Set(selectedConditions().flatMap(c=>getPairs(c).flatMap(p=>[p.source,p.target]))),all=baseNodes();if(!ids.size){for(const id of [...state.sources,...state.targets])ids.add(id);}return all.filter(n=>!ids.size||ids.has(n.id)).map(n=>({...n,label:label(n),color:state.nodeColors[n.id]||PALETTE[nodeOrder(n.id)%PALETTE.length]}));}
function selectedConditions(){return state.compare?['WT','cKO']:[state.condition];}
function selected(cond){return filteredCache[cond]||(filteredCache[cond]=M.filter(data.datasets[state.scope+'|'+cond].rows,state));}
function getPairs(cond){return M.eligible(data.datasets[state.scope+'|'+cond].eligible_pairs,state);}
function sync(){
 for(const k of Object.keys(DEFAULT)){const el=$(k);if(!el)continue;if(el.multiple){for(const option of el.options)option.selected=state[k].includes(option.value);}else if(el.type==='checkbox')el.checked=state[k];else el.value=state[k];}
 for(const k of Object.keys(NUMBERS))if($(k+'-value'))$(k+'-value').textContent=state[k];
 $('condition').disabled=state.compare;
 $('tier-note').textContent=state.tier==='replicated_specific'?'Adds specificity rank ≤0.05, CellPhoneDB P≤0.05, ≥3 supporting WT and cKO mice, and reference-supported atomic IN identities. Magnitude cutoff capped at0.05.':state.tier==='all_candidates'?'Includes identity-review groups and candidates without matching endpoint DEGs. Native expression gates remain fixed.':'At least one original DESeq2 DEG must belong to the ligand’s source or receptor’s target. This is gene-level evidence.';
 document.querySelectorAll('[data-view]').forEach(b=>b.classList.toggle('active',b.dataset.view===state.view));
 $('table-condition').value=tableCondition;
}
function populationOptions(){
 const all=baseNodes(),html=all.map(n=>'<option value="'+esc(n.id)+'">'+esc(n.id+' · '+n.annotation+(n.needs_review?' * review':''))+'</option>').join('');
 $('sources').innerHTML=html;$('targets').innerHTML=html;$('node-target').innerHTML=html;sync();updateNodeEditor();
}
function updateNodeEditor(){const n=data.nodes.find(n=>n.id===$('node-target').value);if(n)$('population-color').value=state.nodeColors[n.id]||PALETTE[nodeOrder(n.id)%PALETTE.length];}
function change(key,value){state[key]=value;state=valid(state);pairPage=0;tablePage=0;filteredCache={};if(key==='condition')tableCondition=state.condition;if(key==='compare')tableCondition=state.compare?'WT':state.condition;if(['lowColor','midColor','highColor'].includes(key))state.palette='custom';if(key==='edgeColor')state.edgeColorMode='custom';if(key==='arcColor')state.arcColorMode='custom';if(key==='nodeColor')state.nodeColorMode='custom';if(key==='scope'){state.sources=[];state.targets=[];populationOptions();}sync();save();schedule();}
function schedule(){clearTimeout(timer);timer=setTimeout(render,100);}
function options(){return Object.fromEntries(Object.entries(state).filter(([k])=>!['sources','targets','nodeColors'].includes(k)));}
function dotOrdering(rows){return P.helpers.selectDotplot(rows,{options:{...options(),topN:120,pairLimit:100}});}
function commonContext(rowsByCondition){
 const o=options(),conds=Object.keys(rowsByCondition),union=conds.flatMap(c=>rowsByCondition[c]),ns=nodes(),ordering=dotOrdering(union);
 const allPairs=ordering.pairs.map(p=>p.key),maxPage=Math.max(0,Math.ceil(allPairs.length/state.pairLimit)-1);pairPage=Math.min(pairPage,maxPage);
 const pairOrder=allPairs.slice(pairPage*state.pairLimit,(pairPage+1)*state.pairLimit),lrOrder=ordering.lrs.slice(0,state.topN).map(l=>l.key);
 const stats=conds.map(c=>M.counts(rowsByCondition[c]));
 const ctx={nodes:ns,options:o,scope:state.scope,direction:state.direction,lrOrder,pairOrder,totalLrs:ordering.totalLrs,totalPairs:ordering.totalPairs,edgeMax:Math.max(1,...stats.map(s=>s.edgeMax)),degreeMax:Math.max(1,...stats.map(s=>s.degreeMax))};
 const keys=new Set(lrOrder),pairs=new Set(pairOrder);
 if(state.view==='dotplot')ctx.colorDomain=P.helpers.extent(union.filter(r=>keys.has(P.helpers.lrKey(r))&&pairs.has(P.helpers.pairKey(r))).map(r=>r[state.dotMetric]),state.dotMetric!=='lr_logfc');
 if(state.view==='heatmap'){const values=conds.flatMap(c=>P.helpers.aggregateHeatmap(rowsByCondition[c],{nodes:ns,eligiblePairs:getPairs(c)}).map(v=>v[state.heatMetric]));ctx.colorDomain=P.helpers.extent(values,true);}
 if(state.view==='circle'||state.view==='circos')ctx.colorDomain=[0,ctx.edgeMax];
 const sorted=M.sortRows(union),seen=new Set();ctx.geneKeys=[];for(const r of sorted){const k=N.geneKey(r);if(!seen.has(k)){seen.add(k);ctx.geneKeys.push(k);if(ctx.geneKeys.length>=state.geneEdges)break;}}
 if(state.view==='geneCircle'){const geneRows=union.filter(r=>ctx.geneKeys.includes(N.geneKey(r)));ctx.colorDomain=P.helpers.extent(geneRows.map(r=>r.lr_gmeans),true);ctx.edgeMax=Math.max(1e-12,...geneRows.map(r=>Number(r.lr_gmeans)||0));}
 if(state.scaleMode==='manual')ctx.colorDomain=[state.colorMin,state.colorMax];
 $('pair-page').textContent='Pairs '+(allPairs.length?pairPage*state.pairLimit+1:0)+'–'+Math.min(allPairs.length,(pairPage+1)*state.pairLimit)+' of '+allPairs.length+' · top '+lrOrder.length+' LR identities';
 $('pair-prev').disabled=pairPage===0;$('pair-next').disabled=pairPage>=maxPage;
 return ctx;
}
const GUIDES={dotplot:['Ligand–receptor dotplot','Color: selected native score. Dot area: minimum detected fraction across every complex subunit. Shared LR rows, population-pair columns and numeric scales in WT/cKO comparisons.'],heatmap:['Interaction heatmap','Rows are senders; columns are receivers. Gray NA cells are untested directions. Zero counts only apply to eligible pairs with no candidates passing the selected filters.'],circos:['Directed circos ribbons','Arcs identify populations. Ribbons show qualifying ordered LR counts as each sector’s fraction; arrow tips identify the receiver. Circle and heatmap show absolute counts.'],circle:['Directed population circle','Arrow width encodes distinct qualifying LR identities; node area encodes incoming + outgoing counts. Stable population positions and shared size scales support WT/cKO inspection.'],geneCircle:['Gene / complex circle','Nodes combine population, complex and ligand/receptor role. Edge width is native LR geometric-mean expression; node area is detected fraction. Selected interaction keys and positions are shared across conditions.']};
function render(){
 try{
  renderToken++;filteredCache={};$('selected-detail').hidden=true;const conds=selectedConditions(),rowsByCondition=Object.fromEntries(conds.map(c=>[c,selected(c)]));
  const primary=rowsByCondition[conds[0]],summary=M.counts(primary),stats=conds.map(c=>({c,...M.counts(rowsByCondition[c])}));
  $('stats').innerHTML=[['Filtered records',stats.map(s=>s.c+': '+format(s.records)).join(' · ')],['Distinct molecular LR pairs',stats.map(s=>s.c+': '+format(s.molecules)).join(' · ')],['Directed population pairs',stats.map(s=>s.c+': '+format(s.pairs)).join(' · ')],['Independent mice','6 WT + 6 cKO']].map(([name,v])=>'<div class="stat"><b style="font-size:'+(v.length>21?'1.1rem':'1.6rem')+'">'+esc(v)+'</b><span>'+esc(name)+'</span></div>').join('');
  $('selection-note').textContent=(state.scope==='IN_total'?'Inhibitory total is a separate fit; constituent clusters are not mixed into this dataset. ':'All 20 inhibitory clusters remain available, including groups flagged for identity review. ')+(state.compare?'WT/cKO counts and scores are descriptive, with shared display axes. ':'The same molecular LR in different population pairs is a separate record. ');
  ['chart','table','findings','downloads'].forEach(v=>$(v+'-view').hidden=v==='chart'?!GUIDES[state.view]:state.view!==v);
  if(GUIDES[state.view]){
   $('chart-title').textContent=GUIDES[state.view][0];$('chart-guide').textContent=GUIDES[state.view][1];$('dot-pagination').hidden=state.view!=='dotplot';
   if(state.scaleMode==='manual'&&state.colorMax<=state.colorMin)throw Error('The color maximum must be greater than the color minimum.');
   const ctx=commonContext(rowsByCondition);displayed=[];$('plots').classList.toggle('compare',state.compare);
   $('plots').innerHTML=conds.map((c,i)=>{const rows=rowsByCondition[c],context={...ctx,condition:c,eligiblePairs:getPairs(c)},svg=state.view==='dotplot'?P.dotplot(rows,context):state.view==='heatmap'?P.heatmap(rows,context):N[state.view](rows,context);displayed.push({condition:c,svg,rows,context,settings:JSON.parse(JSON.stringify(state))});return '<article class="plot"><h3>'+esc(c==='pooled'?'Pooled · all mice':c)+'</h3><p class="plot-subtitle small">'+format(rows.length)+' filtered records · '+esc(state.scope==='IN_clusters'?'individual IN clusters':'inhibitory total')+'</p><div class="svg-wrap">'+svg+'</div><div class="plot-actions"><button data-plot="'+i+'" data-format="pdf">PDF</button><button data-plot="'+i+'" data-format="svg">SVG</button><button data-plot="'+i+'" data-format="png">PNG</button><button data-plot="'+i+'" data-format="json">Plot state &amp; data</button></div></article>';}).join('');
   // Each SVG is standalone. Namespace its CSS and IDs so side-by-side SVGs cannot affect each other.
   $('plots').querySelectorAll('.svg-wrap>svg').forEach((svg,i)=>scopeSVG(svg,'ccc-live-'+i));
   $('render-status').textContent=state.view==='geneCircle'?'Gene-network limits apply to population × LR records, not necessarily distinct molecular LR identities. Missing condition records retain their fixed nodes with NA evidence.':state.tier==='replicated_specific'&&state.scope==='IN_total'?'This stricter tier is empty for inhibitory total because its heterogeneous identity fails the atomic-cluster reference-supported policy. Other evidence tiers remain available.':'Hover over a dot, tile, node, ribbon or arrow for native scores and support. Downloads use the current appearance. * marks a retained identity-review population.';
  }
  if(state.view==='table')renderTable();
 }catch(e){displayed=[];$('plots').innerHTML='<p class="empty">'+esc(e.message)+'</p>';$('render-status').textContent=e.message;console.error(e);}
}
function scopeSVG(svg,prefix){
 // Renderer styles only refer to SVG children. Scope selectors to this plot, not the surrounding HTML.
 svg.setAttribute('id',prefix);
 svg.querySelectorAll('style').forEach(el=>{el.textContent=el.textContent.replace(/([^{}]+)\{/g,(m,selectors)=>selectors.split(',').map(s=>'#'+prefix+' '+s.trim()).join(',')+'{');});
 const ids=new Map();svg.querySelectorAll('[id]').forEach(el=>{const old=el.id;ids.set(old,prefix+'-'+old);el.id=ids.get(old);});
 svg.querySelectorAll('*').forEach(el=>{for(const a of [...el.attributes]){let v=a.value;for(const [old,n] of ids)v=v.replaceAll('url(#'+old+')','url(#'+n+')');if(v!==a.value)el.setAttribute(a.name,v);}});
}
function componentRows(population,complex){return data.components[population]&&data.components[population][complex]||[];}
function renderTable(){
 const rows=M.sortRows(selected(tableCondition),$('table-sort').value,$('table-sort').value!=='lr_gmeans');currentTable=rows;const size=50,max=Math.max(0,Math.ceil(rows.length/size)-1);tablePage=Math.min(tablePage,max);
 $('table-note').textContent=format(rows.length)+' filtered '+tableCondition+' records. Native lr_logfc is group-versus-rest specificity. Component DESeq2 cKO/WT effects appear only after selecting Inspect.';
 const columns=[['source','Source'],['target','Target'],['ligand_complex','Ligand complex'],['receptor_complex','Receptor complex'],['lr_logfc','LR logFC¹'],['lr_gmeans','LR gmean'],['magnitude_rank','Magnitude rank'],['specificity_rank','Specificity rank'],['cellphone_pvals','CPDB P²']];
 $('results-table').innerHTML='<table><thead><tr>'+columns.map(([,l])=>'<th>'+l+'</th>').join('')+'<th>WT support</th><th>cKO support</th><th>Details</th></tr></thead><tbody>'+rows.slice(tablePage*size,(tablePage+1)*size).map((r,i)=>'<tr>'+columns.map(([k])=>'<td>'+esc(['source','target','ligand_complex','receptor_complex'].includes(k)?r[k]:format(r[k]))+'</td>').join('')+'<td>'+r.sample_support_WT+'/'+r.sample_eligible_WT+'</td><td>'+r.sample_support_cKO+'/'+r.sample_eligible_cKO+'</td><td><button data-inspect="'+(tablePage*size+i)+'">Inspect'+(r.both_sides_deg_supported?' · both DEG':r.deg_supported?' · DEG':'')+'</button></td></tr>').join('')+'</tbody></table>'+(rows.length?'':'<p class="empty">No records pass these filters. Try another evidence tier or clear the search.</p>');
 $('table-page').textContent='Page '+(tablePage+1)+' / '+(max+1);$('table-prev').disabled=tablePage===0;$('table-next').disabled=tablePage>=max;
}
function renderDetail(r){
 const geneTable=(pop,complex,role)=>{const parts=componentRows(pop,complex);return '<section><h3>'+role+' · '+esc(complex)+'</h3><p>'+esc(pop)+' <span class="tag">'+esc(r[role==='Ligand source'?'source_review':'target_review'])+'</span></p><div class="table-wrap"><table><thead><tr><th>Subunit</th><th>Stable gene ID</th><th>log2(cKO / WT)</th><th>Adjusted P</th><th>Evidence</th></tr></thead><tbody>'+parts.flatMap(p=>{const results=p.stable_id_results&&p.stable_id_results.length?p.stable_id_results:[{}];return results.map(g=>'<tr><td>'+esc(p.symbol)+'</td><td>'+esc(g.gene_id||'NA')+'</td><td>'+format(g.log2FoldChange)+'</td><td>'+format(g.padj)+'</td><td>'+esc(!p.tested_endpoint?'NOT TESTED':p.ambiguous_symbol?'Ambiguous symbol · excluded from primary support':p.primary_deg_support?'DEG · matching endpoint':'No primary DEG support')+'</td></tr>');}).join('')+'</tbody></table></div>'+(parts.length?'':'<p class="small">No matching component evidence was packaged; this is unavailable, not zero.</p>')+'</section>';};
 $('selected-detail').hidden=false;$('selected-detail').innerHTML='<div class="card-head"><h2>'+esc(r.source+' '+r.ligand_complex+' → '+r.target+' '+r.receptor_complex)+'</h2><button id="close-detail">Close</button></div><p class="small">'+esc(r.condition)+' native LR result · group-versus-rest lr_logfc '+format(r.lr_logfc)+' · geometric mean '+format(r.lr_gmeans)+' · magnitude rank '+format(r.magnitude_rank)+'</p><p>WT expression support <b>'+r.sample_support_WT+' / '+r.sample_eligible_WT+'</b> eligible mice; cKO <b>'+r.sample_support_cKO+' / '+r.sample_eligible_cKO+'</b>. Native source / target review flags remain attached. Gene-level DE evidence refers to the original six-WT versus six-cKO fit regardless of the selected scoring dataset.</p><div class="detail-grid">'+geneTable(r.source,r.ligand_complex,'Ligand source')+geneTable(r.target,r.receptor_complex,'Receptor target')+'</div><p class="small">¹ LR logFC is native group-versus-rest specificity. ² CPDB P is cell-label specificity, not a genotype test. Complex support requires each subunit’s expression gate; a DEG in one subunit does not mean every subunit changes. Opposite gene effects do not establish a net activity change.</p><button id="download-component">Download exact component evidence JSON</button>';
 $('download-component').onclick=()=>download(JSON.stringify({scope:r.scope,condition:r.condition,record:r,ligand_components:componentRows(r.source,r.ligand_complex),receptor_components:componentRows(r.target,r.receptor_complex),score_meanings:'lr_logfc=group-versus-rest, component log2FoldChange=cKO/WT; no genotype interaction test'},null,2),'components_'+r.scope+'_'+r.condition+'_'+r.source+'_'+r.ligand_complex+'_'+r.target+'_'+r.receptor_complex+'.json','application/json');
 $('close-detail').onclick=()=>{$('selected-detail').hidden=true;};$('selected-detail').scrollIntoView({behavior:'smooth',block:'start'});
}
function download(content,name,type){const blob=content instanceof Blob?content:new Blob([content],{type:type||'text/plain;charset=utf-8'}),url=URL.createObjectURL(blob),a=document.createElement('a');a.href=url;a.download=name.replace(/[^A-Za-z0-9_.-]/g,'_');document.body.append(a);a.click();a.remove();setTimeout(()=>URL.revokeObjectURL(url),20000);}
function filename(c,suffix,snapshot=state){return 'SATB2_LR_'+snapshot.scope+'_'+c+'_'+snapshot.direction+'_'+snapshot.view+'.'+suffix;}
function tableDownload(extended,condition){const c=condition||(!state.compare?state.condition:tableCondition),rows=M.sortRows(selected(c));download(M.tsv(rows,extended?['scope','condition',...data.columns]:M.SEVEN),filename(c,extended?'extended.tsv':'tsv'),'text/tab-separated-values;charset=utf-8');}
async function exportPlot(i,type){
 const plot=displayed[i];if(!plot)return;const exportState=plot.settings,exportView=exportState.view;const button=document.querySelector('[data-plot="'+i+'"][data-format="'+type+'"]');if(button)button.disabled=true;
 try{
  if(type==='svg')return download(plot.svg,filename(plot.condition,'svg',exportState),'image/svg+xml;charset=utf-8');
  if(type==='json')return download(JSON.stringify({schema:'satb2-communication-plot/v1',settings:exportState,condition:plot.condition,scope:exportState.scope,context:plot.context,records:plot.rows,provenance:data.provenance},null,2),filename(plot.condition,'json',exportState),'application/json');
  const parsed=new DOMParser().parseFromString(plot.svg,'image/svg+xml').documentElement,w=Number(parsed.getAttribute('width')),h=Number(parsed.getAttribute('height'));
  const url=URL.createObjectURL(new Blob([plot.svg],{type:'image/svg+xml;charset=utf-8'}));
  const img=new Image();try{img.src=url;await img.decode();}finally{URL.revokeObjectURL(url);}
  const scale=Math.min(3,6500/Math.max(w,h)),canvas=document.createElement('canvas');canvas.width=Math.ceil(w*scale);canvas.height=Math.ceil(h*scale);const context=canvas.getContext('2d');context.fillStyle='white';context.fillRect(0,0,canvas.width,canvas.height);context.drawImage(img,0,0,canvas.width,canvas.height);
  if(type==='png'){const blob=await new Promise(resolve=>canvas.toBlob(resolve,'image/png'));if(!blob)throw Error('PNG export did not produce a file.');download(blob,filename(plot.condition,'png',exportState));}
  else{if(!window.jspdf)throw Error('The local PDF library was not loaded.');const pdf=new window.jspdf.jsPDF({orientation:w>h?'landscape':'portrait',unit:'pt',format:[w,h],compress:true});pdf.setProperties({title:'SATB2 '+GUIDES[exportView][0]+' · '+plot.condition,subject:'Saved Seurat LIANA candidates; display-only settings. No genotype interaction test.',author:'SATB2 atlas'});pdf.addImage(canvas.toDataURL('image/png'),'PNG',0,0,w,h,undefined,'FAST');pdf.save(filename(plot.condition,'pdf',exportState));}
  $('render-status').textContent=type.toUpperCase()+' downloaded with the current plot appearance.';
 }catch(e){$('render-status').textContent='Download failed: '+e.message;console.error(e);}finally{if(button)button.disabled=false;}
}
function bind(){
 for(const k of Object.keys(DEFAULT)){const el=$(k);if(!el)continue;el.addEventListener(el.type==='search'||el.type==='range'||el.type==='color'?'input':'change',()=>change(k,el.multiple?[...el.selectedOptions].map(o=>o.value):el.type==='checkbox'?el.checked:NUMBERS[k]?Number(el.value):el.value));}
 document.querySelectorAll('[data-view]').forEach(b=>b.addEventListener('click',()=>change('view',b.dataset.view)));
 $('plots').addEventListener('click',e=>{const b=e.target.closest('[data-format]');if(b)exportPlot(Number(b.dataset.plot),b.dataset.format);});
 $('results-table').addEventListener('click',e=>{const b=e.target.closest('[data-inspect]');if(b)renderDetail(currentTable[Number(b.dataset.inspect)]);});
 $('clear-populations').onclick=()=>{state.sources=[];state.targets=[];change('sources',[]);};
 $('node-target').onchange=updateNodeEditor;$('apply-node-color').onclick=()=>{state.nodeColors[$('node-target').value]=$('population-color').value;state.nodeColorMode='palette';sync();save();schedule();};$('reset-node-color').onclick=()=>{delete state.nodeColors[$('node-target').value];save();updateNodeEditor();schedule();};
 $('pair-prev').onclick=()=>{pairPage=Math.max(0,pairPage-1);render();};$('pair-next').onclick=()=>{pairPage++;render();};
 $('table-condition').onchange=()=>{tableCondition=$('table-condition').value;tablePage=0;$('selected-detail').hidden=true;renderTable();};$('table-sort').onchange=()=>{tablePage=0;renderTable();};$('table-prev').onclick=()=>{tablePage--;renderTable();};$('table-next').onclick=()=>{tablePage++;renderTable();};
 $('download-values').onclick=()=>selectedConditions().forEach(c=>tableDownload(false,c));$('download-extended').onclick=()=>selectedConditions().forEach(c=>tableDownload(true,c));$('table-download').onclick=()=>tableDownload(false,tableCondition);$('table-extended').onclick=()=>tableDownload(true,tableCondition);
 const settings=()=>download(JSON.stringify({schema:'satb2-communication-settings/v1',settings:state},null,2),'SATB2_LR_settings.json','application/json');$('export-settings').onclick=settings;$('download-settings-copy').onclick=settings;
 $('import-settings').onclick=()=>$('settings-file').click();$('settings-file').onchange=async()=>{try{const file=$('settings-file').files[0];if(!file)return;const imported=JSON.parse(await file.text());if(imported.schema!=='satb2-communication-settings/v1')throw Error('Choose a SATB2 communication settings file.');state=valid(imported.settings);tableCondition=state.condition;pairPage=0;tablePage=0;filteredCache={};populationOptions();save();render();}catch(e){$('render-status').textContent=e.message;}finally{$('settings-file').value='';}};
 $('reset-settings').onclick=()=>{state=valid(DEFAULT);tableCondition='pooled';tablePage=0;pairPage=0;populationOptions();save();render();};
 $('show-lead').onclick=()=>{state=valid({...DEFAULT,tier:'replicated_specific',bothDEG:true,search:'Cntn3',view:'table'});tableCondition='pooled';tablePage=0;populationOptions();save();render();};$('show-strong').onclick=()=>{state=valid({...DEFAULT,tier:'replicated_specific'});tableCondition='pooled';populationOptions();save();render();};
}
async function start(){
 try{
  await new Promise(r=>setTimeout(r,30));const bytes=Uint8Array.from(atob(window.CCC_DATA_GZIP),c=>c.charCodeAt(0));data=M.decode(JSON.parse(new TextDecoder().decode(window.fflate.gunzipSync(bytes))));
  // Scope/condition are transport context, never inferred from gene effects.
  for(const [key,dataset] of Object.entries(data.datasets)){const [scope,condition]=key.split('|');for(const row of dataset.rows){row.scope=scope;row.condition=condition;}}
  let saved=null;try{saved=JSON.parse(localStorage.getItem(STORE)||'null');}catch(_){}state=valid(saved);tableCondition=state.condition;
  populationOptions();bind();$('workspace').hidden=false;$('load-status').hidden=true;$('provenance').textContent='Interactive payload:64,936 native records across six fits; matching component DESeq2 evidence is embedded.';
  const genes=new Set();for(const d of Object.values(data.datasets))for(const r of d.rows){genes.add(r.ligand_complex);genes.add(r.receptor_complex);} $('gene-list').innerHTML=[...genes].sort().map(g=>'<option value="'+esc(g)+'"></option>').join('');$('search').setAttribute('list','gene-list');
  // Exposed read-only access is useful for reproducible browser acceptance checks.
  window.CCCAtlas={get state(){return JSON.parse(JSON.stringify(state));},get data(){return data;},get plots(){return displayed;},filter:M.filter,render,exportPlot};render();
 }catch(e){$('load-status').className='error';$('load-status').textContent='Could not load the embedded communication data: '+e.message;console.error(e);}
}
start();
})();
