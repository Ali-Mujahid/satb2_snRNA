'use strict';
const $ = id => document.getElementById(id);
const sources = {
  connectome: ['Connectome-seq marker overlap', 'Compare the 30 published Purkinje-cell markers with your CA1 and interneuron DEGs. Color shows the local RNA effect for significant genes. The default view includes genes with at least one overlap.'],
  syndive: ['CA1 stratum proteomics', 'Inspect each ligand or receptor across CA1 strata. The FASS fraction contains sorted synaptosomes. Values are median stratum log₂ fold changes across unambiguous protein groups. A star marks support from at least one group passing the authors’ significance criterion.'],
  synaptome: ['Synaptic compartment evidence', 'Count the publications reporting each protein in a synaptic fraction. Select a species and region to see which experiments contribute. The table includes the supporting PubMed IDs.'],
  g2c: ['Postsynaptic proteins and complexes', 'Trace ligand and receptor gene products into published cortical PSD and receptor-complex tables. Mouse and human calls are shown separately. Each cell has a source status in the table below.'],
  overlap: ['Unique and shared pairs across references', 'Each row is a directed, population-specific LIANA candidate. A filled cell means both components meet that reference’s coverage rule. The table below groups identical combinations of resources.'],
  circos: ['Directed CA1–interneuron Circos', 'Follow each ribbon from the ligand gene to the arrow tip at the receptor gene. Select a molecular direction, reference or gene to isolate its connections. Arrow direction follows the LIANA annotation.'],
  pairs: ['Directed ligand–receptor candidates', 'Read the cKO RNA change at each end of the interaction. Each row identifies the ligand-producing and receptor-expressing populations. The table includes native pooled LIANA ranks, mouse support and the reference evidence.']
};
const findings = {
  connectome: 'Deep CA1 shares six of 25 testable panel genes (one up, five down, BH FDR 0.023). Superficial CA1 shares three of 25 (BH FDR 1.00). The enriched C5 result rests on Arhgap31, its sole DEG.',
  syndive: 'CNTN3 and PTPRG each have unambiguous sorted-synaptosome estimates across all four CA1 strata. In your mice, Cntn3 falls in inhibitory C6 and Ptprg rises in both pyramidal sublayers.',
  synaptome: 'Both proteins have mouse hippocampal study records for 12 of 14 directed candidates. The separate interaction lookup returns three publications for NLGN1–NRXN1.',
  g2c: 'PTN and PTPRZ1 both appear in the mouse cortical PSD preparation. Your RNA data place lower Ptn in both CA1 sublayers and higher Ptprz1 in inhibitory C1. Their physical site of action in CA1 remains open.',
  overlap: 'Only Synaptome.db has exclusive two-end matches in this set. SynDive and G2Cdb each overlap Synaptome.db on different pairs. No candidate matches three or four references.',
  circos: 'The left half shows pyramidal genes by sublayer. The right half shows inhibitory C1 (PV basket-like) and parent C6 (SST). Each ribbon represents one candidate in one layer and cluster. Click a ribbon to isolate the pair, then choose “All matching pairs” to return.',
  pairs: 'Fourteen directed rows connect the two CA1 sublayers with C1 or C6. NLGN1–NRXN1 has reciprocal resource entries and three indexed interaction publications. These entries describe one molecular pairing.'
};
const downloads = {
  'CA1_IN_directed_pair_roles.tsv':'CA1 and interneuron roles, effects and references',
  'exact_reference_overlap.tsv':'All 16 exact reference combinations',
  'lr_pair_evidence.tsv':'Directed LR evidence',
  'all_DEGs_reference_lookup.tsv.gz':'Every DEG with four reference lookups',
  'Connectome_seq_overlap_tests.tsv':'Connectome-seq enrichment tests',
  'Connectome_seq_30_genes_by_endpoint.tsv':'All 30 markers × 22 populations',
  'SynDive_gene_strata.tsv.gz':'SynDive gene and stratum estimates',
  'SynDive_protein_groups.tsv.gz':'SynDive protein groups and q-values',
  'Synaptome_study_hits.tsv.gz':'Synaptome.db studies and PMIDs',
  'Synaptome_pair_PPI.tsv':'Synaptome.db interaction assays',
  'G2C_memberships.tsv.gz':'G2C source-table memberships',
  'provenance.json':'Input checksums'
};
const inputs = [...document.querySelectorAll('#controls input,#controls select,#appearance input,#appearance select')];
const value = id => $(id).value;
const pop = s => ({CA1_Superficial:'CA1 superficial',CA1_Deep:'CA1 deep'}[s] || s.replace('IN_','IN '));
const fmt = n => n == null ? '—' : n===0 ? '0' : Math.abs(n)<0.001 ? n.toExponential(2) : Number(n.toPrecision(3)).toString();
const plus = n => (n>0 ? '+' : '') + n.toFixed(2);
const esc = s => String(s??'').replace(/[&<>"']/g,c=>({'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;',"'":'&#39;'}[c]));
const fields = ['Presynaptic','Postsynaptic','Synaptosome','Synaptic_Vesicle'];
const defaults = {};
let data, view, figure, plottedRows, columns;
function option(id, val, title) { $(id).add(new Option(title,val)); }
function table(id,rows,cols) {
  if(!rows.length){$(id).textContent='No matching records.';return;}
  const row = (r,tag) => '<tr>'+r.map(x=>'<'+tag+'>'+esc(x)+'</'+tag+'>').join('')+'</tr>';
  $(id).innerHTML='<table><thead>'+row(cols.map(c=>c[1]),'th')+'</thead><tbody>'+rows.map(r=>row(cols.map(([k])=>r[k]??'—'),'td')).join('')+'</tbody></table>';
}
function filteredPairs(){
  return (isPairView(view)?data.pairRoles:data.pairs).filter(p=>(!value('layer')||[p.source,p.target].includes(value('layer'))) &&
    (!value('partner')||[p.source,p.target].includes(value('partner'))) &&
    (!value('direction')||p.source.startsWith(value('direction'))));
}
function selectedPopulations(){
  const ep=value('endpoint'),hit=new Set(data.stats.filter(s=>s.panel_DEG).map(s=>s.endpoint));
  return Object.keys(data.status).filter(k=> ep==='all'||ep===k|| (ep==='layers'&&k.startsWith('CA1')) ||
    (ep==='in'&&k.startsWith('IN')) || (ep==='hits'&&hit.has(k)));
}
function matrix() {
  const pairRows=filteredPairs(), query=value('filter').trim().toLowerCase();
  const match=s=>s.toLowerCase().includes(query);
  const genes=[...new Set(pairRows.flatMap(p=>[p.ligand,p.receptor]))].sort().filter(match);
  let cells=[], rows, cols, legend, scale=1, tableRows, headings;
  function cell(row,col,n,display,details={},type='value'){cells.push({row,col,n,display,type,...details});}
  if(view==='connectome'){
    const eps=selectedPopulations(), subset=data.con.filter(x=>eps.includes(x.endpoint));
    rows=data.markers.filter(g=>match(g)&&(value('panel')==='all'||subset.some(x=>x.gene===g&&x.status==='DEG')));
    cols=eps.map(pop);scale=Math.max(.01,...subset.filter(x=>x.status==='DEG'&&rows.includes(x.gene)).map(x=>Math.abs(x.log2FC)));
    for(const x of subset.filter(x=>rows.includes(x.gene))){
      const sig=x.status==='DEG',tested=x.status==='tested';
      cell(x.gene,pop(x.endpoint),sig?x.log2FC:null,sig?plus(x.log2FC):tested?'·':'—',{status:x.status,log2FC:x.log2FC,FDR:x.FDR},sig?'value':tested?'tested':'missing');
    }
    legend='Local DEG log₂(cKO/WT)    ·  tested, FDR ≥ 0.05     —  unavailable';
    headings=[['row','Gene'],['col','Population'],['status','DESeq2 status'],['log2FC','cKO/WT log₂FC'],['FDR','Adjusted P']];
    const stats=data.stats.filter(x=>eps.includes(x.endpoint)).map(x=>({...x,pop:pop(x.endpoint),p:fmt(x.p),BH_FDR:fmt(x.BH_FDR)}));
    if(eps.includes('IN_C18')) stats.push({pop:'IN C18',genes:'No eligible model'});
    table('stats',stats,[['pop','Population'],['tested_unique_symbols','Tested genes'],['DEG_unique_symbols','DEG genes'],['panel_tested','Markers tested'],['panel_DEG','Overlap'],['up','Up'],['down','Down'],['p','Enrichment P'],['BH_FDR','BH FDR'],['genes','Overlapping genes']]);
  } else if(view==='syndive'){
    rows=genes;cols=data.strata;const fr=value('fraction');
    for(const g of rows) for(const a of cols){
      const x=data.syn.find(t=>t.gene===g&&t.stratum===a&&t.fraction===fr);
      cell(g,a,x?.log2FC??null,x?plus(x.log2FC)+(x.significant_groups?'*':''):'—',{fraction:fr,protein_groups:x?.groups??0,significant_groups:x?.significant_groups??0},x?'value':'missing');
    }
    scale=Math.max(.01,...cells.map(x=>Math.abs(x.n||0)));
    legend=(fr==='FASS'?'FASS log₂FC vs other strata':'Tissue log₂FC vs four-stratum mean')+'    *  author-significant group     —  no unambiguous estimate';
    headings=[['row','Gene'],['col','Stratum'],['fraction','Fraction'],['n','Median log₂FC'],['protein_groups','Groups'],['significant_groups','Significant groups']];
  } else if(view==='synaptome'){
    rows=genes;cols=fields.map(c=>c.replace('_',' '));const region=value('region'),species=value('species');
    const eligible=data.db.filter(x=>(!region||x.region===region)&&(!species||String(x.species_taxid)===species));
    for(const g of rows) for(const [i,c] of fields.entries()){
      const found=eligible.filter(x=>x.gene===g&&x.compartment===c);
      const pmids=[...new Set(found.map(x=>x.pmid))].sort();
      cell(g,cols[i],pmids.length,String(pmids.length),{region:region||'All regions',species:species||'All species',PMIDs:pmids.join(', ')});
    }
    scale=Math.max(1,...cells.map(x=>x.n));
    legend=($('species').selectedOptions[0].text)+' · '+(region||'All regions')+'    Color = unique PMIDs    0 = no indexed hit';
    headings=[['row','Gene'],['col','Fraction'],['n','Unique PMIDs'],['region','Region'],['species','Species (taxonomy ID)'],['PMIDs','PubMed IDs']];
  } else if(view==='g2c'){
    rows=genes;cols=data.collections;
    for(const x of data.g2.filter(x=>rows.includes(x.gene))){
      const n={listed:1,ambiguous:.45,unmarked:0}[x.status]??null;
      cell(x.gene,x.collection,n,({listed:'+',ambiguous:'?',unmarked:'·','no eligible row':'—'}[x.status]),{status:x.status},n===null?'missing':'value');
    }
    legend='+  listed     ?  ambiguous molecular assignment     ·  eligible, unmarked     —  no eligible row';
    headings=[['row','Gene'],['col','Reference collection'],['status','Source status']];
  } else {
    tableRows=selectedPairRows(pairRows).filter(x=>match(x.ligand)||match(x.receptor));
    if(view==='circos')tableRows=isolateRibbon(tableRows);
    rows=[];cols=view==='overlap'?referenceNames:['Ligand RNA','Receptor RNA'];
    for(const p of tableRows){
      const r=p.ligand+' → '+p.receptor+'   ·   '+pop(p.source)+' → '+pop(p.target);
      rows.push(r);
      if(view==='overlap'){
        for(const name of referenceNames)cell(r,name,p[name]?1:0,p[name]?'+':'·');
        continue;
      }
      cell(r,cols[0],p.ligand_log2FC,plus(p.ligand_log2FC));cell(r,cols[1],p.receptor_log2FC,plus(p.receptor_log2FC));
    }
    scale=Math.max(.01,...cells.map(x=>Math.abs(x.n)));
    legend='cKO/WT RNA log₂FC at the source and target    |    Pooled LIANA candidates';
    if(view==='overlap'){scale=1;legend='+  both components meet the reference criterion     ·  no both-end match';}
    headings=pairColumns;
  }
  return {rows,cols,cells,legend,scale,headings,tableRows:tableRows||cells};
}
function blend(hex,t){
  const ch=[1,3,5].map(p=>parseInt(hex.slice(p,p+2),16));
  return 'rgb('+ch.map(n=>Math.round(247+(n-247)*t)).join(',')+')';
}
function svg(m){
  const fs=Number(value('fontsize')),dy=Number(value('spacing')),dot=value('style')==='dot';
  const wide=['pairs','connectome'].includes(view);
  const left=Math.max(105,Math.max(0,...m.rows.map(s=>s.length))*.60*fs+18);
  const dx=view==='g2c'?109:view==='overlap'?136:wide?116:115;
  const tilt=view==='g2c',top=tilt?176:105;
  const W=Math.ceil(Math.max(740,left+dx*m.cols.length+30)), H=top+dy*m.rows.length+32;
  const txt=(x,y,t,attrs='')=>'<text x="'+x+'" y="'+y+'" '+attrs+'>'+esc(t)+'</text>';
  let parts=['<svg xmlns="http://www.w3.org/2000/svg" width="'+W+'" height="'+H+'" viewBox="0 0 '+W+' '+H+'" role="img" aria-label="'+esc(sources[view][0])+'"><rect width="100%" height="100%" fill="white"/><g font-family="Arial, sans-serif" fill="#1c333d" font-size="'+fs+'">',
    txt(13,27,sources[view][0],'font-size="19" font-weight="700"'),txt(13,51,m.legend,'font-size="13"')];
  const signed=['syndive','pairs','connectome'].includes(view);
  if(signed) {
    const high=value('high'),low=value('low'); let x=13;
    for(const [c,s] of [[low,'−'+fmt(m.scale)],['#f7f7f7','0'],[high,'+'+fmt(m.scale)]]){
      parts.push('<rect x="'+x+'" y="62" width="12" height="10" fill="'+c+'" stroke="#acbab8"/>',txt(x+17,72,s,'font-size="11"'));x+=78;
    }
  }
  m.cols.forEach((col,j)=>{const x=left+(j+.5)*dx,y=top-10;
    parts.push(txt(x,y,col,tilt?'transform="rotate(-36 '+x+' '+y+')" text-anchor="start" font-size="'+fs+'"':'text-anchor="middle" font-weight="600"'));});
  m.rows.forEach((row,i)=>{parts.push(txt(left-10,top+(i+.5)*dy+4,row,'text-anchor="end"'));});
  for(const c of m.cells){
    const x=left+dx*m.cols.indexOf(c.col),y=top+dy*m.rows.indexOf(c.row);
    const t=c.n===null?0:Math.min(1,Math.abs(c.n)/m.scale);
    const color=c.type==='missing'?'#e1e6e8':c.type==='tested'?'#f2f5f4':blend(c.n<0?value('low'):value('high'),t);
    parts.push('<g><title>'+esc(c.row+' · '+c.col+' · '+(c.status??c.display)+(typeof c.n==='number'?' ('+fmt(c.n)+')':''))+'</title>');
    if(dot){
      if(c.n!==null)parts.push('<circle cx="'+(x+dx*.5)+'" cy="'+(y+dy*.5)+'" r="'+Math.min(dy*.46,Number(value('dotsize'))*(.4+.6*Math.sqrt(t)))+'" fill="'+color+'" stroke="#8b9da2" stroke-width=".5"/>');
      if($('values').checked||c.n===null)parts.push(txt(x+dx*.5+($('values').checked&&c.n!==null?Number(value('dotsize'))+2:0),y+dy*.5+4,c.display,c.n===null?'text-anchor="middle"':'font-size="11"'));
    }else{
      parts.push('<rect x="'+(x+1)+'" y="'+(y+1)+'" width="'+(dx-2)+'" height="'+(dy-2)+'" fill="'+color+'"/>');
      if($('values').checked){
        const rgb=(color.match(/\d+/g)||[]).map(Number);
        const lightText=color.startsWith('rgb')&&(0.2126*rgb[0]+0.7152*rgb[1]+0.0722*rgb[2]<137);
        parts.push(txt(x+dx*.5,y+dy*.5+4,c.display,'text-anchor="middle" fill="'+(lightText?'white':'#263c46')+'"'));
      }
    }
    parts.push('</g>');
  }
  if(!m.rows.length)parts.push(txt(13,H-9,'No matching genes or pairs. Change the filters or display all panel genes.'));
  return parts.join('')+'</g></svg>';
}
function draw() {
  document.querySelectorAll('[data-for]').forEach(el=>{el.hidden=!(el.dataset.for===view||(el.dataset.for==='other'&&view!=='connectome')||(el.dataset.for==='pairViews'&&isPairView(view))||(el.dataset.for==='matrix'&&view!=='circos'));});
  document.querySelectorAll('[data-view]').forEach(el=>el.setAttribute('aria-pressed',String(el.dataset.view===view)));
  $('statsBox').hidden=view!=='connectome'; $('title').textContent=sources[view][0];$('guide').textContent=sources[view][1];$('message').textContent='';
  $('finding').textContent=findings[view];
  const m=matrix();
  $('toPairs').hidden=isPairView(view);
  if(!isPairView(view))$('toPairs').textContent='View pairs matching '+resourceForView[view];
  $('referenceBox').hidden=!isPairView(view);
  if(isPairView(view))renderOverlapTable(m.tableRows);
  figure=view==='circos'?circosFigure(m.tableRows):svg(m);$('plot').innerHTML=figure;plottedRows=m.tableRows;columns=m.headings;
  table('table',plottedRows,columns);$('tableTitle').textContent=isPairView(view)?'CA1 and interneuron genes and roles':'Plotted values';
}
function saveBlob(blob,name){
  const a=document.createElement('a');a.href=URL.createObjectURL(blob);a.download=name;a.click();setTimeout(()=>URL.revokeObjectURL(a.href),3000);
}
function downloadRoot(){
  if(location.protocol==='file:')return new URL('./results/',location.href);
  if(location.pathname.includes('/site/synaptic_evidence/'))return new URL('../../analysis/synaptic_evidence_20261125/results/',location.href);
  return 'https://raw.githubusercontent.com/Ali-Mujahid/satb2_snRNA/ATLAS_RELEASE_TAG/analysis/synaptic_evidence_20261125/results/';
}
function settings(){return Object.fromEntries(inputs.map(e=>[e.id,e.type==='checkbox'?e.checked:e.value]));}
async function save(format){
  const base='SATB2_'+view,blob=new Blob([figure],{type:'image/svg+xml;charset=utf-8'});
  if(format==='svg')return saveBlob(blob,base+'.svg');
  if(format==='tsv'){
    const keys=columns.map(c=>c[0]);const clean=s=>String(s??'').replace(/[\t\n\r]/g,' ');
    const lines=[keys.join('\t'),...plottedRows.map(r=>keys.map(k=>clean(r[k])).join('\t'))];
    return saveBlob(new Blob([lines.join('\n')+'\n'],{type:'text/tab-separated-values'}),base+'.tsv');
  }
  if(format==='json')return saveBlob(new Blob([JSON.stringify({view,settings:settings(),data:plottedRows},null,2)],{type:'application/json'}),base+'.json');
  const im=new Image(),url=URL.createObjectURL(blob);
  try {
    im.src=url;await im.decode();const canvas=document.createElement('canvas');canvas.width=im.width*2;canvas.height=im.height*2;
    canvas.getContext('2d').drawImage(im,0,0,canvas.width,canvas.height);
    if(format==='png'){const b=await new Promise(res=>canvas.toBlob(res,'image/png'));if(!b)throw Error('PNG conversion failed');saveBlob(b,base+'.png');}
    if(format==='pdf'){const pdf=new window.jspdf.jsPDF({unit:'pt',orientation:im.width>=im.height?'landscape':'portrait',format:[im.width,im.height],compress:true});
      pdf.addImage(canvas.toDataURL('image/png'),'PNG',0,0,im.width,im.height,undefined,'FAST');pdf.save(base+'.pdf');}
  }finally{URL.revokeObjectURL(url);}
}
async function init(){
  data=JSON.parse($('sourceData').textContent);
  $('toPairs').onclick=()=>{
    const ref=resourceForView[view];
    $('resource').value=ref;$('matchCount').value='';$('assay').checked=false;
    $('filter').value='';$('layer').value='';$('partner').value='';$('direction').value='';
    view='overlap';history.replaceState(null,'','#overlap');draw();
  };
  for(const ep of Object.keys(data.status))option('endpoint',ep,pop(ep));
  for(const area of [...new Set(data.db.map(r=>r.region))].sort())if(area!=='Hippocampus')option('region',area,area);
  Object.assign(defaults,settings());view=Object.hasOwn(sources,location.hash.slice(1))?location.hash.slice(1):'connectome';
  for(const e of inputs)e.addEventListener('input',draw);
  function pickRibbon(event){
    const ribbon=event.target.closest('.lr-ribbon');
    if(!ribbon || (event.type==='keydown' && !['Enter',' '].includes(event.key)))return;
    event.preventDefault();
    $('oneRibbon').value=ribbon.dataset.key;
    draw();
    $('oneRibbon').focus();
  }
  $('plot').addEventListener('click',pickRibbon);$('plot').addEventListener('keydown',pickRibbon);
  document.querySelectorAll('[data-view]').forEach(e=>e.addEventListener('click',()=>{view=e.dataset.view;history.replaceState(null,'','#'+view);draw();}));
  $('reset').onclick=()=>{for(const el of inputs){if(el.type==='checkbox')el.checked=defaults[el.id];else el.value=defaults[el.id];}draw();};
  document.querySelectorAll('[data-save]').forEach(e=>e.addEventListener('click',()=>save(e.dataset.save).catch(err=>{$('message').textContent='Export failed. '+err.message;})));
  for(const [url,title] of Object.entries(downloads)){const a=document.createElement('a');a.href=new URL(url,downloadRoot()).href;a.download=url;a.textContent=title;$('files').append(a);}
  draw();
}
init().catch(err=>{$('plot').textContent=err.message;});
