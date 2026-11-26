'use strict';
const referenceNames = ['SynDive','Synaptome.db','G2Cdb','Connectome-seq'];
const resourceForView = {syndive:'SynDive',synaptome:'Synaptome.db',g2c:'G2Cdb',connectome:'Connectome-seq'};
const isPairView = name => ['pairs','overlap','circos','binding'].includes(name);
const pairColumns = [
  ['binding_verdict','Reviewed binding evidence'],['corrected_pair','Reviewed molecular partners'],
  ['CA1_compartment','CA1 protein compartment'],['IN_compartment','IN protein compartment'],
  ['synaptic_direction','Predicted synaptic direction'],['synaptic_effect','Predicted synaptic effect'],
  ['binding_assay','Binding assay'],['binding_context','Experimental context'],['molecular_form','Fragment / isoform'],
  ['interpretation','Biological interpretation'],['contact_validation','P60 contact evidence'],['binding_sources','Primary sources'],
  ['CA1_sublayer','CA1 sublayer'],['CA1_gene','CA1 gene'],['CA1_role','Archived CA1 resource role'],['CA1_change','cKO RNA'],
  ['IN_cluster','Inhibitory cluster'],['IN_gene','Interneuron gene'],['IN_role','Archived interneuron resource role'],['IN_change','cKO RNA'],
  ['exact_reference_set','Exact resource matches'],['reference_count','Resources'],['PPI_pmids','Pair assay PMIDs'],
  ['source','Archived LIANA source'],['ligand','Archived ligand'],['target','Archived LIANA target'],['receptor','Archived receptor'],
  ['CA1_log2FC','CA1 log₂FC'],['CA1_FDR','CA1 FDR'],['IN_log2FC','IN log₂FC'],['IN_FDR','IN FDR'],
  ...referenceNames.map(n=>[n,n]),['FASS_common_strata','Shared FASS strata'],['G2C_shared_lists','Shared G2C lists'],
  ['PPI_types','Indexed PPI types'],['connectome_marker_ends','Connectome-seq marker ends'],
  ['magnitude_rank','LIANA magnitude rank'],['specificity_rank','LIANA specificity rank'],
  ['WT_mouse_support','WT mice'],['cKO_mouse_support','cKO mice']
];
function selectedPairRows(rows) {
  const ref=value('resource'), size=value('matchCount');
  return rows.filter(p=>(!ref||p[ref])&&(!size||(size==='2+'?p.reference_count>=2:p.reference_count===Number(size)))&&
    (!$('assay').checked||Boolean(p.PPI_pmids))&&(!value('bindingStatus')||p.binding_status===value('bindingStatus')));
}
function ribbonKey(p) { return [p.source,p.ligand,p.target,p.receptor].join('|'); }
function isolateRibbon(rows) {
  const select=$('oneRibbon'), chosen=select.value;
  select.replaceChildren(new Option('All matching pairs',''));
  for(const p of rows) select.add(new Option(p.corrected_pair+' · '+pop(p.corrected_source)+' / '+pop(p.corrected_target)+' · '+p.binding_verdict,ribbonKey(p)));
  if(rows.some(p=>ribbonKey(p)===chosen))select.value=chosen;
  return select.value?rows.filter(p=>ribbonKey(p)===select.value):rows;
}
function pairCounts(rows) {
  const ordered=new Set(rows.map(p=>p.ligand+' → '+p.receptor));
  const unordered=new Set(rows.map(p=>[p.ligand,p.receptor].sort().join(' / ')));
  return {rows:rows.length,ordered:ordered.size,unordered:unordered.size,genes:[...ordered].sort().join(', ')};
}
function renderOverlapTable(rows) {
  const counts=pairCounts(rows),groups=new Map();
  for(const p of rows){
    if(!groups.has(p.exact_reference_set))groups.set(p.exact_reference_set,[]);
    groups.get(p.exact_reference_set).push(p);
  }
  $('countNote').textContent=counts.rows+' population-specific rows · '+counts.ordered+' ordered ligand → receptor pairs · '+counts.unordered+
    ' pairs after also collapsing reciprocal roles. Each row below names exactly the references with a both-end match in the current selection.';
  const summary=[...groups].map(([name,items])=>({name,...pairCounts(items)})).sort((a,b)=>b.rows-a.rows||a.name.localeCompare(b.name));
  table('referenceTable',summary,[['name','Exact reference set'],['rows','Population rows'],['ordered','Ordered gene pairs'],['unordered','Unordered pairs'],['genes','Ligand → receptor']]);
}
const circosPopulations=[
  {id:'CA1_Superficial',label:'Superficial CA1 · C2 + C3',from:-169,to:-98,left:true,top:true},
  {id:'CA1_Deep',label:'Deep CA1 · C1 + C5',from:98,to:169,left:true,top:false},
  {id:'IN_C1',label:'IN C1 · PV basket-like',from:-82,to:-11,left:false,top:true},
  {id:'IN_C6',label:'IN C6 · SST',from:11,to:82,left:false,top:false}
];
function circosFigure(rows) {
  const width=960,height=633,cx=480,cy=319,r=190,outer=r+15;
  const labelSize=Number(value('fontsize')),alpha=Number(value('ribbonOpacity'))/100;
  const caColor=value('caColor'),inColor=value('inColor'),labels=$('geneLabels').checked;
  const point=(a,radius)=>[cx+radius*Math.cos(a),cy+radius*Math.sin(a)];
  const xy=p=>p.map(n=>n.toFixed(2)).join(',');
  const key=(i,role)=>i+':'+role;
  const angles=new Map(),arcs=[],names=[];
  const text=(x,y,str,extra='')=>'<text x="'+x+'" y="'+y+'" '+extra+'>'+esc(str)+'</text>';
  for(const group of circosPopulations){
    const ends=[];
    rows.forEach((p,i)=>{
      if(p.source===group.id)ends.push({gene:p.ligand,i,role:'ligand',lfc:p.ligand_log2FC});
      if(p.target===group.id)ends.push({gene:p.receptor,i,role:'receptor',lfc:p.receptor_log2FC});
    });
    ends.sort((a,b)=>a.gene.localeCompare(b.gene)||a.role.localeCompare(b.role)||a.i-b.i);
    const from=group.from*Math.PI/180,to=group.to*Math.PI/180,span=to-from;
    const side=group.left?caColor:inColor,groupX=group.left?223:737,anchor=group.left?'end':'start';
    arcs.push('<path d="M'+xy(point(from,outer))+' A'+outer+' '+outer+' 0 0 1 '+xy(point(to,outer))+
      '" stroke="'+side+'" opacity="'+(ends.length?1:.18)+'" stroke-width="8" stroke-linecap="round" fill="none"/>');
    names.push(text(groupX,group.top?97:369,group.label,'font-size="17" font-weight="700" text-anchor="'+anchor+'"'));
    const groups=new Map(),step=span/Math.max(1,ends.length);
    ends.forEach((e,j)=>{
      const angle=from+(j+.5)*step;
      const half=Math.min(.061,step*.44)*Number(value('ribbonWidth'))/100;
      angles.set(key(e.i,e.role),{angle,half});
      if(!groups.has(e.gene))groups.set(e.gene,{gene:e.gene,angles:[],lfc:e.lfc});
      groups.get(e.gene).angles.push(angle);
    });
    if(!labels)continue;
    const ports=[...groups.values()].map(g=>{
      const mid=g.angles.reduce((a,b)=>a+b,0)/g.angles.length;
      return {...g,pt:point(mid,outer+5)};
    }).sort((a,b)=>a.pt[1]-b.pt[1]);
    ports.forEach((g,i)=>{
      const minY=group.top?119:393,maxY=group.top?256:527;
      const y=ports.length===1?(minY+maxY)/2:minY+(maxY-minY)*i/(ports.length-1);
      const x=groupX+(group.left?9:-9), elbow=group.left?265:695;
      names.push('<path d="M'+xy(g.pt)+' Q'+elbow+','+g.pt[1].toFixed(2)+' '+x+','+y.toFixed(2)+'" fill="none" stroke="#a8b8ba" stroke-width=".9"/>',
        text(groupX,y+4,g.gene+(g.lfc>0?' ↑':' ↓'),'text-anchor="'+anchor+'" font-size="'+labelSize+'" fill="#233941"'));
    });
  }
  const ribbons=rows.map((p,i)=>{
    const s=angles.get(key(i,'ligand')),t=angles.get(key(i,'receptor'));
    const a=s.angle-s.half,b=s.angle+s.half,c=t.angle-t.half,d=t.angle+t.half;
    const stroke=reviewColor(p,caColor,inColor);
    const arrow = connectionMode(p);
    const targetTip = arrow==='forward'||arrow==='both';
    const sourceTip = arrow==='reverse'||arrow==='both';
    const path='M'+xy(point(a,r))+' A'+r+' '+r+' 0 0 1 '+xy(point(b,r))+
      ' C'+xy(point(b,r*.36))+' '+xy(point(c,r*.36))+' '+xy(point(c,r))+
      (targetTip?' L'+xy(point(t.angle,r+12)):'')+' L'+xy(point(d,r))+
      ' C'+xy(point(d,r*.36))+' '+xy(point(a,r*.36))+' '+xy(point(a,r))+(sourceTip?' L'+xy(point(s.angle,r+12)):'')+' Z';
    const desc=p.corrected_pair+' · '+p.binding_verdict+' · '+p.synaptic_direction+' · '+p.CA1_compartment+' / '+p.IN_compartment;
    return '<path class="lr-ribbon" tabindex="0" role="button" data-key="'+esc(ribbonKey(p))+'" aria-label="Isolate '+esc(desc)+'" data-source="'+esc(p.source)+'" data-target="'+esc(p.target)+'" data-ligand="'+esc(p.ligand)+'" data-receptor="'+esc(p.receptor)+'" d="'+path+'" fill="'+stroke+'" fill-opacity="'+(arrow==='none'?alpha*.38:alpha)+'" stroke="'+stroke+'" stroke-width="1"'+(p.binding_status==='direct'?'':' stroke-dasharray="5 4"')+'><title>'+esc(desc)+'</title></path>';
  });
  return '<svg xmlns="http://www.w3.org/2000/svg" width="'+width+'" height="'+height+'" viewBox="0 0 '+width+' '+height+'" role="img" aria-label="Directed CA1 and interneuron ligand–receptor ribbons"><rect width="100%" height="100%" fill="white"/><g font-family="Arial, sans-serif" fill="#213b43">'+
    text(26,32,'CA1–interneuron molecular partners','font-size="22" font-weight="700"')+
    text(26,57,rows.length+' archived rows  ·  '+arrowCaption(),'font-size="14"')+
    arcs.join('')+ribbons.join('')+names.join('')+(!rows.length?text(cx,cy,'No candidates match the selected filters.','text-anchor="middle" font-size="16"'):'')+
    text(30,570,'Gray = unverified / indirect binding; amber = screen. Pale tip-free ribbons have unresolved synaptic orientation.','font-size="12"')+
    text(30,594,arrowMeaning(),'font-size="12"')+
    text(30,618,'Gene ↑ / ↓ = cKO RNA versus WT. Ribbon width is adjustable. The reciprocal NRXN1–NLGN1 rows describe one predicted orientation.','font-size="11"')+'</g></svg>';
}
