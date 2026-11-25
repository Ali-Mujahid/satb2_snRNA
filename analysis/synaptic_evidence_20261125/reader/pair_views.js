'use strict';
const referenceNames = ['SynDive','Synaptome.db','G2Cdb','Connectome-seq'];
const resourceForView = {syndive:'SynDive',synaptome:'Synaptome.db',g2c:'G2Cdb',connectome:'Connectome-seq'};
const isPairView = name => ['pairs','overlap','circos'].includes(name);
const pairColumns = [
  ['CA1_sublayer','CA1 sublayer'],['CA1_gene','CA1 gene'],['CA1_role','CA1 role'],['CA1_change','cKO RNA'],
  ['IN_cluster','Inhibitory cluster'],['IN_gene','Interneuron gene'],['IN_role','Interneuron role'],['IN_change','cKO RNA'],
  ['exact_reference_set','Exact resource matches'],['reference_count','Resources'],['PPI_pmids','Pair assay PMIDs'],
  ['source','Ligand source'],['ligand','Ligand'],['target','Receptor target'],['receptor','Receptor'],
  ['CA1_log2FC','CA1 log₂FC'],['CA1_FDR','CA1 FDR'],['IN_log2FC','IN log₂FC'],['IN_FDR','IN FDR'],
  ...referenceNames.map(n=>[n,n]),['FASS_common_strata','Shared FASS strata'],['G2C_shared_lists','Shared G2C lists'],
  ['PPI_types','Indexed PPI types'],['connectome_marker_ends','Connectome-seq marker ends'],
  ['magnitude_rank','LIANA magnitude rank'],['specificity_rank','LIANA specificity rank'],
  ['WT_mouse_support','WT mice'],['cKO_mouse_support','cKO mice']
];
function selectedPairRows(rows) {
  const ref=value('resource'), size=value('matchCount');
  return rows.filter(p=>(!ref||p[ref])&&(!size||(size==='2+'?p.reference_count>=2:p.reference_count===Number(size)))&&
    (!$('assay').checked||Boolean(p.PPI_pmids)));
}
function ribbonKey(p) { return [p.source,p.ligand,p.target,p.receptor].join('|'); }
function isolateRibbon(rows) {
  const select=$('oneRibbon'), chosen=select.value;
  select.replaceChildren(new Option('All matching pairs',''));
  for(const p of rows) select.add(new Option(pop(p.source)+' '+p.ligand+' → '+pop(p.target)+' '+p.receptor,ribbonKey(p)));
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
    const stroke=p.CA1_role==='ligand'?caColor:inColor;
    const path='M'+xy(point(a,r))+' A'+r+' '+r+' 0 0 1 '+xy(point(b,r))+
      ' C'+xy(point(b,r*.36))+' '+xy(point(c,r*.36))+' '+xy(point(c,r))+
      ' L'+xy(point(t.angle,r+12))+' L'+xy(point(d,r))+
      ' C'+xy(point(d,r*.36))+' '+xy(point(a,r*.36))+' '+xy(point(a,r))+' Z';
    const desc=pop(p.source)+' '+p.ligand+' (ligand) → '+pop(p.target)+' '+p.receptor+' (receptor) · '+p.exact_reference_set;
    return '<path class="lr-ribbon" tabindex="0" role="button" data-key="'+esc(ribbonKey(p))+'" aria-label="Isolate '+esc(desc)+'" data-source="'+esc(p.source)+'" data-target="'+esc(p.target)+'" data-ligand="'+esc(p.ligand)+'" data-receptor="'+esc(p.receptor)+'" d="'+path+'" fill="'+stroke+'" fill-opacity="'+alpha+'" stroke="'+stroke+'" stroke-width=".8"><title>'+esc(desc)+'</title></path>';
  });
  const nCA=rows.filter(r=>r.CA1_role==='ligand').length;
  const legend=(x,y,color,caption)=>'<path d="M'+x+' '+(y-9)+'h25l8 5 -8 5h-25z" fill="'+color+'"/>'+text(x+40,y,caption,'font-size="14"');
  return '<svg xmlns="http://www.w3.org/2000/svg" width="'+width+'" height="'+height+'" viewBox="0 0 '+width+' '+height+'" role="img" aria-label="Directed CA1 and interneuron ligand–receptor ribbons"><rect width="100%" height="100%" fill="white"/><g font-family="Arial, sans-serif" fill="#213b43">'+
    text(26,32,'CA1–interneuron ligand–receptor pairs','font-size="22" font-weight="700"')+
    text(26,57,rows.length+' population-specific candidates  ·  '+(value('resource')||'All references')+($('assay').checked?'  ·  indexed pair assay':''),'font-size="14"')+
    arcs.join('')+ribbons.join('')+names.join('')+(!rows.length?text(cx,cy,'No candidates match the selected filters.','text-anchor="middle" font-size="16"'):'')+
    legend(120,570,caColor,'CA1 ligand → interneuron receptor  ('+nCA+')')+
    legend(495,570,inColor,'Interneuron ligand → CA1 receptor  ('+(rows.length-nCA)+')')+
    text(480,599,'Arrow tip = receptor end.  Gene ↑ / ↓ = cKO RNA relative to WT.','text-anchor="middle" font-size="13"')+
    text(480,619,'One ribbon per row. Width has no quantitative meaning.','text-anchor="middle" font-size="12"')+'</g></svg>';
}
