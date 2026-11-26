'use strict';
const bindingStatusLabels={direct:'Direct binding',screen:'Binding screen',indirect:'Indirect crosstalk',unverified:'UNVERIFIED binding'};
function arrowCaption(){return {molecular:'Reviewed molecular interaction',synaptic:'Predicted synaptic transmission',resource:'Archived LIANA convention'}[value('arrowMode')];}
function arrowMeaning(){
  return {
    molecular:'Single tip = molecular receptor end; double tips = adhesion partners. Synaptic transmission direction is a separate annotation.',
    synaptic:'Tips show a predicted presynaptic-to-postsynaptic direction. Tip-free links have unresolved synaptic orientation.',
    resource:'Tips retain the archived ligand-to-receptor convention. This convention alone does not assign pre- or postsynaptic location.'
  }[value('arrowMode')];
}
function connectionMode(p){
  const mode=value('arrowMode');
  if(mode==='resource')return 'forward';
  if(mode==='synaptic'){
    if(!p.synaptic_source)return 'none';
    return p.synaptic_source===p.source?'forward':'reverse';
  }
  return p.molecular_mode==='adhesion'?'both':p.molecular_mode==='signal'?'forward':'none';
}
function reviewColor(p,caColor,inColor){
  if(['unverified','indirect'].includes(p.binding_status))return '#7f8b90';
  if(p.binding_status==='screen')return '#9b742c';
  return (value('arrowMode')==='resource'?p.source:p.corrected_source).startsWith('CA1')?caColor:inColor;
}
function renderBindingAudit(rows,audit){
  $('auditBox').hidden=!isPairView(view);
  if(!isPairView(view))return;
  const keys=new Set(rows.map(p=>[p.ligand,p.receptor].sort().join('|')));
  const selected=audit.pairs.filter(a=>keys.has([...a.genes].sort().join('|')));
  const counts=Object.keys(bindingStatusLabels).map(k=>rows.filter(p=>p.binding_status===k).length);
  $('auditCounts').textContent=rows.length+' archived population rows · '+selected.length+' molecular pairs · '+counts[0]+' direct-binding rows · '+counts[1]+' screen-supported · '+counts[2]+' indirect · '+counts[3]+' unverified. Reviewed '+audit.review_date+'.';
  $('auditCards').innerHTML=selected.map(a=>{
    const examples=rows.filter(p=>[p.ligand,p.receptor].sort().join('|')===[...a.genes].sort().join('|'));
    const p=examples[0];
    return '<article class="audit-card '+esc(a.status)+'"><h3>'+esc(a.genes.join(' / '))+'</h3><p class="audit-status">'+esc(a.title)+'</p><p>'+esc(a.assay)+'</p><dl><dt>CA1 '+esc(p.CA1_gene)+'</dt><dd>'+esc(p.CA1_compartment)+'</dd><dt>'+esc(pop(p.IN_cluster))+' '+esc(p.IN_gene)+'</dt><dd>'+esc(p.IN_compartment)+'</dd><dt>Synaptic direction</dt><dd>'+esc(p.synaptic_direction)+'</dd><dt>Fragment and isoform</dt><dd>'+esc(a.molecular_form)+'</dd></dl><p>'+esc(a.interpretation)+'</p><details><summary>Experimental context and sources</summary><p>'+esc(a.experimental_context)+'</p>'+a.references.map(r=>'<p><a href="'+esc(r.url)+'" target="_blank" rel="noopener">'+esc(r.label)+'</a></p>').join('')+'</details></article>';
  }).join('');
}
const shortLocations={
  Nrxn1:'Presynaptic axon membrane',Nlgn1:'Postsynaptic excitatory membrane',
  Grm7:'Presynaptic release-site receptor',Penk:'Secreted, processed enkephalins',
  App:'Axonal / dendritic membrane; cleavage',Cntn3:'GPI cell surface; pre/post unresolved',Ptprg:'Cell surface; pre/post unresolved',
  Slit1:'Secreted extracellular cue',Slit2:'Secreted cue; cleavage matters',
  Robo1:'Cell surface; local compartment unresolved',Dcc:'Netrin receptor; local side unresolved',
  Plxna1:'Cell-surface receptor for SLIT-C',Ptn:'Secreted / extracellular matrix',
  Ptprz1:'Membrane receptor / soluble phosphacan',Fstl5:'Secreted extracellular protein',Cntn5:'GPI cell-surface adhesion protein'
};
function bindingFigure(rows){
  const fs=Number(value('fontsize')),width=1180,rowHeight=110,height=160+rows.length*rowHeight;
  const text=(x,y,s,attrs='')=>'<text x="'+x+'" y="'+y+'" '+attrs+'>'+esc(s)+'</text>';
  const caColor=value('caColor'),inColor=value('inColor');
  const parts=['<svg xmlns="http://www.w3.org/2000/svg" width="'+width+'" height="'+height+'" viewBox="0 0 '+width+' '+height+'" role="img" aria-label="Binding evidence and reviewed synaptic roles"><rect width="100%" height="100%" fill="white"/><defs><marker id="tip" markerWidth="8" markerHeight="8" refX="7" refY="4" orient="auto-start-reverse"><path d="M0 0 L8 4 L0 8z" fill="context-stroke"/></marker></defs><g font-family="Arial, sans-serif" fill="#213b43">',text(20,30,'CA1–interneuron binding and compartment audit','font-size="21" font-weight="700"'),text(20,57,rows.length+' archived rows · '+arrowCaption()+' · protein locations from published studies','font-size="13"'),text(20,82,arrowMeaning(),'font-size="12"'),text(20,110,'CA1 sublayer and protein','font-size="14" font-weight="700"'),text(875,110,'Interneuron and protein','font-size="14" font-weight="700"')];
  rows.forEach((p,i)=>{
    const y=131+i*rowHeight;
    const caSource=p.source.startsWith('CA1');
    const ca=caSource?p.ligand:p.receptor,inh=caSource?p.receptor:p.ligand;
    const leftPop=caSource?p.source:p.target,rightPop=caSource?p.target:p.source;
    const c=reviewColor(p,caColor,inColor),mode=connectionMode(p);
    const rightTip=mode==='both'||(mode==='forward'&&caSource)||(mode==='reverse'&&!caSource);
    const leftTip=mode==='both'||(mode==='forward'&&!caSource)||(mode==='reverse'&&caSource);
    const dashed=p.binding_status!=='direct'||mode==='none';
    if(i%2===0)parts.push('<rect x="10" y="'+(y-9)+'" width="1160" height="106" fill="#f3f6f5"/>');
    parts.push(text(20,y+12,pop(leftPop),'font-size="12"'),text(20,y+35,ca+((caSource?p.ligand_log2FC:p.receptor_log2FC)>0?' ↑':' ↓'),'font-size="'+(fs+2)+'" font-weight="700"'),text(20,y+60,shortLocations[ca],'font-size="11"'));
    parts.push(text(875,y+12,pop(rightPop),'font-size="12"'),text(875,y+35,inh+((caSource?p.receptor_log2FC:p.ligand_log2FC)>0?' ↑':' ↓'),'font-size="'+(fs+2)+'" font-weight="700"'),text(875,y+60,shortLocations[inh],'font-size="11"'));
    parts.push(text(590,y+11,bindingStatusLabels[p.binding_status],'text-anchor="middle" font-size="12" font-weight="700" fill="'+c+'"'));
    const x1=342,x2=838,yy=y+35;
    parts.push('<path d="M'+x1+' '+yy+'H'+x2+'" fill="none" stroke="'+c+'" stroke-width="'+(Number(value('ribbonWidth'))/35)+'"'+(dashed?' stroke-dasharray="7 5"':'')+(leftTip?' marker-start="url(#tip)"':'')+(rightTip?' marker-end="url(#tip)"':'')+'/>');
    parts.push(text(590,y+62,value('arrowMode')==='synaptic'?p.synaptic_direction:value('arrowMode')==='resource'?p.ligand+' → '+p.receptor+' (archived)':p.corrected_pair,'text-anchor="middle" font-size="'+fs+'"'));
    if(p.binding_status==='unverified')parts.push(text(590,y+82,'Direct partner relationship requires validation','text-anchor="middle" font-size="11" fill="#765b34"'));
    if(p.binding_status==='indirect')parts.push(text(590,y+82,'SLIT/ROBO and netrin/DCC pathway crosstalk','text-anchor="middle" font-size="11"'));
  });
  parts.push(text(20,height-12,'↑ / ↓ = cKO RNA versus WT. Reciprocal NRXN1–NLGN1 entries describe the same proposed excitatory orientation. All local contacts are predictions.','font-size="11"'));
  if(!rows.length)parts.push(text(20,145,'No rows match the current filters.'));
  return parts.join('')+'</g></svg>';
}
