(function(root){
'use strict';
const SEVEN=['source','target','ligand_complex','receptor_complex','lr_logfc','lr_gmeans','magnitude_rank'];
const finite=v=>v!==null&&v!==undefined&&v!==''&&Number.isFinite(Number(v));
const pairKey=(s,t)=>JSON.stringify([s,t]);
const lrKey=r=>JSON.stringify([r.ligand_complex,r.receptor_complex]);
const interactionKey=r=>JSON.stringify([r.source,r.target,r.ligand_complex,r.receptor_complex]);
function directionPass(s,t,d){return d==='both'||(d==='CA1_to_IN'?s.startsWith('CA1_')&&t.startsWith('IN_'):s.startsWith('IN_')&&t.startsWith('CA1_'));}
function filter(rows,state){
 const query=(state.search||'').trim().toLowerCase(),rank=Number(state.maxRank);
 return rows.filter(r=>{
  if(!finite(r.magnitude_rank)||Number(r.magnitude_rank)>rank)return false;
  if(!directionPass(r.source,r.target,state.direction))return false;
  if(state.sources.length&&!state.sources.includes(r.source))return false;
  if(state.targets.length&&!state.targets.includes(r.target))return false;
  if(query&&!([r.source,r.target,r.source_label,r.target_label,r.ligand_complex,r.receptor_complex].join(' ').toLowerCase().includes(query)))return false;
  if(state.tier!=='all_candidates'&&!r.deg_supported)return false;
  if(state.bothDEG&&!r.both_sides_deg_supported)return false;
  if(state.tier==='replicated_specific'){
   if(Number(r.magnitude_rank)>Math.min(rank,.05)||!finite(r.specificity_rank)||Number(r.specificity_rank)>.05||!finite(r.cellphone_pvals)||Number(r.cellphone_pvals)>.05)return false;
   if(Number(r.sample_support_WT)<3||Number(r.sample_support_cKO)<3)return false;
   if(r.source.startsWith('IN_')&&r.source_review!=='REFERENCE_SUPPORTED')return false;
   if(r.target.startsWith('IN_')&&r.target_review!=='REFERENCE_SUPPORTED')return false;
  }
  return true;
 });
}
function eligible(pairs,state){return pairs.filter(p=>directionPass(p.source,p.target,state.direction)&&(!state.sources.length||state.sources.includes(p.source))&&(!state.targets.length||state.targets.includes(p.target)));}
function counts(rows){
 const pairs=new Map(),seen=new Set(),molecules=new Set(),degrees=new Map();
 for(const r of rows){const ik=interactionKey(r);if(seen.has(ik))continue;seen.add(ik);molecules.add(lrKey(r));const pk=pairKey(r.source,r.target);pairs.set(pk,(pairs.get(pk)||0)+1);degrees.set(r.source,(degrees.get(r.source)||0)+1);degrees.set(r.target,(degrees.get(r.target)||0)+1);}
 return {records:seen.size,molecules:molecules.size,pairs:pairs.size,edgeMax:Math.max(1,...pairs.values()),degreeMax:Math.max(1,...degrees.values())};
}
function sortRows(rows,key='magnitude_rank',ascending=true){return [...rows].sort((a,b)=>{const av=a[key],bv=b[key];if(!finite(av)&&!finite(bv))return interactionKey(a).localeCompare(interactionKey(b));if(!finite(av))return 1;if(!finite(bv))return -1;return (Number(av)-Number(bv))*(ascending?1:-1)||interactionKey(a).localeCompare(interactionKey(b));});}
function tsv(rows,columns=SEVEN){const cell=v=>v===null||v===undefined?'NA':typeof v==='object'?JSON.stringify(v):String(v);return columns.join('\t')+'\n'+rows.map(r=>columns.map(c=>cell(r[c]).replace(/[\t\r\n]/g,' ')).join('\t')).join('\n')+'\n';}
function decode(data){const datasets={};for(const [k,v] of Object.entries(data.datasets)){datasets[k]={...v,rows:v.rows.map(a=>Object.fromEntries(data.columns.map((c,i)=>[c,a[i]])))};}return {...data,datasets};}
root.CCCModel={SEVEN,finite,pairKey,lrKey,interactionKey,directionPass,filter,eligible,counts,sortRows,tsv,decode};
})(typeof window!=='undefined'?window:globalThis);
