(function(root){'use strict';
const schema='satb2-ca1-in-expanded/v1';
const finite=x=>x!==null&&x!==undefined&&x!==''&&Number.isFinite(Number(x));
const yes=x=>x===true||x===1||x==='true'||x==='True'||x==='TRUE'||x==='1';
const no=x=>x===false||x===0||x==='false'||x==='False'||x==='FALSE'||x==='0';
function status(v){return yes(v)?'DEG':no(v)?'not-DEG':'unavailable';}
function key(r){return JSON.stringify([r.source,r.target,r.ligand_complex,r.receptor_complex]);}
function filter(rows,state,nodes){const side=new Map(nodes.map(n=>[n.id,n.side])),term=String(state.search||'').trim().toLowerCase();return rows.filter(r=>{
 const ss=side.get(r.source),ts=side.get(r.target);if(!ss||!ts||ss===ts)return false;
 if(state.direction==='ca1_in'&&!(ss==='CA1'&&ts==='IN'))return false;if(state.direction==='in_ca1'&&!(ss==='IN'&&ts==='CA1'))return false;
 const ca=ss==='CA1'?r.source:r.target, i=ss==='IN'?r.source:r.target;
 if(state.ca1.length&&!state.ca1.includes(ca))return false;if(state.inhibitory.length&&!state.inhibitory.includes(i))return false;
 if(state.threshold==='5'&&!yes(r.passes_5pct))return false;if(state.threshold==='10'&&!yes(r.passes_10pct_projection))return false;
 const l=yes(r.ligand_any_deg),rec=yes(r.receptor_any_deg),la=yes(r.ligand_all_subunits_deg),ra=yes(r.receptor_all_subunits_deg);
 if(state.rule==='both'&&!(la&&ra))return false;if(state.rule==='both_any'&&!(l&&rec))return false;if(state.rule==='ligand_all'&&!la)return false;if(state.rule==='receptor_all'&&!ra)return false;if(state.rule==='ligand'&&!l)return false;if(state.rule==='receptor'&&!rec)return false;if(state.rule==='either'&&!(l||rec))return false;
 if(state.ligandSign!=='any'&&r.ligand_deg_direction!==state.ligandSign)return false;if(state.receptorSign!=='any'&&r.receptor_deg_direction!==state.receptorSign)return false;
 if(finite(state.maxRank)&&(!finite(r.magnitude_rank)||Number(r.magnitude_rank)>Number(state.maxRank)))return false;
 return !term||[r.source,r.target,r.ligand_complex,r.receptor_complex].some(x=>String(x||'').toLowerCase().includes(term));
 });}
function eligible(pairs,state,nodes){let side=new Map(nodes.map(n=>[n.id,n.side]));return pairs.filter(p=>{let ss=side.get(p.source),ts=side.get(p.target);if(!ss||!ts||ss===ts)return false;if(state.direction==='ca1_in'&&ss!=='CA1')return false;if(state.direction==='in_ca1'&&ss!=='IN')return false;let c=ss==='CA1'?p.source:p.target,i=ss==='IN'?p.source:p.target;return (!state.ca1.length||state.ca1.includes(c))&&(!state.inhibitory.length||state.inhibitory.includes(i));});}
function tsv(rows,columns){let quote=x=>x==null?'':String(typeof x==='object'?JSON.stringify(x):x).replace(/[\t\r\n]+/g,' ');return columns.join('\t')+'\n'+rows.map(r=>columns.map(k=>quote(r[k])).join('\t')).join('\n')+'\n';}
const api={schema,finite,yes,no,status,key,filter,eligible,tsv};root.SATB2ExpandedModel=api;if(typeof module!=='undefined'&&module.exports)module.exports=api;
})(typeof window!=='undefined'?window:globalThis);
