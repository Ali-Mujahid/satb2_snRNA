#!/usr/bin/env node
'use strict';
/* Regression for the 23 October user correction. This independently reconstructs
 * the strict, population-matched, BOTH-end DEG subset from the unchanged saved
 * component DESeq2 evidence. No analysis model or browser is launched. */
const fs=require('node:fs'),path=require('node:path'),vm=require('node:vm'),zlib=require('node:zlib'),assert=require('node:assert/strict');
const base=path.join(__dirname,'../site/communication'),sandbox={window:{}};
vm.createContext(sandbox);
for(const name of ['data/data.js','model.js'])vm.runInContext(fs.readFileSync(path.join(base,name),'utf8'),sandbox,{filename:name});
const payload=JSON.parse(zlib.gunzipSync(Buffer.from(sandbox.window.CCC_DATA_GZIP,'base64'))),M=sandbox.window.CCCModel,data=M.decode(payload);
const ALLOWED=new Set(['CA1_Deep','CA1_Superficial','IN_C1','IN_C6','IN_C8']);
const expectedAll={pooled:16,WT:15,cKO:16},expectedRank05={pooled:8,WT:9,cKO:7};
const state=extra=>({scope:'IN_clusters',condition:'pooled',direction:'both',sources:[],targets:[],search:'',tier:'DEG_supported',maxRank:1,bothDEG:true,...extra});
const key=r=>JSON.stringify([r.source,r.target,r.ligand_complex,r.receptor_complex]);
const keySet=rows=>new Set(rows.map(key));
const equalKeys=(actual,expected,note)=>{assert.equal(actual.size,expected.size,note);for(const k of expected)assert(actual.has(k),note+' '+k);};
function parts(population,complex){return payload.components[population]?.[complex]||[];}
function allSignificant(population,complex){
 const p=parts(population,complex);
 return p.length>0&&p.every(g=>g.tested_endpoint===true&&g.ambiguous_symbol===false&&g.primary_deg_support===true&&g.stable_id_results.some(s=>Number.isFinite(s.padj)&&s.padj<.05&&Number.isFinite(s.log2FoldChange)&&s.log2FoldChange!==0));
}
function independentlyStrict(r){return ALLOWED.has(r.source)&&ALLOWED.has(r.target)&&allSignificant(r.source,r.ligand_complex)&&allSignificant(r.target,r.receptor_complex);}
const checks=[];
function check(name,fn){fn();checks.push(name);}

check('Independent component evidence identifies the five-population strict subset',()=>{
 for(const condition of ['pooled','WT','cKO']){
  const rows=data.datasets['IN_clusters|'+condition].rows,strict=rows.filter(independentlyStrict);
  assert.equal(strict.length,expectedAll[condition],condition+' all ranks');
  assert.equal(strict.filter(r=>r.magnitude_rank<=.05).length,expectedRank05[condition],condition+' rank <= .05');
  assert(strict.every(r=>r.ligand_all_subunits_deg&&r.receptor_all_subunits_deg));
  assert(!strict.some(r=>r.source==='IN_C8'||r.target==='IN_C8'),'C8 retains context, but no LR qualifies under BOTH-end DEG support.');
 }
});
check('Active filters cannot re-enable OR-only support or other populations through legacy settings',()=>{
 for(const condition of ['pooled','WT','cKO']){
  const rows=data.datasets['IN_clusters|'+condition].rows,strict=rows.filter(independentlyStrict),expected=keySet(strict);
  for(const legacy of [{tier:'all_candidates',bothDEG:false},{tier:'DEG_supported',bothDEG:false},{scope:'IN_total',tier:'all_candidates',bothDEG:false}])equalKeys(keySet(M.filter(rows,state(legacy))),expected,condition+' legacy option cannot loosen the invariant');
  equalKeys(keySet(M.filter(rows,state({maxRank:.05}))),keySet(strict.filter(r=>r.magnitude_rank<=.05)),condition+' default rank');
  assert.equal(M.filter(data.datasets['IN_total|'+condition].rows,state({tier:'all_candidates',bothDEG:false})).length,0,'IN total has no place in the active strict cohort');
 }
});
check('All-subunit and ambiguity gates are independent mandatory invariants, not optional checkboxes',()=>{
 const rows=data.datasets['IN_clusters|pooled'].rows,lead=rows.find(r=>r.source==='IN_C6'&&r.target==='CA1_Deep'&&r.ligand_complex==='Cntn3'&&r.receptor_complex==='Ptprg');
 assert(lead&&independentlyStrict(lead));
 for(const patch of [
  {ligand_all_subunits_deg:false},
  {receptor_all_subunits_deg:false},
  {both_sides_deg_supported:false},
  {ligand_de_symbol_ambiguous:true},
  {receptor_de_symbol_ambiguous:true},
  {source:'IN_C2'},
  {target:'IN_Total'}
 ])assert.equal(M.filter([{...lead,...patch}],state({tier:'all_candidates',bothDEG:false})).length,0,JSON.stringify(patch));
});
check('Direction and population filters restrict the fixed cohort without adding rows',()=>{
 for(const condition of ['pooled','WT','cKO']){
  const rows=data.datasets['IN_clusters|'+condition].rows,strict=rows.filter(independentlyStrict);
  const f=M.filter(rows,state({direction:'CA1_to_IN'})),r=M.filter(rows,state({direction:'IN_to_CA1'}));
  equalKeys(keySet([...f,...r]),keySet(strict),condition+' union');
  assert(f.every(x=>x.source.startsWith('CA1_')&&x.target.startsWith('IN_')));
  assert(r.every(x=>x.source.startsWith('IN_')&&x.target.startsWith('CA1_')));
  assert.equal(M.filter(rows,state({sources:['IN_C2']})).length,0);
  assert.equal(M.filter(rows,state({targets:['IN_C20']})).length,0);
 }
});
check('Strict active payload contains only 47 preserved native records and exact matching signed DESeq2 evidence',()=>{
 const script=fs.readFileSync(path.join(base,'data/deg_only_data.js'),'utf8');
 vm.runInContext(script,sandbox,{filename:'data/deg_only_data.js'});
 const strictPayload=JSON.parse(zlib.gunzipSync(Buffer.from(sandbox.window.CCC_DATA_GZIP,'base64'))),strictData=M.decode(strictPayload);
 assert.equal(Object.keys(strictData.datasets).length,3);assert.equal(strictData.nodes.length,5);assert(strictData.nodes.every(n=>ALLOWED.has(n.id)));
 assert.equal(Object.values(strictData.datasets).reduce((n,d)=>n+d.rows.length,0),47);
 const directionExpected={pooled:{'up→up':3,'up→down':2,'down→up':7,'down→down':4},WT:{'up→up':2,'up→down':2,'down→up':7,'down→down':4},cKO:{'up→up':3,'up→down':2,'down→up':7,'down→down':4}};
 for(const condition of ['pooled','WT','cKO']){
  const nativeRows=data.datasets['IN_clusters|'+condition].rows,expected=nativeRows.filter(independentlyStrict),strict=strictData.datasets['IN_clusters|'+condition];
  assert.equal(strict.rows.length,expectedAll[condition]);assert.equal(strict.eligible_pairs.length,12);equalKeys(keySet(strict.rows),keySet(expected),condition+' active payload exact keys');
  const native=new Map(nativeRows.map(r=>[key(r),r]));
  for(const r of strict.rows){
   for(const col of payload.columns)assert.deepEqual(r[col],native.get(key(r))[col],condition+' preserved native column '+col);
   for(const [endpoint,complex,field]of [[r.source,r.ligand_complex,'ligand'],[r.target,r.receptor_complex,'receptor']]){
    const genes=parts(endpoint,complex);assert.equal(genes.length,1,'This saved strict subset contains monomeric pairs.');
    const previous=genes[0].stable_id_results.find(g=>g.padj<.05),current=strictPayload.components[endpoint][complex][0].stable_id_results.find(g=>g.padj<.05),direction=previous.log2FoldChange>0?'up':'down';
    // The earlier JSON copy rounded a few final digits. New component values
    // come directly from original seven-column DESeq2 exports, so assert their
    // exact row/component match and numerical continuity with the legacy copy.
    assert(Math.abs(current.log2FoldChange-previous.log2FoldChange)<1e-12);assert(Math.abs(current.padj-previous.padj)<1e-12);
    assert.equal(r[field+'_deg_direction'],direction);assert.equal(r[field+'_deseq2_log2fc'],current.log2FoldChange);assert.equal(r[field+'_deseq2_padj'],current.padj);
   }
   assert.equal(r.deg_direction_pair,r.ligand_deg_direction+'→'+r.receptor_deg_direction);
  }
  for(const [combo,count]of Object.entries(directionExpected[condition])){
   const [ligandDirection,receptorDirection]=combo.split('→');
   assert.equal(M.filter(strict.rows,state({directionPair:combo})).length,count,condition+' '+combo+' compound selector');
   assert.equal(M.filter(strict.rows,state({ligandDirection,receptorDirection})).length,count,condition+' '+combo+' independent selectors');
  }
  assert.equal(M.filter(strict.rows,state({directionPair:'mixed'})).length,0);
  assert.equal(M.filter(strict.rows,state({directionPair:'up→up',ligandDirection:'down'})).length,0);
 }
});
check('Twelve eligible strict directed pairs retain C8 zero-result cells, not its old broad candidates',()=>{
 for(const condition of ['pooled','WT','cKO']){
  const pairs=data.datasets['IN_clusters|'+condition].eligible_pairs,eligible=M.eligible(pairs,state());
  assert.equal(eligible.length,12);
  assert(eligible.every(p=>ALLOWED.has(p.source)&&ALLOWED.has(p.target)));
  assert.equal(eligible.filter(p=>p.source==='IN_C8'||p.target==='IN_C8').length,4);
  assert.equal(M.eligible(pairs,state({direction:'CA1_to_IN'})).length,6);
  assert.equal(M.eligible(pairs,state({direction:'IN_to_CA1'})).length,6);
 }
});
check('Seven-column native export contract remains exact while component gene effects stay signed',()=>{
 assert.equal(Array.from(M.SEVEN).join('\t'),'source\ttarget\tligand_complex\treceptor_complex\tlr_logfc\tlr_gmeans\tmagnitude_rank');
 const cntn=parts('IN_C6','Cntn3')[0].stable_id_results.find(g=>g.padj<.05),deep=parts('CA1_Deep','Ptprg')[0].stable_id_results.find(g=>g.padj<.05),sup=parts('CA1_Superficial','Ptprg')[0].stable_id_results.find(g=>g.padj<.05);
 assert(cntn.log2FoldChange<0);assert(deep.log2FoldChange>0);assert(sup.log2FoldChange>0);
 assert.match(payload.provenance.score_semantics.lr_logfc,/not cKO\/WT/i,'Native lr_logfc cannot become the gene cKO/WT effect.');
});
console.log(JSON.stringify({status:'PASS',checks:checks.length,names:checks,expectedAll,expectedRank05,mode:'PURE_DATA_AND_VM_NO_BROWSER_NO_FITS'},null,2));
