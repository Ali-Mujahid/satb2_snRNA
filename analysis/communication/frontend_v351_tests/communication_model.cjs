#!/usr/bin/env node
/* Independent acceptance of the saved data, filtering and native TSV contract.
 * No browser automation or scientific tools are launched by this test.
 */
'use strict';
const fs=require('fs'),path=require('path'),vm=require('vm'),zlib=require('zlib'),crypto=require('crypto'),assert=require('assert/strict');
const repo=path.resolve(__dirname,'..'),base=path.join(repo,'site/communication'),dataDir=path.join(base,'data');
const read=p=>fs.readFileSync(p,'utf8'),sha=b=>crypto.createHash('sha256').update(b).digest('hex');
const checks=[];function check(name,fn){fn();checks.push(name);}
const context={window:{},Uint8Array,TextDecoder,TextEncoder};vm.createContext(context);
vm.runInContext(read(path.join(repo,'site/vendor/fflate.min.js')),context);
vm.runInContext(read(path.join(dataDir,'data.js')),context);
vm.runInContext(read(path.join(base,'model.js')),context);
const compressed=Buffer.from(context.window.CCC_DATA_GZIP,'base64');
const raw=Buffer.from(context.fflate.gunzipSync(compressed));
const payload=JSON.parse(raw.toString('utf8')),manifest=JSON.parse(read(path.join(dataDir,'manifest.json'))),M=context.window.CCCModel;
const data=M.decode(payload);
// Match the browser's explicit context stamping; gene effects do not set context.
for(const [key,d] of Object.entries(data.datasets))for(const r of d.rows){[r.scope,r.condition]=key.split('|');}
const expected={'IN_clusters|pooled':20471,'IN_clusters|WT':19467,'IN_clusters|cKO':21680,'IN_total|pooled':1108,'IN_total|WT':1040,'IN_total|cKO':1170};
const baseState={direction:'both',sources:[],targets:[],search:'',tier:'DEG_supported',maxRank:.05,bothDEG:false};
const state=patch=>({...baseState,...patch});
const rows=data.datasets['IN_clusters|pooled'].rows;
const keys=values=>new Set(values.map(M.interactionKey));
function equalSets(a,b){assert.equal(a.size,b.size);for(const value of a)assert(b.has(value),value);}
function walk(directory){return fs.readdirSync(directory,{withFileTypes:true}).flatMap(e=>e.isDirectory()?walk(path.join(directory,e.name)):[path.join(directory,e.name)]);}

check('actual bundled fflate transport and all cryptographic payload hashes',()=>{
 assert.equal(payload.schema,'satb2-liana-atlas/v1');
 assert.equal(sha(raw),manifest.transport.json_sha256);
 assert.equal(sha(compressed),manifest.transport.gzip_sha256);
 assert.equal(sha(fs.readFileSync(path.join(dataDir,'data.js'))),manifest.transport.script_sha256);
 assert.equal(payload.columns.length,48);
 assert.equal(payload.columns.includes('scope'),false);
 assert.equal(payload.columns.includes('condition'),false);
});
check('all six saved fit counts, independent contexts and native interaction identities',()=>{
 assert.equal(Object.keys(data.datasets).length,6);
 let total=0;
 for(const [key,expectedCount] of Object.entries(expected)){
  const dataset=data.datasets[key],observed=dataset.rows;
  assert.equal(observed.length,expectedCount,key);
  assert.equal(keys(observed).size,expectedCount,key+' keys must be unique');
  const [scope,condition]=key.split('|');
  assert.equal(dataset.eligible_pairs.length,scope==='IN_clusters'?80:4);
  for(const r of observed){assert.equal(r.scope,scope);assert.equal(r.condition,condition);assert.equal(r.direction,r.source.startsWith('CA1_')?'CA1_to_IN':'IN_to_CA1');}
  assert.equal(M.filter(observed,state({tier:'all_candidates',maxRank:1})).length,expectedCount);
  if(scope==='IN_total')assert(observed.every(r=>!/^IN_C\d+$/.test(r.source)&&!/^IN_C\d+$/.test(r.target)));
  else assert(observed.every(r=>r.source!=='IN_Total'&&r.target!=='IN_Total'));
  total+=observed.length;
 }
 assert.equal(total,64936);
});
check('default 1589, strict 324, two exact both-end component candidates',()=>{
 const defaults=M.filter(rows,state()),strict=M.filter(rows,state({tier:'replicated_specific'}));
 assert.equal(defaults.length,1589);
 assert.equal(M.filter(rows,state({direction:'CA1_to_IN'})).length,691);
 assert.equal(M.filter(rows,state({direction:'IN_to_CA1'})).length,898);
 assert.equal(strict.length,324);
 assert.equal(strict.filter(r=>r.direction==='CA1_to_IN').length,112);
 assert.equal(strict.filter(r=>r.direction==='IN_to_CA1').length,212);
 const both=M.filter(rows,state({tier:'replicated_specific',bothDEG:true}));
 assert.equal(both.length,2);
 assert(both.every(r=>r.source==='IN_C6'&&r.ligand_complex==='Cntn3'&&r.receptor_complex==='Ptprg'&&r.sample_support_WT===6&&r.sample_support_cKO===6));
 equalSets(new Set(both.map(r=>r.target)),new Set(['CA1_Deep','CA1_Superficial']));
 equalSets(keys(M.filter(rows,state({tier:'replicated_specific',maxRank:1}))),keys(strict));
 for(const condition of ['pooled','WT','cKO'])assert.equal(M.filter(data.datasets['IN_total|'+condition].rows,state({tier:'replicated_specific'})).length,0);
});
check('direction, source, target, rank, query and both-DEG restrictions combine with AND',()=>{
 const patch={direction:'IN_to_CA1',sources:['IN_C6'],targets:['CA1_Deep'],search:'cntn3',tier:'replicated_specific',bothDEG:true};
 const actual=M.filter(rows,state(patch));
 assert.equal(actual.length,1);assert.equal(actual[0].receptor_complex,'Ptprg');
 equalSets(keys(actual),keys(rows.filter(r=>r.magnitude_rank<=.05&&r.deg_supported&&r.both_sides_deg_supported&&r.source==='IN_C6'&&r.target==='CA1_Deep'&&r.ligand_complex.toLowerCase().includes('cntn3')&&r.specificity_rank<=.05&&r.cellphone_pvals<=.05&&r.sample_support_WT>=3&&r.sample_support_cKO>=3&&r.source_review==='REFERENCE_SUPPORTED')));
 assert.equal(M.filter(rows,state({...patch,targets:['CA1_Superficial'],direction:'CA1_to_IN'})).length,0);
 equalSets(keys(M.filter(rows,state({search:'pTpRg'}))),keys(M.filter(rows,state({search:'ptprg'}))));
 const forward=M.filter(rows,state({direction:'CA1_to_IN'})),reverse=M.filter(rows,state({direction:'IN_to_CA1'}));
 equalSets(keys([...forward,...reverse]),keys(M.filter(rows,state())));
 for(const direction of ['both','CA1_to_IN','IN_to_CA1'])for(const rank of [.001,.01,.05,1]){
  const all=M.filter(rows,state({direction,maxRank:rank,tier:'all_candidates'}));
  const deg=M.filter(rows,state({direction,maxRank:rank}));
  const both=M.filter(rows,state({direction,maxRank:rank,bothDEG:true}));
  const allKeys=keys(all),degKeys=keys(deg);
  assert(deg.every(r=>allKeys.has(M.interactionKey(r))));
  assert(both.every(r=>degKeys.has(M.interactionKey(r))));
 }
});
check('eligible pairs retain zero-count possibilities and all source-only/target-only endpoints',()=>{
 const pairs=data.datasets['IN_clusters|pooled'].eligible_pairs;
 assert.equal(M.eligible(pairs,state()).length,80);
 assert.equal(M.eligible(pairs,state({direction:'CA1_to_IN'})).length,40);
 const sourceOnly=M.eligible(pairs,state({sources:['CA1_Deep']}));
 assert.equal(sourceOnly.length,20);
 assert.equal(new Set(sourceOnly.flatMap(p=>[p.source,p.target])).size,21);
 const targetOnly=M.eligible(pairs,state({targets:['CA1_Deep']}));
 assert.equal(targetOnly.length,20);
 assert.equal(new Set(targetOnly.flatMap(p=>[p.source,p.target])).size,21);
 assert.equal(M.eligible(pairs,state({sources:['CA1_Deep'],targets:['CA1_Superficial']})).length,0);
 assert(pairs.every(p=>!(p.source.startsWith('IN_')&&p.target.startsWith('IN_'))));
 assert(pairs.every(p=>!(p.source.startsWith('CA1_')&&p.target.startsWith('CA1_'))));
});
check('exact seven-column exports reproduce every native saved score and retain full precision',()=>{
 assert.equal(Array.from(M.SEVEN).join('\t'),'source\ttarget\tligand_complex\treceptor_complex\tlr_logfc\tlr_gmeans\tmagnitude_rank');
 for(const [key,dataset] of Object.entries(data.datasets)){
  const [scope,condition]=key.split('|');
  const text=read(path.join(repo,'communication_data/full_result_tables',scope,condition,'bidirectional_rank_aggregate.tsv')).trimEnd();
  const [header,...nativeLines]=text.split(/\r?\n/);
  assert.equal(header,Array.from(M.SEVEN).join('\t'));
  assert.equal(nativeLines.length,dataset.rows.length);
  const byKey=new Map(dataset.rows.map(r=>[M.interactionKey(r),r]));
  for(const line of nativeLines){const a=line.split('\t');assert.equal(a.length,7);const r=byKey.get(JSON.stringify(a.slice(0,4)));assert(r);for(let i=4;i<7;i++)assert.equal(Number(a[i]),r[M.SEVEN[i]]);}
  const filtered=M.filter(dataset.rows,state()),exported=M.tsv(M.sortRows(filtered)),lines=exported.trimEnd().split('\n');
  assert.equal(lines[0],header);assert.equal(lines.length,filtered.length+1);
  const filteredKeys=keys(filtered);
  for(const line of lines.slice(1)){const a=line.split('\t');assert.equal(a.length,7);assert(filteredKeys.has(JSON.stringify(a.slice(0,4))));}
  const extended=M.tsv(dataset.rows.slice(0,2),['scope','condition',...payload.columns]).trimEnd().split('\n');
  assert.equal(extended[1].split('\t')[0],scope);assert.equal(extended[1].split('\t')[1],condition);
 }
});
check('component gene evidence uses the exact population, stable IDs and gene-level signed effects',()=>{
 const cntn=data.components.IN_C6.Cntn3.find(p=>p.symbol==='Cntn3');
 const deep=data.components.CA1_Deep.Ptprg.find(p=>p.symbol==='Ptprg');
 const superficial=data.components.CA1_Superficial.Ptprg.find(p=>p.symbol==='Ptprg');
 assert(cntn.primary_deg_support&&deep.primary_deg_support&&superficial.primary_deg_support);
 assert(cntn.stable_id_results.some(r=>r.log2FoldChange<0&&r.padj<.05));
 assert(deep.stable_id_results.some(r=>r.log2FoldChange>0&&r.padj<.05));
 assert(superficial.stable_id_results.some(r=>r.log2FoldChange>0&&r.padj<.05));
 const untested=data.nodes.find(n=>n.id==='IN_C18');assert.equal(untested.de_tested,false);assert.equal(untested.native_DEG_FDR05,null);
 for(const group of Object.values(data.components))for(const complex of Object.values(group))for(const part of complex){
  assert(Array.isArray(part.stable_id_results));
  if(part.ambiguous_symbol||!part.tested_endpoint)assert.equal(part.primary_deg_support,false);
 }
 for(const d of Object.values(data.datasets))for(const r of d.rows){
  const ligand=data.components[r.source][r.ligand_complex],receptor=data.components[r.target][r.receptor_complex];
  assert(ligand&&receptor);
  assert.equal(r.ligand_deg_any,ligand.some(p=>p.primary_deg_support));
  assert.equal(r.receptor_deg_any,receptor.some(p=>p.primary_deg_support));
  assert.equal(r.deg_supported,r.ligand_deg_any||r.receptor_deg_any);
  assert.equal(r.both_sides_deg_supported,r.ligand_deg_any&&r.receptor_deg_any);
  assert(r.lr_prop>=-1e-9&&r.lr_prop<=1+1e-9);
 }
});
check('explicit interpretation guardrails avoid genotype interaction significance claims',()=>{
 assert.equal(data.provenance.genotype_interaction_test,false);
 assert.equal(data.provenance.no_new_scientific_fits,true);
 assert.match(data.provenance.score_semantics.lr_logfc,/not cKO\/WT/i);
 assert.match(data.provenance.score_semantics.magnitude_rank,/not a genotype/i);
 assert.match(data.provenance.score_semantics.cellphone_pvals,/cell-label specificity/i);
 const index=read(path.join(base,'index.html')),app=read(path.join(base,'app.js')),interpretation=read(path.join(base,'RESULTS_AND_INTERPRETATION.md'));
 assert.match(index,/no replicate-aware differential-interaction test was fitted/i);
 assert.match(index,/not a new genotype FDR test/i);
 assert.match(app,/No genotype interaction test/i);
 assert.match(interpretation,/no interaction is called significantly altered between WT and cKO/i);
 assert(!payload.columns.some(c=>/genotype.*(?:pval|padj|fdr)/i.test(c)));
});
let verifiedNative=0,verifiedOriginal=0;
check('all copied native scientific bytes and local original source hashes remain unchanged',()=>{
 const receipt=JSON.parse(read(path.join(dataDir,'TRANSPORT_RECEIPT.json')));
 const originalRoot=process.env.CCC_ORIGINAL_REPORT||'/Users/mujahid357/Documents/codex_new/SATB2_CCC_20261022/deliverables/SATB2_communication_report';
 for(const item of receipt.files){
  const transport=fs.readFileSync(path.join(repo,item.path));
  assert.equal(sha(transport),item.transport_sha256,item.path+' transport SHA');
  if(item.native_bytes_preserved){const native=item.compression==='gzip'?zlib.gunzipSync(transport):transport;assert.equal(sha(native),item.source_sha256,item.path+' native SHA');verifiedNative++;}
  const original=path.isAbsolute(item.source)?item.source:path.join(originalRoot,item.source);
  if(fs.existsSync(original)){assert.equal(sha(fs.readFileSync(original)),item.source_sha256,original+' original SHA');verifiedOriginal++;}
 }
 const evidence=fs.readFileSync(path.join(repo,'communication_data/full_result_tables/interaction_evidence.tsv.gz'));
 assert.equal(sha(evidence),data.provenance.source_evidence_sha256);
 const sevenHeader='gene\tbaseMean\tlog2FoldChange\tlfcSE\tstat\tpvalue\tpadj';
 for(const file of walk(path.join(repo,'communication_data/existing_DESeq2_exports')).filter(p=>/(?:ALL_RESULTS|DEG_FDR05)\.txt$/.test(p)))assert.equal(read(file).split(/\r?\n/,1)[0],sevenHeader,file);
});
console.log(JSON.stringify({passed:checks.length,checks,datasets:expected,nativeRecords:64936,defaultDEG:1589,strict:324,bothEndpointDEG:2,nativeFilesHashVerified:verifiedNative,availableOriginalSourceHashesVerified:verifiedOriginal},null,2));
