#!/usr/bin/env node
/* Independent, read-only verification of strict DEG-only transport and exports. */
'use strict';
const fs=require('fs'),path=require('path'),vm=require('vm'),zlib=require('zlib'),crypto=require('crypto'),assert=require('assert/strict');
const repo=path.resolve(__dirname,'..'),base=path.join(repo,'site/communication/data'),download=path.join(repo,'communication_data/deg_only_v36');
const read=p=>fs.readFileSync(p,'utf8'),hash=b=>crypto.createHash('sha256').update(b).digest('hex');
const checks=[];function check(name,fn){fn();checks.push(name);}
function payload(name){const context={window:{}};vm.createContext(context);vm.runInContext(read(path.join(base,name)),context);const gzip=Buffer.from(context.window.CCC_DATA_GZIP,'base64'),raw=zlib.gunzipSync(gzip);return {value:JSON.parse(raw),gzip,raw};}
function decode(data){return Object.fromEntries(Object.entries(data.datasets).map(([key,d])=>[key,{...d,rows:d.rows.map(values=>Object.fromEntries(data.columns.map((c,i)=>[c,values[i]])))}]));}
const current=payload('deg_only_data.js'),old=payload('data.js'),manifest=JSON.parse(read(path.join(base,'DEG_ONLY_MANIFEST.json'))),audit=JSON.parse(read(path.join(repo,'docs/DEG_ONLY_ENDPOINT_AUDIT.json'))),datasets=decode(current.value),legacy=decode(old.value);
const expectedNodes=['CA1_Deep','CA1_Superficial','IN_C1','IN_C6','IN_C8'];
const ik=r=>JSON.stringify([r.source,r.target,r.ligand_complex,r.receptor_complex]);
const all=Object.values(datasets).flatMap(d=>d.rows),seven=['source','target','ligand_complex','receptor_complex','lr_logfc','lr_gmeans','magnitude_rank'];
function tsv(p){const [header,...lines]=read(p).trimEnd().split('\n'),columns=header.split('\t');return {columns,rows:lines.map(l=>Object.fromEntries(l.split('\t').map((v,i)=>[columns[i],v])))};}
check('new payload cryptographic transport hashes and unchanged legacy files',()=>{
 assert.equal(hash(current.raw),manifest.transport.json_sha256);assert.equal(hash(current.gzip),manifest.transport.gzip_sha256);assert.equal(hash(fs.readFileSync(path.join(base,'deg_only_data.js'))),manifest.transport.script_sha256);
 for(const [p,sha] of Object.entries(manifest.legacy_preserved_hashes))assert.equal(hash(fs.readFileSync(path.resolve(repo,p))),sha);
 assert.equal(audit.status,'PASS');assert.equal(audit.checks.original_DESeq2_exports_byte_identical,true);
});
check('only three corrected fits and five exact populations; C8 remains explicitly eligible',()=>{
 assert.deepEqual(current.value.nodes.map(n=>n.id),expectedNodes);
 assert.equal(Object.keys(datasets).length,3);assert.equal(all.length,47);
 for(const [condition,count] of Object.entries({pooled:16,WT:15,cKO:16})){
  const d=datasets['IN_clusters|'+condition];assert.equal(d.rows.length,count);assert.equal(new Set(d.rows.map(ik)).size,count);assert.equal(d.eligible_pairs.length,12);
  assert(d.rows.every(r=>expectedNodes.includes(r.source)&&expectedNodes.includes(r.target)));
  const zeros=d.eligible_pairs.filter(p=>p.source==='IN_C8'||p.target==='IN_C8');assert.equal(zeros.length,4);assert(zeros.every(p=>!d.rows.some(r=>r.source===p.source&&r.target===p.target)));
 }
 assert.equal(current.value.nodes.find(n=>n.id==='IN_C8').needs_review,true);
});
check('both complete endpoint complexes are significant and signed from exact original DESeq2 evidence',()=>{
 const used={};for(const r of all){
  const signs=[];
  assert(r.ligand_all_subunits_deg&&r.receptor_all_subunits_deg&&r.both_sides_deg_supported);
  for(const [side,role] of [['source','ligand'],['target','receptor']]){
   const pop=r[side],complex=r[role+'_complex'],parts=current.value.components[pop][complex];assert(parts.length>0);used[pop]??=new Set();used[pop].add(complex);
   for(const p of parts){assert(p.primary_deg_support&&p.tested_endpoint&&!p.ambiguous_symbol);assert(p.stable_id_results.every(g=>g.padj!==null&&g.padj<.05));}
   assert.equal(parts.length,1,'all currently observed strict complexes are single-subunit');const g=parts[0].stable_id_results[0];
   assert.equal(r[role+'_deseq2_log2fc'],g.log2FoldChange);assert.equal(r[role+'_deseq2_padj'],g.padj);const sign=g.log2FoldChange>0?'up':'down';assert.equal(r[role+'_deg_direction'],sign);signs.push(sign);
  }
  assert.equal(r.deg_direction_pair,signs.join('→'));
 }
 for(const [population,complexes] of Object.entries(current.value.components))for(const complex of Object.keys(complexes))assert(used[population]?.has(complex)||(population==='IN_C8'&&complex==='Calb1'),'no non-selected non-DEG partner in payload');
 const calb=current.value.components.IN_C8.Calb1[0];assert.equal(calb.in_mouse_resource,false);assert.equal(calb.stable_id_results[0].log2FoldChange,1.93206760368337);
});
check('all 47 saved scalar fields and identity flags equal the untouched original full-background payload',()=>{
 for(const [key,d] of Object.entries(datasets)){
  const lookup=new Map(legacy[key].rows.map(r=>[ik(r),r]));
  for(const row of d.rows){const original=lookup.get(ik(row));assert(original);for(const column of old.value.columns)assert.deepEqual(row[column],original[column],key+':'+column);}
 }
});
check('all-rank and native rank <=0.05 directional counts and ordered gene-sign counts',()=>{
 for(const [condition,expected] of Object.entries({pooled:[6,10,3,5],WT:[5,10,3,6],cKO:[6,10,3,4]})){
  const rows=datasets['IN_clusters|'+condition].rows;let i=0;
  for(const rank of [1,.05])for(const direction of ['CA1_to_IN','IN_to_CA1'])assert.equal(rows.filter(r=>r.direction===direction&&r.magnitude_rank<=rank).length,expected[i++]);
 }
 for(const count of manifest.direction_pair_counts){const limit=count.rank_selection==='all_saved_ranks'?1:.05;assert.equal(datasets['IN_clusters|'+count.condition].rows.filter(r=>(count.direction==='bidirectional'||r.direction===count.direction)&&r.magnitude_rank<=limit&&r.deg_direction_pair===count.deg_direction_pair).length,count.records);}
});
check('nine native seven-column projections and nine signed exports match strict payload values',()=>{
 for(const condition of ['pooled','WT','cKO'])for(const direction of ['CA1_to_IN','IN_to_CA1','bidirectional']){
  const base=path.join(download,condition,direction),table=tsv(base+'_rank_aggregate.tsv'),extended=tsv(base+'_signed_extended.tsv'),expected=datasets['IN_clusters|'+condition].rows.filter(r=>direction==='bidirectional'||r.direction===direction),lookup=new Map(expected.map(r=>[ik(r),r]));
  assert.deepEqual(table.columns,seven);assert.equal(table.rows.length,expected.length);assert.equal(extended.rows.length,expected.length);
  for(const row of table.rows){const original=lookup.get(ik(row));assert(original);for(const c of seven)assert.equal(typeof original[c]==='number'?Number(row[c]):row[c],original[c]);}
  for(const row of extended.rows){const original=lookup.get(ik(row));assert(original);for(const c of ['ligand_deg_direction','receptor_deg_direction','deg_direction_pair'])assert.equal(row[c],original[c]);}
 }
});
check('five original DEG copies retain exact source bytes and standard seven columns, signed subsets partition them',()=>{
 for(const group of expectedNodes){const original=path.join(repo,'communication_data/existing_DESeq2_exports',group,'DEG_FDR05.txt'),folder=path.join(download,'DESeq2',group),copy=path.join(folder,'DEG_FDR05.txt');assert.equal(hash(fs.readFileSync(original)),hash(fs.readFileSync(copy)));const all=tsv(copy),up=tsv(path.join(folder,'UP_in_cKO.txt')),down=tsv(path.join(folder,'DOWN_in_cKO.txt'));assert.deepEqual(all.columns,['gene','baseMean','log2FoldChange','lfcSE','stat','pvalue','padj']);assert.deepEqual(up.columns,all.columns);assert.deepEqual(down.columns,all.columns);assert.equal(up.rows.length+down.rows.length,all.rows.length);assert(up.rows.every(r=>Number(r.log2FoldChange)>0));assert(down.rows.every(r=>Number(r.log2FoldChange)<0));assert.equal(all.rows.length,manifest.endpoint_summaries[group].native_DEG_FDR05);}
});
check('every documented immutable download exists with exact bytes and hash',()=>{
 for(const file of manifest.export_files){const p=path.join(download,file.path),bytes=fs.readFileSync(p);assert.equal(bytes.length,file.bytes);assert.equal(hash(bytes),file.sha256);assert(file.url.startsWith('https://raw.githubusercontent.com/Ali-Mujahid/satb2_snRNA/atlas-v3.6.0/communication_data/deg_only_v36/'));}
});
const result={status:'PASS',checks:checks.length,passed:checks,strict_records:all.length,source_scores_unchanged:true,no_refits:true};
fs.writeFileSync(path.join(repo,'tests/communication_deg_only_acceptance.json'),JSON.stringify(result,null,2)+'\n');console.log(JSON.stringify(result,null,2));
