'use strict';
const assert=require('node:assert/strict'),fs=require('node:fs'),path=require('node:path'),vm=require('node:vm');
const root=path.resolve(__dirname,'../analysis/synaptic_evidence_20261125');
const frozen=JSON.parse(fs.readFileSync(path.join(root,'results/data.json'),'utf8'));
const html=fs.readFileSync(path.join(root,'offline.html'),'utf8');
const payload=JSON.parse(html.match(/<script type="application\/json" id="sourceData">([\s\S]*?)<\/script>/)[1]);
assert.equal(payload.pairs.length,14);
assert.equal(payload.pairRoles.length,14);
assert.equal(new Set(payload.pairs.map(p=>[p.ligand,p.receptor].sort().join('|'))).size,10);
// Published annotations must never modify the source RNA values or resource roles.
for(let i=0;i<frozen.pairs.length;i++)for(const [key,val] of Object.entries(frozen.pairs[i]))assert.deepEqual(payload.pairs[i][key],val,key);
assert.deepEqual(payload.pairRoles.reduce((a,p)=>(a[p.binding_status]=(a[p.binding_status]||0)+1,a),{}),{unverified:2,direct:10,screen:1,indirect:1});
let mode='molecular';
const ctx=vm.createContext({value:id=>id==='arrowMode'?mode:'',console});
vm.runInContext(fs.readFileSync(path.join(root,'reader/binding_views.js'),'utf8'),ctx);
const connection=p=>{ctx.p=p;return vm.runInContext('connectionMode(p)',ctx)};
for(const p of payload.pairRoles){
 assert(p.CA1_compartment&&p.IN_compartment&&p.binding_sources);
 assert.equal(p.contact_validation,'Predicted P60 population relationship; cell-specific contact untested');
 if(['App','Penk'].includes(p.CA1_gene)&&p.IN_gene==='Grm7'){
  assert.equal(p.binding_status,'unverified');assert.match(p.IN_compartment,/presynaptic/);assert.equal(connection(p),'none');
 }
 if(p.ligand==='Slit2'&&p.receptor==='Dcc'){assert.equal(p.binding_status,'indirect');assert.equal(connection(p),'none');}
 if(p.receptor==='Plxna1')assert.match(p.molecular_form,/SLIT-C/);
 if(p.receptor==='Ptprz1')assert.match(p.molecular_form,/phosphacan/);
}
const nr=payload.pairRoles.filter(p=>[p.ligand,p.receptor].includes('Nrxn1'));
assert.equal(nr.length,2);assert(nr.every(p=>connection(p)==='both'));
mode='synaptic';assert.equal(payload.pairRoles.filter(p=>connection(p)!=='none').length,2);
assert.deepEqual(nr.map(connection).sort(),['forward','reverse']);
assert(nr.every(p=>p.synaptic_source==='CA1_Superficial'&&p.synaptic_target==='IN_C1'));
mode='resource';assert(payload.pairRoles.every(p=>connection(p)==='forward'));
const lines=fs.readFileSync(path.join(root,'results/binding_audit_14_rows.tsv'),'utf8').trim().split('\n');
assert.equal(lines.length,15);assert.match(lines[0],/CA1_gene/);assert.match(lines[0],/IN_compartment/);
console.log('Synaptic binding audit: frozen RNA preserved; all 14 rows covered; unsupported arrows suppressed; reciprocal synaptic roles corrected.');
