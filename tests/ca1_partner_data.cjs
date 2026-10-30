process.chdir(require('path').join(__dirname, '../site/communication/ca1_partner'));
/* Structural data checks. Browser rendering and actual PDF saving are separate QA. */
const assert=require('assert'),fs=require('fs'),vm=require('vm');
const M=require('../site/communication/ca1_partner/model.js');let w={};vm.runInNewContext(fs.readFileSync('data.js','utf8'),{window:w});const D=w.SATB2_CA1_IN_DATA;
const html=fs.readFileSync('index.html','utf8'), appSource=fs.readFileSync('app.js');
const expectedHash=require('crypto').createHash('sha256').update(appSource).digest('hex').slice(0,12);
const appScript=html.match(/<script src="(app\.js\?v=[a-f0-9]+)"><\/script>/);
assert.ok(appScript);assert.equal(appScript[1],`app.js?v=${expectedHash}`);assert.ok(fs.existsSync(appScript[1].split('?')[0]));
assert.equal(D.metadata.state,'reviewed_complete');assert.equal(D.pairs.length,1364);assert.equal(D.subunits.length,2875);assert.equal(D.mouse_expression.length,34500);
const state={threshold:'5',direction:'',ca1:'',inpop:'',category:'',majority:false,replicated:false,query:''};
assert.equal(M.filter(D.pairs,state).length,1364);assert.equal(M.filter(D.pairs,{...state,threshold:'10'}).length,937);
assert.equal(M.filter(D.pairs,{...state,replicated:true}).length,1330);assert.equal(M.filter(D.pairs,{...state,threshold:'10',replicated:true}).length,891);
assert.equal(M.filter(D.pairs,{...state,direction:'CA1 to inhibitory'}).length,675);assert.equal(M.filter(D.pairs,{...state,direction:'inhibitory to CA1'}).length,689);
assert.equal(M.filter(D.pairs,{...state,category:'IN differential status unavailable'}).length,363);
const idset=new Set(D.pairs.map(r=>r.record_id));assert.equal(idset.size,1364);assert.ok(D.subunits.every(r=>idset.has(r.record_id)));assert.ok(D.mouse_expression.every(r=>idset.has(r.record_id)));
assert.ok(D.subunits.some(r=>r.padj===''));assert.ok(D.pairs.every(r=>['CA1_Deep','CA1_Superficial'].includes(r.ca1_population)&&['IN_C1','IN_C6','IN_C8'].includes(r.in_population)));
console.log('actual data integrity and filter counts passed; no real-browser/PDF claim');
