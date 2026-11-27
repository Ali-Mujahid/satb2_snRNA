const assert=require('assert'),fs=require('fs'),path=require('path'),zlib=require('zlib'),crypto=require('crypto');
const root=path.join(__dirname,'..'),dir=path.join(root,'analysis/ca1_volcano_lr_size_20261127'),L=JSON.parse(fs.readFileSync(path.join(dir,'lr_data.json'),'utf8')),A=require('../site/ca1_layer_volcano/app.js');
assert.equal(L.rows.length,14);assert.equal(crypto.createHash('sha256').update(fs.readFileSync(path.join(root,L.source_path))).digest('hex'),L.source_sha256);
for(const [key,n] of Object.entries({superficial:7,deep:5})){
 const fit=JSON.parse(zlib.gunzipSync(fs.readFileSync(path.join(root,'analysis/ca1_layer_volcano_20261127',key+'.json.gz')))),before=JSON.stringify(fit),state={...A.defaults,fit:key},c=A.build(fit,state,L.rows);
 assert.equal(c.lrByGene.size,n);assert.equal(c.ops.filter(o=>o.role==='lr-outline').length,n);assert.equal(c.pairs.length,key==='superficial'?9:5);assert(c.lrVisible.every(r=>r.padj<.05));assert(c.pairs.every(r=>Number(r.IN_FDR)<.05));
 const targets=L.rows.filter(r=>r.layer===key);for(const r of targets){const native=c.rows.find(x=>x.gene_id===r.gene_id);assert.equal(native.gene_symbol,r.CA1_gene);assert.equal(native.log2FoldChange,r.CA1_log2FC_exact);assert.equal(native.padj,r.CA1_FDR_exact);}
 const small=A.build(fit,{...state,width:480,height:380,'show-title':false,'show-notes':false,labels:false},L.rows),large=A.build(fit,{...state,width:1200,height:850},L.rows);
 for(const chart of [small,large]){assert.equal(chart.visible.length,c.visible.length);assert(chart.ops.filter(o=>o.type==='circle'&&o.status).every(o=>o.r===state['dot-size']));assert(chart.ops.filter(o=>o.role==='gene-label').every(o=>o.size===state['label-size']));assert.equal(chart.ops.filter(o=>o.role==='lr-outline').length,n);}
 const smallSvg=A.svg(small);assert.match(smallSvg,/width="127mm"/);assert.match(smallSvg,/data-role="lr-outline"/);
 const id=targets[0].gene_id,chosen=A.build(fit,{...state,labels:false,'lr-label-genes':[id]},L.rows);assert.equal(chosen.labelRows.length,1);assert.equal(chosen.labelRows[0].gene_id,id);assert.equal(chosen.lrVisible.length,n);
 const clear=A.build(fit,{...state,labels:false,'lr-label-genes':[]},L.rows);assert.equal(clear.labelRows.length,0);assert.equal(clear.lrVisible.length,n);
 const direct=A.build(fit,{...state,'lr-evidence':'direct'},L.rows);assert(direct.pairs.every(r=>r.binding_status==='direct'));assert(!direct.pairs.some(r=>r.receptor==='Grm7'));
 const noRows=A.build(fit,{...state,'lr-rank':0},L.rows);assert.equal(noRows.lrByGene.size,0);assert.equal(noRows.pairs.length,0);
 const inOnly=A.build(fit,{...state,'lr-direction':'in_to_ca1'},L.rows);assert(inOnly.pairs.every(r=>r.source.startsWith('IN_')));
 assert.equal(JSON.stringify(fit),before);console.log(key,`${n} CA1 LR genes, ${c.labelRows.length}/${c.requested} default labels; small ${small.labelRows.length}/${small.requested}`);
}
assert.equal(L.rows.filter(r=>r.binding_status==='unverified').length,2);assert(L.rows.filter(r=>r.receptor==='Grm7').every(r=>r.IN_gene==='Grm7'&&r.binding_status==='unverified'));assert(L.rows.filter(r=>r.CA1_gene==='Nrxn1').every(r=>r.synaptic_direction.includes('CA1 superficial')));
console.log('PASS archived LR mapping, CA1-only highlights, evidence flags, filters, gene selection, fixed dot/font sizes, physical SVG sizing, source immutability.');
if(!process.argv.includes('--browser'))process.exit(0);
(async()=>{
 const {chromium}=require(process.env.PLAYWRIGHT_MODULE||'/Users/mujahid357/.cache/codex-runtimes/codex-primary-runtime/dependencies/node/node_modules/playwright'),http=require('http');
 const model=path.join(root,'analysis/ca1_layer_volcano_20261127'),M=JSON.parse(fs.readFileSync(path.join(model,'PROVENANCE.json'),'utf8'));
 const server=http.createServer((req,res)=>{const p=path.join(root,'site',decodeURIComponent(req.url.split('?')[0]));try{const f=fs.statSync(p).isDirectory()?path.join(p,'index.html'):p;res.setHeader('Content-Type',f.endsWith('.js')?'application/javascript':f.endsWith('.css')?'text/css':'text/html');res.end(fs.readFileSync(f));}catch(e){res.statusCode=404;res.end('missing');}});await new Promise(ok=>server.listen(0,'127.0.0.1',ok));
 const browser=await chromium.launch({headless:true,executablePath:'/Applications/Google Chrome.app/Contents/MacOS/Google Chrome'}),page=await browser.newPage({viewport:{width:1500,height:1000},acceptDownloads:true}),errors=[];page.on('pageerror',e=>errors.push(e.message));
 await page.route('https://raw.githubusercontent.com/Ali-Mujahid/satb2_snRNA/'+M.tag+'/**',route=>route.fulfill({body:fs.readFileSync(path.join(model,path.basename(new URL(route.request().url()).pathname))),contentType:'application/gzip'}));
 const out=path.join(dir,'browser_checks');fs.mkdirSync(out,{recursive:true});const download=async(kind,suffix)=>{const w=page.waitForEvent('download');await page.locator(`[data-export="${kind}"]`).click();const d=await w;await d.saveAs(path.join(out,suffix||d.suggestedFilename()));};
 try{
 await page.goto(`http://127.0.0.1:${server.address().port}/ca1_layer_volcano/index.html`);await page.waitForFunction(()=>window.CA1_LAYER_VOLCANO?.ready&&window.CA1_LAYER_VOLCANO.current);
 for(const [key,n]of Object.entries({superficial:7,deep:5})){
  await page.locator('#fit').selectOption(key);await page.waitForFunction(k=>window.CA1_LAYER_VOLCANO?.state.fit===k&&window.CA1_LAYER_VOLCANO.current,key);
  assert.equal(await page.locator('[data-role="lr-outline"]').count(),n);assert.equal(await page.locator('input[data-lr-gene]').count(),n);
  await page.locator('#figure').screenshot({path:path.join(out,key+'_LR_default.png')});await download('svg',key+'_LR_default.svg');await download('pairs',key+'_LR_pairs.tsv');
  if(key==='superficial'){assert.match(await page.locator('#lr-pair-table').innerText(),/Direct binding unverified/);await page.locator('#lr-evidence').selectOption('direct');await page.waitForTimeout(200);assert(!/Direct binding unverified/.test(await page.locator('#lr-pair-table').innerText()));await page.locator('#lr-evidence').selectOption('all');await page.waitForTimeout(200);}
  await page.locator('#lr-only').click();await page.locator('#lr-select-none').click();assert.equal(await page.locator('[data-role="gene-label"]').count(),0);await page.locator('input[data-lr-gene]').first().check();assert.equal(await page.locator('[data-role="gene-label"]').count(),1);
  const radiiBefore=await page.locator('circle[data-status]').evaluateAll(nodes=>[...new Set(nodes.map(n=>n.getAttribute('r')))]);
  await page.locator('[data-size="small"]').click();assert.equal(await page.locator('#dot-size').inputValue(),'2.4');assert.deepEqual(await page.locator('circle[data-status]').evaluateAll(nodes=>[...new Set(nodes.map(n=>n.getAttribute('r')))]),radiiBefore);assert.equal(await page.locator('#label-size').inputValue(),'15');
  await page.locator('#lr-select-all').click();await page.locator('#figure').screenshot({path:path.join(out,key+'_small.png')});for(const kind of ['pdf','svg','png','state'])await download(kind,key+'_small.'+({pdf:'pdf',svg:'svg',png:'png',state:'json'}[kind]));
  await page.locator('#lr-focus').click();assert.equal(await page.locator('[data-role="lr-outline"]').count(),n);await page.locator('#figure').screenshot({path:path.join(out,key+'_LR_focus.png')});await download('pdf',key+'_LR_focus.pdf');await download('svg',key+'_LR_focus.svg');assert.deepEqual(await page.locator('circle[data-status]').evaluateAll(nodes=>[...new Set(nodes.map(n=>n.getAttribute('r')))]),radiiBefore);
  await page.locator('#width-mm').fill('150');await page.waitForTimeout(250);assert.equal(await page.locator('#width').inputValue(),'567');assert.deepEqual(await page.locator('circle[data-status]').evaluateAll(nodes=>[...new Set(nodes.map(n=>n.getAttribute('r')))]),radiiBefore);
  await page.locator('#state-file').setInputFiles(path.join(out,key+'_small.json'));await page.waitForFunction(()=>document.getElementById('width').value==='480');assert.equal(await page.locator('[data-role="lr-outline"]').count(),n);
  await page.locator('[data-size="poster"]').click();assert.equal(await page.locator('#dot-size').inputValue(),'2.4');assert.equal(await page.locator('#label-size').inputValue(),'15');await download('pdf',key+'_poster.pdf');
  await page.locator('#reset').click();
 }
 assert.deepEqual(errors,[]);fs.writeFileSync(path.join(out,'ACCEPTANCE.json'),JSON.stringify({status:'PASS',checks:['CA1-only exact LR highlights','unverified partner flags','direct-evidence filter','label selection','small/presentation/poster reflow','physical size input','dot/font invariance','SVG/PDF/PNG/pairs/state exports','selection and size state roundtrip'],page_errors:errors},null,2));console.log('PASS browser LR choices, flags, dimensions, invariant marks, exports and settings.');
 }finally{await browser.close();server.close();}
})().catch(e=>{console.error(e);process.exit(1)});
