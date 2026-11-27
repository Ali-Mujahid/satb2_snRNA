const assert=require('assert'),fs=require('fs'),path=require('path'),zlib=require('zlib'),crypto=require('crypto');
const root=path.join(__dirname,'..'),base=path.join(root,'analysis/ca1_layer_volcano_20261127'),A=require('../site/ca1_layer_volcano/app.js');
const manifest=JSON.parse(fs.readFileSync(path.join(base,'PROVENANCE.json'),'utf8')),fits={};
for(const [key,counts]of Object.entries({superficial:[1775,938,837,14206],deep:[742,360,382,13439]})){
 const packed=fs.readFileSync(path.join(base,key+'.json.gz')),raw=zlib.gunzipSync(packed),fit=JSON.parse(raw);fits[key]=fit;
 const hash=b=>crypto.createHash('sha256').update(b).digest('hex');assert.equal(hash(packed),manifest.fits[key].sha256);assert.equal(hash(raw),manifest.fits[key].payload_sha256);
 assert.equal(hash(zlib.gunzipSync(fs.readFileSync(path.join(base,key+'_source.tsv.gz')))),manifest.fits[key].source_sha256);
 const before=JSON.stringify(fit),c=A.build(fit,{...A.defaults,fit:key});assert.deepEqual([c.degs.length,c.up,c.down,c.rows.length],counts);assert.equal(c.labelRows.length,10);assert.equal(new Set(c.labelRows.map(r=>r.gene_id)).size,10);assert(!c.omitted);assert.equal(c.ops.filter(o=>o.role==='fc-guide').length,2);assert.equal(c.ops.filter(o=>o.role==='plot-boundary').length,2);
 assert.deepEqual(fit.metadata.clusters,key==='superficial'?[2,3]:[1,5]);assert.equal(fit.metadata.C8_included,false);assert.equal(fit.mice.length,12);
 const changed=A.build(fit,{...A.defaults,up:'#00aa11',down:'#1100aa',fc:1,ymax:50,xmin:-2,xmax:2,'guide-fc':.5,boundary:false});assert(changed.visible.every(r=>Math.abs(r.x)<=2&&r.y<=50));assert(changed.ops.filter(o=>o.gene_id&&o.type==='circle'&&o.status==='up').every(o=>o.fill==='#00aa11'));assert(changed.outside>0);assert(!changed.ops.some(o=>o.role==='plot-boundary'));assert.equal(changed.ops.filter(o=>o.role==='fc-guide').length,2);
 assert.throws(()=>A.build(fit,{...A.defaults,xmin:2,xmax:1}),/minimum/);assert.throws(()=>A.build(fit,{...A.defaults,xtick:.00001}),/80 marks/);assert.equal(JSON.stringify(fit),before);
 const hidden=A.build(fit,{...A.defaults,labels:false});assert.equal(hidden.labelRows.length,0);
 const svg=A.svg(c);assert(svg.includes('stroke-dasharray="2 4"'));assert.equal((svg.match(/data-role="gene-label"/g)||[]).length,10);
}
const synthetic={metadata:{label:'Synthetic'},columns:fits.deep.columns,rows:[['g1','Zero',1,2,1,1,0,0],['g2','Missing',1,-1,1,1,.01,null],['g3','Absent',1,1,1,1,null,null]]};
assert.equal(A.build(synthetic,{...A.defaults}).zeros,1);assert.equal(A.build(synthetic,{...A.defaults}).missing,2);assert.equal(A.build(synthetic,{...A.defaults,metric:'pvalue'}).visible.length,2);
console.log('PASS exact saved CA1 fits, 1775/742 DEGs, C8 excluded, top 10, axes, borders, guides, colors, missing/zero handling; no source mutation.');
if(!process.argv.includes('--browser'))process.exit(0);
(async()=>{
 const {chromium}=require(process.env.PLAYWRIGHT_MODULE||'/Users/mujahid357/.cache/codex-runtimes/codex-primary-runtime/dependencies/node/node_modules/playwright'),http=require('http');
 const server=http.createServer((req,res)=>{const p=path.join(root,'site',decodeURIComponent(req.url.split('?')[0]));try{const f=fs.statSync(p).isDirectory()?path.join(p,'index.html'):p;const ext=path.extname(f);res.setHeader('Content-Type',ext==='.js'?'application/javascript':ext==='.css'?'text/css':'text/html');res.end(fs.readFileSync(f));}catch(e){res.statusCode=404;res.end('missing');}});await new Promise(ok=>server.listen(0,'127.0.0.1',ok));
 const browser=await chromium.launch({headless:true,executablePath:'/Applications/Google Chrome.app/Contents/MacOS/Google Chrome'});const page=await browser.newPage({viewport:{width:1450,height:1000},acceptDownloads:true}),errors=[];page.on('pageerror',e=>errors.push(e.message));
 await page.route('https://raw.githubusercontent.com/Ali-Mujahid/satb2_snRNA/'+manifest.tag+'/**',route=>route.fulfill({body:fs.readFileSync(path.join(base,path.basename(new URL(route.request().url()).pathname))),contentType:'application/gzip'}));
 const out=path.join(base,'browser_checks');fs.mkdirSync(out,{recursive:true});
 try{
  await page.goto(`http://127.0.0.1:${server.address().port}/ca1_layer_volcano/index.html`);await page.waitForFunction(()=>window.CA1_LAYER_VOLCANO?.ready&&window.CA1_LAYER_VOLCANO.current);
  for(const key of ['superficial','deep']){
   await page.locator('#fit').selectOption(key);await page.waitForFunction(k=>window.CA1_LAYER_VOLCANO?.state.fit===k&&window.CA1_LAYER_VOLCANO.ready&&window.CA1_LAYER_VOLCANO.current,key);
   await page.locator('#lr-labels').uncheck();await page.waitForTimeout(200);assert.equal(await page.locator('[data-role="gene-label"]').count(),10);
   await page.locator('#figure').screenshot({path:path.join(out,key+'_default.png')});
   for(const ext of ['svg','pdf','png','degs','all','labels','state']){const waiting=page.waitForEvent('download');await page.locator(`[data-export="${ext}"]`).click();const d=await waiting;await d.saveAs(path.join(out,d.suggestedFilename()));}
   await page.locator('#fc').fill('1');await page.locator('#guide-fc').fill('1.5');await page.locator('#ymax').fill('50');await page.locator('#up').fill('#8844aa');await page.locator('#boundary').uncheck();await page.waitForTimeout(250);assert.equal(await page.locator('[data-role="plot-boundary"]').count(),0);assert.equal((await page.evaluate(()=>window.CA1_LAYER_VOLCANO.current.axes.ymax)),50);assert.equal(await page.locator('[data-role="fc-guide"]').count(),2);
   await page.locator('#xmin').fill('20');await page.waitForTimeout(200);assert(await page.locator('[data-export="pdf"]').isDisabled());assert.match(await page.locator('#error').textContent(),/minimum/);
   await page.locator('#reset').click();await page.locator('#lr-labels').uncheck();await page.waitForTimeout(200);assert.equal(await page.locator('[data-role="gene-label"]').count(),10);
   const setting=path.join(out,`Seurat_CA1_${key==='superficial'?'Superficial_C2_C3':'Deep_C1_C5'}_volcano_settings.json`);await page.locator('#label-size').fill('20');await page.waitForTimeout(200);await page.locator('#state-file').setInputFiles(setting);await page.waitForFunction(()=>document.getElementById('label-size').value==='15');
  }
  assert.deepEqual(errors,[]);fs.writeFileSync(path.join(out,'ACCEPTANCE.json'),JSON.stringify({status:'PASS',checks:['two exact fits','default top 10','all seven downloads','axis limits and validation','fold change cutoff and guides','dot colour','plot boundary','state roundtrip'],page_errors:errors},null,2));console.log('PASS browser: both fits, labels, all exports, colors, axes, dotted guides, boundary, settings import.');
 }finally{await browser.close();server.close();}
})().catch(e=>{console.error(e);process.exit(1)});
