const assert = require('node:assert/strict');
const fs = require('node:fs');
const zlib = require('node:zlib');

global.window = global;
global.SATB2_NO_AUTO_INIT = true;
require('../site/v2_features.js');
require('../site/presentation_layouts.js');
require('../site/data/c6_subclusters.js');
require('../site/app.js');

async function run() {
  const manifest = JSON.parse(fs.readFileSync('site/data/manifest.json'));
  manifest.defaults = {...manifest.defaults, population:'Inhibitory', branch:'original_full', annotation:'c6_program', tab:'umap'};
  const app = await SATB2_APP.init({
    manifest, headless:true, ignoreSavedState:true,
    shardLoader:async (id, desc)=>({encoding:'json', data:JSON.parse(zlib.gunzipSync(fs.readFileSync('site/'+desc.path)))})
  });
  const {coords, c6Groups} = app.active.seurat;
  assert.equal(coords.nuclei.length, 8911);
  assert.equal(c6Groups.size, 501);
  assert.equal(app.active.scanpy.c6Groups, null);
  const original = app.presentationCoordinates('seurat');
  assert.equal(original.clusters.length, 22);
  assert.ok(original.points.every((p,i)=>p.x===coords.nuclei[i].x_harmony && p.y===coords.nuclei[i].y_harmony));
  for(const [id,wt,cko] of [['6.0',147,144],['6.1',41,69],['6.2',26,74]]) {
    app.state.c6Focus = id;
    for(const [genotype,n] of [['WT',wt],['cKO',cko]]) {
      app.state.umapCondition = genotype;
      assert.equal(app.umapIndices('seurat').length,n);
    }
  }
  const oneMouse = coords.nuclei.find(n=>c6Groups.get(n.nucleus_id)==='6.2' && n.condition==='cKO').sample;
  app.state.c6Mouse = oneMouse;
  const selected = app.umapIndices('seurat');
  assert.ok(selected.length>0 && selected.every(i=>coords.nuclei[i].sample===oneMouse));
  app.state.layoutAdjustments[JSON.stringify(['seurat','original_full','Inhibitory/C6-split'])] = {'6.2':{dx:.1,dy:0,scale:1}};
  const moved = app.presentationCoordinates('seurat').points.filter(p=>p.x!==p.source_x);
  assert.equal(moved.length,100);
  assert.ok(moved.every(p=>p.cluster_id==='6.2'));
  app.state.annotation = 'cluster_id';
  const parent = app.presentationCoordinates('seurat');
  assert.equal(parent.clusters.length,20);
  assert.ok(parent.points.every(p=>p.x===p.source_x && p.y===p.source_y));
  assert.equal(coords.nuclei.filter(n=>String(n.cluster_id)==='6').length,501);
  console.log('PASS: verified C6 join, source positions, genotype and mouse filters, independent layouts and parent recovery');
}
run().catch(err=>{console.error(err);process.exitCode=1;});
