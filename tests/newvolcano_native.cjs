/* Exact existing CA1 total rows, independently counted; no browser/network/fitting. */
'use strict';
const assert = require('node:assert/strict'), fs = require('node:fs'), path = require('node:path'), vm = require('node:vm'), zlib = require('node:zlib');
const root = path.resolve(__dirname, '..'), site = path.join(root, 'site');
const catalog = JSON.parse(fs.readFileSync(path.join(site, 'data/DE_CATALOG.json')));
const remote = 'https://raw.githubusercontent.com/Ali-Mujahid/satb2_snRNA/atlas-v2.1.0/data_downloads/';
function payload(id) {const desc = catalog.shards[id], file = desc.path.startsWith(remote) ? path.join(root, 'data_downloads', desc.path.slice(remote.length)) : path.join(site, desc.path); return JSON.parse(zlib.gunzipSync(fs.readFileSync(file)));}
const sandbox = {window: {}}; vm.runInNewContext(fs.readFileSync(path.join(site, 'volcano_controls.js'), 'utf8'), sandbox);
class Canvas {
  constructor() {this.font = '10px system-ui'; this.colors = {}; this.alpha = {}; this.stack = [];}
  save() {this.stack.push({fillStyle: this.fillStyle, globalAlpha: this.globalAlpha, font: this.font});}
  restore() {Object.assign(this, this.stack.pop());}
  measureText(s) {return {width: String(s).length * parseFloat(this.font) * .55};}
  fill() {this.colors[this.fillStyle] = (this.colors[this.fillStyle] || 0) + 1; this.alpha[this.fillStyle] = this.globalAlpha;}
}
for (const key of ['clearRect','fillRect','strokeRect','beginPath','moveTo','lineTo','stroke','arc','rect','clip','translate','rotate','setLineDash','fillText']) Canvas.prototype[key] = function () {};
const checks = [];
for (const tool of ['scanpy', 'seurat']) {
  const dataset = catalog.datasets.find(x => x.tool === tool && x.branch === 'original_full' && x.population === 'CA1_ProS');
  const fit = payload(dataset.fits_shard).fits.find(x => x.analysis_kind === 'total' && (x.endpoint || 'PRIMARY') === 'PRIMARY');
  const rows = payload(fit.rows_shard).rows, before = JSON.stringify(rows);
  const finite = rows.filter(x => Number.isFinite(x.padj) && x.padj >= 0 && x.padj <= 1 && Number.isFinite(x.log2FoldChange));
  const nonDE = finite.filter(x => x.padj > .05);
  const above = finite.filter(x => -Math.log10(Math.max(x.padj, Number.MIN_VALUE)) > 5);
  const c = new Canvas(), meta = sandbox.window.SATB2_VOLCANO.draw({c, width: 448, height: 800}, rows, {deFdr: .05, deLfc: 0, volcanoView: 'both', volcanoLowYMax: 5, volcanoNS: '#cfd8dc'});
  assert.equal(finite.length, 16255); assert.equal(nonDE.length, 14230); assert.equal(meta.eligibleRows, 2025);
  assert.equal(meta.counts.nonsignificant, 14230); assert.equal(meta.counts.up + meta.counts.down, 2025); assert.equal(meta.omittedRows, 0);
  assert.equal(meta.axisBounds.y[1], 200); assert.equal(meta.lowPanel.axisBounds.y[1], 5);
  assert.equal(meta.lowPanel.visibleCounts.nonsignificant, 14230, 'Every finite non-DEG is retained in expanded view');
  assert.equal(meta.lowPanel.aboveViewRows, above.length); assert.equal(meta.lowPanel.visibleRows + above.length, finite.length);
  assert.equal(meta.labelCount, 20); assert.equal(c.colors['#cfd8dc'], nonDE.length * 2); assert.equal(c.alpha['#cfd8dc'], .9);
  assert.equal(JSON.stringify(rows), before);
  checks.push({tool, total_finite: finite.length, DEGs: meta.eligibleRows, non_DEGs: nonDE.length, low_view_visible: meta.lowPanel.visibleRows, above_low_view_only: above.length, NA_omitted: meta.omittedRows, non_DEG_opacity: meta.nonDEGOpacity, non_DEG_radius: meta.nonDEGSize});
}
console.log(JSON.stringify({status: 'PASS', scope: 'UNCHANGED_NATIVE_CA1_PRIMARY_TOTAL_ROWS_RENDERER_ONLY', checks}));
