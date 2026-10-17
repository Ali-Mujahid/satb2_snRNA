/* Publication display controls: no leader lines and optional clean PDF text. */
'use strict';
const assert = require('node:assert/strict');
const fs = require('node:fs');
const path = require('node:path');
const vm = require('node:vm');
const sandbox = {window: {}};
vm.runInNewContext(fs.readFileSync(path.join(__dirname, '../site/volcano_controls.js'), 'utf8'), sandbox);
const {draw} = sandbox.window.SATB2_VOLCANO;
class Canvas {
  constructor() {this.calls = []; this.font = '10px system-ui'; this.stack = [];}
  save() {this.stack.push({font: this.font});}
  restore() {Object.assign(this, this.stack.pop());}
  measureText(text) {return {width: String(text).length * parseFloat(this.font) * .55};}
}
for (const method of ['clearRect', 'fillRect', 'strokeRect', 'beginPath', 'moveTo', 'lineTo', 'stroke', 'fill', 'arc', 'rect', 'clip', 'translate', 'rotate', 'setLineDash', 'fillText']) Canvas.prototype[method] = function (...args) {this.calls.push({method, args});};
const rows = [
  {gene: 'HighUp', gene_id: 'G1', padj: 1e-100, log2FoldChange: 3},
  {gene: 'HighDown', gene_id: 'G2', padj: 1e-80, log2FoldChange: -3},
  {gene: 'LowUp', gene_id: 'G3', padj: .01, log2FoldChange: 2},
  {gene: 'NotDE', gene_id: 'G4', padj: .8, log2FoldChange: .1},
  {gene: 'Unavailable', gene_id: 'G5', padj: null, log2FoldChange: 2}
];
const before = JSON.stringify(rows);
const base = {deFdr: .05, deLfc: 0, volcanoLabels: true, volcanoLabelCount: 20, volcanoView: 'full'};
for (const view of ['full', 'low', 'both']) {
  const screen = new Canvas(), clean = new Canvas(), unlabeled = new Canvas();
  const settings = {...base, volcanoView: view};
  const shape = {width: 700, height: view === 'both' ? 800 : 450};
  const screenResult = draw({c: screen, ...shape}, rows, settings);
  const cleanResult = draw({c: clean, ...shape}, rows, settings, {summaryText: false});
  draw({c: unlabeled, ...shape}, rows, {...settings, volcanoLabels: false}, {summaryText: false});
  assert.equal(cleanResult.summaryText, false);
  assert.equal(screenResult.summaryText, true);
  assert.equal(cleanResult.labelLeaderLines, false);
  assert.equal(clean.calls.filter(call => call.method === 'stroke').length, unlabeled.calls.filter(call => call.method === 'stroke').length, 'Gene annotations never add connecting strokes');
  const texts = clean.calls.filter(call => call.method === 'fillText').map(call => String(call.args[0]));
  assert.ok(texts.includes('log2 fold change (cKO / WT)'));
  assert.ok(texts.includes('−log10 adjusted P'));
  assert.ok(texts.includes(view === 'low' ? 'LowUp' : 'HighUp'));
  assert.ok(!texts.some(text => /Full range|Low-y zoom|shown|gene labels|Higher in cKO|Lower in cKO|Below thresholds|Non-DEG|points retained|omitted|NA\/invalid|above view/.test(text)), 'Clean PDF canvas excludes status, legend and provenance captions');
  assert.equal(clean.calls.filter(call => call.method === 'arc').length, screen.calls.filter(call => call.method === 'arc').length, 'Clean rendering retains all visible source points');
  for (const field of ['totalRows', 'validRows', 'omittedRows', 'eligibleGenes', 'labelCount']) assert.equal(cleanResult[field], screenResult[field], field + ' remains unchanged');
  assert.equal(JSON.stringify(cleanResult.counts), JSON.stringify(screenResult.counts));
  assert.ok(cleanResult.note.includes('unavailable FDR'), 'Audit information remains available as metadata');
  for (const label of cleanResult.labels) {
    assert.ok(label.box.left >= cleanResult.plotBounds.left && label.box.right <= cleanResult.plotBounds.right);
    assert.ok(label.box.top >= cleanResult.plotBounds.top && label.box.bottom <= cleanResult.plotBounds.bottom);
  }
  const bySide = [0, 1].map(side => cleanResult.labels.filter(label => label.side === side).sort((a, b) => a.y - b.y));
  for (const lane of bySide) for (let i = 1; i < lane.length; i++) assert.ok(lane[i].box.top > lane[i - 1].box.bottom);
}
assert.equal(JSON.stringify(rows), before, 'Presentation options never alter source rows');
console.log('Volcano publication display passed: full/low/both views, no leader strokes, clean canvas text, scientific axes, gene labels, points, audit metadata, collision spacing and source preservation.');
