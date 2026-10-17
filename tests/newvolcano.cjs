/* No dependencies: record native Canvas 2D calls and inspect geometry and ranking. */
'use strict';
const assert = require('node:assert/strict');
const fs = require('node:fs');
const vm = require('node:vm');
const path = require('node:path');
const context = {window: {}};
vm.runInNewContext(fs.readFileSync(path.join(__dirname, '../site/volcano_controls.js'), 'utf8'), context);
const {draw, prepare} = context.window.SATB2_VOLCANO;
class RecordingCanvas {
  constructor() {this.calls = []; this.fillStyle = ''; this.font = '10px system-ui'; this.textAlign = 'left'; this.globalAlpha = 1; this.stack = [];}
  save() {this.stack.push({fillStyle: this.fillStyle, font: this.font, textAlign: this.textAlign, globalAlpha: this.globalAlpha});}
  restore() {Object.assign(this, this.stack.pop());}
  measureText(text) {return {width: String(text).length * parseFloat(this.font) * 0.55};}
  record(method, args) {for (const value of args) if (typeof value === 'number') assert.ok(Number.isFinite(value), `${method} received ${value}`); this.calls.push({method, args: [...args], fillStyle: this.fillStyle, align: this.textAlign, font: this.font, alpha: this.globalAlpha});}
}
for (const name of ['clearRect', 'fillRect', 'strokeRect', 'beginPath', 'moveTo', 'lineTo', 'stroke', 'fill', 'arc', 'rect', 'clip', 'translate', 'rotate', 'setLineDash', 'fillText']) RecordingCanvas.prototype[name] = function (...args) {this.record(name, args);};
const row = (id, p, lfc, gene = id) => ({gene_id: id, gene, padj: p, log2FoldChange: lfc});
const rows = Array.from({length: 26}, (_, i) => row('ID' + String(i).padStart(2, '0'), 10 ** (-(26 - i)), i % 2 ? -2 : 2));
rows.push(row('ZERO', 0, 4), row('NONSIG', 0.8, 10), row('NA', null, 7), row('BADP', 2, 1), row('BADNEG', -0.1, 1), row('BADFC', 0.01, Infinity));
const settings = {deFdr: 0.05, deLfc: 1, volcanoLabels: true, volcanoLabelCount: 20, volcanoUp: '#ee1122', volcanoDown: '#22aaee', volcanoNS: '#777777', volcanoThresholds: true, volcanoView: 'full'};
const originalRows = JSON.stringify(rows);
const canvas = new RecordingCanvas(), result = draw({c: canvas, width: 900, height: 450}, rows, settings);
assert.equal(result.labelCount, 20); assert.equal(result.eligibleGenes, 26); assert.equal(result.validRows, 28); assert.equal(result.omittedRows, 4);
assert.equal(result.unavailableFdrRows, 1); assert.equal(result.invalidFdrRows, 2); assert.equal(result.invalidEffectRows, 1); assert.equal(result.zeroAdjustedPRows, 1); assert.equal(result.cappedRows, 1);
assert.equal(result.axisBounds.y[1], 320); assert.equal(result.omittedLabelsForSpace, 0);
assert.deepEqual(Array.from(result.labels, x => x.gene_id).sort(), ['ZERO', ...Array.from({length: 19}, (_, i) => 'ID' + String(i).padStart(2, '0'))].sort());
for (const label of result.labels) {assert.ok(label.box.left >= result.plotBounds.left); assert.ok(label.box.right <= result.plotBounds.right); assert.ok(label.box.top >= result.plotBounds.top); assert.ok(label.box.bottom <= result.plotBounds.bottom);}
for (const side of [0, 1]) {const labels = result.labels.filter(x => x.side === side).sort((a, b) => a.y - b.y); for (let i = 1; i < labels.length; i++) assert.ok(labels[i].box.top > labels[i - 1].box.bottom, 'Labels in the same lane do not collide');}
assert.ok(canvas.calls.some(x => x.method === 'fillText' && x.args[0] === 'log2 fold change (cKO / WT)'));
assert.ok(canvas.calls.some(x => x.method === 'fillText' && x.args[0] === '−log10 adjusted P'));
for (const color of ['#ee1122', '#22aaee', '#777777']) assert.ok(canvas.calls.some(x => x.method === 'fill' && x.fillStyle === color));
assert.ok(result.ticks.x.length >= 5 && result.ticks.y.length >= 4);
const tieRows = [row('B', 0.001, 2, 'Same'), row('A', 0.001, -3, 'Same'), row('C', 0.001, 2), row('A', 0.001, -3, 'Same')];
const tied = prepare(tieRows, {...settings, volcanoLabelCount: 20});
assert.deepEqual(Array.from(tied.labels, x => x.id), ['A', 'B', 'C']);
assert.equal(tied.labels[0].label, 'Same [A]'); assert.equal(tied.labels[1].label, 'Same [B]');
assert.equal(prepare([row('LOW', 0.01, 0.2)], settings).labels.length, 0, 'Never pad top labels with nonsignificant rows');
assert.equal(draw({c: new RecordingCanvas(), width: 900, height: 450}, rows, {...settings, volcanoLabels: false}).labelCount, 0);
const positive = Array.from({length: 20}, (_, i) => row('P' + i, 10 ** (-(i + 2)), 3));
assert.equal(draw({c: new RecordingCanvas(), width: 900, height: 340}, positive, settings).labelCount, 20, 'All-positive top20 fit by spilling to the other lane');
const crowded = draw({c: new RecordingCanvas(), width: 300, height: 240}, positive, {...settings, volcanoLabelCount: 100});
assert.ok(crowded.omittedLabelsForSpace > 0); assert.equal(crowded.labelCount + crowded.omittedLabelsForSpace, 20);
const empty = draw({c: new RecordingCanvas(), width: 900, height: 450}, [row('NA', null, 1)], settings);
assert.equal(empty.validRows, 0); assert.equal(empty.labelCount, 0);
assert.equal(draw({c: new RecordingCanvas(), width: 900, height: 450}, rows, {...settings, volcanoUp: 'garbage'}).colors.up, '#d18a36');
assert.equal(prepare([row('ZERO', 0, 1), row('SMALL', 1e-40, 1)], {...settings, deFdr: 0}).labels.length, 1, 'Zero FDR threshold includes only zero adjusted P');
assert.equal(prepare([row('EXACT', 0.05, 1)], settings).labels.length, 1, 'Displayed threshold boundaries are inclusive');
assert.equal(JSON.stringify(rows), originalRows, 'All source rows remain unchanged');
// A full-strength near-baseline fixture models the published CA1 issue: many
// finite non-DEGs, extreme significant points and unavailable FDR are distinct.
const lowRows = [row('NS_ZERO', 1, 0), row('NS_SMALL', 0.8, 0.2), row('NS_NEAR', 0.051, -1), row('UP_LOW', 0.001, 2), row('DOWN_LOW', 0.01, -2), row('UP_HIGH', 1e-198, 3), row('ZERO_P', 0, -4), row('NA_FDR', null, 2)];
const lowBefore = JSON.stringify(lowRows), lowCanvas = new RecordingCanvas();
const both = draw({c: lowCanvas, width: 448, height: 800}, lowRows, {...settings, volcanoView: 'both', volcanoNS: '#cfd8dc'});
assert.equal(both.view, 'both'); assert.equal(both.lowYMax, 5); assert.equal(both.validRows, 7); assert.equal(both.omittedRows, 1);
assert.deepEqual(JSON.parse(JSON.stringify(both.counts)), {up: 2, down: 2, nonsignificant: 3});
assert.equal(both.lowPanel.visibleRows, 5); assert.equal(both.lowPanel.aboveViewRows, 2);
assert.deepEqual(Array.from(both.axisBounds.y), [0, 320]); assert.deepEqual(Array.from(both.lowPanel.axisBounds.y), [0, 5]);
assert.deepEqual(Array.from(both.axisBounds.x), Array.from(both.lowPanel.axisBounds.x));
assert.ok(both.plotBounds.bottom < both.lowPanel.plotBounds.top, 'Full and zoom plots occupy separate canvas regions');
assert.equal(both.labelCount, 4, 'Overview retains every eligible top label'); assert.equal(both.lowPanel.labelCount, 0);
assert.equal(both.nonDEGOpacity, 0.9); assert.equal(both.nonDEGSize, 1.8);
assert.equal(lowCanvas.calls.filter(x => x.method === 'arc').length, 12, 'Seven full points plus five visible zoom points, without clamped high-y points');
assert.equal(lowCanvas.calls.filter(x => x.method === 'fill' && x.fillStyle === '#cfd8dc').length, 6, 'All three non-DEGs appear in both panels');
assert.ok(lowCanvas.calls.filter(x => x.method === 'fill' && x.fillStyle === '#cfd8dc').every(x => x.alpha === 0.9));
assert.ok(lowCanvas.calls.some(x => x.method === 'arc' && x.args[1] === both.plotBounds.bottom && x.args[2] === 1.8), 'True y=0 centers stay exactly on baseline');
assert.ok(both.note.includes('2 finite points above 5') && both.note.includes('not NA'));
const lowOnlyCanvas = new RecordingCanvas(), lowOnly = draw({c: lowOnlyCanvas, width: 448, height: 450}, lowRows, {...settings, volcanoView: 'low'});
assert.equal(lowOnly.visibleRows, 5); assert.equal(lowOnly.aboveViewRows, 2); assert.equal(lowOnlyCanvas.calls.filter(x => x.method === 'arc').length, 5);
assert.ok(!lowOnlyCanvas.calls.some(x => x.method === 'arc' && x.args[1] === lowOnly.plotBounds.top), 'Above-cutoff points are omitted instead of moved to the top border');
assert.equal(lowOnly.labelCount, 2, 'Only visible members of the global top-label set appear in the zoom');
for (const [request, expected] of [[-10, .1], [0, .1], [400, 320]]) assert.equal(draw({c: new RecordingCanvas(), width: 448, height: 450}, lowRows, {...settings, volcanoView: 'low', volcanoLowYMax: request}).lowYMax, expected);
const styled = draw({c: new RecordingCanvas(), width: 448, height: 450}, lowRows, {...settings, volcanoView: 'low', volcanoNSOpacity: 7, volcanoNSSize: .1});
assert.equal(styled.nonDEGOpacity, 1); assert.equal(styled.nonDEGSize, .6);
assert.equal(draw({c: new RecordingCanvas(), width: 448, height: 800}, lowRows, {...settings, volcanoView: undefined}).view, 'both', 'Combined view is the default');
assert.equal(JSON.stringify(lowRows), lowBefore, 'Zoom and stronger non-DEG styling preserve every input record');
const pdfCanvas = new RecordingCanvas(), repeated = draw({c: pdfCanvas, width: 448, height: 800}, lowRows, {...settings, volcanoView: 'both', volcanoNS: '#cfd8dc'});
assert.equal(JSON.stringify(repeated), JSON.stringify(both), 'The same single canvas draw used for PDF reproduces overview/zoom metadata exactly');
assert.equal(JSON.stringify(pdfCanvas.calls), JSON.stringify(lowCanvas.calls), 'Screen and PDF share the same complete two-panel drawing');
console.log('Volcano renderer: numeric axes, 20-label ranking, stable-ID ties, duplicate symbols, NA/invalid/zero adjusted P, clipping, collision spacing, colors and crowded-canvas handling passed.');
console.log('Full/low/combined views: true zero coordinates, stronger visible non-DEGs, exact low cutoff, above-view counts, shared x range, independent panels and identical PDF drawing passed.');
