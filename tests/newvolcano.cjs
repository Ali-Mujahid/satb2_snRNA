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
  record(method, args) {for (const value of args) if (typeof value === 'number') assert.ok(Number.isFinite(value), `${method} received ${value}`); this.calls.push({method, args: [...args], fillStyle: this.fillStyle, align: this.textAlign, font: this.font});}
}
for (const name of ['clearRect', 'fillRect', 'strokeRect', 'beginPath', 'moveTo', 'lineTo', 'stroke', 'fill', 'arc', 'rect', 'clip', 'translate', 'rotate', 'setLineDash', 'fillText']) RecordingCanvas.prototype[name] = function (...args) {this.record(name, args);};
const row = (id, p, lfc, gene = id) => ({gene_id: id, gene, padj: p, log2FoldChange: lfc});
const rows = Array.from({length: 26}, (_, i) => row('ID' + String(i).padStart(2, '0'), 10 ** (-(26 - i)), i % 2 ? -2 : 2));
rows.push(row('ZERO', 0, 4), row('NONSIG', 0.8, 10), row('NA', null, 7), row('BADP', 2, 1), row('BADNEG', -0.1, 1), row('BADFC', 0.01, Infinity));
const settings = {deFdr: 0.05, deLfc: 1, volcanoLabels: true, volcanoLabelCount: 20, volcanoUp: '#ee1122', volcanoDown: '#22aaee', volcanoNS: '#777777', volcanoThresholds: true};
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
console.log('Volcano renderer: numeric axes, 20-label ranking, stable-ID ties, duplicate symbols, NA/invalid/zero adjusted P, clipping, collision spacing, colors and crowded-canvas handling passed.');
