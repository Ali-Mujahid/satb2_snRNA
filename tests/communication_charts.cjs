'use strict';
const assert = require('node:assert/strict');
const fs = require('node:fs');
const vm = require('node:vm');
const path = require('node:path');
const sandbox = { window: {} };
vm.runInNewContext(fs.readFileSync(path.join(__dirname, '../site/communication/charts.js'), 'utf8'), sandbox);
const plots = sandbox.window.CCCPlots;
const ctx = { nodes: [{ id: 'CA1_Deep', label: 'CA1 Deep', color: '#123456' }, { id: 'IN_C6', label: 'Sst C6', color: '#654321' }], eligiblePairs: [{ source: 'CA1_Deep', target: 'IN_C6' }, { source: 'IN_C6', target: 'CA1_Deep' }], condition: 'WT', scope: 'IN_clusters', options: { topN: 24, pairLimit: 20, heatMetric: 'count', dotMetric: 'magnitude_rank' } };
const rows = [
  { source: 'CA1_Deep', target: 'IN_C6', ligand_complex: 'Slit1', receptor_complex: 'Robo1', magnitude_rank: .01, lr_gmeans: 2, lr_prop: .25, lr_logfc: -.2, cellphone_pvals: 0, deg_supported: true, sample_support_WT: 6, sample_support_cKO: 6, sample_eligible_WT: 6, sample_eligible_cKO: 6 },
  { source: 'CA1_Deep', target: 'IN_C6', ligand_complex: 'A_B', receptor_complex: 'C_D', magnitude_rank: .02, lr_gmeans: 4, lr_prop: 1 + 1.2e-13, lr_logfc: .4, cellphone_pvals: .005, deg_supported: true }
];
const original = JSON.stringify(rows);
const aggregate = plots.helpers.aggregateHeatmap(rows.concat(rows[0]), ctx);
assert.equal(aggregate.find(c => c.source === 'CA1_Deep' && c.target === 'IN_C6').count, 2, 'Molecular identities must be counted once per population pair.');
assert.equal(aggregate.find(c => c.source === 'CA1_Deep' && c.target === 'IN_C6').mean_gmeans, 3, 'Mean geometric means use unique molecular identities.');
assert.equal(aggregate.find(c => c.source === 'IN_C6' && c.target === 'CA1_Deep').count, 0, 'Tested empty population pairs have zero candidate counts.');
assert.equal(aggregate.find(c => c.source === 'CA1_Deep' && c.target === 'CA1_Deep').count, null, 'Untested directions remain NA.');
assert.equal(aggregate.find(c => c.source === 'IN_C6' && c.target === 'CA1_Deep').mean_gmeans, null, 'A tested empty pair has no mean, not a zero mean.');
const heat = plots.heatmap(rows, ctx);
assert.match(heat, /data-tested="false" data-value=""/);
assert.match(heat, /data-tested="true" data-value="0"/);
assert.match(heat, /Receiver \/ target/);
assert.match(heat, /no interaction-level genotype test/);
const dot = plots.dotplot(rows, ctx);
assert.match(dot, /data-fraction="0.25"[^>]+r="4.5000"/, 'Radius sqrt(fraction) makes dot area proportional to the true fraction.');
assert.match(dot, /data-fraction="1"[^>]+r="9.0000"/, 'Floating detection values just above 1 are clamped for display only.');
assert.match(dot, /not a genotype FDR/);
assert.match(dot, /WT: 6 \/ 6; cKO: 6 \/ 6/);
assert.equal(JSON.stringify(rows), original, 'Rendering must never mutate saved evidence.');
const comparison = Object.assign({}, ctx, { lrOrder: ['MissingLigand→MissingReceptor', 'Slit1→Robo1'], pairOrder: ['IN_C6→CA1_Deep', 'CA1_Deep→IN_C6'], colorDomain: [0, .05] });
const selection = plots.helpers.selectDotplot(rows, comparison);
assert.equal(selection.lrs[0].key, 'MissingLigand→MissingReceptor', 'Compared plots retain absent LR identities for aligned rows.');
assert.equal(selection.pairs[0].key, 'IN_C6→CA1_Deep', 'Compared plots retain eligible empty pair columns.');
assert.match(plots.dotplot(rows, comparison), /MissingLigand/);
assert.match(plots.dotplot(rows, comparison), /shared across compared contexts/);
const labelContext = Object.assign({}, ctx, { nodes: ctx.nodes.map(n => Object.assign({}, n, { label: n.id === 'CA1_Deep' ? 'CA1 DL' : 'SST label edit' })) });
assert.match(plots.heatmap(rows, labelContext), />CA1 DL<\//, 'Heatmap axes respect the selected label mode.');
assert.match(plots.dotplot(rows, labelContext), />CA1 DL → SST label edit<\//, 'Dotplot population labels respect the selected label mode.');
assert.match(plots.dotplot(rows, Object.assign({}, ctx, { options: { fontSize: 8, dotRadius: 16 } })), /data-fraction="1"[^>]+r="16.0000"/, 'Spacing grows to honor the complete dot-size slider range.');
assert.match(plots.dotplot(rows, Object.assign({}, ctx, { options: { fontSize: 8 } })), /font-size:8px/, 'The minimum font slider value is honored.');
assert.equal(plots.color(0, 0, 1, { palette: 'custom', lowColor: '#000000', midColor: '#888888', highColor: '#ffffff' }), '#000000');
assert.equal(plots.color(0, 0, 1, { palette: 'custom', lowColor: '#000000', midColor: '#888888', highColor: '#ffffff', reverseScale: true }), '#ffffff');
assert.match(plots.dotplot([], ctx), /No qualifying interactions/);
assert.match(plots.dotplot([Object.assign({}, rows[0], { ligand_complex: '<script>&"' })], ctx), /&lt;script&gt;&amp;&quot;/);
for (const svg of [heat, dot, plots.heatmap([], ctx), plots.dotplot(rows, comparison), plots.dotplot(rows, Object.assign({}, ctx, { options: { dotMetric: 'lr_logfc', scaleMode: 'manual', colorMin: -1, colorMax: 1, fontSize: 16, dotRadius: 15 } }))]) {
  assert.ok(svg.startsWith('<svg ')); assert.ok(svg.endsWith('</svg>'));
  assert.doesNotMatch(svg, /(?:NaN|Infinity|undefined)/, 'Exports contain finite geometry and valid label text.');
}
console.log('Communication heatmap/dotplot tests passed: unique LR counts, eligible zero versus NA, true fraction areas, shared orders/scales, immutable evidence and SVG export semantics.');
