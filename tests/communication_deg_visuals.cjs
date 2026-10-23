'use strict';
const assert=require('node:assert/strict');
const fs=require('node:fs');
const path=require('node:path');
const vm=require('node:vm');
const sandbox={window:{}};
vm.runInNewContext(fs.readFileSync(path.join(__dirname,'../site/communication/charts.js'),'utf8'),sandbox);
const plots=sandbox.window.CCCPlots;
const networks=require('../site/communication/networks.js');
const rows=[
  {source:'IN_C6',target:'CA1_Deep',ligand_complex:'Cntn3',receptor_complex:'Ptprg',magnitude_rank:.01,lr_gmeans:2,lr_logfc:.8,lr_prop:.25,ligand_all_subunit_min_prop:.5,receptor_all_subunit_min_prop:.25,deg_supported:true,ligand_deg_direction:'down',receptor_deg_direction:'up',ligand_deseq2_log2fc:-.5383457365,receptor_deseq2_log2fc:.3966628983,ligand_deseq2_padj:.0037005639,receptor_deseq2_padj:.0164890846},
  {source:'CA1_Deep',target:'IN_C6',ligand_complex:'Slit1',receptor_complex:'Robo1',magnitude_rank:.02,lr_gmeans:1,lr_logfc:1.6,lr_prop:.4,ligand_all_subunit_min_prop:.4,receptor_all_subunit_min_prop:.7,deg_supported:true,ligand_deg_direction:'down',receptor_deg_direction:'down',ligand_deseq2_log2fc:-.30242,receptor_deseq2_log2fc:-.45822656,ligand_deseq2_padj:.0191792,receptor_deseq2_padj:.0000338377},
  {source:'CA1_Deep',target:'IN_C6',ligand_complex:'X',receptor_complex:'Y',magnitude_rank:.3,lr_gmeans:.7,lr_logfc:-2,lr_prop:.8,ligand_all_subunit_min_prop:.8,receptor_all_subunit_min_prop:.9,deg_supported:true,ligand_deg_direction:'up',receptor_deg_direction:'up',ligand_deseq2_log2fc:1.2,receptor_deseq2_log2fc:.6,ligand_deseq2_padj:.001,receptor_deseq2_padj:.01}
];
const ctx={condition:'pooled',scope:'IN_clusters',direction:'both',nodes:[{id:'IN_C6',label:'Sst C6',color:'#123456'},{id:'CA1_Deep',label:'CA1 Deep',color:'#654321'}],eligiblePairs:[{source:'IN_C6',target:'CA1_Deep'},{source:'CA1_Deep',target:'IN_C6'}],options:{showDEGDirections:true,degUpColor:'#aabbcc',degDownColor:'#ddeeff',geneEdges:40}};
const before=JSON.stringify({rows,ctx});
const json=x=>JSON.parse(JSON.stringify(x));
const counts=plots.helpers.directionCounts(rows.concat(rows[0]));
assert.deepEqual(json(counts),{'down→up':1,'down→down':1,'up→up':1},'Counts deduplicate identical population × LR records, without pooling endpoint signs.');
const heat=plots.heatmap(rows,ctx);
assert.match(heat,/Endpoint DEG directions, cKO vs WT/);
assert.match(heat,/L↓\/R↑: 1/);
assert.match(heat,/Directions describe endpoint genes, not net interaction activity/);
const cells=plots.helpers.aggregateHeatmap(rows,ctx);
assert.deepEqual(json(cells.find(c=>c.source==='CA1_Deep'&&c.target==='IN_C6').deg_direction_counts),{'down→down':1,'up→up':1});
const dot=plots.dotplot(rows,ctx);
assert.equal((dot.match(/class="ccc-deg-dot-rim"/g)||[]).length,6,'Each displayed LR has distinct ligand and receptor signed rims.');
assert.match(dot,/data-role="ligand" data-deg-direction="down"[^>]*stroke="#ddeeff"/);
assert.match(dot,/data-role="receptor" data-deg-direction="up"[^>]*stroke="#aabbcc"/);
assert.match(dot,/data-fraction="0.25"[^>]*r="4.5000"/,'The signed rim must preserve the true detection-area interior.');
assert.match(dot,/Sender ligand DESeq2: down in cKO vs WT; log2FC = -0.538; padj = 0.0037/);
assert.match(dot,/Receiver receptor DESeq2: up in cKO vs WT; log2FC = 0.397; padj = 0.0165/);
assert.match(dot,/LR group-versus-rest logFC: 0.8/,'Native LR specificity remains distinct from DESeq2 genotype effects.');
assert.match(dot,/Dot left rim = sender ligand DEG; right rim = receiver receptor DEG/);
const without={...ctx,options:{...ctx.options,showDEGDirections:false}};
assert.doesNotMatch(plots.dotplot(rows,without),/class="ccc-deg-dot-rim"/,'Direction overlay can be turned off without changing records.');
for(const fn of [networks.circle,networks.circos]){
  const svg=fn(rows.concat(rows[0]),ctx);
  assert.match(svg,/Endpoint DEGs, cKO vs WT: L↓\/R↓: 1 · L↓\/R↑: 1 · L↑\/R↑: 1/);
  assert.match(svg,/data-deg-counts=/,'Cell-group aggregates expose separate signed endpoint combination counts.');
  assert.match(svg,/not .*interaction|do not .*interaction|no.*interaction/);
  assert.doesNotMatch(svg,/(?:NaN|Infinity|undefined)/);
}
const gene=networks.geneCircle(rows,ctx);
assert.match(gene,/data-role="ligand"[^>]*data-deg-direction="down"[^>]*>[\s\S]*?fill="#123456" stroke="#ddeeff" stroke-width="3"/,'Population fill is retained while the ligand outline shows its DEG direction.');
assert.match(gene,/data-role="receptor"[^>]*data-deg-direction="up"[^>]*>[\s\S]*?stroke="#aabbcc" stroke-width="3"/);
assert.match(gene,/>Cntn3 ↓<\/text>/);
assert.match(gene,/>Ptprg ↑<\/text>/);
assert.match(gene,/Ligand DESeq2 cKO vs WT: down; log2FC -0.5383; padj 0.003701/);
assert.match(gene,/Gene directions do not define net interaction activity/);
assert.doesNotMatch(networks.geneCircle(rows,without),/class="ccc-deg-legend"/);
for(const fn of [networks.circle,networks.circos,networks.geneCircle])for(const fontSize of [8,12,20,30]){
  const svg=fn(rows,{...ctx,options:{...ctx.options,fontSize}});
  const labels=[...svg.matchAll(/<text class="ccc-(?:node-label|sector-label|gene-population-label)"[^>]* y="([^"]+)"[^>]*font-size="([^"]+)"/g)];
  assert(labels.length>0);
  assert(labels.every(m=>Number(m[1])-Number(m[2])>=105.9),'Population labels reserve vertical space below the endpoint DEG summary.');
}
const old=rows.map(r=>Object.fromEntries(Object.entries(r).filter(([k])=>!k.includes('_deg_direction')&&!k.includes('_deseq2_'))));
for(const fn of [plots.heatmap,plots.dotplot,networks.circle,networks.circos,networks.geneCircle]){
  const oldSvg=fn(old,ctx),signedSvg=fn(rows,ctx);
  assert.doesNotMatch(oldSvg,/data-deg-summary="true"/,'Legacy records without signed evidence do not acquire fabricated signs.');
  for(const svg of [oldSvg,signedSvg]){assert.ok(svg.startsWith('<svg '));assert.ok(svg.endsWith('</svg>'));assert.doesNotMatch(svg,/(?:NaN|Infinity|undefined)/);}
}
assert.equal(JSON.stringify({rows,ctx}),before,'Signed rendering cannot mutate saved biological records or controls.');
console.log('DEG visual tests passed: separate signed endpoint evidence, exact DESeq2 tooltip fields, deduplicated aggregate counts, customizable outlines/rims, immutable detection areas, old-data compatibility and export-safe SVGs.');
