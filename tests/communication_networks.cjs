'use strict';
const assert=require('assert');
const N=require('../site/communication/networks.js');
const checks=[];
function check(name,fn){fn();checks.push(name);console.log('PASS '+name);}
const decode=s=>s.replace(/&quot;/g,'"').replace(/&apos;/g,"'").replace(/&lt;/g,'<').replace(/&gt;/g,'>').replace(/&amp;/g,'&');
const metadata=svg=>JSON.parse(decode(svg.match(/<metadata>([\s\S]*?)<\/metadata>/)[1]));
const attrs=(svg,cls)=>[...svg.matchAll(new RegExp('<(?:g|path) class="'+cls+'" ([^>]+)>','g'))].map(m=>Object.fromEntries([...m[1].matchAll(/([\w-]+)="([^"]*)"/g)].map(a=>[a[1],decode(a[2])])));
const row=(source,target,ligand,receptor,rank,gm,lp=.5,rp=.25)=>({source,target,ligand_complex:ligand,receptor_complex:receptor,magnitude_rank:rank,lr_gmeans:gm,ligand_all_subunit_min_prop:lp,receptor_all_subunit_min_prop:rp});
const rows=[row('CA1_Deep','IN_C1','Slit1','Robo1',.01,1.2),row('CA1_Deep','IN_C1','Slit3','Robo1',.02,.9),row('IN_C1','CA1_Deep','Cntn3','Ptprg',.005,1.7),row('CA1_Superficial','IN_C1','Ntng1','Lrrc4c',.03,1.9),row('CA1_Deep','IN_C1','Slit1','Robo1',.01,1.2)];
const nodes=[{id:'CA1_Deep',label:'CA1 deep',color:'#2166ac'},{id:'CA1_Superficial',label:'CA1 superficial',color:'#b35806'},{id:'IN_C1',label:'C1 Pvalb',color:'#5e3c99'},{id:'IN_C20',label:'C20 retained review group',color:'#1b9e77'}];
const eligiblePairs=[{source:'CA1_Deep',target:'IN_C1'},{source:'CA1_Superficial',target:'IN_C1'},{source:'IN_C1',target:'CA1_Deep'},{source:'IN_C20',target:'CA1_Deep'}];
const ctx={nodes,eligiblePairs,condition:'WT',scope:'IN_clusters',direction:'bidirectional',options:{geneEdges:12}};
const before=JSON.stringify({rows,ctx});
check('Ordered pairs count distinct LR identities and count incident degree correctly',()=>{
  const s=N.summary(rows,ctx);assert.equal(s.total,4);assert.equal(s.pairs.length,3);assert.equal(s.pairs.find(p=>p.source==='CA1_Deep').count,2);assert.deepStrictEqual(s.degree,{CA1_Deep:3,IN_C1:4,CA1_Superficial:1});
});
check('Eligible pairs constrain the summarized native rows without inventing untested links',()=>{
  const s=N.summary(rows,{eligiblePairs:[{source:'IN_C1',target:'CA1_Deep'}]});assert.equal(s.total,1);assert.equal(s.pairs.length,1);assert.equal(s.pairs[0].source,'IN_C1');
  assert.equal(N.summary(rows,{eligiblePairs:[]}).total,0);
});
check('Circle keeps all selected groups and uses node area proportional to incident degree',()=>{
  const svg=N.circle(rows,ctx),ns=attrs(svg,'ccc-node');assert.equal(ns.length,nodes.length);assert.equal(ns.find(n=>n['data-node']==='IN_C20')['data-degree'],'0');assert(svg.includes('Hollow node = no qualifying edges'));
  const a=Number(ns.find(n=>n['data-node']==='CA1_Superficial')['data-radius']),b=Number(ns.find(n=>n['data-node']==='IN_C1')['data-radius']);assert(Math.abs(a*a/(b*b)-1/4)<.0002);
  assert.equal(attrs(svg,'ccc-edge').length,3);assert.equal((svg.match(/class="ccc-arrow"/g)||[]).length,3);assert.equal(metadata(svg).qualifying_population_lr_records,4);
});
check('Shared node order/edge/degree maxima preserve comparability across conditions',()=>{
  const a=N.circle(rows,{...ctx,edgeMax:5,degreeMax:10}),b=N.circle(rows.slice(0,1),{...ctx,condition:'cKO',edgeMax:5,degreeMax:10});
  const positions=svg=>[...svg.matchAll(/class="ccc-node"[^>]*data-node="([^"]*)"[^>]*>[\s\S]*?<circle cx="([^"]*)" cy="([^"]*)"/g)].map(m=>m.slice(1,4));
  assert.deepStrictEqual(positions(a),positions(b));assert.equal(metadata(a).edge_max,5);assert.equal(metadata(b).degree_max,10);
});
check('Circle styles change geometry and colors while preserving counts',()=>{
  const base=N.circle(rows,ctx),styled=N.circle(rows,{...ctx,options:{arrowSize:20,curvature:.8,edgeWidth:13,nodeScale:1.5,edgeColorMode:'custom',edgeColor:'#123456',nodeColorMode:'custom',nodeColor:'#654321'}});
  assert(styled.includes('#123456'));assert(styled.includes('#654321'));assert.notEqual(base,styled);assert.equal(metadata(base).qualifying_population_lr_records,metadata(styled).qualifying_population_lr_records);
  assert.deepStrictEqual(attrs(base,'ccc-edge').map(n=>n['data-count']),attrs(styled,'ccc-edge').map(n=>n['data-count']));
  assert.equal((N.circle(rows,{...ctx,options:{arrowSize:0}}).match(/class="ccc-arrow"/g)||[]).length,0);
});
check('Circos uses closed finite-width ribbons and retains zero-edge sectors',()=>{
  const svg=N.circos(rows,ctx),ribbons=attrs(svg,'ccc-ribbon'),arcs=attrs(svg,'ccc-sector');assert.equal(ribbons.length,3);assert.equal(arcs.length,nodes.length);assert(arcs.some(n=>n['data-node']==='IN_C20'&&n['data-degree']==='0'));
  for(const r of ribbons){assert(/ Z$/.test(r.d));assert(r.d.includes(' A'));assert(r.d.includes(' L'));assert(Number(r['data-source-fraction'])>0);assert(Number(r['data-target-fraction'])>0);}
  const sourceFractions=ribbons.filter(n=>n['data-source']==='CA1_Deep').reduce((s,n)=>s+Number(n['data-source-fraction']),0);assert(Math.abs(sourceFractions-2/3)<1e-12);assert(svg.includes('not an absolute comparison across sectors'));
});
check('Circos arc thickness/gap, ribbon scale and receiver arrow controls are independent display changes',()=>{
  const base=N.circos(rows,ctx),styled=N.circos(rows,{...ctx,options:{arcWidth:40,arcGap:8,ribbonScale:.3,arrowSize:22,curvature:.7,arcColorMode:'custom',arcColor:'#aa11bb',edgeColorMode:'target'}});
  assert(styled.includes('#aa11bb'));assert.notEqual(attrs(base,'ccc-sector')[0].d,attrs(styled,'ccc-sector')[0].d);assert.notEqual(attrs(base,'ccc-ribbon')[0].d,attrs(styled,'ccc-ribbon')[0].d);
  assert.equal(metadata(styled).qualifying_population_lr_records,4);assert.equal(attrs(styled,'ccc-sector').length,nodes.length);
});
check('Metric colors support custom, reversed, manual and shared numeric scales',()=>{
  const custom=N.circle(rows,{...ctx,options:{edgeColorMode:'metric',palette:'custom',lowColor:'#000000',midColor:'#777777',highColor:'#ffffff',scaleMode:'manual',colorMin:1,colorMax:2}});
  assert(custom.includes('stroke="#000000"'));assert(custom.includes('stroke="#ffffff"'));
  const reversed=N.circle(rows,{...ctx,options:{edgeColorMode:'metric',palette:'custom',lowColor:'#000000',midColor:'#777777',highColor:'#ffffff',reverseScale:true},colorDomain:[1,2]});
  assert.notEqual(custom.match(/class="ccc-edge"[^>]*>[\s\S]*?stroke="([^"]*)"/)[1],reversed.match(/class="ccc-edge"[^>]*>[\s\S]*?stroke="([^"]*)"/)[1]);
  assert(custom.includes('Edge color: LR count'));
});
check('Gene selection counts population×LR keys and deduplicates repeated native records',()=>{
  const svg=N.geneCircle(rows,ctx),m=metadata(svg);assert.equal(m.selected_keys,4);assert.equal(m.present_edges,4);assert.equal(attrs(svg,'ccc-gene-edge').length,4);assert(svg.includes('repeated molecular pairs in different groups are separate edges'));
  assert.equal(N.geneKey(rows[0]),JSON.stringify(['CA1_Deep','IN_C1','Slit1','Robo1']));
  const top=N.geneCircle(rows,{...ctx,options:{geneEdges:1}});assert.equal(metadata(top).selected_keys,1);assert.equal(attrs(top,'ccc-gene-edge')[0]['data-key'],N.geneKey(rows[2]));
});
check('Gene nodes carry true all-subunit fractions with separate ligand and receptor roles',()=>{
  const svg=N.geneCircle(rows,ctx),ns=attrs(svg,'ccc-gene-node');assert(ns.some(n=>n['data-role']==='ligand'&&Number(n['data-prop'])===.5));assert(ns.some(n=>n['data-role']==='receptor'&&Number(n['data-prop'])===.25));
  assert(svg.includes('equal-area square = receptor'));assert(svg.includes('lr_gmeans'));assert(!svg.includes('lr_logfc (cKO'));
  const almost=[row('CA1_Deep','IN_C1','A','B',.01,.8,1+1e-13,0)];assert(!N.geneCircle(almost,ctx).includes('data-prop="1.0000000000001"'));
});
check('Fixed gene keys preserve nodes and positions when a condition lacks a qualifying edge',()=>{
  const keys=rows.slice(0,4).map(N.geneKey),a=N.geneCircle(rows,{...ctx,geneKeys:keys}),b=N.geneCircle([],{...ctx,condition:'cKO',geneKeys:keys});
  const positions=svg=>attrs(svg,'ccc-gene-node').map(n=>[n['data-node'],n['data-x'],n['data-y']]);assert.deepStrictEqual(positions(a),positions(b));assert.equal(metadata(b).present_edges,0);assert.equal(metadata(b).selected_keys,4);assert(attrs(b,'ccc-gene-node').every(n=>n['data-prop']==='NA'));assert(b.includes('not zero expression'));assert(b.includes('Hollow dashed node = NA'));
});
check('Empty, singleton, degenerate, self-edge, escaped labels and extreme controls remain finite',()=>{
  const variants=[[],[row('X','X','A','B',.01,0,0,1)],[row('X','Y','A&B','<C>',.01,.8,1,.9)]];
  for(const fn of [N.circle,N.circos,N.geneCircle])for(const input of variants){const svg=fn(input,{condition:'WT & cKO',nodes:[{id:'X',label:'X <test> & stable'},{id:'Y',label:'Y "quoted"'}],options:{arcGap:45,arcWidth:60,ribbonScale:2,nodeRadius:48,nodeScale:4,arrowSize:35,curvature:1,fontSize:30}});assert(!/(?:NaN|Infinity|undefined)/.test(svg));assert(svg.startsWith('<svg '));assert(svg.endsWith('</svg>'));assert(svg.includes('WT &amp; cKO'));if(!input.length)assert(svg.includes('No qualifying ligand–receptor records'));}
});
check('Native input records, selected groups and caller settings are never mutated',()=>{assert.equal(JSON.stringify({rows,ctx}),before);});
const fs=require('fs'),path=require('path'),vm=require('vm'),zlib=require('zlib');
const packedPath=path.join(__dirname,'../site/communication/data/data.js');
if(fs.existsSync(packedPath)){
  const sandbox={window:{}};vm.runInNewContext(fs.readFileSync(packedPath,'utf8'),sandbox);
  const packed=JSON.parse(zlib.gunzipSync(Buffer.from(sandbox.window.CCC_DATA_GZIP,'base64')));
  require('../site/communication/model.js');const model=global.CCCModel,data=model.decode(packed),ds=data.datasets['IN_clusters|pooled'];
  const state={maxRank:.05,direction:'both',sources:[],targets:[],search:'',tier:'replicated_specific',bothDEG:false};
  const selected=model.filter(ds.rows,state),nativeNodes=data.nodes.filter(n=>n.id!=='IN_Total');
  const nativeCtx={nodes:nativeNodes,eligiblePairs:model.eligible(ds.eligible_pairs,state),condition:'pooled',scope:'IN_clusters',direction:'both'};
  const saved=JSON.stringify(selected);
  check('Packaged native strong tier reproduces 324 population-specific records and 112/212 directions',()=>{
    assert.equal(selected.length,324);assert.equal(selected.filter(r=>r.source.startsWith('CA1_')).length,112);assert.equal(selected.filter(r=>r.source.startsWith('IN_')).length,212);
    assert.equal(N.summary(selected,nativeCtx).total,324);
  });
  check('Actual 22-group SVGs render finite geometry with exact native counts and immutable rows',()=>{
    const a=N.circle(selected,nativeCtx),b=N.circos(selected,nativeCtx),c=N.geneCircle(selected,nativeCtx);
    assert.equal(attrs(a,'ccc-node').length,22);assert.equal(attrs(b,'ccc-sector').length,22);assert.equal(metadata(a).qualifying_population_lr_records,324);assert.equal(metadata(b).qualifying_population_lr_records,324);
    for(const svg of [a,b,c]){assert(!/(?:NaN|Infinity|undefined)/.test(svg));assert(svg.includes('Transcript')||svg.includes('transcript'));}
    assert.equal(JSON.stringify(selected),saved);assert.equal(metadata(c).present_edges,12);
  });
  check('Native total-IN strict empty result remains an honest empty chart with retained selected groups',()=>{
    const total=data.datasets['IN_total|pooled'],filtered=model.filter(total.rows,state),context={...nativeCtx,nodes:data.nodes.filter(n=>n.id.startsWith('CA1_')||n.id==='IN_Total'),eligiblePairs:total.eligible_pairs};assert.equal(filtered.length,0);
    assert.equal(attrs(N.circle(filtered,context),'ccc-node').length,3);assert.equal(attrs(N.circos(filtered,context),'ccc-sector').length,3);assert(N.geneCircle(filtered,context).includes('No qualifying ligand–receptor records'));
  });
}
console.log(JSON.stringify({status:'PASS',checks:checks.length,scope:'DISPLAY_ONLY_COMMUNICATION_NETWORK_RENDERERS',check_names:checks},null,2));
