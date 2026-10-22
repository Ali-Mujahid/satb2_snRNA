/* Portable SVG communication views. Styling never changes the native LIANA records. */
(function (global) {
  'use strict';
  const TAU = Math.PI * 2, W = 1280, H = 1080, CX = W / 2, CY = 470;
  const DEFAULTS = Object.freeze({palette:'viridis', lowColor:'#440154', midColor:'#21918c', highColor:'#fde725', reverseScale:false, scaleMode:'auto', colorMin:0, colorMax:1, fontSize:12, opacity:.7, nodeRadius:16, nodeScale:1, edgeWidth:6, arrowSize:9, curvature:.35, arcWidth:16, arcGap:3, ribbonScale:1, edgeColorMode:'source', edgeColor:'#687b9c', arcColorMode:'node', arcColor:'#687b9c', nodeColorMode:'palette', nodeColor:'#4a84a7', geneEdges:12});
  const COLORS = ['#2166ac','#b35806','#5e3c99','#1b9e77','#d95f02','#7570b3','#e7298a','#66a61e','#a6761d','#1f78b4','#b2df8a','#e31a1c','#6a3d9a','#ff7f00','#33a02c','#fb9a99','#cab2d6','#b15928','#8dd3c7','#bc80bd','#80b1d3','#fdb462'];
  const esc = x => String(x == null ? '' : x).replace(/[&<>"']/g, c => ({'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;',"'":'&apos;'}[c]));
  const finite = x => typeof x === 'number' && Number.isFinite(x);
  const number = (x, fallback) => (x !== null && x !== '' && Number.isFinite(Number(x))) ? Number(x) : fallback;
  const clamp = (x, a, b) => Math.min(b, Math.max(a, x));
  const round = x => Number(x.toFixed(3));
  const point = (angle, radius) => ({x:CX + radius * Math.cos(angle), y:CY + radius * Math.sin(angle)});
  const xy = p => round(p.x) + ',' + round(p.y);
  const pairKey = (s,t) => JSON.stringify([String(s),String(t)]);
  const lrKey = row => JSON.stringify([String(row.ligand_complex),String(row.receptor_complex)]);
  const geneKey = row => JSON.stringify([String(row.source),String(row.target),String(row.ligand_complex),String(row.receptor_complex)]);
  const nodeKey = (population,complex,role) => JSON.stringify([population,complex,role]);
  function opt(ctx) {
    const o = {...DEFAULTS,...((ctx || {}).options || {})};
    o.fontSize = clamp(number(o.fontSize,12),8,30); o.opacity = clamp(number(o.opacity,.7),.02,1);
    o.nodeRadius = clamp(number(o.nodeRadius,16),2,48); o.nodeScale = clamp(number(o.nodeScale,1),.1,4);
    o.edgeWidth = clamp(number(o.edgeWidth,6),.1,35); o.arrowSize = clamp(number(o.arrowSize,9),0,35);
    o.curvature = clamp(number(o.curvature,.35),0,1); o.arcWidth = clamp(number(o.arcWidth,16),1,60);
    o.arcGap = clamp(number(o.arcGap,3),.1,45); o.ribbonScale = clamp(number(o.ribbonScale,1),.05,2);
    o.geneEdges = Math.round(clamp(number(o.geneEdges,12),1,100));
    return o;
  }
  function safeColor(value,fallback) { return typeof value === 'string' && /^(#[\da-f]{3,8}|[a-z]{3,25})$/i.test(value) ? value : fallback; }
  function colorDomain(values,ctx,o) {
    const domain = ctx && Array.isArray(ctx.colorDomain) ? ctx.colorDomain : (o.scaleMode === 'manual' ? [o.colorMin,o.colorMax] : null);
    let lo = domain ? number(domain[0],0) : Math.min(0,...values.filter(finite));
    let hi = domain ? number(domain[1],1) : Math.max(0,...values.filter(finite));
    if (lo > hi) [lo,hi] = [hi,lo];
    if (!(hi > lo)) hi = lo + 1;
    return [lo,hi];
  }
  function hexRgb(value) {
    const s = /^#[\da-f]{6}$/i.test(value) ? value : '#888888';
    return [parseInt(s.slice(1,3),16),parseInt(s.slice(3,5),16),parseInt(s.slice(5,7),16)];
  }
  function color(value,domain,o) {
    if (!finite(value)) return '#c7cdd2';
    if (global.CCCPlots && typeof global.CCCPlots.color === 'function') {
      try { const c = global.CCCPlots.color(value,domain[0],domain[1],o); if (typeof c === 'string') return safeColor(c,'#888888'); } catch (_) { /* Standalone fallback for offline exports. */ }
    }
    const palettes = {viridis:['#440154','#21918c','#fde725'],blue:['#eff3ff','#6baed6','#08519c'],rdBu:['#2166ac','#f7f7f7','#b2182b']};
    const colors = o.palette === 'custom' ? [safeColor(o.lowColor,DEFAULTS.lowColor),safeColor(o.midColor,DEFAULTS.midColor),safeColor(o.highColor,DEFAULTS.highColor)] : (palettes[o.palette] || palettes.viridis);
    let t = clamp((value-domain[0])/(domain[1]-domain[0]),0,1); if (o.reverseScale) t=1-t;
    const i = t <= .5 ? 0 : 1, f = i === 0 ? t*2 : t*2-1, a = hexRgb(colors[i]), b = hexRgb(colors[i+1]);
    return '#' + a.map((v,j)=>Math.round(v+(b[j]-v)*f).toString(16).padStart(2,'0')).join('');
  }
  function populations(rows,ctx,o) {
    const out=[], seen=new Set();
    const add = (id,label,c) => {id=String(id); if (seen.has(id)) return; seen.add(id); out.push({id,label:String(label || id),color:o.nodeColorMode==='custom' ? safeColor(o.nodeColor,DEFAULTS.nodeColor) : safeColor(c,COLORS[out.length % COLORS.length])});};
    for (const n of ((ctx || {}).nodes || [])) if (n && n.id != null) add(n.id,n.label,n.color);
    for (const row of rows) { if (row.source != null) add(row.source,row.source_label); if (row.target != null) add(row.target,row.target_label); }
    return out;
  }
  function aggregate(rows,ctx) {
    const eligible = ctx && Array.isArray(ctx.eligiblePairs) ? new Set(ctx.eligiblePairs.map(p=>pairKey(p.source,p.target))) : null;
    const map = new Map();
    for (const row of rows) {
      if (!row || row.source == null || row.target == null || row.ligand_complex == null || row.receptor_complex == null) continue;
      const key=pairKey(row.source,row.target); if (eligible && !eligible.has(key)) continue;
      if (!map.has(key)) map.set(key,{source:String(row.source),target:String(row.target),identities:new Set(),records:0});
      const p=map.get(key); p.identities.add(lrKey(row)); p.records++;
    }
    return Array.from(map.values()).map(p=>({source:p.source,target:p.target,count:p.identities.size,records:p.records})).sort((a,b)=>pairKey(a.source,a.target).localeCompare(pairKey(b.source,b.target)));
  }
  function summary(rows,ctx) {
    const pairs=aggregate(rows,ctx), degree={};
    for (const p of pairs) {degree[p.source]=(degree[p.source] || 0)+p.count; degree[p.target]=(degree[p.target] || 0)+p.count;}
    return {pairs,degree,total:pairs.reduce((s,p)=>s+p.count,0),maxEdge:Math.max(0,...pairs.map(p=>p.count)),maxDegree:Math.max(0,...Object.values(degree))};
  }
  function edgeColor(pair,lookup,domain,o,value) {
    if (o.edgeColorMode === 'custom') return safeColor(o.edgeColor,DEFAULTS.edgeColor);
    if (o.edgeColorMode === 'metric') return color(value,domain,o);
    const node=lookup.get(o.edgeColorMode === 'target' ? pair.target : pair.source);
    return node ? node.color : safeColor(o.edgeColor,DEFAULTS.edgeColor);
  }
  function metadata(kind,ctx,rows,o,extra) {
    return {view:kind,condition:ctx.condition || '',scope:ctx.scope || '',direction:ctx.direction || '',input_records:rows.length,display_only:true,options:o,...extra};
  }
  function start(kind,ctx,rows,o,extra={}) {
    const titles={circle:'Directed cell-group circle',circos:'Directed cell-group circos',geneCircle:'Ligand–receptor complex circle'};
    const subtitle=[ctx.condition,ctx.scope,ctx.direction].filter(Boolean).join(' · ');
    return '<svg xmlns="http://www.w3.org/2000/svg" xmlns:xlink="http://www.w3.org/1999/xlink" width="'+W+'" height="'+H+'" viewBox="0 0 '+W+' '+H+'" role="img" aria-label="'+esc(titles[kind]+' '+subtitle)+'" data-view="'+kind+'" data-records="'+rows.length+'">'+
      '<title>'+esc(titles[kind]+' — '+subtitle)+'</title><desc>Transcript and prior compatible candidates. Arrowheads indicate the listed ligand sender to receptor receiver. No genotype interaction significance test is implied.</desc>'+
      '<metadata>'+esc(JSON.stringify(metadata(kind,ctx,rows,o,extra)))+'</metadata><rect width="100%" height="100%" fill="#ffffff"/>'+
      '<g font-family="Arial,Helvetica,sans-serif" fill="#243b4b"><text x="32" y="35" font-size="22" font-weight="700">'+esc(titles[kind])+'</text>'+
      '<text x="32" y="61" font-size="13" fill="#526977">'+esc(subtitle)+'</text></g>';
  }
  function end(lines,o,domain,label) {
    let s='<g font-family="Arial,Helvetica,sans-serif" fill="#526977" font-size="'+Math.min(14,o.fontSize)+'">';
    lines.forEach((line,i)=>{s+='<text x="34" y="'+(945+i*22)+'">'+esc(line)+'</text>';});
    if (o.edgeColorMode==='metric') {
      s+='<g transform="translate(910,944)">';
      for (let i=0;i<80;i++) s+='<rect x="'+i*3+'" y="0" width="3.2" height="10" fill="'+esc(color(domain[0]+(domain[1]-domain[0])*i/79,domain,o))+'"/>';
      s+='<text x="0" y="27" font-size="11">'+esc(format(domain[0]))+'</text><text x="240" y="27" font-size="11" text-anchor="end">'+esc(format(domain[1]))+'</text><text x="0" y="45" font-size="11">'+esc(label)+'</text></g>';
    }
    return s+'</g></svg>';
  }
  function format(value) { return finite(value) ? (Math.abs(value)<.001 && value!==0 ? value.toExponential(2) : Number(value.toPrecision(4)).toString()) : 'NA'; }
  function textLabel(x,y,label,anchor,o,cls,tooltip) {
    // Keep the label footprint approximately constant when the font grows. The full
    // annotation remains in the native SVG tooltip even when its display is shortened.
    const max=Math.round(clamp(432/o.fontSize,12,36)), words=String(label).split(/\s+/), lines=[]; let line='';
    for (const word of words) { if (line.length+word.length+1>max && line) {lines.push(line);line=word;} else line+=(line?' ':'')+word; }
    if (line) lines.push(line);
    const display=lines.slice(0,2).map((s,i)=>(i===1 && lines.length>2)||s.length>max ? s.slice(0,max-1)+'…' : s);
    return '<text class="'+cls+'" x="'+round(x)+'" y="'+round(y)+'" text-anchor="'+anchor+'" font-family="Arial,Helvetica,sans-serif" font-size="'+o.fontSize+'" fill="#243b4b" stroke="#ffffff" stroke-width="3" paint-order="stroke" stroke-linejoin="round"><title>'+esc(tooltip || label)+'</title>'+display.map((s,i)=>'<tspan x="'+round(x)+'" dy="'+(i?o.fontSize*1.12:0)+'">'+esc(s)+'</tspan>').join('')+'</text>';
  }
  function emptyMessage(o,extra) {
    return '<text x="'+CX+'" y="'+CY+'" text-anchor="middle" font-family="Arial,Helvetica,sans-serif" font-size="'+Math.max(14,o.fontSize)+'" fill="#526977">No qualifying ligand–receptor records'+(extra?': '+esc(extra):'')+'.</text>';
  }
  function arrow(end,tangent,size,fill,opacity,attrs) {
    if (!(size>0)) return '';
    const length=Math.hypot(tangent.x,tangent.y) || 1, ux=tangent.x/length, uy=tangent.y/length, half=size*.48;
    const a={x:end.x-size*ux+half*uy,y:end.y-size*uy-half*ux},b={x:end.x-size*ux-half*uy,y:end.y-size*uy+half*ux};
    return '<polygon class="ccc-arrow" points="'+xy(end)+' '+xy(a)+' '+xy(b)+'" fill="'+esc(fill)+'" fill-opacity="'+opacity+'" '+(attrs||'')+'/>';
  }
  function curve(from,to,sourceRadius,targetRadius,o) {
    const dx=to.x-from.x,dy=to.y-from.y,d=Math.hypot(dx,dy) || 1,ux=dx/d,uy=dy/d;
    if (d<1e-7) {
      const r=Math.max(30,sourceRadius*2),a={x:from.x-sourceRadius*.65,y:from.y-sourceRadius*.65},b={x:from.x+sourceRadius*.65,y:from.y-sourceRadius*.65};
      const c1={x:from.x-r,y:from.y-r*2},c2={x:from.x+r,y:from.y-r*2};
      return {d:'M'+xy(a)+' C'+xy(c1)+' '+xy(c2)+' '+xy(b),end:b,tangent:{x:b.x-c2.x,y:b.y-c2.y}};
    }
    const a={x:from.x+sourceRadius*ux,y:from.y+sourceRadius*uy},b={x:to.x-targetRadius*ux,y:to.y-targetRadius*uy};
    const c={x:(a.x+b.x)/2-dy*o.curvature*.55,y:(a.y+b.y)/2+dx*o.curvature*.55};
    return {d:'M'+xy(a)+' Q'+xy(c)+' '+xy(b),end:b,tangent:{x:b.x-c.x,y:b.y-c.y}};
  }
  function circle(rows,ctx={}) {
    rows=Array.isArray(rows)?rows:[]; const o=opt(ctx),nodes=populations(rows,ctx,o),s=summary(rows,ctx),lookup=new Map(nodes.map(n=>[n.id,n]));
    const edgeMax=Math.max(1,number(ctx.edgeMax,s.maxEdge)),degreeMax=Math.max(1,number(ctx.degreeMax,s.maxDegree)),domain=colorDomain(s.pairs.map(p=>p.count),ctx,o);
    nodes.forEach((n,i)=>{n.angle=-Math.PI/2+TAU*i/Math.max(1,nodes.length);n.position=point(n.angle,325);n.degree=s.degree[n.id] || 0;n.radius=n.degree>0?o.nodeRadius*o.nodeScale*Math.sqrt(n.degree/degreeMax):3;});
    let svg=start('circle',ctx,rows,o,{qualifying_population_lr_records:s.total,pair_count:s.pairs.length,edge_max:edgeMax,degree_max:degreeMax,node_order:nodes.map(n=>n.id)});
    svg+='<g class="ccc-edges" data-edge-count="'+s.pairs.length+'">';
    for (const p of s.pairs.slice().sort((a,b)=>a.count-b.count)) {
      const from=lookup.get(p.source),to=lookup.get(p.target); if (!from || !to) continue;
      const fill=edgeColor(p,lookup,domain,o,p.count),width=o.edgeWidth*Math.sqrt(p.count/edgeMax),path=curve(from.position,to.position,from.radius+2,to.radius+3,o);
      const desc=p.source+' → '+p.target+'; '+p.count+' distinct ligand→receptor identities; listed sender→receiver, not inferred activation or inhibition.';
      svg+='<g class="ccc-edge" data-source="'+esc(p.source)+'" data-target="'+esc(p.target)+'" data-count="'+p.count+'"><title>'+esc(desc)+'</title><path d="'+path.d+'" fill="none" stroke="'+esc(fill)+'" stroke-width="'+round(width)+'" stroke-opacity="'+o.opacity+'"/>'+arrow(path.end,path.tangent,o.arrowSize,fill,o.opacity)+'</g>';
    }
    svg+='</g><g class="ccc-nodes">';
    for (const n of nodes) {
      svg+='<g class="ccc-node" data-node="'+esc(n.id)+'" data-degree="'+n.degree+'" data-radius="'+round(n.radius)+'"><title>'+esc(n.label+'; incident qualifying LR count '+n.degree+(n.degree?'':'; no qualifying edges in this view; group retained'))+'</title><circle cx="'+round(n.position.x)+'" cy="'+round(n.position.y)+'" r="'+round(n.radius)+'" fill="'+esc(n.degree?n.color:'#ffffff')+'" stroke="'+esc(n.color)+'" stroke-width="1.5"'+(n.degree?'':' stroke-dasharray="2 2"')+'/></g>';
      const p=point(n.angle,360+o.nodeRadius*o.nodeScale*.3),anchor=Math.abs(Math.cos(n.angle))<.15?'middle':Math.cos(n.angle)>0?'start':'end';
      svg+=textLabel(p.x,p.y,n.label,anchor,o,'ccc-node-label');
    }
    svg+='</g>'; if (!s.total) svg+=emptyMessage(o,'the selected cutoffs or endpoints');
    return svg+end(['Arrow: listed ligand sender → receptor receiver; width = distinct qualifying LR count (scale max '+edgeMax+').','Node area = incident LR count (scale max '+degreeMax+'); all selected groups retained. Hollow node = no qualifying edges.','Candidate scores describe transcript/prior compatibility. They do not test genotype interaction changes.'],o,domain,'Edge color: LR count');
  }
  function ringSegment(r0,r1,a,b) {
    const p1=point(a,r1),p2=point(b,r1),p3=point(b,r0),p4=point(a,r0),large=b-a>Math.PI?1:0;
    return 'M'+xy(p1)+' A'+round(r1)+','+round(r1)+' 0 '+large+' 1 '+xy(p2)+' L'+xy(p3)+' A'+round(r0)+','+round(r0)+' 0 '+large+' 0 '+xy(p4)+' Z';
  }
  function ribbon(sa,sb,ta,tb,r,o) {
    const a=point(sa,r),b=point(sb,r),c=point(ta,r),d=point(tb,r),f=clamp(1-o.curvature,.02,.98),inside=p=>({x:CX+(p.x-CX)*f,y:CY+(p.y-CY)*f});
    let path='M'+xy(a)+' A'+r+','+r+' 0 '+(sb-sa>Math.PI?1:0)+' 1 '+xy(b)+' C'+xy(inside(b))+' '+xy(inside(c))+' '+xy(c);
    if (o.arrowSize>0) path+=' L'+xy(point((ta+tb)/2,r+o.arrowSize))+' L'+xy(d);
    else path+=' A'+r+','+r+' 0 '+(tb-ta>Math.PI?1:0)+' 1 '+xy(d);
    return path+' C'+xy(inside(d))+' '+xy(inside(a))+' '+xy(a)+' Z';
  }
  function circos(rows,ctx={}) {
    rows=Array.isArray(rows)?rows:[];const o=opt(ctx),nodes=populations(rows,ctx,o),s=summary(rows,ctx),lookup=new Map(nodes.map(n=>[n.id,n])),r=298,step=TAU/Math.max(1,nodes.length),gap=Math.min(o.arcGap*Math.PI/180,step*.75),span=step-gap,domain=colorDomain(s.pairs.map(p=>p.count),ctx,o),arcInner=r+Math.max(12,o.arrowSize+4);
    nodes.forEach((n,i)=>{n.start=-Math.PI/2+i*step+gap/2;n.end=n.start+span;n.degree=s.degree[n.id] || 0;n.cursor=n.start+span*.05;});
    const links=s.pairs.map(p=>({...p}));
    for (const p of links) {
      for (const [field,id] of [['source',p.source],['target',p.target]]) {
        const n=lookup.get(id),block=span*.9*p.count/Math.max(1,n.degree),middle=n.cursor+block/2,half=block*o.ribbonScale/2;
        p[field+'Angles']=[middle-half,middle+half];n.cursor+=block;
      }
    }
    let svg=start('circos',ctx,rows,o,{qualifying_population_lr_records:s.total,pair_count:s.pairs.length,node_order:nodes.map(n=>n.id),sector_policy:'equal sectors; endpoint width is within-sector incident fraction multiplied by display ribbonScale'});
    svg+='<g class="ccc-ribbons" data-edge-count="'+links.length+'">';
    for (const p of links.slice().sort((a,b)=>a.count-b.count)) {
      const fill=edgeColor(p,lookup,domain,o,p.count),from=lookup.get(p.source),to=lookup.get(p.target),sf=p.count/from.degree,tf=p.count/to.degree;
      svg+='<path class="ccc-ribbon" data-source="'+esc(p.source)+'" data-target="'+esc(p.target)+'" data-count="'+p.count+'" data-source-fraction="'+sf+'" data-target-fraction="'+tf+'" d="'+ribbon(...p.sourceAngles,...p.targetAngles,r,o)+'" fill="'+esc(fill)+'" fill-opacity="'+o.opacity+'" stroke="'+esc(fill)+'" stroke-width=".4"><title>'+esc(p.source+' → '+p.target+'; '+p.count+' distinct qualifying LR identities; sender sector fraction '+format(sf)+'; receiver sector fraction '+format(tf)+'. Equal-sized sectors mean ribbon width is not an absolute comparison across sectors.')+'</title></path>';
    }
    svg+='</g><g class="ccc-arcs">';
    for (const n of nodes) {
      const fill=o.arcColorMode==='custom'?safeColor(o.arcColor,DEFAULTS.arcColor):n.color;
      svg+='<path class="ccc-sector" data-node="'+esc(n.id)+'" data-degree="'+n.degree+'" d="'+ringSegment(arcInner,arcInner+o.arcWidth,n.start,n.end)+'" fill="'+esc(fill)+'" fill-opacity="'+(n.degree?1:.3)+'" stroke="#ffffff" stroke-width="1"><title>'+esc(n.label+'; incident qualifying LR count '+n.degree+(n.degree?'':'; no qualifying edges; group retained'))+'</title></path>';
      const a=(n.start+n.end)/2,p=point(a,arcInner+o.arcWidth+24),anchor=Math.abs(Math.cos(a))<.15?'middle':Math.cos(a)>0?'start':'end';
      svg+=textLabel(p.x,p.y,n.label,anchor,o,'ccc-sector-label');
    }
    svg+='</g>';if (!s.total) svg+=emptyMessage(o,'the selected cutoffs or endpoints');
    return svg+end(['Arrow tip points to the receptor receiver. Sectors are equal sized; all selected groups retained.','Ribbon width = fraction of incident qualifying LR records within each endpoint sector × '+format(o.ribbonScale)+' display scale.','Use circle/heatmap counts for absolute comparisons across groups.'+(o.ribbonScale>1?' Width scale >1 can overlap adjacent ribbons.':''),'Transcript/prior candidates; no measured communication or genotype interaction significance is implied.'],o,domain,'Ribbon color: LR count');
  }
  function selectedGenes(rows,ctx,o) {
    const byKey=new Map();
    for (const row of rows) {if (!row || row.source==null || row.target==null || row.ligand_complex==null || row.receptor_complex==null) continue;const key=geneKey(row),old=byKey.get(key);if (!old || number(row.magnitude_rank,Infinity)<number(old.magnitude_rank,Infinity)) byKey.set(key,row);}
    if (Array.isArray(ctx.geneKeys)) return ctx.geneKeys.slice(0,o.geneEdges).map(key=>{try {const a=JSON.parse(key);return Array.isArray(a)&&a.length===4 ? {key,source:String(a[0]),target:String(a[1]),ligand_complex:String(a[2]),receptor_complex:String(a[3]),row:byKey.get(key) || null} : null;} catch (_) {return null;}}).filter(Boolean);
    return Array.from(byKey.entries()).sort((a,b)=>number(a[1].magnitude_rank,Infinity)-number(b[1].magnitude_rank,Infinity)||a[0].localeCompare(b[0])).slice(0,o.geneEdges).map(([key,row])=>({key,...row,row}));
  }
  function geneCircle(rows,ctx={}) {
    rows=Array.isArray(rows)?rows:[];const o=opt(ctx),pops=populations(rows,ctx,o),popMap=new Map(pops.map((p,i)=>[p.id,{...p,order:i}])),selected=selectedGenes(rows,ctx,o),nodeMap=new Map();
    function add(p,c,role,row) {
      const key=nodeKey(p,c,role);if (!nodeMap.has(key)) nodeMap.set(key,{id:key,population:p,complex:c,role,props:[],color:popMap.has(p)?popMap.get(p).color:COLORS[nodeMap.size%COLORS.length]});
      const n=nodeMap.get(key),prop=row && number(row[role==='ligand'?'ligand_all_subunit_min_prop':'receptor_all_subunit_min_prop'],null);if (finite(prop)) n.props.push(clamp(prop,0,1));
    }
    for (const e of selected) {add(e.source,e.ligand_complex,'ligand',e.row);add(e.target,e.receptor_complex,'receptor',e.row);}
    const nodes=Array.from(nodeMap.values()).sort((a,b)=>((popMap.get(a.population)||{}).order??999)-((popMap.get(b.population)||{}).order??999)||a.population.localeCompare(b.population)||a.role.localeCompare(b.role)||a.complex.localeCompare(b.complex));
    const groups=[];for (const n of nodes) {let g=groups[groups.length-1];if (!g||g.population!==n.population) {g={population:n.population,nodes:[],color:n.color};groups.push(g);}g.nodes.push(n);}
    const gap=Math.min(o.arcGap*Math.PI/180,.4),available=TAU-groups.length*gap,step=Math.max(.001,available/Math.max(1,nodes.length));let cursor=-Math.PI/2;
    for (const g of groups) {g.start=cursor;for (const n of g.nodes) {n.angle=cursor+step/2;n.position=point(n.angle,292);n.prop=n.props.length?Math.min(...n.props):null;n.radius=n.prop===null?4:o.nodeRadius*o.nodeScale*Math.sqrt(n.prop);n.radius=Math.max(n.prop>0?.6:0,n.radius);cursor+=step;}g.end=cursor;cursor+=gap;}
    const present=selected.filter(e=>e.row),weights=present.map(e=>number(e.row.lr_gmeans,null)).filter(finite),weightMax=Math.max(1e-12,number(ctx.edgeMax,Math.max(0,...weights))),domain=colorDomain(weights,ctx,o),lookup=new Map(pops.map(p=>[p.id,p]));
    let svg=start('geneCircle',ctx,rows,o,{selected_population_lr_keys:selected.map(e=>e.key),selected_keys:selected.length,present_edges:present.length,node_order:nodes.map(n=>n.id),node_detection_policy:'minimum native all-subunit detection fraction; NA remains missing',edge_max:weightMax});
    svg+='<g class="ccc-gene-edges" data-edge-count="'+present.length+'">';
    for (const e of present.slice().sort((a,b)=>number(a.row.lr_gmeans,0)-number(b.row.lr_gmeans,0))) {
      const from=nodeMap.get(nodeKey(e.source,e.ligand_complex,'ligand')),to=nodeMap.get(nodeKey(e.target,e.receptor_complex,'receptor')),gm=number(e.row.lr_gmeans,null);if (!finite(gm) || gm<0) continue;
      const fill=edgeColor(e,lookup,domain,o,gm),path=curve(from.position,to.position,from.radius+1,to.radius+2,o),width=o.edgeWidth*Math.sqrt(gm/weightMax);
      svg+='<g class="ccc-gene-edge" data-key="'+esc(e.key)+'" data-gmeans="'+gm+'"><title>'+esc(e.source+': '+e.ligand_complex+' → '+e.target+': '+e.receptor_complex+'; native lr_gmeans '+format(gm)+'; magnitude_rank '+format(number(e.row.magnitude_rank,null))+'. Molecular priors and RNA expression do not establish localization, contact, or signaling activity; condition contrasts are descriptive.')+'</title><path d="'+path.d+'" fill="none" stroke="'+esc(fill)+'" stroke-opacity="'+o.opacity+'" stroke-width="'+round(width)+'"/>'+arrow(path.end,path.tangent,o.arrowSize,fill,o.opacity)+'</g>';
    }
    svg+='</g><g class="ccc-gene-group-arcs">';
    const r0=329,r1=r0+o.arcWidth;
    for (const g of groups) {
      const fill=o.arcColorMode==='custom'?safeColor(o.arcColor,DEFAULTS.arcColor):g.color,label=(popMap.get(g.population)||{}).label || g.population,a=(g.start+g.end)/2,p=point(a,r1+65),anchor=Math.abs(Math.cos(a))<.15?'middle':Math.cos(a)>0?'start':'end';
      svg+='<path class="ccc-gene-group-arc" data-population="'+esc(g.population)+'" d="'+ringSegment(r0,r1,g.start,g.end)+'" fill="'+esc(fill)+'"><title>'+esc(label+'; '+g.nodes.length+' population × complex × role nodes')+'</title></path>'+textLabel(p.x,p.y,label,anchor,{...o,fontSize:Math.max(9,o.fontSize-1)},'ccc-gene-population-label');
    }
    svg+='</g><g class="ccc-gene-nodes">';
    for (const n of nodes) {
      const p=n.position,r=Math.max(.6,n.radius),fill=n.prop===null?'#ffffff':n.color,description=n.population+'; '+n.role+' complex '+n.complex+'; all-subunit minimum detection fraction '+format(n.prop)+(n.prop===null?'; no matching qualifying record in this condition (NA), not zero expression.':'.');
      svg+='<g class="ccc-gene-node" data-node="'+esc(n.id)+'" data-population="'+esc(n.population)+'" data-role="'+n.role+'" data-prop="'+(n.prop===null?'NA':n.prop)+'" data-x="'+round(p.x)+'" data-y="'+round(p.y)+'"><title>'+esc(description)+'</title>';
      if (n.role==='ligand') svg+='<circle cx="'+round(p.x)+'" cy="'+round(p.y)+'" r="'+round(r)+'" fill="'+esc(fill)+'" stroke="'+esc(n.color)+'" stroke-width="1.4"'+(n.prop===null?' stroke-dasharray="2 2"':'')+'/>';
      else {const half=r*Math.sqrt(Math.PI)/2;svg+='<rect x="'+round(p.x-half)+'" y="'+round(p.y-half)+'" width="'+round(half*2)+'" height="'+round(half*2)+'" fill="'+esc(fill)+'" stroke="'+esc(n.color)+'" stroke-width="1.4"'+(n.prop===null?' stroke-dasharray="2 2"':'')+'/>';}
      svg+='</g>';const label=point(n.angle,309),anchor=Math.cos(n.angle)>0?'start':'end',angle=n.angle*180/Math.PI+(Math.cos(n.angle)>0?0:180),complexMax=Math.round(clamp(370/Math.max(8,o.fontSize-1),12,35));
      svg+='<text class="ccc-gene-complex-label" x="'+round(label.x)+'" y="'+round(label.y)+'" transform="rotate('+round(angle)+' '+round(label.x)+' '+round(label.y)+')" text-anchor="'+anchor+'" font-family="Arial,Helvetica,sans-serif" font-size="'+Math.max(8,o.fontSize-1)+'" fill="#243b4b" stroke="#ffffff" stroke-width="2.5" paint-order="stroke"><title>'+esc(description)+'</title>'+esc(n.complex.length>complexMax?n.complex.slice(0,complexMax-1)+'…':n.complex)+'</text>';
    }
    svg+='</g>';if (!present.length) svg+=emptyMessage(o,selected.length?'the fixed selection in this condition':'the selected cutoffs');
    return svg+end(['Shown '+present.length+' / '+selected.length+' selected population × LR records; repeated molecular pairs in different groups are separate edges.','Circle = ligand; equal-area square = receptor. Node area = true minimum detection fraction over all complex subunits.','Hollow dashed node = NA in this condition. Edge width = native lr_gmeans (max '+format(weightMax)+'); arrow = listed receiver.','Candidates require localization/contact validation. Native rank/p-values are not genotype interaction tests.'],o,domain,'Edge color: native lr_gmeans');
  }
  const api=Object.freeze({circle,circos,geneCircle,geneKey,summary,defaults:DEFAULTS});
  global.CCCNetworks=api;
  if (typeof module!=='undefined' && module.exports) module.exports=api;
})(typeof window!=='undefined'?window:globalThis);
