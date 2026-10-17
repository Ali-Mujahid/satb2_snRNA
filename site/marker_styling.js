/* Reversible dot appearance only. Expression, detection fractions and cluster membership are never edited. */
(function (global) {
  'use strict';
  const own=(value,key)=>value!=null&&Object.prototype.hasOwnProperty.call(value,key);
  const TOOLS=Object.freeze(['scanpy','seurat']);
  const BASE=Object.freeze({areaMultiplier:1,opacityMultiplier:1,hidden:false});
  const CAUTION='Marker styling changes dot appearance only. Hidden, smaller or lighter dots retain their saved expression and detected fractions; edited dots are identified in settings and exported marker tables.';
  function defaults(){return {};}
  function identity(value){return (typeof value==='string'||typeof value==='number')&&String(value).trim()!==''?String(value):null;}
  function cluster(value){return value==null||value==='ALL'||value==='all'||value==='*'?null:identity(value);}
  function context(value,allowBoth=false){
    if(!value||typeof value!=='object')return null;
    const tool=value.tool,branch=identity(value.branch),population=identity(value.population),geneId=identity(value.geneId);
    if(!(TOOLS.includes(tool)||(allowBoth&&tool==='both'))||!branch||!population||!geneId)return null;
    if(value.clusterId!=null&&!['string','number'].includes(typeof value.clusterId))return null;
    // An empty cluster choice is a whole-gene selection, matching a blank select.
    const clusterId=cluster(value.clusterId);
    return {tool,branch,population,clusterId,geneId};
  }
  function key(value){const c=context(value);return c?JSON.stringify([c.tool,c.branch,c.population,c.clusterId,c.geneId]):null;}
  function keyContext(value){
    try{const parts=JSON.parse(value);if(!Array.isArray(parts)||parts.length!==5)return null;const c=context({tool:parts[0],branch:parts[1],population:parts[2],clusterId:parts[3],geneId:parts[4]});return c&&key(c)===value?c:null;}catch(_){return null;}
  }
  function multiplier(value,fallback=1){return typeof value==='number'&&Number.isFinite(value)?Math.min(1,Math.max(0,value)):fallback;}
  function sanitizeStyle(value){
    const s=value&&typeof value==='object'?value:{};
    return {areaMultiplier:multiplier(s.areaMultiplier),opacityMultiplier:multiplier(s.opacityMultiplier),hidden:s.hidden===true};
  }
  function normalize(map){
    const out={};if(!map||typeof map!=='object'||Array.isArray(map))return out;
    for(const k of Object.keys(map))if(keyContext(k)&&map[k]&&typeof map[k]==='object'&&!Array.isArray(map[k]))out[k]=sanitizeStyle(map[k]);
    return out;
  }
  function selections(value){const c=context(value,true);return c?(c.tool==='both'?TOOLS.map(tool=>({...c,tool})):[c]):[];}
  function apply(map,value,style){const out=normalize(map),s=sanitizeStyle(style);for(const c of selections(value))out[key(c)]={...s};return out;}
  function resetSelected(map,value){const out=normalize(map);for(const c of selections(value))delete out[key(c)];return out;}
  function resetPanel(map,value){
    const out=normalize(map),scope=value&&typeof value==='object'?value:{},branch=identity(scope.branch),population=identity(scope.population);
    if(!branch||!population||!(TOOLS.includes(scope.tool)||scope.tool==='both'))return out;
    for(const k of Object.keys(out)){const c=keyContext(k);if(c.branch===branch&&c.population===population&&(scope.tool==='both'||c.tool===scope.tool))delete out[k];}
    return out;
  }
  function resolveStyle(map,value){
    const c=context(value);if(!c)return {...BASE,source:'default'};
    const specific=key(c),whole=key({...c,clusterId:null});
    if(c.clusterId!==null&&own(map,specific))return {...sanitizeStyle(map[specific]),source:'dot'};
    if(own(map,whole))return {...sanitizeStyle(map[whole]),source:'gene'};
    return {...BASE,source:'default'};
  }
  function styleDot(visual,override){
    const original=visual&&typeof visual==='object'?visual:{},s=sanitizeStyle(override),out={...original};
    if(typeof original.area==='number'&&Number.isFinite(original.area))out.area=original.area*s.areaMultiplier;
    if(typeof original.alpha==='number'&&Number.isFinite(original.alpha))out.alpha=original.alpha*s.opacityMultiplier;
    out.hidden=s.hidden;
    return out;
  }
  // Sorted, source-keyed records let CSV/settings exports disclose presentation edits.
  function records(map,value){
    const scope=value&&typeof value==='object'?value:{};
    return Object.entries(normalize(map)).map(([k,s])=>({...keyContext(k),...s})).filter(c=>(!scope.tool||scope.tool==='both'||scope.tool===c.tool)&&(!scope.branch||String(scope.branch)===c.branch)&&(!scope.population||String(scope.population)===c.population)).sort((a,b)=>key(a).localeCompare(key(b)));
  }
  const api=Object.freeze({defaults,normalize,key,keyContext,apply,resetSelected,resetPanel,resolveStyle,styleDot,records,sanitizeStyle,caution:CAUTION});
  global.SATB2_MARKER_STYLE=api;
  if(typeof module!=='undefined'&&module.exports)module.exports=api;
})(typeof window!=='undefined'?window:globalThis);
