/* Display-only label and categorical-color helpers. No annotation, grouping or fit changes. */
(function (global) {
  'use strict';
  const own = (object, key) => object && Object.prototype.hasOwnProperty.call(object, key);
  const PALETTES = Object.freeze({
    tableau: Object.freeze(['#4e79a7','#f28e2b','#e15759','#76b7b2','#59a14f','#edc949','#af7aa1','#ff9da7','#9c755f','#bab0ab']),
    bright: Object.freeze(['#4477aa','#ee6677','#228833','#ccbb44','#66ccee','#aa3377','#bb5588','#44aa99','#882255','#999933']),
    colorblind: Object.freeze(['#0072b2','#e69f00','#009e73','#cc79a7','#56b4e9','#d55e00','#f0e442','#000000'])
  });
  const LABEL_CAUTION = 'SL-like and DL-like describe superficial- and deep-associated expression, not a physical layer assignment. Short labels and colors affect display only; original annotations, grouping keys, nucleus membership and DE targets are preserved.';
  function hash(value) { let h=2166136261;for(const c of String(value))h=Math.imul(h^c.charCodeAt(0),16777619);return h>>>0; }
  function validateHex(value) { return typeof value==='string' && /^#[0-9a-f]{6}$/i.test(value); }
  function nativeColor(value) {
    if(value==='WT')return '#377eb8';if(value==='cKO')return '#e09133';
    const s=String(value).toLowerCase();
    if(/unavailable|unresolved|mixed|weak/.test(s))return '#9ca6ad';
    if(s.includes('conflict'))return '#8b6b9f';
    if(s.includes('superficial-associated'))return '#458ab3';
    if(s.includes('deep-associated'))return '#b66942';
    return `hsl(${hash(value)%360} 47% 47%)`;
  }
  // categoryColors is an exact ORIGINAL category-value -> #rrggbb map for the
  // current field/scope. The caller namespaces maps by field; shortened labels
  // must never be used as color keys, because several original labels may collide.
  function colorForCategory(value, options={}) {
    const key=String(value), colors=options.categoryColors;
    if(own(colors,key)&&validateHex(colors[key]))return colors[key].toLowerCase();
    const palette=PALETTES[options.categoricalPalette];
    return palette ? palette[hash(key)%palette.length] : nativeColor(value);
  }
  function splitPrefix(text) {
    // Native cluster/reference identifiers remain visible in every automatic style.
    const match=text.match(/^((?:C\d+|Cluster\s+\S+)\s*[·|:]\s*)(.*)$/s);
    return match ? {prefix:match[1],body:match[2]} : {prefix:'',body:text};
  }
  function uncertainty(text) {
    const result=[];
    if(/conflict/i.test(text))result.push('conflict');
    if(/mixed/i.test(text))result.push('mixed');
    if(/weak/i.test(text))result.push('weak');
    if(/provisional/i.test(text))result.push('provisional');
    if(/unresolved/i.test(text))result.push('unresolved');
    if(/signature review|review required|evidence review|evidence-review/i.test(text))result.push('review');
    if(/reference[- ]like/i.test(text))result.push('reference-like');
    return result;
  }
  function retainUncertainty(label, full) {
    const missing=uncertainty(full).filter(flag=>!label.toLowerCase().includes(flag));
    return missing.length ? label+' ('+missing.join('; ')+')' : label;
  }
  function automaticLabel(full, style) {
    const {prefix,body}=splitPrefix(full), short=style==='short';let result;
    if(/^Superficial-associated expression program enriched$/i.test(body))result=short?'SL-like':'Sup.assoc.exp.';
    else if(/^Deep-associated expression program enriched$/i.test(body))result=short?'DL-like':'Deep.assoc.exp.';
    else if(/^Mixed or weak layer-associated expression programs?$/i.test(body))result=short?'Mixed (weak)':'Mixed/weak.exp.';
    else if(/^Layer program not assigned:\s*conflicting identity signatures$/i.test(body))result=short?'Unassigned (conflict)':'Unassigned/conflict';
    else if(/^Not assigned$/i.test(body))result='Unassigned';
    else {
      result=body
        .replace(/expression program enriched/gi,short?'exp.':'exp.enriched')
        .replace(/expression programs?/gi,'exp.')
        .replace(/superficial-associated/gi,'Sup.assoc.')
        .replace(/deep-associated/gi,'Deep.assoc.')
        .replace(/layer-associated/gi,'layer.assoc.')
        .replace(/associated signature review/gi,'assoc.sig.review')
        .replace(/conflicting identity signatures/gi,'identity conflict')
        .replace(/sparse WT marker evidence/gi,'sparse WT evidence');
    }
    // No arbitrary character-count clipping: it can hide mixed/review qualifiers.
    return prefix+retainUncertainty(result,body);
  }
  function formatLabel(full, options={}) {
    const source=full==null?'Unavailable':String(full), style=options.labelStyle||'full';
    if(style==='full'||!['compact','short','custom'].includes(style))return source;
    const {prefix,body}=splitPrefix(source);
    // A category alias also applies when a renderer adds a native cluster ID.
    // An explicit whole-label alias always takes precedence over a body alias.
    const aliasKey=own(options.labelAliases,source)?source:prefix&&own(options.labelAliases,body)?body:null;
    if(style==='custom'&&aliasKey!==null) {
      const raw=options.labelAliases[aliasKey];
      if(typeof raw==='string'&&raw.trim()) {
        const alias=raw.trim();
        const text=prefix&&!alias.startsWith(prefix.trim())?prefix+alias:alias;
        return retainUncertainty(text,body);
      }
    }
    return automaticLabel(source,style==='short'?'short':'compact');
  }
  // Feed ORIGINAL category values to this function. Never feed already-shortened
  // labels back into category(), grouping, resolveFit(), or marker source rows.
  const api=Object.freeze({formatLabel,validateHex,colorForCategory,nativeColor,labelCaution:LABEL_CAUTION,palettes:PALETTES});
  global.SATB2_DISPLAY=api;
  if(typeof module!=='undefined'&&module.exports)module.exports=api;
})(typeof window!=='undefined'?window:globalThis);
