/* Display helpers over produced arrays only. No models, embeddings or RNA extraction. */
(function(global) {
  'use strict';
  const finite=value=>typeof value==='number'&&Number.isFinite(value);
  function clamp(value,min,max,fallback=min) {
    if(!finite(min)||!finite(max)||min>max)throw Error('Invalid display limits');
    return Math.max(min,Math.min(max,finite(value)?value:fallback));
  }
  function range(values,override) {
    let lo=null,hi=null;for(const value of values||[])if(finite(value)){lo=lo===null?value:Math.min(lo,value);hi=hi===null?value:Math.max(hi,value);}
    if(override){const lower=Array.isArray(override)?override[0]:override.min,upper=Array.isArray(override)?override[1]:override.max;if(finite(lower)&&finite(upper)){if(lower>=upper)throw Error('Display minimum must be below maximum');lo=lower;hi=upper;}}
    return lo===null?null:[lo,hi];
  }
  function dotGeometry(options={}) {
    const {fraction,z,metric='scaled_expression_z',negativeZ=true,maxArea=100,negativeArea=.25,negativeAlpha=.18,alpha=.9,glyphScale=1}=options;
    const value=Object.prototype.hasOwnProperty.call(options,'value')?options.value:metric==='scaled_expression_z'?z:metric==='fraction_detected'?fraction:options[metric];
    const negative=!!negativeZ&&metric==='scaled_expression_z'&&finite(z)&&z<0;
    const fractionValid=finite(fraction)&&fraction>=0&&fraction<=1,scale=clamp(glyphScale,.05,8,1),base=clamp(maxArea,0,100000,100),attenuation=negative?clamp(negativeArea,0,1,.25):1;
    const area=fractionValid?base*fraction*attenuation*scale*scale:null;
    return {value:finite(value)?value:null,fraction:finite(fraction)?fraction:null,negative,area,radius:area===null?null:Math.sqrt(area/Math.PI),alpha:negative?clamp(negativeAlpha,0,1,.18):clamp(alpha,0,1,.9),fraction_valid:fractionValid};
  }
  function transformCoordinates(nuclei,basis='harmony',spacing=1) {
    if(!['harmony','unintegrated'].includes(basis))throw Error('Unsupported original coordinate basis');
    const factor=clamp(spacing,.25,2.5,1),xKey='x_'+basis,yKey='y_'+basis,centroids=new Map();let sx=0,sy=0,n=0;
    for(const point of nuclei||[]){if(!finite(point[xKey])||!finite(point[yKey]))continue;const id=String(point.cluster_id),c=centroids.get(id)||{x:0,y:0,n:0};c.x+=point[xKey];c.y+=point[yKey];c.n++;centroids.set(id,c);sx+=point[xKey];sy+=point[yKey];n++;}
    const gx=n?sx/n:0,gy=n?sy/n:0;for(const c of centroids.values()){c.x/=c.n;c.y/=c.n;}
    return (nuclei||[]).map(point=>{const x=point[xKey],y=point[yKey],valid=finite(x)&&finite(y),c=centroids.get(String(point.cluster_id));return {x:!valid?null:factor===1?x:x+(factor-1)*(c.x-gx),y:!valid?null:factor===1?y:y+(factor-1)*(c.y-gy),cluster_id:String(point.cluster_id),nucleus_id:point.nucleus_id,source_x:finite(x)?x:null,source_y:finite(y)?y:null,spacing:factor,display_only:factor!==1,source_coordinates_modified:false};});
  }
  function quantile(sorted,p) {
    if(!sorted.length)return null;const q=clamp(p,0,1,.5),index=(sorted.length-1)*q,lower=Math.floor(index),upper=Math.ceil(index);return sorted[lower]+(sorted[upper]-sorted[lower])*(index-lower);
  }
  function boxStats(values) {
    const sorted=(values||[]).filter(finite).slice().sort((a,b)=>a-b),n=sorted.length;
    if(!n)return {n:0,q1:null,median:null,q3:null,iqr:null,lowerFence:null,upperFence:null,whiskerLow:null,whiskerHigh:null,outliers:[],minimum:null,maximum:null,quantile_definition:'linear_interpolation_(n-1)*p'};
    const q1=quantile(sorted,.25),median=quantile(sorted,.5),q3=quantile(sorted,.75),iqr=q3-q1,lowerFence=q1-1.5*iqr,upperFence=q3+1.5*iqr,inliers=sorted.filter(v=>v>=lowerFence&&v<=upperFence);
    return {n,q1,median,q3,iqr,lowerFence,upperFence,whiskerLow:inliers[0]??null,whiskerHigh:inliers[inliers.length-1]??null,outliers:sorted.filter(v=>v<lowerFence||v>upperFence),minimum:sorted[0],maximum:sorted[n-1],quantile_definition:'linear_interpolation_(n-1)*p'};
  }
  function mousePoints(rows,options={}) {
    const countKey=options.countKey||'count',factorKey=options.factorKey||'size_factor',sampleKey=options.sampleKey||'sample',conditionKey=options.conditionKey||'condition',seen=new Set();
    return (rows||[]).map(row=>{
      const sample=row[sampleKey];if(sample===null||sample===undefined||sample==='')throw Error('Mouse identity is missing');if(seen.has(String(sample)))throw Error('Duplicate mouse row: '+sample);seen.add(String(sample));
      const count=row[countKey],factor=row[factorKey],countValid=finite(count)&&count>=0,positiveSF=finite(factor)&&factor>0;let missing_reason=null;
      if(!countValid)missing_reason='SOURCE_COUNT_UNAVAILABLE';else if(!positiveSF)missing_reason='RECORDED_POSITIVE_SIZE_FACTOR_UNAVAILABLE';
      const normalized=!missing_reason?count/factor:null,condition=row[conditionKey];
      return {...row,sample,condition,source_count:countValid?count:null,recorded_size_factor:positiveSF?factor:null,normalized_value:finite(normalized)?normalized:null,value:finite(normalized)?normalized:null,missing_reason,normalization:'source_integer_pseudobulk_count_divided_by_recorded_source_fit_size_factor',biological_unit:'independent_mouse',statistical_test_performed:false};
    });
  }
  function mouseRows(countPayload,fit,geneIndex,samples=[]) {
    if(countPayload?.gene_axis!=='global'||!Number.isInteger(countPayload.n_genes)||!Number.isInteger(geneIndex)||geneIndex<0||geneIndex>=countPayload.n_genes)throw Error('Unknown source gene index');
    const order=countPayload.sample_order,fitOrder=fit?.sample_order;if(!Array.isArray(order)||!Array.isArray(fitOrder)||order.length!==12||fitOrder.length!==12||order.some((id,i)=>id!==fitOrder[i])||new Set(order).size!==12)throw Error('Recorded mouse source axes differ');
    for(const key of ['size_factors','eligible','n_nuclei'])if(!Array.isArray(fit[key])||fit[key].length!==12)throw Error('Missing all12 source audit: '+key);
    const hits=(countPayload.rows||[]).filter(row=>row[0]===geneIndex);if(hits.length>1)throw Error('Duplicate source gene index');
    let counts=hits[0]?.[1];if(!counts){if(countPayload.omitted_rows!=='EXPLICIT_ALL_ZERO')throw Error('Missing counts do not declare explicit zeros');counts=Array(12).fill(0);}
    if(!Array.isArray(counts)||counts.length!==12||counts.some(v=>!Number.isSafeInteger(v)||v<0))throw Error('Source integer mouse counts are invalid');
    const conditions=new Map((samples||[]).map(s=>[s.sample,s.condition]));return order.map((sample,i)=>({sample,condition:conditions.get(sample)??null,count:counts[i],size_factor:fit.size_factors[i],included_in_model:fit.eligible[i],n_nuclei:fit.n_nuclei[i],gene_index:geneIndex,source_fit_id:fit.id??null,source_count_sha256:fit.counts_sha256??null,source_counts_scope:'EXISTING_FULL_GENE_INTEGER_MOUSE_PSEUDOBULK'}));
  }
  function labelSuggestions(row={},options={}) {
    const descriptors=[['reviewed_label','Existing reviewed molecular label'],['identity_label','Existing identity review'],['layer_program_label','WT-associated expression program; physical layer undetermined'],['expression_program_label','Existing expression program; morphology unvalidated'],['reference_subclass_name','Source MapMyCells subclass'],['reference_supertype_name','Source MapMyCells supertype; provisional'],['reference_cluster_name','Source MapMyCells reference cluster; provisional']],result=[],seen=new Set();
    for(const [key,note]of descriptors){const value=row[key]??options[key];if(typeof value!=='string'||!value.trim()||seen.has(value.trim())||/^(NA|NONE|UNAVAILABLE|UNDETERMINED)$/i.test(value.trim()))continue;seen.add(value.trim());result.push({text:value.trim(),source:key,note,source_confidence:row.label_confidence??null,provisional:true,physical_layer_assignment:'UNDETERMINED',scientific_annotation_modified:false});}
    return result;
  }
  function markerLayout(nGenes,nRows,options={}) {
    const count=Math.max(0,Math.floor(finite(nGenes)?nGenes:0)),rows=Math.max(0,Math.floor(finite(nRows)?nRows:0)),left=clamp(options.labelWidth,80,520,260),column=clamp(options.columnWidth,8,60,18),rowHeight=clamp(options.rowHeight,12,80,26);
    return {gene_count:count,row_count:rows,label_width:left,column_width:column,row_height:rowHeight,width:Math.max(clamp(options.minWidth,320,5000,700),left+count*column+35),height:Math.max(340,rows*rowHeight+115),requires_horizontal_scroll:count>30,page_capacity_options:[60,90,120,180],values_modified:false};
  }

  /* MARKER_ORDERING_HELPERS_START: descriptive display order, never marker/DE evidence. */
  function markerGroupKey(row) {return (row.source_cluster_ids||[String(row.cluster_id)]).map(String).join('+');}
  function naturalCompare(a,b) {return String(a).localeCompare(String(b),undefined,{numeric:true});}
  function markerGroups(rows,options={}) {
    const unique=new Map();for(const row of rows||[]){const key=markerGroupKey(row);if(!unique.has(key))unique.set(key,{key,label:row.display_label||String(row.cluster_id),ids:(row.source_cluster_ids||[String(row.cluster_id)]).map(String)});}
    const groups=[...unique.values()].sort((a,b)=>naturalCompare(a.ids[0],b.ids[0])||naturalCompare(a.key,b.key));
    if(options.mode==='descending')groups.reverse();
    if(options.mode==='custom') {const tokens=String(options.custom||'').split(/[\s,;]+/).filter(Boolean).map(t=>t.replace(/^C(?=\d)/i,'')),ranks=new Map(tokens.map((t,i)=>[t,i]));groups.sort((a,b)=>(ranks.get(a.key)??Infinity)-(ranks.get(b.key)??Infinity)||naturalCompare(a.ids[0],b.ids[0])||naturalCompare(a.key,b.key));}
    return groups;
  }
  function markerOrdering(geneIds,rows,options={}) {
    const ids=[...new Set(geneIds||[])],mode=['selection','cluster','focus'].includes(options.mode)?options.mode:'selection',groups=markerGroups(rows),byGene=new Map(ids.map(id=>[id,[]])),groupRank=new Map(groups.map((g,i)=>[g.key,i])),meta={};
    for(const row of rows||[])if(byGene.has(row.gene_id))byGene.get(row.gene_id).push(row);
    const score=(a,b)=>(b.contrast-a.contrast)||(b.detected_contrast-a.detected_contrast)||(b.mean-a.mean)||(groupRank.get(a.group)-groupRank.get(b.group));
    for(const id of ids) {
      const records=byGene.get(id),valid=records.length===groups.length&&groups.length>1&&records.every(r=>finite(r.mean_log1p_normalized)&&finite(r.n_nuclei)&&r.n_nuclei>0),candidates=[];
      if(valid)for(const r of records){const other=records.filter(x=>markerGroupKey(x)!==markerGroupKey(r)),denom=other.reduce((n,x)=>n+x.n_nuclei,0),rest=other.reduce((n,x)=>n+x.mean_log1p_normalized*x.n_nuclei,0)/denom,detectValid=records.every(x=>finite(x.fraction_detected)),restDetected=detectValid?other.reduce((n,x)=>n+x.fraction_detected*x.n_nuclei,0)/denom:null;candidates.push({group:markerGroupKey(r),contrast:r.mean_log1p_normalized-rest,detected_contrast:detectValid?r.fraction_detected-restDetected:0,mean:r.mean_log1p_normalized});}
      candidates.sort(score);const best=candidates[0],tied=best&&best.contrast>1e-12?candidates.filter(c=>Math.abs(c.contrast-best.contrast)<=1e-12&&Math.abs(c.detected_contrast-best.detected_contrast)<=1e-12):[],assigned=best&&best.mean>0&&best.contrast>1e-12&&tied.length===1,focus=candidates.find(c=>c.group===options.focus);
      meta[id]={gene_id:id,peak_group:assigned?best.group:null,peak_contrast:best?.contrast??null,peak_detected_contrast:best?.detected_contrast??null,focus_contrast:focus?.contrast??null,focus_detected_contrast:focus?.detected_contrast??null,peak_tied:tied.length>1,tied_groups:tied.map(c=>c.group),ordering_status:!valid?'INCOMPLETE_OR_SINGLE_GROUP':assigned?'DESCRIPTIVE_UNIQUE_POSITIVE_PATTERN':tied.length>1?'TIED_PATTERN_UNASSIGNED':'FLAT_OR_ZERO_UNASSIGNED'};
    }
    const original=new Map(ids.map((id,i)=>[id,i]));let ordered=ids.slice();
    if(mode==='cluster')ordered.sort((a,b)=>{const x=meta[a],y=meta[b];return (groupRank.get(x.peak_group)??Infinity)-(groupRank.get(y.peak_group)??Infinity)||((x.peak_group!==null&&y.peak_group!==null)?y.peak_contrast-x.peak_contrast||y.peak_detected_contrast-x.peak_detected_contrast:0)||original.get(a)-original.get(b);});
    if(mode==='focus')ordered.sort((a,b)=>{const x=meta[a],y=meta[b];return (finite(x.focus_contrast)?0:1)-(finite(y.focus_contrast)?0:1)||(finite(x.focus_contrast)&&finite(y.focus_contrast)?y.focus_contrast-x.focus_contrast||y.focus_detected_contrast-x.focus_detected_contrast:0)||original.get(a)-original.get(b);});
    return {gene_ids:ordered,metadata:meta,groups,mode,focus:options.focus||null,method:'DESCRIPTIVE_MEAN_LOG1P_MINUS_NUCLEUS_WEIGHTED_OTHER_RETAINED_ROWS',statistical_test_performed:false,validated_representative_marker_selection:false};
  }
  /* MARKER_ORDERING_HELPERS_END */
  global.SATB2_V2={version:'display-helpers-2.0.0',finite,clamp,range,dotGeometry,transformCoordinates,quantile,boxStats,mousePoints,mouseRows,labelSuggestions,markerLayout,markerGroupKey,markerGroups,markerOrdering};
})(typeof window==='undefined'?globalThis:window);
