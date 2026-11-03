(function () {
  'use strict';
  const M=window.Satb2ExonModel, $=id=>document.getElementById(id), ns='http://www.w3.org/2000/svg';
  const state={table:null,source:'Reviewed source table'};
  function optionize(el,values){for(const [value,label] of Object.entries(values)){const o=document.createElement('option');o.value=value;o.textContent=label;el.appendChild(o);}}
  optionize($('metric'),M.metrics);optionize($('denominator'),M.denominators);$('denominator').value='per_1000_nuclei';
  function node(tag,attrs,text){const e=document.createElementNS(ns,tag);for(const [k,v] of Object.entries(attrs||{}))e.setAttribute(k,String(v));if(text!=null)e.textContent=text;return e;}
  function svgText(svg,x,y,text,attrs={}){svg.appendChild(node('text',{x,y,...attrs},text));}
  function choose(){return M.visible(state.table,$('cohort').value,[$('showWT').checked?'WT':'', $('showCKO').checked?'cKO':''], $('sampleSet').value);}
  function load(tsv,label){try{const parsed=M.parse(tsv);state.table=parsed;state.source=label;$('status').textContent='Showing '+label+'. Descriptive RNA evidence.';render();}catch(e){$('status').textContent='File was not loaded: '+e.message;}}
  function render(){
    const rows=choose(),metric=$('metric').value,den=$('denominator').value;
    const yLabel=M.metrics[metric]+' · '+M.denominators[den];$('plotTitle').textContent=yLabel;
    const values=rows.map(r=>M.plotValue(r,metric,den)).filter(x=>x!==null);
    const ymax=Math.max(1,...values)*1.16, left=105,right=825,top=55,bottom=430,width=900,height=530;
    const svg=node('svg',{xmlns:ns,viewBox:`0 0 ${width} ${height}`,role:'img','aria-label':yLabel+': per-mouse values'});
    svg.appendChild(node('rect',{x:0,y:0,width,height,fill:'#fff'}));
    svgText(svg,26,27,'Satb2 exon and junction evidence · '+state.source,{fill:'#835200','font-size':14,'font-family':'Arial, sans-serif','font-weight':'bold'});
    const cohortLabel=$('cohort').value==='all_CA1'?'All CA1/ProS':'CA1/ProS excluding C11';
    svgText(svg,left,49,yLabel+' · '+cohortLabel,{fill:'#173546','font-size':15,'font-family':'Arial, sans-serif','font-weight':'bold'});
    for(let i=0;i<=5;i++){const yy=bottom-(bottom-top)*i/5, val=ymax*i/5;svg.appendChild(node('line',{x1:left,y1:yy,x2:right,y2:yy,stroke:i?'#e2ebef':'#7893a0'}));svgText(svg,left-10,yy+4,Number(val.toPrecision(3)).toString(),{'text-anchor':'end',fill:'#506e7e','font-family':'Arial, sans-serif','font-size':12});}
    svg.appendChild(node('line',{x1:left,y1:top,x2:left,y2:bottom,stroke:'#7893a0'}));
    const colors={WT:$('wtColor').value,cKO:$('ckoColor').value}, groupX={WT:285,cKO:645}, radius=Number($('pointSize').value);
    for(const group of ['WT','cKO']){
      const rr=rows.filter(r=>r.condition===group).sort((a,b)=>a.sample.localeCompare(b.sample,undefined,{numeric:true}));
      svgText(svg,groupX[group],bottom+31,group,{'text-anchor':'middle','font-size':16,'font-family':'Arial, sans-serif','font-weight':'bold',fill:'#173546'});
      rr.forEach((r,i)=>{
        const val=M.plotValue(r,metric,den); if(val===null)return;
        const x=groupX[group]+(i-(rr.length-1)/2)*34,y=bottom-(bottom-top)*val/ymax;
        const shape=group==='WT'?node('circle',{cx:x,cy:y,r:radius,fill:colors[group]}):node('path',{d:`M ${x} ${y-radius-1} L ${x-radius-1} ${y+radius} L ${x+radius+1} ${y+radius} Z`,fill:colors[group]});
        shape.appendChild(node('title',{},`${r.sample}: ${val.toPrecision(5)}; ${r[metric+'_distinct_cb_ub']} pairs; ${r.n_matched_nuclei} matched nuclei; ${r.locus_distinct_cb_ub} locus pairs`));svg.appendChild(shape);
        if($('labels').checked)svgText(svg,x,y-radius-7,r.sample.replace(/^(wt|cKO)_/,''),{'text-anchor':'middle','font-size':11,'font-family':'Arial, sans-serif',fill:colors[group]});
      });
    }
    const sampleLabel=M.selectionLabel(rows,$('sampleSet').value);
    svgText(svg,left,bottom+72,'Each point = one mouse · '+sampleLabel+' · no cell-level inference',{'font-size':12,'font-family':'Arial, sans-serif',fill:'#506e7e'});
    $('chart').replaceChildren(svg);
    const missing=rows.filter(r=>M.plotValue(r,metric,den)===null).map(r=>r.sample);$('missing').textContent=missing.length?'Unavailable or insufficient coverage: '+missing.join(', '):'No unavailable values in the displayed mice.';
    const body=$('rows');body.replaceChildren();for(const r of rows){const tr=document.createElement('tr'),v=M.plotValue(r,metric,den);for(const text of [r.sample,r.condition,r.n_matched_nuclei,r.locus_distinct_cb_ub,r[metric+'_distinct_cb_ub'],v===null?'Unavailable':Number(v.toPrecision(5)),r.coverage_status]){const td=document.createElement('td');td.textContent=String(text);tr.appendChild(td);}body.appendChild(tr);}
  }
  function filename(ext){return 'SATB2_'+$('cohort').value+'_'+$('sampleSet').value+'_'+$('metric').value+'_'+$('denominator').value+'.'+ext;}
  function save(blob,name){const url=URL.createObjectURL(blob),a=document.createElement('a');a.href=url;a.download=name;a.click();setTimeout(()=>URL.revokeObjectURL(url),10000);}
  function svgString(){const el=$('chart').querySelector('svg');if(!el)throw Error('No plot to export');return new XMLSerializer().serializeToString(el);}
  function raster(){return new Promise((resolve,reject)=>{const blob=new Blob([svgString()],{type:'image/svg+xml;charset=utf-8'}),url=URL.createObjectURL(blob),im=new Image();im.onload=()=>{URL.revokeObjectURL(url);const canvas=document.createElement('canvas');canvas.width=1800;canvas.height=1060;const c=canvas.getContext('2d');c.fillStyle='#fff';c.fillRect(0,0,canvas.width,canvas.height);c.drawImage(im,0,0,canvas.width,canvas.height);resolve(canvas);};im.onerror=()=>{URL.revokeObjectURL(url);reject(Error('Could not render plot'));};im.src=url;});}
  $('saveSvg').addEventListener('click',()=>save(new Blob([svgString()],{type:'image/svg+xml;charset=utf-8'}),filename('svg')));
  $('saveTsv').addEventListener('click',()=>save(new Blob([M.toTsv(state.table,choose())],{type:'text/tab-separated-values;charset=utf-8'}),filename('tsv')));
  $('savePng').addEventListener('click',async()=>{try{const canvas=await raster();canvas.toBlob(blob=>{if(blob)save(blob,filename('png'));else $('status').textContent='PNG conversion failed.';},'image/png');}catch(e){$('status').textContent=e.message;}});
  $('savePdf').addEventListener('click',async()=>{try{if(!window.jspdf||!window.jspdf.jsPDF)throw Error('Bundled PDF library not available');const canvas=await raster();const pdf=new window.jspdf.jsPDF({orientation:'landscape',unit:'pt',format:[900,530]});pdf.addImage(canvas.toDataURL('image/png'),'PNG',0,0,900,530);pdf.save(filename('pdf'));}catch(e){$('status').textContent=e.message;}});
  $('reset').addEventListener('click',()=>{$('cohort').value='all_CA1';$('sampleSet').value='all_12';$('metric').value='candidate_exon';$('denominator').value='per_1000_nuclei';$('showWT').checked=$('showCKO').checked=true;$('labels').checked=false;$('pointSize').value=7;$('wtColor').value='#0878b7';$('ckoColor').value='#dc7200';render();});
  for(const id of ['cohort','sampleSet','metric','denominator','showWT','showCKO','labels','pointSize','wtColor','ckoColor'])$(id).addEventListener('input',render);
  load(window.SATB2_REVIEWED_TSV,'reviewed mouse-level source');
})();
