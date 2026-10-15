/* Static release inventory and saved-figure navigation; no data extraction or fits. */
(function(){
 'use strict';
 const data=window.SATB2_PUBLISHED_DOWNLOADS||{assets:[],figures:[]},$=id=>document.getElementById(id);
 function node(tag,text,cls){const e=document.createElement(tag);if(text!==undefined)e.textContent=text;if(cls)e.className=cls;return e;}
 function link(text,url){const a=node('a',text);a.href=url;a.rel='noopener';return a;}
 function size(bytes){return typeof bytes==='number'?(bytes/1048576).toFixed(bytes>10485760?0:1)+' MiB':'See versioned release';}
 if($('release-link')){$('release-link').href=data.release_url;$('release-link').textContent='Open release '+(data.version||'');}
 if($('asset-list')){
  if(data.assets.some(r=>r.kind==='processed_object')){const select=$('asset-kind');if(![...(select.options||[])].some(o=>o.value==='processed_object')){const o=node('option','Processed AnnData / Seurat objects');o.value='processed_object';select.append(o);}}
  const render=()=>{const search=$('asset-search').value.toLowerCase(),kind=$('asset-kind').value,rows=data.assets.filter(r=>(kind==='ALL'||r.kind===kind)&&(!search||[r.label,r.scope,r.release_filename].join(' ').toLowerCase().includes(search))),box=$('asset-list');box.replaceChildren();$('asset-count').textContent=rows.length+' source assets';
   for(const r of rows){const card=node('article',undefined,'asset-card');card.append(link(r.label,r.url),node('p',r.scope||'Saved source file','small'),node('p',size(r.bytes)+' · '+r.release_filename,'small'));const details=node('details');details.append(node('summary','Integrity and version'),node('pre',r.sha256?'SHA256 '+r.sha256:'The archive integrity manifest is supplied by the versioned release.'));card.append(details);box.append(card);}
   if(!rows.length)box.append(node('p','No assets match these filters.','empty-state'));
  };$('asset-search').oninput=render;$('asset-kind').onchange=render;render();
 }
 if($('gallery-grid'))for(const r of data.figures){const card=node('article',undefined,'gallery-card'),a=link('',r.pdf_url),img=node('img');img.src=r.png_url;img.alt=r.title;img.loading='lazy';a.append(img);card.append(node('h2',r.title),a,node('p',r.caption),link('Download the illustrated PDF',r.pdf_url));$('gallery-grid').append(card);}
})();
