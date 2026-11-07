/* Combine existing marker panels without altering expression summaries. */
(() => {
 'use strict';
 const anchor=document.getElementById('marker-preset'); if(!anchor)return;
 const container=document.createElement('div');container.className='controls-card';
 container.innerHTML='<h3>Combine marker panels</h3><p>Select multiple panels (on a Mac, hold Command while clicking), then combine their genes. Repeated genes are included once.</p><label for="multi-marker-presets">Marker panels</label><select id="multi-marker-presets" multiple size="10" style="width:100%;max-width:580px"></select><div class="inline-controls"><button type="button" id="combine-marker-presets">Use selected panels</button><button type="button" id="add-marker-presets">Add to current genes</button></div><p id="multi-marker-note" aria-live="polite"></p>';
 anchor.closest('.marker-controls').after(container);
 const select=container.querySelector('select'),note=container.querySelector('#multi-marker-note');
 function refresh(){const prior=new Set([...select.selectedOptions].map(o=>o.value));select.replaceChildren(...[...anchor.options].filter(o=>o.value).map(o=>{const n=document.createElement('option');n.value=o.value;n.textContent=o.textContent;n.selected=prior.has(o.value);return n}));}
 function combine(add){const app=window.SATB2_APP,ids=[...select.selectedOptions].map(o=>o.value);if(!app?.manifest || !ids.length){note.textContent='Select at least one available marker panel.';return}const sets=ids.map(id=>app.manifest.gene_sets.find(s=>s.id===id)).filter(Boolean);if(sets.length!==ids.length){note.textContent='A selected panel is unavailable; reopen the list.';return}const valid=new Set(app.manifest.genes.map(g=>g.gene_id)),union=new Set(add?app.state.markers:[]);for(const s of sets)for(const id of s.gene_ids)if(valid.has(id))union.add(id);app.state.markers=[...union];app.state.markerPage=0;app.state.markerOrder='selection';app.state.markerPreset='';app.applyMarkerOrder('');note.textContent=`${sets.length} panels combined; ${union.size} unique available genes selected. Plot values and cell selection are unchanged.`;}
 select.addEventListener('focus',refresh);container.querySelector('#combine-marker-presets').addEventListener('click',()=>combine(false));container.querySelector('#add-marker-presets').addEventListener('click',()=>combine(true));refresh();
})();
