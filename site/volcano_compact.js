/* Display-only controls for the existing main atlas volcano renderer. */
(function () {
  'use strict';
  const app=window.SATB2_APP, grid=document.querySelector('.volcano-control-grid');
  if (!app || !grid) return;
  const container=document.createElement('div');
  container.className='volcano-compact-controls';
  container.innerHTML='<label for="volcano-compact-width">Volcano plot width (pixels)</label><input id="volcano-compact-width" type="number" min="360" max="1100" step="20" value="520"><label for="volcano-label-size">Gene label size (pixels)</label><input id="volcano-label-size" type="number" min="8" max="18" step="1" value="10"><p class="small">The figure is redrawn at the chosen width, preserving label size. At smaller widths some labels may be omitted if they overlap; see the note under each plot. The PDF uses its existing page dimensions but keeps the chosen label size.</p>';
  grid.append(container);
  const width=document.getElementById('volcano-compact-width');
  const font=document.getElementById('volcano-label-size');
  font.value=Number.isFinite(Number(app.state.volcanoLabelSize))?app.state.volcanoLabelSize:10;
  function update(redraw) {
    const w=Math.max(360,Math.min(1100,Number(width.value)||520));
    const f=Math.max(8,Math.min(18,Number(font.value)||10));
    for (const id of ['volcano-scanpy','volcano-seurat']) {
      const canvas=document.getElementById(id);
      if (canvas) {canvas.style.width=w+'px';canvas.style.maxWidth='100%';}
    }
    app.state.volcanoLabelSize=f;
    if (redraw) app.render();
  }
  width.addEventListener('input',()=>update(true));
  font.addEventListener('input',()=>update(true));
  update(false);
})();
