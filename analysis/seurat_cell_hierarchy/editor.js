'use strict';

const byId = id => document.getElementById(id);
const storageKey = 'satb2:cell-hierarchy:1';
const rootId = 'Hippocampal nuclei';
const widthOK = n => Number.isInteger(n) && n >= 50 && n <= 150 && n % 5 === 0;
const colorOK = c => typeof c === 'string' && /^#[0-9a-f]{6}$/i.test(c);

function mixWhite(hex, fraction) {
  const channels = hex.slice(1).match(/../g).map(c => {
    const n = parseInt(c, 16);
    return Math.round(n + (255 - n) * fraction).toString(16).padStart(2, '0');
  });
  return '#' + channels.join('');
}

// The source SVG uses paired cubic curves with horizontal tangents at both ends.
function readRibbon(path) {
  const original = path.getAttribute('d');
  const n = original.match(/-?[\d.]+/g).map(Number);
  if (n.length !== 16) throw new Error('Unsupported ribbon geometry');
  return {
    path, original, fill: path.getAttribute('fill'), from: path.dataset.from, to: path.dataset.to,
    x: [n[0], n[2], n[4], n[6]], start: [n[1], n[15]], end: [n[7], n[9]]
  };
}

function writeRibbon(edge) {
  const [x, a, b, z] = edge.x;
  const [s, t] = edge.placedStart;
  const [u, v] = edge.placedEnd;
  return 'M '+[x,s]+' C '+[a,s,b,u,z,u]+' L '+[z,v]+' C '+[b,v,a,t,x,t]+' Z';
}

function branchId(path) {
  return path.dataset.kind === 'CA1-assignment' ? path.dataset.from : path.dataset.to;
}

function mount(svg) {
  const nodes = new Map();
  for (const bar of svg.querySelectorAll('#node-bars rect[data-node]')) {
    const id = bar.dataset.node;
    nodes.set(id, {id, label:bar.dataset.label, fill:bar.getAttribute('fill'), incoming:[], outgoing:[]});
  }
  const edges = [...svg.querySelectorAll('path[data-kind]')].map(readRibbon);
  for (const e of edges) {
    nodes.get(e.from).outgoing.push(e);
    nodes.get(e.to).incoming.push(e);
  }
  const rectangles = [...svg.querySelectorAll('rect[data-node]')].map(rect => ({
    rect, id:rect.dataset.node, y:Number(rect.getAttribute('y')), h:Number(rect.getAttribute('height')),
    fill:rect.getAttribute('fill'), bridge:rect.parentElement.id === 'ribbons'
  }));

  const state = {overall:100, selected:'Excitatory neurons', edits:{}};
  try {
    const saved = JSON.parse(localStorage.getItem(storageKey));
    if (saved) {
      if (widthOK(saved.overall)) state.overall = saved.overall;
      if (nodes.has(saved.selected)) state.selected = saved.selected;
      for (const [id, edit] of Object.entries(saved.edits || {})) {
        if (!nodes.has(id) || !edit) continue;
        const clean = {};
        if (widthOK(edit.width) && edit.width !== 100) clean.width = edit.width;
        if (colorOK(edit.color)) clean.color = edit.color.toLowerCase();
        if (Object.keys(clean).length) state.edits[id] = clean;
      }
    }
  } catch { /* Storage can be unavailable in private browsing. */ }
  function persist() {
    try { localStorage.setItem(storageKey, JSON.stringify(state)); } catch { /* Editor still works. */ }
  }
  const size = id => state.edits[id]?.width || 100;
  const color = id => state.edits[id]?.color;
  function scale(id) {
    return state.overall * size(id) / 10000;
  }

  const ports = [];
  for (const n of nodes.values()) {
    for (const [edgesAtNode, side, placed, other] of [
      [n.outgoing,'start','placedStart','to'],[n.incoming,'end','placedEnd','from']
    ]) {
      if (!edgesAtNode.length) continue;
      const ordered = edgesAtNode.slice().sort((a,b) => a[side][0] - b[side][0]);
      const lo = ordered[0][side][0], hi = ordered[ordered.length-1][side][1];
      ports.push({id:n.id, ordered, side, placed, other, span:hi-lo, center:(hi+lo)/2});
    }
  }
  const glialId = 'Oligodendrocytes';
  const glialLabel = svg.querySelector('#editable-labels [data-node="Oligodendrocytes"]');
  const labelY = Number(glialLabel.getAttribute('y'));
  const geometry = id => rectangles.find(r => r.id === id && !r.bridge);
  const glia = geometry(glialId), excitatory = geometry('Excitatory neurons');
  function render() {
    // Lift the adjacent terminal when enlarged major branches would otherwise touch.
    const excitatoryTop = excitatory.y + excitatory.h * (1-scale(excitatory.id))/2;
    const glialBottom = glia.y + glia.h * (1+scale(glialId))/2;
    const glialShift = Math.min(0, excitatoryTop - glialBottom - 12);
    glialLabel.setAttribute('y', labelY + glialShift);
    for (const {id,ordered,side,placed,other,span,center} of ports) {
      const factor = scale(id);
      if (factor === 1 && (ordered.length === 1 || ordered.every(e => size(e[other]) === 100))) {
        for (const e of ordered) e[placed] = e[side];
        continue;
      }
      // Repack adjacent ribbons into each junction so they remain attached as widths change.
      const weights = ordered.map(e => (e[side][1] - e[side][0]) * size(e[other]));
      const sum = weights.reduce((a,b) => a+b, 0);
      let y = center - span * factor / 2;
      ordered.forEach((e,i) => {
        const next = y + span * factor * weights[i] / sum;
        e[placed] = [Number(y.toFixed(3)), Number(next.toFixed(3))];
        y = next;
      });
    }
    for (const e of edges) {
      if (glialShift && e.to === glialId) e.placedEnd = e.placedEnd.map(y => y+glialShift);
      const originalShape = e.placedStart === e.start && e.placedEnd === e.end;
      e.path.setAttribute('d', originalShape ? e.original : writeRibbon(e));
      const assignment = e.path.dataset.kind === 'CA1-assignment';
      const fill = assignment ? color(e.from) || color(e.to) : color(e.to);
      e.path.setAttribute('fill', fill || e.fill);
    }
    for (const r of rectangles) {
      const h = r.h * scale(r.id);
      r.rect.setAttribute('height', h);
      r.rect.setAttribute('y', r.y + (r.h-h)/2 + (r.id === glialId ? glialShift : 0));
      const c = color(r.id);
      r.rect.setAttribute('fill', c ? (r.bridge ? mixWhite(c,0.70) : c) : r.fill);
    }
  }

  const select = byId('branch');
  const groups = new Map();
  for (const n of nodes.values()) {
    const group = n.id.startsWith('CA1-C') ? 'CA1 clusters' : n.id.startsWith('IN-C') ? 'Inhibitory clusters'
      : [rootId,'Oligodendrocytes','Excitatory neurons','Inhibitory neurons','Astrocytes','Microglia'].includes(n.id) ? 'Stem and cell types' : 'CA1 layer or region';
    if (!groups.has(group)) groups.set(group, []);
    groups.get(group).push(n);
  }
  for (const [name, group] of groups) {
    const og = document.createElement('optgroup');
    og.label = name;
    if (name.endsWith('clusters')) group.sort((a,b) => a.id.localeCompare(b.id,'en',{numeric:true}));
    for (const n of group) {
      const label = n.label ? n.id.replace(/^(CA1|IN)-/,'')+' · '+n.label : n.id;
      og.append(new Option(label,n.id));
    }
    select.append(og);
  }

  function syncControls() {
    const id = state.selected;
    select.value = id;
    byId('width').value = size(id);
    byId('width-value').value = size(id)+'%';
    byId('overall').value = state.overall;
    byId('overall-value').value = state.overall+'%';
    byId('color').value = color(id) || nodes.get(id).fill;
    byId('hex').value = byId('color').value.toUpperCase();
    byId('thinner').disabled = size(id) === 50;
    byId('thicker').disabled = size(id) === 150;
  }
  function update() { render(); syncControls(); persist(); }
  function choose(id) {
    if (!nodes.has(id)) return;
    state.selected = id; syncControls(); persist();
  }
  function editSelected(field, value) {
    const id = state.selected;
    const edit = {...state.edits[id], [field]:value};
    if (edit.width === 100) delete edit.width;
    if (edit.color === nodes.get(id).fill) delete edit.color;
    if (Object.keys(edit).length) state.edits[id] = edit;
    else delete state.edits[id];
    update();
  }
  select.addEventListener('change', () => choose(select.value));
  svg.addEventListener('click', event => {
    const hit = event.target.closest('[data-node],path[data-kind]');
    if (hit) choose(hit.dataset.node || branchId(hit));
  });
  byId('width').addEventListener('input', e => editSelected('width',Number(e.target.value)));
  byId('thinner').addEventListener('click', () => editSelected('width', Math.max(50,size(state.selected)-5)));
  byId('thicker').addEventListener('click', () => editSelected('width', Math.min(150,size(state.selected)+5)));
  for (const event of ['input','change']) byId('color').addEventListener(event, e => editSelected('color',e.target.value));
  byId('overall').addEventListener('input', e => { state.overall = Number(e.target.value); update(); });
  byId('reset-branch').addEventListener('click', () => { delete state.edits[state.selected]; update(); });
  byId('reset-all').addEventListener('click', () => { state.edits = {}; state.overall = 100; update(); });
  byId('download').addEventListener('click', event => {
    event.preventDefault();
    const xml = new XMLSerializer().serializeToString(svg);
    const url = URL.createObjectURL(new Blob([xml],{type:'image/svg+xml;charset=utf-8'}));
    const a = document.createElement('a');
    a.href = url; a.download = 'SATB2_Seurat_cell_hierarchy.svg';
    document.body.append(a); a.click(); a.remove();
    setTimeout(() => URL.revokeObjectURL(url), 60000);
  });
  render(); syncControls();
  byId('plot').replaceChildren(svg);
  byId('controls').disabled = false;
  byId('status').textContent = '';
}

(async () => {
  try {
    const response = await fetch('cell-hierarchy.svg');
    if (!response.ok) throw new Error('Figure unavailable');
    const doc = new DOMParser().parseFromString(await response.text(),'image/svg+xml');
    if (doc.documentElement.localName !== 'svg') throw new Error('Invalid figure');
    mount(document.importNode(doc.documentElement,true));
  } catch {
    byId('status').textContent = 'The editor could not load. Reload to try again, or download the original SVG.';
  }
})();
