/* Presentation geometry only. Source UMAP coordinates and scientific objects remain untouched. */
(function (global) {
  'use strict';
  const finite = value => typeof value === 'number' && Number.isFinite(value);
  const own = (object, key) => object && Object.prototype.hasOwnProperty.call(object, key);
  const clamp = (value, lo, hi, fallback) => Math.max(lo, Math.min(hi, finite(value) ? value : fallback));
  const MODES = Object.freeze(['original', 'compact', 'packed', 'rounded', 'hippocampus']);
  function normalizeOptions(options = {}) {
    if (!options || typeof options !== 'object' || Array.isArray(options)) options = {};
    const mode = MODES.includes(options.mode) ? options.mode : 'original';
    const basis = options.basis === 'unintegrated' ? 'unintegrated' : 'harmony';
    const rotation = finite(options.rotation) ? ((options.rotation % 360) + 360) % 360 : 0;
    const adjustments = Object.create(null);
    for (const key of Object.keys(options.adjustments || {})) {
      const source = options.adjustments[key];
      if (!source || typeof source !== 'object' || Array.isArray(source)) continue;
      const value = {dx: clamp(source.dx, -2, 2, 0), dy: clamp(source.dy, -2, 2, 0), scale: clamp(source.scale, .1, 3, 1)};
      if (value.dx || value.dy || value.scale !== 1) adjustments[key] = value;
    }
    return {mode, basis, spacing: clamp(options.spacing, .25, 2.5, 1), packCoverage: clamp(options.packCoverage, .8, 1, .95), rotation, flipX: options.flipX === true, adjustments};
  }
  function bounds(points) {
    let xmin = Infinity, xmax = -Infinity, ymin = Infinity, ymax = -Infinity, n = 0;
    for (const point of points) {
      if (!finite(point.x) || !finite(point.y)) continue;
      xmin = Math.min(xmin, point.x); xmax = Math.max(xmax, point.x);
      ymin = Math.min(ymin, point.y); ymax = Math.max(ymax, point.y); n++;
    }
    return n ? {xmin, xmax, ymin, ymax} : null;
  }
  function idCompare(a, b) {
    const aa = Number(a), bb = Number(b);
    if (String(a).trim() && String(b).trim() && finite(aa) && finite(bb) && aa !== bb) return aa - bb;
    return String(a).localeCompare(String(b), 'en', {numeric: true});
  }
  function noticeFor(options, changed) {
    if (options.mode === 'rounded') return 'SCHEMATIC · Within- and between-cluster geometry replaced by packed rounded discs. Nuclei are arranged by stable ID; disc area approximates nucleus count with a minimum readable radius. This is not original UMAP geometry, measured anatomy or spatial mapping. Original coordinates, membership, expression and DE results are preserved.';
    if (options.mode === 'packed') return `DISPLAY LAYOUT · Cluster islands were packed using ${Math.round(options.packCoverage * 100)}% core outlines; peripheral nuclei are retained and can overlap other islands. Intercluster distances and neighbor relationships are artificial. Internal shapes stay unchanged unless resized manually. Original coordinates, membership, expression and DE results are preserved.`;
    if (options.mode === 'hippocampus') return 'SCHEMATIC · Hippocampus-like presentation only. Cluster positions follow cluster ID order, not measured anatomy. Shapes may be uniformly resized; no spatial mapping or UMAP refit was performed.';
    if (options.mode === 'compact') return options.spacing === 1 && Object.keys(options.adjustments).length === 0 ? 'DISPLAY LAYOUT · Cluster spacing 1× retains original distances. Original UMAP coordinates, nucleus membership, expression and DE results are preserved.' : 'DISPLAY LAYOUT · Cluster center spacing or placement changed for presentation; relative distances are altered and clusters may overlap. Original UMAP coordinates, nucleus membership, expression and DE results are preserved.';
    if (changed) return 'DISPLAY LAYOUT · Rotation, reflection or manual cluster placement applied for presentation. Original UMAP coordinates, nucleus membership, expression and DE results are preserved.';
    return 'Original saved UMAP coordinates.';
  }
  // Deterministic greedy packing of robust core discs. Dispersed peripheral
  // nuclei otherwise make a small isolated body occupy the entire figure.
  // Every nucleus is retained and translated with its cluster, including those
  // outside the core outline; peripheral points may overlap other islands.
  function packClusters(clusters, span, center, spacing, coverage) {
    if (!clusters.length) return {gap: 0, order: []};
    const gap = span * .012 * spacing, minimumRadius = span * .002;
    for (const cluster of clusters) {
      cluster.distances.sort((a, b) => a - b);
      cluster.core_radius = finite(cluster.layout_radius) ? cluster.layout_radius : cluster.distances[Math.max(0, Math.ceil(cluster.distances.length * coverage) - 1)] || 0;
    }
    const sorted = clusters.slice().sort((a, b) => b.core_radius - a.core_radius || idCompare(a.cluster_id, b.cluster_id));
    const placed = [], eps = span * 1e-9;
    for (const cluster of sorted) {
      const radius = Math.max(cluster.core_radius, minimumRadius), candidates = [];
      if (!placed.length) candidates.push({x: 0, y: 0});
      // Cardinal directions give a valid fallback even for one existing disc.
      for (const a of placed) {
        const distance = a.pack_radius + radius + gap;
        for (let k = 0; k < 24; k++) {
          const theta = k * Math.PI / 12;
          candidates.push({x: a.x + distance * Math.cos(theta), y: a.y + distance * Math.sin(theta)});
        }
      }
      for (let i = 0; i < placed.length; i++) for (let j = i + 1; j < placed.length; j++) {
        const a = placed[i], b = placed[j], dx = b.x - a.x, dy = b.y - a.y, distance = Math.hypot(dx, dy);
        const ra = a.pack_radius + radius + gap, rb = b.pack_radius + radius + gap;
        if (distance < eps || distance > ra + rb + eps || distance < Math.abs(ra - rb) - eps) continue;
        const along = (ra * ra - rb * rb + distance * distance) / (2 * distance);
        const height = Math.sqrt(Math.max(0, ra * ra - along * along)), ux = dx / distance, uy = dy / distance;
        candidates.push({x: a.x + along * ux - height * uy, y: a.y + along * uy + height * ux}, {x: a.x + along * ux + height * uy, y: a.y + along * uy - height * ux});
      }
      let best = null, bestScore = Infinity;
      for (const candidate of candidates) {
        if (placed.some(a => Math.hypot(candidate.x - a.x, candidate.y - a.y) < radius + a.pack_radius + gap - eps)) continue;
        let xmin = candidate.x - radius, xmax = candidate.x + radius, ymin = candidate.y - radius, ymax = candidate.y + radius;
        for (const a of placed) {xmin = Math.min(xmin, a.x - a.pack_radius); xmax = Math.max(xmax, a.x + a.pack_radius); ymin = Math.min(ymin, a.y - a.pack_radius); ymax = Math.max(ymax, a.y + a.pack_radius);}
        const width = xmax - xmin, height = ymax - ymin;
        // Prioritize a compact, nearly square figure, then deterministic tie breaks.
        const score = Math.max(width, height) ** 2 + .1 * width * height + .00001 * (candidate.x ** 2 + candidate.y ** 2);
        if (score < bestScore - eps * eps) {best = candidate; bestScore = score;}
      }
      // Finite candidate construction guarantees a solution at the outer edge;
      // keep an explicit safe fallback for degenerate floating point geometry.
      if (!best) best = {x: placed.reduce((v, a) => Math.max(v, a.x + a.pack_radius), 0) + radius + gap, y: 0};
      cluster.x = best.x; cluster.y = best.y; cluster.pack_radius = radius; placed.push(cluster);
    }
    const xmin = Math.min(...placed.map(c => c.x - c.pack_radius)), xmax = Math.max(...placed.map(c => c.x + c.pack_radius));
    const ymin = Math.min(...placed.map(c => c.y - c.pack_radius)), ymax = Math.max(...placed.map(c => c.y + c.pack_radius));
    const shiftX = center.x - (xmin + xmax) / 2, shiftY = center.y - (ymin + ymax) / 2;
    for (const c of placed) {c.x += shiftX; c.y += shiftY;}
    return {gap, coverage, order: placed.map(c => c.cluster_id), peripheral_nuclei_retained: true};
  }
  // Convert a requested display-space move into the stored source-span units.
  // This keeps a drag or directional nudge intuitive after rotate/mirror.
  function displayDeltaToSource(layout, dx, dy) {
    if (!layout || !finite(dx) || !finite(dy) || !finite(layout.source_span) || layout.source_span <= 0) throw Error('A finite layout and display delta are required');
    const options = layout.options || {}, radians = (options.rotation || 0) * Math.PI / 180, cos = Math.cos(radians), sin = Math.sin(radians);
    const sourceX = (dx * cos + dy * sin) * (options.flipX ? -1 : 1), sourceY = -dx * sin + dy * cos;
    return {dx: sourceX / layout.source_span, dy: sourceY / layout.source_span};
  }
  function transform(nuclei, options = {}) {
    if (!Array.isArray(nuclei)) throw Error('Presentation layout requires a nucleus array');
    const opts = normalizeOptions(options), xKey = 'x_' + opts.basis, yKey = 'y_' + opts.basis;
    const source = nuclei.map(n => ({x: finite(n[xKey]) ? n[xKey] : null, y: finite(n[yKey]) ? n[yKey] : null}));
    const sourceBounds = bounds(source), sourceSpan = sourceBounds ? Math.max(sourceBounds.xmax - sourceBounds.xmin, sourceBounds.ymax - sourceBounds.ymin, 1e-9) : 1;
    const span = sourceSpan > 1e-8 ? sourceSpan : 1;
    const center = sourceBounds ? {x: (sourceBounds.xmin + sourceBounds.xmax) / 2, y: (sourceBounds.ymin + sourceBounds.ymax) / 2} : {x: 0, y: 0};
    const groups = new Map(); let totalX = 0, totalY = 0, validN = 0;
    for (let i = 0; i < nuclei.length; i++) {
      const point = source[i]; if (!finite(point.x) || !finite(point.y)) continue;
      const id = String(nuclei[i].cluster_id), c = groups.get(id) || {cluster_id: id, source_x: 0, source_y: 0, n: 0, radius: 0, shape_scale: 1, distances: [], indices: []};
      c.source_x += point.x; c.source_y += point.y; c.n++;
      if (opts.mode === 'rounded') c.indices.push(i);
      totalX += point.x; totalY += point.y; validN++; groups.set(id, c);
    }
    for (const c of groups.values()) {c.source_x /= c.n; c.source_y /= c.n; c.x = c.source_x; c.y = c.source_y;}
    for (let i = 0; i < nuclei.length; i++) {
      const point = source[i]; if (!finite(point.x) || !finite(point.y)) continue;
      const c = groups.get(String(nuclei[i].cluster_id));
      const distance = Math.hypot(point.x - c.source_x, point.y - c.source_y);
      c.radius = Math.max(c.radius, distance); if (opts.mode === 'packed') c.distances.push(distance);
    }
    const clusters = [...groups.values()].sort((a, b) => idCompare(a.cluster_id, b.cluster_id));
    let guide = null, packing = null, schematic = null, schematicOffsets = null;
    if (opts.mode === 'compact') {
      const gx = validN ? totalX / validN : 0, gy = validN ? totalY / validN : 0;
      if (opts.spacing !== 1) for (const c of clusters) {c.x = gx + (c.source_x - gx) * opts.spacing; c.y = gy + (c.source_y - gy) * opts.spacing;}
    } else if (opts.mode === 'packed') {
      packing = packClusters(clusters, span, center, opts.spacing, opts.packCoverage);
    } else if (opts.mode === 'rounded') {
      const maxN = Math.max(1, ...clusters.map(c => c.n)), goldenAngle = Math.PI * (3 - Math.sqrt(5));
      schematicOffsets = new Array(nuclei.length);
      for (const c of clusters) {
        c.layout_radius = span * Math.max(.025, .12 * Math.sqrt(c.n / maxN));
        c.indices.sort((a, b) => {const aa = String(nuclei[a].nucleus_id), bb = String(nuclei[b].nucleus_id); return aa < bb ? -1 : aa > bb ? 1 : a - b;});
        c.indices.forEach((index, rank) => {
          const radius = c.layout_radius * Math.sqrt((rank + .5) / c.n), theta = rank * goldenAngle;
          schematicOffsets[index] = {x: radius * Math.cos(theta), y: radius * Math.sin(theta)};
        });
      }
      packing = packClusters(clusters, span, center, opts.spacing, 1);
      schematic = {kind: 'rounded_cluster_discs', anatomical_mapping: false, original_umap_geometry: false, within_cluster_geometry_preserved: false, nucleus_order: 'stable nucleus ID ascending', disc_area: 'proportional to nucleus count with a minimum readable radius'};
    } else if (opts.mode === 'hippocampus') {
      const start = 50 * Math.PI / 180, end = 310 * Math.PI / 180, radius = span * .45;
      const angleStep = clusters.length > 1 ? (end - start) / (clusters.length - 1) : 0;
      const targetRadius = clusters.length > 1 ? .38 * 2 * radius * Math.sin(angleStep / 2) : radius * .35;
      clusters.forEach((c, index) => {
        const theta = clusters.length > 1 ? start + index * angleStep : Math.PI;
        c.x = clusters.length > 1 ? center.x + radius * Math.cos(theta) : center.x;
        c.y = clusters.length > 1 ? center.y + radius * Math.sin(theta) : center.y;
        c.shape_scale = c.radius ? Math.min(1, targetRadius / c.radius) : 1;
      });
      guide = {kind: 'hippocampus_like_schematic_arc', anatomical_mapping: false, order: clusters.map(c => c.cluster_id), points: Array.from({length: 121}, (_, i) => {const theta = start + (end - start) * i / 120; return {x: center.x + radius * Math.cos(theta), y: center.y + radius * Math.sin(theta)};})};
    }
    for (const c of clusters) {
      const adjustment = own(opts.adjustments, c.cluster_id) ? opts.adjustments[c.cluster_id] : {dx: 0, dy: 0, scale: 1};
      c.x += adjustment.dx * span; c.y += adjustment.dy * span; c.shape_scale *= adjustment.scale;
    }
    const radians = opts.rotation * Math.PI / 180, cos = Math.cos(radians), sin = Math.sin(radians);
    function orient(x, y) {
      if (!opts.rotation && !opts.flipX) return {x, y};
      const dx = (x - center.x) * (opts.flipX ? -1 : 1), dy = y - center.y;
      return {x: center.x + dx * cos - dy * sin, y: center.y + dx * sin + dy * cos};
    }
    const displayOnly = opts.mode !== 'original' || !!opts.rotation || opts.flipX || Object.keys(opts.adjustments).length > 0;
    const points = nuclei.map((n, i) => {
      const original = source[i], c = groups.get(String(n.cluster_id)); let value = {x: null, y: null};
      if (finite(original.x) && finite(original.y)) {
        // The unchanged path must not perform centroid arithmetic: retain exact source doubles.
        const localUnchanged = c.x === c.source_x && c.y === c.source_y && c.shape_scale === 1;
        const local = schematicOffsets?.[i];
        value = local ? orient(c.x + local.x * c.shape_scale, c.y + local.y * c.shape_scale) : localUnchanged ? orient(original.x, original.y) : orient(c.x + (original.x - c.source_x) * c.shape_scale, c.y + (original.y - c.source_y) * c.shape_scale);
      }
      return {...value, nucleus_id: n.nucleus_id, cluster_id: String(n.cluster_id), source_x: original.x, source_y: original.y, display_only: displayOnly, source_coordinates_modified: false};
    });
    const outputClusters = clusters.map(c => {const {distances, indices, ...metadata} = c; return {...metadata, ...orient(c.x, c.y)};});
    if (guide) guide.points = guide.points.map(point => orient(point.x, point.y));
    return {version: 'presentation-layouts-1.2.0', mode: opts.mode, basis: opts.basis, options: opts, points, clusters: outputClusters, bounds: bounds(points), source_bounds: sourceBounds, source_span: span, source_center: center, guide, packing, schematic, display_only: displayOnly, source_coordinates_modified: false, notice: noticeFor(opts, displayOnly)};
  }
  const api = Object.freeze({version: 'presentation-layouts-1.2.0', modes: MODES, transform, normalizeOptions, bounds, displayDeltaToSource});
  global.SATB2_LAYOUTS = api;
  if (typeof module !== 'undefined' && module.exports) module.exports = api;
})(typeof window !== 'undefined' ? window : globalThis);
