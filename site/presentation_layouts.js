/* Presentation geometry only. Source UMAP coordinates and scientific objects remain untouched. */
(function (global) {
  'use strict';
  const finite = value => typeof value === 'number' && Number.isFinite(value);
  const own = (object, key) => object && Object.prototype.hasOwnProperty.call(object, key);
  const clamp = (value, lo, hi, fallback) => Math.max(lo, Math.min(hi, finite(value) ? value : fallback));
  const MODES = Object.freeze(['original', 'compact', 'hippocampus']);
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
    return {mode, basis, spacing: clamp(options.spacing, .25, 2.5, 1), rotation, flipX: options.flipX === true, adjustments};
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
    if (options.mode === 'hippocampus') return 'SCHEMATIC · Hippocampus-like presentation only. Cluster positions follow cluster ID order, not measured anatomy. Shapes may be uniformly resized; no spatial mapping or UMAP refit was performed.';
    if (options.mode === 'compact') return options.spacing === 1 && Object.keys(options.adjustments).length === 0 ? 'DISPLAY LAYOUT · Cluster spacing 1× retains original distances. Original UMAP coordinates, nucleus membership, expression and DE results are preserved.' : 'DISPLAY LAYOUT · Cluster center spacing or placement changed for presentation; relative distances are altered and clusters may overlap. Original UMAP coordinates, nucleus membership, expression and DE results are preserved.';
    if (changed) return 'DISPLAY LAYOUT · Rotation, reflection or manual cluster placement applied for presentation. Original UMAP coordinates, nucleus membership, expression and DE results are preserved.';
    return 'Original saved UMAP coordinates.';
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
      const id = String(nuclei[i].cluster_id), c = groups.get(id) || {cluster_id: id, source_x: 0, source_y: 0, n: 0, radius: 0, shape_scale: 1};
      c.source_x += point.x; c.source_y += point.y; c.n++;
      totalX += point.x; totalY += point.y; validN++; groups.set(id, c);
    }
    for (const c of groups.values()) {c.source_x /= c.n; c.source_y /= c.n; c.x = c.source_x; c.y = c.source_y;}
    for (let i = 0; i < nuclei.length; i++) {
      const point = source[i]; if (!finite(point.x) || !finite(point.y)) continue;
      const c = groups.get(String(nuclei[i].cluster_id));
      c.radius = Math.max(c.radius, Math.hypot(point.x - c.source_x, point.y - c.source_y));
    }
    const clusters = [...groups.values()].sort((a, b) => idCompare(a.cluster_id, b.cluster_id));
    let guide = null;
    if (opts.mode === 'compact') {
      const gx = validN ? totalX / validN : 0, gy = validN ? totalY / validN : 0;
      if (opts.spacing !== 1) for (const c of clusters) {c.x = gx + (c.source_x - gx) * opts.spacing; c.y = gy + (c.source_y - gy) * opts.spacing;}
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
        value = localUnchanged ? orient(original.x, original.y) : orient(c.x + (original.x - c.source_x) * c.shape_scale, c.y + (original.y - c.source_y) * c.shape_scale);
      }
      return {...value, nucleus_id: n.nucleus_id, cluster_id: String(n.cluster_id), source_x: original.x, source_y: original.y, display_only: displayOnly, source_coordinates_modified: false};
    });
    const outputClusters = clusters.map(c => ({...c, ...orient(c.x, c.y)}));
    if (guide) guide.points = guide.points.map(point => orient(point.x, point.y));
    return {version: 'presentation-layouts-1.0.0', mode: opts.mode, basis: opts.basis, options: opts, points, clusters: outputClusters, bounds: bounds(points), source_bounds: sourceBounds, source_span: span, source_center: center, guide, display_only: displayOnly, source_coordinates_modified: false, notice: noticeFor(opts, displayOnly)};
  }
  const api = Object.freeze({version: 'presentation-layouts-1.0.0', modes: MODES, transform, normalizeOptions, bounds});
  global.SATB2_LAYOUTS = api;
  if (typeof module !== 'undefined' && module.exports) module.exports = api;
})(typeof window !== 'undefined' ? window : globalThis);
