/* Dynamic SVG plots for the saved SATB2 LIANA results. No biological fits are run here. */
(function (root) {
  'use strict';
  var ARROW = '\u2192';
  var palettes = {
    viridis: ['#440154', '#3b528b', '#21918c', '#5ec962', '#fde725'],
    blue: ['#f1f6fb', '#aacbe4', '#4b93c3', '#08519c'],
    rdBu: ['#2166ac', '#92c5de', '#f7f7f7', '#f4a582', '#b2182b']
  };
  function escape(value) {
    return String(value == null ? '' : value).replace(/[&<>"']/g, function (c) {
      return { '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;', "'": '&#39;' }[c];
    });
  }
  function num(value) {
    if (value === null || value === undefined || value === '') return null;
    var n = Number(value);
    return Number.isFinite(n) ? n : null;
  }
  function clamp(value, min, max) { return Math.max(min, Math.min(max, value)); }
  function setting(options, key, fallback, min, max) {
    var n = num(options[key]);
    return n == null ? fallback : clamp(n, min, max);
  }
  function pretty(value) {
    if (value == null || !Number.isFinite(Number(value))) return 'NA';
    var n = Number(value), a = Math.abs(n);
    if (n === 0) return '0';
    if (a < .001 || a >= 100000) return n.toExponential(2);
    if (Number.isInteger(n)) return n.toLocaleString('en-US');
    return n.toPrecision(3).replace(/\.0+$/, '').replace(/(\.\d*?)0+$/, '$1');
  }
  function extent(values, includeZero) {
    var valid = values.map(num).filter(function (v) { return v !== null; });
    if (!valid.length) return [0, 1];
    var lo = valid.reduce(function (a, b) { return Math.min(a, b); }, Infinity), hi = valid.reduce(function (a, b) { return Math.max(a, b); }, -Infinity);
    if (includeZero) lo = Math.min(0, lo);
    if (hi === lo) {
      if (hi === 0) hi = 1;
      else if (hi > 0) lo = 0;
      else { lo -= .5; hi += .5; }
    }
    return [lo, hi];
  }
  function domain(values, ctx, includeZero) {
    var o = ctx.options || {};
    if (o.scaleMode === 'manual' && num(o.colorMin) !== null && num(o.colorMax) !== null && Number(o.colorMax) > Number(o.colorMin)) return [Number(o.colorMin), Number(o.colorMax)];
    if (Array.isArray(ctx.colorDomain) && num(ctx.colorDomain[0]) !== null && num(ctx.colorDomain[1]) !== null && Number(ctx.colorDomain[1]) > Number(ctx.colorDomain[0])) return [Number(ctx.colorDomain[0]), Number(ctx.colorDomain[1])];
    return extent(values, includeZero);
  }
  function rgb(hex, fallback) {
    var h = /^#[0-9a-f]{6}$/i.test(String(hex)) ? String(hex) : fallback;
    return [parseInt(h.slice(1, 3), 16), parseInt(h.slice(3, 5), 16), parseInt(h.slice(5, 7), 16)];
  }
  function color(value, min, max, options) {
    var n = num(value); if (n === null) return '#e5e9ee';
    var o = options || {}, colors = palettes[o.palette] || palettes.viridis;
    if (o.palette === 'custom') colors = [o.lowColor || '#2166ac', o.midColor || '#f7f7f7', o.highColor || '#b2182b'];
    var t = max > min ? clamp((n - min) / (max - min), 0, 1) : .5;
    if (o.reverseScale) t = 1 - t;
    var pos = t * (colors.length - 1), a = Math.min(colors.length - 2, Math.floor(pos)), f = pos - a;
    var c1 = rgb(colors[a], '#2166ac'), c2 = rgb(colors[a + 1], '#b2182b');
    return '#' + c1.map(function (c, i) { return Math.round(c + (c2[i] - c) * f).toString(16).padStart(2, '0'); }).join('');
  }
  function lrKey(row) { return String(row.ligand_complex) + ARROW + String(row.receptor_complex); }
  function pairKey(row) { return String(row.source) + ARROW + String(row.target); }
  function rank(row) { var n = num(row.magnitude_rank); return n == null ? Infinity : n; }
  function baseContext(ctx) {
    return Object.assign({ nodes: [], eligiblePairs: [], options: {}, condition: '', scope: '', direction: '' }, ctx || {});
  }
  function metricInfo(metric) {
    return {
      count: { label: 'Distinct ordered ligand\u2013receptor pairs', note: 'Counts within the current filters; expression-compatible candidates.' },
      mean_gmeans: { label: 'Mean native LR geometric mean', note: 'Mean of lr_gmeans across distinct qualifying molecular pairs.' },
      magnitude_rank: { label: 'LIANA magnitude rank', note: 'Lower = stronger consensus priority; not a genotype FDR.' },
      lr_gmeans: { label: 'Native LR geometric mean', note: 'Geometric mean of endpoint mean log-normalized expression.' },
      lr_logfc: { label: 'Native LR group specificity (logFC)', note: 'Cell group versus rest; not the cKO / WT fold change.' },
      cellphone_pvals: { label: 'CellPhoneDB specificity P value', note: 'Cell-label permutation P value; not a genotype interaction test.' }
    }[metric] || { label: metric, note: '' };
  }
  function startSvg(width, height, title, description, fontSize) {
    return '<svg xmlns="http://www.w3.org/2000/svg" width="' + width + '" height="' + height + '" viewBox="0 0 ' + width + ' ' + height + '" role="img" aria-label="' + escape(title) + '">' +
      '<title>' + escape(title) + '</title><desc>' + escape(description) + '</desc>' +
      '<style>text{font-family:Arial,Helvetica,sans-serif;fill:#263f50;font-size:' + fontSize + 'px}.title{font-size:' + (fontSize + 6) + 'px;font-weight:700}.muted{fill:#5d7381}.small{font-size:' + Math.max(10, fontSize - 1) + 'px}.axis{font-weight:700}.grid{stroke:#dce5eb;stroke-width:1}.cell-number{font-size:' + Math.max(8, fontSize - 2) + 'px}</style><rect width="' + width + '" height="' + height + '" fill="#ffffff"/>';
  }
  function text(x, y, str, attrs) { return '<text x="' + x + '" y="' + y + '"' + (attrs ? ' ' + attrs : '') + '>' + escape(str) + '</text>'; }
  function truncate(value, length) { var s = String(value); return s.length > length ? s.slice(0, length - 1) + '\u2026' : s; }
  function shortNode(node, limit) {
    var label = node.label;
    if (!label) {
      if (node.id === 'CA1_Deep') label = 'CA1 Deep';
      else if (node.id === 'CA1_Superficial') label = 'CA1 Sup.';
      else if (node.id === 'IN_Total') label = 'IN Total';
      else label = node.id.replace(/^IN_/, '');
    }
    return truncate(label, limit || 18);
  }
  function legend(x, y, width, limits, o, label, note) {
    var s = text(x, y, label, 'class="axis"');
    for (var i = 0; i < 80; i++) {
      var v = limits[0] + (limits[1] - limits[0]) * i / 79;
      s += '<rect x="' + (x + width * i / 80).toFixed(2) + '" y="' + (y + 9) + '" width="' + (width / 80 + .2).toFixed(2) + '" height="13" fill="' + color(v, limits[0], limits[1], o) + '"/>';
    }
    [0, .5, 1].forEach(function (t) {
      s += text(x + width * t, y + 39, pretty(limits[0] + (limits[1] - limits[0]) * t), 'text-anchor="' + (t === 0 ? 'start' : t === 1 ? 'end' : 'middle') + '" class="small"');
    });
    if (note) s += text(x, y + 59, note, 'class="small muted"');
    return s;
  }
  function emptySvg(title, note, ctx) {
    var f = setting(ctx.options || {}, 'fontSize', 12, 8, 24);
    return startSvg(1050, 230, title, note, f) + text(35, 42, title, 'class="title"') + text(35, 95, 'No qualifying interactions in this view.', 'class="axis"') + text(35, 124, note, 'class="muted"') + '</svg>';
  }
  function aggregateHeatmap(rows, ctx) {
    ctx = baseContext(ctx);
    var groups = new Map(), eligible = new Set(ctx.eligiblePairs.map(pairKey));
    rows.forEach(function (r) {
      var p = pairKey(r);
      if (!groups.has(p)) groups.set(p, new Map());
      if (!groups.get(p).has(lrKey(r))) groups.get(p).set(lrKey(r), r);
    });
    var cells = [];
    ctx.nodes.forEach(function (source) {
      ctx.nodes.forEach(function (target) {
        var key = source.id + ARROW + target.id, group = groups.get(key), values = group ? Array.from(group.values()) : [];
        var means = values.map(function (r) { return num(r.lr_gmeans); }).filter(function (n) { return n !== null; });
        var tested = eligible.has(key);
        cells.push({ source: source.id, target: target.id, tested: tested, count: tested ? values.length : null, mean_gmeans: tested && means.length ? means.reduce(function (s, v) { return s + v; }, 0) / means.length : null });
      });
    });
    return cells;
  }
  function heatmap(rows, suppliedCtx) {
    var ctx = baseContext(suppliedCtx), o = ctx.options, nodes = ctx.nodes;
    var title = 'Ligand\u2013receptor heatmap \u00b7 ' + (ctx.condition || 'selected context');
    if (!nodes.length) return emptySvg(title, 'Select source and target populations to display the matrix.', ctx);
    var f = setting(o, 'fontSize', 12, 8, 24), metric = o.heatMetric === 'mean_gmeans' ? 'mean_gmeans' : 'count';
    var cells = aggregateHeatmap(rows, ctx), info = metricInfo(metric), limits = domain(cells.map(function (c) { return c[metric]; }), ctx, true);
    var size = Math.max(24, Math.min(43, Math.floor(730 / nodes.length))), left = Math.max(185, Math.ceil(18 * f * .55) + 45), top = 115, side = size * nodes.length;
    var labelDepth = Math.max(142, Math.ceil(18 * f * .55 * .82) + 40);
    var width = Math.max(1050, left + side + 100), legendY = top + side + labelDepth, height = legendY + 117;
    var svg = startSvg(width, height, title, 'Source populations are rows and target populations are columns. ' + info.note, f);
    svg += text(28, 38, title, 'class="title"');
    svg += text(28, 64, info.note, 'class="muted small"');
    svg += text(left + side / 2, top - 20, 'Receiver / target', 'text-anchor="middle" class="axis"');
    svg += '<text transform="translate(28 ' + (top + side / 2) + ') rotate(-90)" text-anchor="middle" class="axis">Sender / source</text>';
    cells.forEach(function (cell) {
      var row = nodes.findIndex(function (n) { return n.id === cell.source; }), col = nodes.findIndex(function (n) { return n.id === cell.target; });
      var x = left + col * size, y = top + row * size, value = cell[metric];
      var fill = !cell.tested ? '#e5e9ee' : value === null ? '#fbfcfd' : color(value, limits[0], limits[1], o);
      var tooltip = cell.source + ' \u2192 ' + cell.target + '\n' + (!cell.tested ? 'NA: this ordered direction was not tested in the selected context.' : metric === 'count' ? 'Distinct qualifying LR identities: ' + cell.count + (cell.count === 0 ? '\nEligible pair; zero candidates passed these filters. This does not establish biological absence.' : '') : 'Mean lr_gmeans: ' + pretty(value) + '\nQualifying LR identities: ' + cell.count + (value === null ? '\nNo qualifying geometric mean is available.' : ''));
      svg += '<g data-source="' + escape(cell.source) + '" data-target="' + escape(cell.target) + '" data-tested="' + cell.tested + '" data-value="' + (value === null ? '' : value) + '"><title>' + escape(tooltip) + '</title><rect x="' + x + '" y="' + y + '" width="' + size + '" height="' + size + '" fill="' + fill + '" fill-opacity="' + setting(o, 'opacity', .95, .1, 1) + '" stroke="#ffffff" stroke-width="1"/>';
      if (size >= 27) {
        var displayed = !cell.tested ? 'NA' : value === null ? '\u2014' : pretty(value);
        var c = rgb(fill, '#ffffff'), luminance = .2126 * c[0] + .7152 * c[1] + .0722 * c[2];
        svg += text(x + size / 2, y + size / 2 + 3, displayed, 'text-anchor="middle" class="cell-number" style="fill:' + (luminance < 125 && cell.tested && value !== null ? '#ffffff' : '#263f50') + '"');
      }
      svg += '</g>';
    });
    nodes.forEach(function (node, i) {
      var cx = left + (i + .5) * size, cy = top + (i + .5) * size;
      svg += '<g><title>' + escape(node.label || node.id) + '</title>';
      svg += '<rect x="' + (left - 11) + '" y="' + (top + i * size + 2) + '" width="5" height="' + (size - 4) + '" fill="' + (/^#[0-9a-f]{6}$/i.test(node.color || '') ? node.color : '#6d8797') + '"/>';
      svg += text(left - 19, cy + 4, shortNode(node), 'text-anchor="end"');
      svg += '<text transform="translate(' + cx + ' ' + (top + side + 12) + ') rotate(-55)" text-anchor="end">' + escape(shortNode(node)) + '</text></g>';
    });
    svg += legend(left, legendY, Math.min(410, side), limits, o, info.label, metric === 'count' ? 'Zero = eligible pair with no candidates passing the current filters.' : 'Empty eligible cells have no mean; they are not zero expression.');
    svg += '<rect x="' + (left + 465) + '" y="' + (legendY + 9) + '" width="17" height="13" fill="#e5e9ee"/>' + text(left + 489, legendY + 20, 'NA: untested direction', 'class="small"');
    svg += text(28, height - 22, 'Display summaries only. WT\u2013cKO interaction differences are descriptive; no interaction-level genotype test was performed.', 'class="muted small"');
    return svg + '</svg>';
  }
  function selectDotplot(rows, ctx) {
    ctx = baseContext(ctx);
    var o = ctx.options, lrs = new Map(), pairs = new Map(), index = new Map();
    rows.forEach(function (row, i) {
      var l = lrKey(row), p = pairKey(row), r = rank(row);
      if (!lrs.has(l)) lrs.set(l, { key: l, ligand: row.ligand_complex, receptor: row.receptor_complex, rank: r, order: i });
      else lrs.get(l).rank = Math.min(lrs.get(l).rank, r);
      if (!pairs.has(p)) pairs.set(p, { key: p, source: row.source, target: row.target, rank: r, order: i });
      else pairs.get(p).rank = Math.min(pairs.get(p).rank, r);
      var k = l + '\u0000' + p;
      if (!index.has(k) || rank(index.get(k)) > r) index.set(k, row);
    });
    function compare(a, b) { return o.order === 'alphabetical' ? a.key.localeCompare(b.key) : a.rank - b.rank || a.order - b.order; }
    var allLrs = Array.from(lrs.values()).sort(compare), allPairs = Array.from(pairs.values()).sort(compare);
    var orderedLrs = Array.isArray(ctx.lrOrder) ? ctx.lrOrder.map(function (key) { var parts = key.split(ARROW); return lrs.get(key) || { key: key, ligand: parts[0], receptor: parts.slice(1).join(ARROW), rank: Infinity, order: 0 }; }) : allLrs;
    var orderedPairs = Array.isArray(ctx.pairOrder) ? ctx.pairOrder.map(function (key) { var parts = key.split(ARROW); return pairs.get(key) || { key: key, source: parts[0], target: parts.slice(1).join(ARROW), rank: Infinity, order: 0 }; }) : allPairs;
    return { lrs: orderedLrs.slice(0, setting(o, 'topN', 24, 1, 120)), pairs: orderedPairs.slice(0, setting(o, 'pairLimit', 20, 1, 100)), totalLrs: orderedLrs.length, totalPairs: orderedPairs.length, index: index };
  }
  function tooltip(row, metric) {
    return row.source + ' \u2192 ' + row.target + '\n' + row.ligand_complex + ' \u2192 ' + row.receptor_complex + '\n' +
      metricInfo(metric).label + ': ' + pretty(num(row[metric])) + '\nMagnitude rank: ' + pretty(num(row.magnitude_rank)) + '\nLR geometric mean: ' + pretty(num(row.lr_gmeans)) + '\nLR group-versus-rest logFC: ' + pretty(num(row.lr_logfc)) + '\nCellPhoneDB specificity P: ' + pretty(num(row.cellphone_pvals)) + '\nMinimum detection across ALL complex subunits: ' + pretty(num(row.lr_prop) == null ? null : clamp(Number(row.lr_prop), 0, 1) * 100) + '%\nMouse expression support WT: ' + pretty(num(row.sample_support_WT)) + ' / ' + pretty(num(row.sample_eligible_WT)) + '; cKO: ' + pretty(num(row.sample_support_cKO)) + ' / ' + pretty(num(row.sample_eligible_cKO)) + '\nMatching endpoint DESeq2 support: ' + (row.deg_supported === true || row.deg_supported === 'True' || row.deg_supported === 1 ? 'yes' : 'no');
  }
  function dotplot(rows, suppliedCtx) {
    var ctx = baseContext(suppliedCtx), o = ctx.options, data = selectDotplot(rows, ctx);
    var title = 'Ligand\u2013receptor dotplot \u00b7 ' + (ctx.condition || 'selected context');
    if (!data.lrs.length || !data.pairs.length) return emptySvg(title, 'No pairs pass the selected population, evidence and rank filters.', ctx);
    var f = setting(o, 'fontSize', 12, 8, 24), requestedRadius = setting(o, 'dotRadius', 9, 1, 24);
    var rowHeight = Math.max(26, f + 13, Math.ceil(requestedRadius * 2.2)), colWidth = Math.max(39, f + 29, Math.ceil(requestedRadius * 2.2));
    var longest = data.lrs.reduce(function (m, l) { return Math.max(m, (String(l.ligand) + ' \u2192 ' + String(l.receptor)).length); }, 0);
    var left = clamp(Math.ceil(longest * f * .55) + 35, 215, 520), top = 133;
    var plotWidth = data.pairs.length * colWidth, plotHeight = data.lrs.length * rowHeight;
    var labelDepth = Math.max(153, Math.ceil(31 * f * .55 * .82) + 44);
    var width = Math.max(1080, left + plotWidth + 60), legendY = top + plotHeight + labelDepth, height = legendY + 135;
    var metric = ['magnitude_rank', 'lr_gmeans', 'lr_logfc', 'cellphone_pvals'].indexOf(o.dotMetric) >= 0 ? o.dotMetric : 'magnitude_rank';
    var displayedRows = [];
    data.lrs.forEach(function (l) { data.pairs.forEach(function (p) { var r = data.index.get(l.key + '\u0000' + p.key); if (r) displayedRows.push(r); }); });
    var limits = domain(displayedRows.map(function (r) { return r[metric]; }), ctx, metric !== 'lr_logfc');
    var info = metricInfo(metric), radius = Math.min(requestedRadius, rowHeight * .46, colWidth * .46);
    var svg = startSvg(width, height, title, info.note + ' Dot area encodes the minimum detected fraction across all complex subunits.', f);
    svg += text(28, 38, title, 'class="title"');
    svg += text(28, 64, info.note, 'class="muted small"');
    svg += text(28, 88, data.lrs.length + ' / ' + data.totalLrs + ' LR identities \u00b7 ' + data.pairs.length + ' / ' + data.totalPairs + ' ordered population pairs shown. Full filtered rows remain in the table export.', 'class="muted small"');
    svg += text(left + plotWidth / 2, top - 18, 'Sender \u2192 receiver', 'text-anchor="middle" class="axis"');
    data.lrs.forEach(function (lr, j) {
      var y = top + j * rowHeight;
      svg += '<rect x="' + left + '" y="' + y + '" width="' + plotWidth + '" height="' + rowHeight + '" fill="' + (j % 2 ? '#f6f9fb' : '#ffffff') + '"/>';
      svg += '<g><title>' + escape(lr.ligand + ' \u2192 ' + lr.receptor) + '</title>' + text(left - 13, y + rowHeight / 2 + 4, truncate(lr.ligand + ' \u2192 ' + lr.receptor, Math.floor((left - 30) / (f * .53))), 'text-anchor="end"') + '</g>';
      data.pairs.forEach(function (pair, i) {
        var row = data.index.get(lr.key + '\u0000' + pair.key);
        if (!row) return;
        var prop = num(row.lr_prop), value = num(row[metric]);
        var rr = prop === null ? 0 : radius * Math.sqrt(clamp(prop, 0, 1));
        if (rr > 0) svg += '<circle data-lr="' + escape(lr.key) + '" data-pair="' + escape(pair.key) + '" data-fraction="' + clamp(prop, 0, 1) + '" cx="' + (left + (i + .5) * colWidth) + '" cy="' + (y + rowHeight / 2) + '" r="' + rr.toFixed(4) + '" fill="' + color(value, limits[0], limits[1], o) + '" opacity="' + setting(o, 'opacity', .85, .05, 1) + '"><title>' + escape(tooltip(row, metric)) + '</title></circle>';
      });
    });
    data.pairs.forEach(function (pair, i) {
      var x = left + (i + .5) * colWidth;
      svg += '<line x1="' + (left + i * colWidth) + '" y1="' + top + '" x2="' + (left + i * colWidth) + '" y2="' + (top + plotHeight) + '" class="grid"/>';
      var source = ctx.nodes.find(function (n) { return n.id === pair.source; }) || { id: pair.source }, target = ctx.nodes.find(function (n) { return n.id === pair.target; }) || { id: pair.target };
      svg += '<g><title>' + escape((source.label || source.id) + ' \u2192 ' + (target.label || target.id)) + '</title><text transform="translate(' + x + ' ' + (top + plotHeight + 16) + ') rotate(-55)" text-anchor="end">' + escape(shortNode(source, 14) + ' \u2192 ' + shortNode(target, 14)) + '</text></g>';
    });
    svg += '<rect x="' + left + '" y="' + top + '" width="' + plotWidth + '" height="' + plotHeight + '" fill="none" stroke="#b8cbd7"/>';
    svg += legend(left, legendY, 330, limits, o, info.label, 'Color limits ' + (o.scaleMode === 'manual' ? 'manually set' : ctx.colorDomain ? 'shared across compared contexts' : 'computed from the displayed interactions') + '.');
    var sizeX = Math.max(left + 430, width - 340), baseline = legendY + 26;
    svg += text(sizeX, legendY, 'Minimum all-subunit detection', 'class="axis"');
    [.1, .25, .5, 1].forEach(function (p, i) { var x = sizeX + 24 + i * 63; svg += '<circle cx="' + x + '" cy="' + baseline + '" r="' + (radius * Math.sqrt(p)).toFixed(4) + '" fill="#617d8f"/>' + text(x, legendY + 53, Math.round(p * 100) + '%', 'text-anchor="middle" class="small"'); });
    svg += text(28, height - 43, 'Dot area = detected fraction (0\u2013100%). An absent dot means no qualifying displayed record, not proof of absent communication.', 'class="small muted"');
    svg += text(28, height - 22, 'Saved RNA/prior-compatible candidates; expression support and component DESeq2 evidence do not establish interaction activity.', 'class="small muted"');
    return svg + '</svg>';
  }
  root.CCCPlots = {
    heatmap: heatmap, dotplot: dotplot, color: color,
    helpers: { aggregateHeatmap: aggregateHeatmap, selectDotplot: selectDotplot, extent: extent, escape: escape, lrKey: lrKey, pairKey: pairKey, metricInfo: metricInfo }
  };
})(typeof window !== 'undefined' ? window : globalThis);
