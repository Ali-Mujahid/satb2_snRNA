/* Shared screen/PDF renderer. Descriptive display of unchanged saved DESeq2 rows. */
(function (root) {
  'use strict';
  const finite = value => typeof value === 'number' && Number.isFinite(value);
  const clamp = (value, low, high) => Math.max(low, Math.min(high, value));
  const DEFAULTS = {up: '#d18a36', down: '#397eb5', nonsignificant: '#bcc5cb'};
  const cleanColor = (value, fallback) => /^#[0-9a-f]{6}$/i.test(value || '') ? value.toLowerCase() : /^#[0-9a-f]{3}$/i.test(value || '') ? '#' + value.slice(1).split('').map(x => x + x).join('').toLowerCase() : fallback;
  const tickText = value => Math.abs(value) >= 10000 || (value !== 0 && Math.abs(value) < 0.001) ? value.toExponential(1) : String(Number(value.toPrecision(4)));
  function niceStep(raw) {
    const power = 10 ** Math.floor(Math.log10(Math.max(raw, Number.MIN_VALUE)));
    const unit = raw / power;
    return (unit <= 1 ? 1 : unit <= 2 ? 2 : unit <= 2.5 ? 2.5 : unit <= 5 ? 5 : 10) * power;
  }
  function ticks(low, high, desired) {
    const step = niceStep((high - low) / desired), result = [];
    for (let value = Math.ceil(low / step) * step, n = 0; value <= high + step * 1e-7 && n < 100; value += step, n++) result.push(Math.abs(value) < step * 1e-8 ? 0 : value);
    return result;
  }
  function trimText(context, text, maxWidth) {
    const width = s => typeof context.measureText === 'function' ? context.measureText(s).width : s.length * 6;
    if (width(text) <= maxWidth) return text;
    let end = text.length;
    while (end > 1 && width(text.slice(0, end) + '…') > maxWidth) end--;
    return text.slice(0, end) + '…';
  }
  function prepare(savedRows, settings) {
    const fdr = finite(settings.deFdr) ? clamp(settings.deFdr, 0, 1) : 0.05;
    const lfc = finite(settings.deLfc) ? Math.max(0, settings.deLfc) : 0;
    const metadata = {totalRows: savedRows.length, validRows: 0, unavailableFdrRows: 0, invalidFdrRows: 0, invalidEffectRows: 0, omittedRows: 0, zeroAdjustedPRows: 0, cappedRows: 0, eligibleRows: 0, eligibleGenes: 0, labelCount: 0, requestedLabelCount: 0, omittedLabelsForSpace: 0, displayYCap: 320, thresholds: {fdr, absoluteLog2FC: lfc}};
    const rows = [];
    for (const row of savedRows) {
      const missingP = !finite(row.padj), invalidP = !missingP && (row.padj < 0 || row.padj > 1), invalidEffect = !finite(row.log2FoldChange);
      if (missingP) metadata.unavailableFdrRows++;
      if (invalidP) metadata.invalidFdrRows++;
      if (invalidEffect) metadata.invalidEffectRows++;
      if (missingP || invalidP || invalidEffect) {metadata.omittedRows++; continue;}
      const rawY = -Math.log10(Math.max(row.padj, Number.MIN_VALUE));
      if (row.padj === 0) metadata.zeroAdjustedPRows++;
      if (rawY > 320) metadata.cappedRows++;
      const hit = row.padj <= fdr && Math.abs(row.log2FoldChange) >= lfc;
      const id = String(row.gene_id || row.gene || 'unnamed');
      rows.push({row, id, symbol: String(row.gene || row.gene_symbol || id), x: row.log2FoldChange, y: Math.min(320, Math.max(0, rawY)), hit});
      if (hit) metadata.eligibleRows++;
    }
    metadata.validRows = rows.length;
    const ranked = rows.filter(row => row.hit).sort((a, b) => a.row.padj - b.row.padj || Math.abs(b.x) - Math.abs(a.x) || (a.id < b.id ? -1 : a.id > b.id ? 1 : 0));
    const seen = new Set(), unique = ranked.filter(row => {if (seen.has(row.id)) return false; seen.add(row.id); return true;});
    metadata.eligibleGenes = unique.length;
    metadata.requestedLabelCount = settings.volcanoLabels === false ? 0 : clamp(finite(settings.volcanoLabelCount) ? Math.floor(settings.volcanoLabelCount) : 20, 0, 100);
    const duplicateSymbols = new Map();
    for (const row of rows) {if (!duplicateSymbols.has(row.symbol)) duplicateSymbols.set(row.symbol, new Set()); duplicateSymbols.get(row.symbol).add(row.id);}
    const labels = unique.slice(0, metadata.requestedLabelCount).map(row => ({...row, label: duplicateSymbols.get(row.symbol).size > 1 ? row.symbol + ' [' + row.id + ']' : row.symbol}));
    return {rows, labels, metadata};
  }
  function draw(plot, savedRows, settings) {
    settings = settings || {};
    savedRows = Array.isArray(savedRows) ? savedRows : [];
    const {c, width, height} = plot;
    if (!c || !finite(width) || !finite(height) || width < 240 || height < 220) throw new Error('Volcano plot needs a canvas context and dimensions of at least 240 × 220.');
    const prepared = prepare(savedRows, settings), {rows, labels, metadata} = prepared;
    const colors = {up: cleanColor(settings.volcanoUp, DEFAULTS.up), down: cleanColor(settings.volcanoDown, DEFAULTS.down), nonsignificant: cleanColor(settings.volcanoNS, DEFAULTS.nonsignificant)};
    metadata.colors = colors;
    const left = 72, top = 32, right = width - 22, bottom = height - 91, w = right - left, h = bottom - top;
    let maxX = 1, maxY = 1;
    for (const point of rows) {maxX = Math.max(maxX, Math.abs(point.x)); maxY = Math.max(maxY, point.y);}
    const xStep = niceStep(maxX / 3), roundedX = Math.ceil(maxX / xStep) * xStep, xBound = finite(roundedX) ? roundedX : maxX;
    const yStep = niceStep(maxY / 5), yBound = Math.min(320, Math.ceil(maxY / yStep) * yStep);
    metadata.axisBounds = {x: [-xBound, xBound], y: [0, yBound]};
    metadata.plotBounds = {left, top, right, bottom};
    metadata.ticks = {x: ticks(-xBound, xBound, 6), y: ticks(0, yBound, 5)};
    const mapX = value => left + (value / xBound + 1) * w / 2;
    const mapY = value => bottom - value / yBound * h;
    c.save();
    c.clearRect(0, 0, width, height);
    c.fillStyle = '#ffffff'; c.fillRect(0, 0, width, height);
    c.font = '11px system-ui, sans-serif'; c.lineWidth = 1;
    c.strokeStyle = '#e5ebef'; c.textBaseline = 'middle';
    for (const tick of metadata.ticks.y) {const y = mapY(tick); c.beginPath(); c.moveTo(left, y); c.lineTo(right, y); c.stroke();}
    c.strokeStyle = '#8ba0ae'; c.strokeRect(left, top, w, h);
    c.fillStyle = '#3d5868';
    for (const tick of metadata.ticks.x) {
      const x = mapX(tick); c.beginPath(); c.moveTo(x, bottom); c.lineTo(x, bottom + 5); c.stroke(); c.textAlign = 'center'; c.fillText(tickText(tick), x, bottom + 16);
    }
    for (const tick of metadata.ticks.y) {
      const y = mapY(tick); c.beginPath(); c.moveTo(left - 5, y); c.lineTo(left, y); c.stroke(); c.textAlign = 'right'; c.fillText(tickText(tick), left - 9, y);
    }
    c.font = '12px system-ui, sans-serif'; c.textAlign = 'center'; c.fillText('log2 fold change (cKO / WT)', left + w / 2, height - 49);
    c.save(); c.translate(19, top + h / 2); c.rotate(-Math.PI / 2); c.fillText('−log10 adjusted P', 0, 0); c.restore();
    if (settings.volcanoThresholds !== false) {
      c.save(); c.strokeStyle = '#7c8d97'; c.lineWidth = 0.8; if (c.setLineDash) c.setLineDash([4, 4]);
      for (const x of [...new Set([-metadata.thresholds.absoluteLog2FC, metadata.thresholds.absoluteLog2FC])]) if (Math.abs(x) <= xBound) {c.beginPath(); c.moveTo(mapX(x), top); c.lineTo(mapX(x), bottom); c.stroke();}
      const thresholdY = metadata.thresholds.fdr === 0 ? Infinity : -Math.log10(metadata.thresholds.fdr);
      if (thresholdY <= yBound) {c.beginPath(); c.moveTo(left, mapY(thresholdY)); c.lineTo(right, mapY(thresholdY)); c.stroke();}
      c.restore();
    }
    c.save(); c.beginPath(); c.rect(left, top, w, h); c.clip();
    for (const hit of [false, true]) {
      c.globalAlpha = hit ? 0.85 : 0.4;
      for (const point of rows) if (point.hit === hit) {
        point.px = clamp(mapX(point.x), left + 1.6, right - 1.6); point.py = clamp(mapY(point.y), top + 1.6, bottom - 1.6);
        c.fillStyle = hit ? point.x >= 0 ? colors.up : colors.down : colors.nonsignificant;
        c.beginPath(); c.arc(point.px, point.py, hit ? 1.7 : 1.2, 0, Math.PI * 2); c.fill();
      }
    }
    c.restore();
    const fontSize = clamp(finite(settings.volcanoLabelSize) ? settings.volcanoLabelSize : 10, 8, 18), gap = fontSize + 4;
    const laneCapacity = Math.max(1, Math.floor((h - 8) / gap)), lanes = [[], []];
    // Keep labels away from dense points; spill into the other lane before omitting a label.
    for (const label of labels) {
      label.px = clamp(mapX(label.x), left + 1.6, right - 1.6); label.py = clamp(mapY(label.y), top + 1.6, bottom - 1.6);
      let lane = label.x >= 0 ? 1 : 0;
      if (lanes[lane].length >= laneCapacity) lane = 1 - lane;
      if (lanes[lane].length >= laneCapacity) {metadata.omittedLabelsForSpace++; continue;}
      lanes[lane].push(label);
    }
    c.font = fontSize + 'px system-ui, sans-serif'; c.textBaseline = 'middle';
    const layout = [];
    for (let side = 0; side < 2; side++) {
      const lane = lanes[side].sort((a, b) => a.py - b.py || (a.id < b.id ? -1 : a.id > b.id ? 1 : 0));
      const minY = top + fontSize / 2 + 4, maxLabelY = bottom - fontSize / 2 - 4;
      for (let i = 0; i < lane.length; i++) lane[i].labelY = Math.max(i ? lane[i - 1].labelY + gap : minY, lane[i].py);
      if (lane.length && lane[lane.length - 1].labelY > maxLabelY) {
        lane[lane.length - 1].labelY = maxLabelY;
        for (let i = lane.length - 2; i >= 0; i--) lane[i].labelY = Math.min(lane[i].labelY, lane[i + 1].labelY - gap);
      }
      const labelX = side === 0 ? left + 8 : right - 8, maxWidth = Math.max(20, w / 2 - 24);
      for (const label of lane) {
        const visibleText = trimText(c, label.label, maxWidth);
        const textWidth = typeof c.measureText === 'function' ? c.measureText(visibleText).width : visibleText.length * fontSize * 0.6;
        const box = {left: side === 0 ? labelX : labelX - textWidth, right: side === 0 ? labelX + textWidth : labelX, top: label.labelY - fontSize / 2 - 1, bottom: label.labelY + fontSize / 2 + 1};
        c.strokeStyle = '#687e8d'; c.lineWidth = 0.6; c.beginPath(); c.moveTo(label.px, label.py); c.lineTo(side === 0 ? box.right + 3 : box.left - 3, label.labelY); c.stroke();
        c.fillStyle = 'rgba(255,255,255,0.9)'; c.fillRect(box.left - 2, box.top, box.right - box.left + 4, box.bottom - box.top);
        c.fillStyle = '#233e50'; c.textAlign = side === 0 ? 'left' : 'right'; c.fillText(visibleText, labelX, label.labelY);
        layout.push({gene_id: label.id, gene: label.symbol, label: label.label, visibleText, adjustedP: label.row.padj, log2FoldChange: label.x, x: labelX, y: label.labelY, box, side});
      }
    }
    metadata.labels = layout; metadata.labelCount = layout.length;
    c.font = '10px system-ui, sans-serif'; c.textAlign = 'left'; c.fillStyle = '#3d5868';
    const status = rows.length ? `${metadata.labelCount} / ${metadata.eligibleGenes} eligible genes labeled; FDR ≤ ${tickText(metadata.thresholds.fdr)}, |log2FC| ≥ ${tickText(metadata.thresholds.absoluteLog2FC)}` : 'No rows with finite log2 fold change and adjusted P in [0, 1]';
    c.fillText(trimText(c, status, w), left, 13);
    const narrow = w < 430;
    const legend = [{text: narrow ? 'cKO ↑' : 'Higher in cKO', color: colors.up}, {text: narrow ? 'cKO ↓' : 'Lower in cKO', color: colors.down}, {text: narrow ? 'Other' : 'Below thresholds', color: colors.nonsignificant}];
    let legendX = left;
    for (const item of legend) {c.fillStyle = item.color; c.fillRect(legendX, height - 32, 8, 8); c.fillStyle = '#3d5868'; c.fillText(item.text, legendX + 12, height - 28); legendX += Math.max(narrow ? 58 : 110, item.text.length * 5.4 + 30);}
    const notes = [`Omitted ${metadata.omittedRows} rows: ${metadata.unavailableFdrRows} unavailable FDR, ${metadata.invalidFdrRows} invalid FDR, ${metadata.invalidEffectRows} unavailable effect.`];
    if (metadata.cappedRows) notes.push(`${metadata.cappedRows} point(s) capped at −log10 adjusted P = 320 (${metadata.zeroAdjustedPRows} zero adjusted P).`);
    if (metadata.omittedLabelsForSpace) notes.push(`${metadata.omittedLabelsForSpace} label(s) omitted for space; enlarge the plot or reduce label count.`);
    metadata.note = notes.join(' ');
    c.fillText(trimText(c, metadata.note, width - left - 10), left, height - 9);
    c.restore();
    return metadata;
  }
  root.SATB2_VOLCANO = {draw, prepare};
})(typeof window === 'undefined' ? globalThis : window);
