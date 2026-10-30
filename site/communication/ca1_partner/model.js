/* Pure filtering helpers; source data and original tests are never changed. */
(function (root) {
  'use strict';
  function yes(v) { return v === true || v === 1 || ['true', '1', 'yes'].includes(String(v).toLowerCase()); }
  function num(v) { if (v === null || v === undefined || String(v).trim() === '') return null; const n = Number(v); return Number.isFinite(n) ? n : null; }
  function qualifies(row, threshold) { return yes(row[threshold === '10' ? 'passes_10pct' : 'passes_5pct']); }
  function filter(rows, state) {
    return rows.filter(r => qualifies(r, state.threshold) &&
      (!state.direction || r.direction === state.direction) &&
      (!state.ca1 || r.ca1_population === state.ca1) &&
      (!state.inpop || r.in_population === state.inpop) &&
      (!state.category || r.in_status === state.category) &&
      (!state.majority || yes(r[state.threshold === '10' ? 'majority_mouse_support_10pct' : 'majority_mouse_support_5pct'])) &&
      (!state.replicated || yes(r[state.threshold === '10' ? 'replicated_majority_10pct' : 'replicated_majority_5pct'])) &&
      (!state.query || [r.ligand_complex, r.receptor_complex, r.ca1_significant_genes, r.in_significant_genes].join(' ').toLowerCase().includes(state.query.toLowerCase())));
  }
  const api = {yes, num, qualifies, filter};
  if (typeof module !== 'undefined' && module.exports) module.exports = api;
  root.CA1INModel = api;
})(typeof window !== 'undefined' ? window : globalThis);
