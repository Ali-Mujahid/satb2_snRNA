(function (root, factory) {
  if (typeof module === 'object' && module.exports) module.exports = factory();
  else root.Satb2ExonModel = factory();
})(typeof globalThis !== 'undefined' ? globalThis : this, function () {
  'use strict';
  const metrics = {
    candidate_exon: 'Candidate exon E4',
    E3_to_E5_skip_E4: 'E3 → E5: skips E4',
    E3_to_E4: 'E3 → E4: includes E4',
    E4_to_E5: 'E4 → E5: includes E4',
    E5_to_E6_control: 'E5 → E6: control junction'
  };
  const denominators = {
    distinct_cb_ub: 'Distinct CB/UB pairs',
    per_1000_nuclei: 'Per 1,000 matched nuclei',
    per_1000_locus: 'Per 1,000 Satb2-locus pairs'
  };
  const expected = ['cKO_933','cKO_991','cKO_1138','cKO_1139','cKO_1160','cKO_1165',
    'wt_1122','wt_1123','wt_1178','wt_1220','wt_1226','wt_1228'];
  function num(s) {
    if (s === '' || s == null || s === 'NA') return null;
    const n = Number(s);
    if (!Number.isFinite(n) || n < 0) throw Error('Non-negative numeric value expected');
    return n;
  }
  function parse(tsv) {
    // Preserve trailing tab-delimited empty fields in the last data row.
    const lines = String(tsv).replace(/\r\n/g, '\n').replace(/\n+$/, '').split('\n');
    if (lines.length < 2) throw Error('The table has no records');
    const headers = lines[0].split('\t');
    if (new Set(headers).size !== headers.length) throw Error('Duplicate table columns');
    const required = ['sample','condition','cohort','n_matched_nuclei','locus_distinct_cb_ub','coverage_status'];
    Object.keys(metrics).forEach(m => Object.keys(denominators).forEach(d => required.push(m+'_'+d)));
    const missing = required.filter(k => !headers.includes(k));
    if (missing.length) throw Error('Missing columns: ' + missing.join(', '));
    const rows = lines.slice(1).filter(Boolean).map((line, ix) => {
      const values = line.split('\t');
      if (values.length !== headers.length) throw Error('Malformed row ' + (ix+2));
      const row = Object.fromEntries(headers.map((h, i) => [h,values[i]]));
      if (!expected.includes(row.sample) || !['all_CA1','exclude_C11'].includes(row.cohort) ||
          row.condition !== (row.sample.startsWith('wt_') ? 'WT' : 'cKO')) throw Error('Unexpected sample, condition or cohort');
      const cells = num(row.n_matched_nuclei), locus = num(row.locus_distinct_cb_ub);
      if (!Number.isInteger(cells) || !Number.isInteger(locus)) throw Error('Invalid denominator');
      row.n_matched_nuclei=cells; row.locus_distinct_cb_ub=locus;
      const available = row.coverage_status === 'evidence_present_requires_review';
      if (!available && row.coverage_status !== 'unavailable_or_insufficient_locus_coverage') throw Error('Unknown evidence status');
      for (const metric of Object.keys(metrics)) {
        const raw = num(row[metric+'_distinct_cb_ub']);
        if (!Number.isInteger(raw) || raw > locus) throw Error('Invalid raw count');
        row[metric+'_distinct_cb_ub']=raw;
        for (const den of ['per_1000_nuclei','per_1000_locus']) {
          const k=metric+'_'+den, v=num(row[k]); row[k]=v;
          if (!available && v !== null) throw Error('Unavailable row must not report a rate');
          if (available) {
            const base=den==='per_1000_nuclei'?cells:locus;
            if (!base || v === null || Math.abs(v-1000*raw/base) > 1e-6*Math.max(1,v)) throw Error('Inconsistent numerator and denominator');
          }
        }
      }
      return row;
    });
    const keys=new Set(rows.map(r=>r.sample+'|'+r.cohort));
    if (rows.length!==24 || keys.size!==24 || expected.some(s=>!keys.has(s+'|all_CA1') || !keys.has(s+'|exclude_C11'))) throw Error('Expected exactly 12 mice in each cohort');
    return {headers,rows};
  }
  function visible(table,cohort,conditions,sampleSet='all_12') {
    if (!['all_12','exclude_933'].includes(sampleSet)) throw Error('Unknown mouse selection');
    return table.rows.filter(r=>r.cohort===cohort && conditions.includes(r.condition) &&
      (sampleSet!=='exclude_933' || r.sample!=='cKO_933'));
  }
  function selectionLabel(rows,sampleSet) {
    if (!['all_12','exclude_933'].includes(sampleSet)) throw Error('Unknown mouse selection');
    const count=condition=>rows.filter(r=>r.condition===condition).length;
    return (sampleSet==='all_12'?'All 12 eligible':'cKO_933 excluded')+' · '+count('WT')+' WT + '+count('cKO')+' cKO selected';
  }
  function plotValue(row,metric,den) {if(row.coverage_status!=='evidence_present_requires_review')return null;return row[metric+'_'+den];}
  function toTsv(table,rows) {
    return table.headers.join('\t')+'\n'+rows.map(r=>table.headers.map(h=>r[h] == null ? '' : String(r[h]).replace(/[\r\n\t]/g,' ')).join('\t')).join('\n')+'\n';
  }
  return {metrics,denominators,expected,parse,visible,selectionLabel,plotValue,toTsv};
});
