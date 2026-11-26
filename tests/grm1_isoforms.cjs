const fs = require('node:fs');
const path = require('node:path');
const assert = require('node:assert/strict');
const root = path.join(__dirname, '..', 'site', 'grm1_isoforms');
const html = fs.readFileSync(path.join(root, 'index.html'), 'utf8');
const match = html.match(/<script id="data" type="application\/json">([\s\S]*?)<\/script>/);
assert(match, 'Portable embedded data must be present');
const packed = JSON.parse(match[1]);
const data = packed.encoding === 'gzip-base64' ? JSON.parse(require('node:zlib').gunzipSync(Buffer.from(packed.payload, 'base64')).toString('utf8')) : packed;
assert.equal(data.samples.length, 12);
assert.equal(new Set(data.samples).size, 12);
assert.equal(data.samples.filter(s => s.startsWith('wt_')).length, 6);
assert.equal(data.populations.length, 35);
assert.equal(data.annotation.gene.strand, '-');
assert.equal(data.annotation.transcripts.length, 6);
const coding = data.annotation.transcripts.filter(t => t.CDS.length);
assert.equal(coding.length, 3);
const beta = coding.filter(t => t.protein_family === 'beta');
assert.equal(beta.length, 2);
assert.equal(new Set(beta.map(t => t.protein_sha256)).size, 1);
assert(beta.every(t => t.protein_length === 906));
assert.equal(coding.find(t => t.protein_family === 'alpha').protein_length, 1199);
assert.equal(data.coverage.length, 420);
assert.equal(data.coverage.filter(r => r.status === 'NOT_SAMPLED').length, 9);
for (const row of data.rows) {
  assert(!('barcode' in row), 'Public page contains aggregated records only');
  assert(Number.isInteger(row.count) && row.count >= 0);
  const expected = row.count / row.total_gene_UMIs * 1e6;
  assert(Math.abs(expected - row.per_million_library_UMIs) < 1e-8);
}
for (const row of data.family_summary) {
  const selected = data.rows.filter(r => r.population === row.population && r.kind === 'family_splice_supported' && r.feature === row.protein_family);
  assert.equal(selected.reduce((n, r) => n + r.count, 0), Number(row.strict_UMIs));
}
assert.equal(data.events.length, 3);
const svgs = fs.readdirSync(path.join(root, 'figures')).filter(n => n.endsWith('.svg'));
assert.equal(svgs.length, 28);
for (const svg of svgs) for (const ext of ['.svg', '.png', '.pdf']) assert(fs.existsSync(path.join(root, 'figures', svg.replace('.svg', ext))));
const validation = JSON.parse(fs.readFileSync(path.join(root, 'VALIDATION.json'), 'utf8'));
assert.equal(validation.status, 'PASS');
assert.equal(validation.full_library_denominators_verified, true);
assert(fs.readFileSync(path.join(__dirname, '..', 'site/index.html'), 'utf8').includes('grm1_isoforms/index.html'));
console.log('Grm1 atlas: twelve mice, six competing models, full-library normalization, family counts and all figure formats PASS.');
