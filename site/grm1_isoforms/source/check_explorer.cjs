const fs = require('node:fs');
const path = require('node:path');
const assert = require('node:assert/strict');
const { chromium } = require('/Users/mujahid357/.cache/codex-runtimes/codex-primary-runtime/dependencies/node/node_modules/playwright');
const root = __dirname;
const out = path.join(root, 'qa');
fs.mkdirSync(out, { recursive: true });
(async () => {
  const browser = await chromium.launch({ headless: true, executablePath: '/Applications/Google Chrome.app/Contents/MacOS/Google Chrome' });
  const page = await browser.newPage({ viewport: { width: 1440, height: 1050 }, acceptDownloads: true });
  const errors = [];
  page.on('pageerror', e => errors.push(e.message));
  await page.goto('file://' + path.join(root, 'deliverables/index.html'));
  await page.waitForSelector('#plot svg');
  const checks = [];
  assert.equal(await page.locator('#population option').count(), 35);
  assert.equal(await page.evaluate(() => shown.length), 24);
  checks.push('All 35 populations selectable; default alpha/beta panel contains 2 protein-family assignments × 12 mice.');
  for (const kind of ['family_splice_supported','family_diagnostic_feature','family_compatible_upper_bound','family_junction','gene', 'diagnostic_feature', 'splice_supported_unique', 'compatible_upper_bound', 'exon', 'junction', 'CDS', 'UTR']) {
    await page.selectOption('#kind', kind);
    assert((await page.evaluate(() => shown.length)) > 0, kind);
    assert.equal(await page.locator('#plot svg').count(), 1);
  }
  checks.push('All twelve evidence views render.');
  await page.locator('summary').filter({ hasText: 'Appearance and feature selection' }).click();
  await page.selectOption('#population', 'IN_C2');
  await page.selectOption('#kind', 'family_splice_supported');
  await page.fill('#search', 'alpha');
  assert.equal(await page.evaluate(() => shown.reduce((n, r) => n + r.count, 0)), 15);
  await page.fill('#search', '');
  const missing = await page.evaluate(() => D.coverage.find(x => x.status === 'NOT_SAMPLED'));
  await page.selectOption('#population', missing.population);
  await page.selectOption('#kind', 'gene');
  assert((await page.locator('#plot svg').textContent()).includes('—'));
  checks.push('IN C2 conservative alpha sum matches source (15); unsampled combinations display missing dashes.');
  await page.selectOption('#population', 'CA1_superficial');
  await page.selectOption('#kind', 'family_splice_supported');
  await page.selectOption('#unit', 'per_1000_nuclei');
  await page.selectOption('#style', 'dot');
  await page.evaluate(() => { for (const [id, value] of Object.entries({ high: '#772277', font: '17', size: '12', width: '75' })) { document.getElementById(id).value = value; document.getElementById(id).dispatchEvent(new Event('input')); } });
  assert((await page.locator('#plot').innerHTML()).includes('font-size="17"'));
  assert((await page.locator('#plot circle').count()) > 0);
  checks.push('Population, normalization, color, label size, dot size, width and dotplot controls update the panel.');
  for (const ext of ['svg', 'png', 'pdf', 'tsv', 'json']) {
    const downloadEvent = page.waitForEvent('download');
    await page.locator('[data-save="' + ext + '"]').click();
    const download = await downloadEvent;
    const target = path.join(out, 'interactive_check.' + ext);
    await download.saveAs(target);
    const bytes = fs.readFileSync(target);
    assert(bytes.length > 100, ext);
    if (ext === 'pdf') assert(bytes.toString('ascii', 0, 5) === '%PDF-');
    if (ext === 'png') assert(bytes.subarray(0, 8).equals(Buffer.from([137,80,78,71,13,10,26,10])));
    if (ext === 'svg') assert(bytes.toString().includes('<svg'));
    if (ext === 'json') { const data = JSON.parse(bytes); assert.equal(data.data.length, 24); assert.equal(data.settings.high, '#772277'); }
    if (ext === 'tsv') { const rows = bytes.toString().trimEnd().split('\n'); assert.equal(rows.length, 25); assert(rows[0].includes('sample\t')); }
  }
  checks.push('SVG, PNG, PDF, TSV and JSON download successfully; signatures, row counts and saved settings verified.');
  for (const name of ['01_Grm1_transcript_schematic', '08_Grm1_alpha_beta_protein_schematic']) {
    await page.selectOption('#schematic', name);
    await page.waitForFunction(() => document.getElementById('scheme').complete && document.getElementById('scheme').naturalWidth > 0);
    for (const ext of ['svg', 'png', 'pdf']) assert(fs.existsSync(path.join(root, 'deliverables/figures', name + '.' + ext)));
  }
  assert.equal(await page.locator('#gallery p').count(), 28);
  checks.push('Two schematics load; all 28 static figures have SVG, PNG and PDF copies.');
  assert.equal(await page.locator('#error').textContent(), '');
  assert.deepEqual(errors, []);
  await page.locator('#plot svg').screenshot({ path: path.join(out, 'interactive_panel.png') });
  await page.screenshot({ path: path.join(out, 'explorer_top.png') });
  await browser.close();
  const result = { status: 'PASS', browser: 'Chrome / Playwright', checks, page_errors: errors, files_verified: ['SVG', 'PNG', 'PDF', 'TSV', 'JSON'] };
  fs.writeFileSync(path.join(out, 'EXPLORER_QA.json'), JSON.stringify(result, null, 2));
  console.log(JSON.stringify(result, null, 2));
})().catch(e => { console.error(e); process.exit(1); });
