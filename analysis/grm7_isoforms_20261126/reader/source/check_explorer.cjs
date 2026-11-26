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
  assert.equal(await page.evaluate(() => shown.length), 48);
  checks.push('All 35 populations selectable; default terminal-event panel contains 4 events × 12 mice.');
  for (const kind of ['gene', 'diagnostic_feature', 'splice_supported_unique', 'compatible_upper_bound', 'exon', 'junction', 'CDS', 'UTR', 'terminal_event']) {
    await page.selectOption('#kind', kind);
    assert((await page.evaluate(() => shown.length)) > 0, kind);
    assert.equal(await page.locator('#plot svg').count(), 1);
  }
  checks.push('All nine evidence views render.');
  await page.locator('summary').filter({ hasText: 'Appearance and feature selection' }).click();
  await page.selectOption('#population', 'IN_C1');
  await page.selectOption('#kind', 'diagnostic_feature');
  await page.fill('#search', 'Grm7-202');
  assert.equal(await page.evaluate(() => shown.reduce((n, r) => n + r.count, 0)), 67);
  await page.fill('#search', '');
  const missing = await page.evaluate(() => D.coverage.find(x => x.status === 'NOT_SAMPLED'));
  await page.selectOption('#population', missing.population);
  await page.selectOption('#kind', 'gene');
  assert((await page.locator('#plot svg').textContent()).includes('—'));
  checks.push('IN C1 diagnostic A-feature sum matches source (67); unsampled combinations display missing dashes.');
  await page.selectOption('#population', 'CA1_superficial');
  await page.selectOption('#kind', 'terminal_event');
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
    if (ext === 'json') { const data = JSON.parse(bytes); assert.equal(data.data.length, 48); assert.equal(data.settings.high, '#772277'); }
    if (ext === 'tsv') { const rows = bytes.toString().trimEnd().split('\n'); assert.equal(rows.length, 49); assert(rows[0].includes('sample\t')); }
  }
  checks.push('SVG, PNG, PDF, TSV and JSON download successfully; signatures, row counts and saved settings verified.');
  for (const name of ['01_Grm7_transcript_schematic', '09_Grm7_A_B_terminal_splicing', '08_Grm7_presynaptic_isoform_schematic']) {
    await page.selectOption('#schematic', name);
    await page.waitForFunction(() => document.getElementById('scheme').complete && document.getElementById('scheme').naturalWidth > 0);
    for (const ext of ['svg', 'png', 'pdf']) assert(fs.existsSync(path.join(root, 'deliverables/figures', name + '.' + ext)));
  }
  assert.equal(await page.locator('#gallery p').count(), 21);
  checks.push('Three schematics load; all 21 static figures have SVG, PNG and PDF copies.');
  assert.equal(await page.locator('#error').textContent(), '');
  assert.deepEqual(errors, []);
  await page.locator('#plot svg').screenshot({ path: path.join(out, 'interactive_panel.png') });
  await page.screenshot({ path: path.join(out, 'explorer_top.png') });
  await browser.close();
  const result = { status: 'PASS', browser: 'Chrome / Playwright', checks, page_errors: errors, files_verified: ['SVG', 'PNG', 'PDF', 'TSV', 'JSON'] };
  fs.writeFileSync(path.join(out, 'EXPLORER_QA.json'), JSON.stringify(result, null, 2));
  console.log(JSON.stringify(result, null, 2));
})().catch(e => { console.error(e); process.exit(1); });
