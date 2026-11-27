/* Drive the real frontend against the live server in jsdom.
 *
 * WebGL does not exist here, so $3Dmol is replaced with a recording stub; the
 * page is written to cope with a missing viewer, and this checks that path as
 * well as the parts that matter: form wiring, submission, polling, and the
 * results table.
 */
import { readFileSync } from 'node:fs';
import { JSDOM, VirtualConsole } from 'jsdom';

const BASE = process.argv[2] || 'http://127.0.0.1:8079';
const root = process.argv[3] || '/home/claude/pydock_kiattawee/pydock_web/static';
const pass = [], fail = [];
const check = (ok, label, extra = '') =>
  (ok ? pass : fail).push(label + (extra ? `  (${extra})` : ''));

// --- a stub that records what the page asks the viewer to draw -----------
const calls = [];
const rec = (name) => (...a) => {
  calls.push({ name, args: a.map((x) => (typeof x === 'string' ? `str:${x.length}` : x)) });
  return { id: calls.length };
};
// Models need selectedAtoms: the page computes pocket residues from real
// coordinates rather than trusting 3Dmol's `within` selector.
let modelId = 0;
const fakeAtoms = (n, chain) => Array.from({ length: n }, (_, i) => ({
  x: i * 0.4, y: 0, z: 0, chain, resi: 100 + i, serial: i + 1,
}));
const makeModel = (n, chain) => ({
  id: ++modelId,
  selectedAtoms: () => fakeAtoms(n, chain),
  getID: () => modelId,
});
const viewerStub = {
  addModel: (...a) => { calls.push({ name: 'addModel', args: [] });
    return makeModel(8, calls.length % 2 ? 'A' : 'B'); },
  removeModel: rec('removeModel'),
  setStyle: rec('setStyle'), addStyle: rec('addStyle'),
  addSurface: rec('addSurface'), removeAllSurfaces: rec('removeAllSurfaces'),
  clear: rec('clear'), zoomTo: rec('zoomTo'), zoom: rec('zoom'),
  render: rec('render'),
};

const vc = new VirtualConsole();
const pageErrors = [];
vc.on('jsdomError', (e) => pageErrors.push(String(e.message || e)));
vc.on('error', (...a) => pageErrors.push(a.map(String).join(' ')));

const html = readFileSync(`${root}/index.html`, 'utf8')
  .replace(/<script src="https:\/\/3Dmol[^>]*><\/script>/, '')
  .replace('<script src="app.js"></script>', '');

const dom = new JSDOM(html, {
  url: BASE + '/', runScripts: 'dangerously', pretendToBeVisual: true,
  virtualConsole: vc,
});
const { window } = dom;
window.$3Dmol = {
  createViewer: () => viewerStub,
  SurfaceType: { VDW: 'vdw' },
};
// jsdom has fetch; point relative URLs at the server.
window.fetch = (url, opts) =>
  fetch(url.startsWith('http') ? url : BASE + url, opts);
window.FormData = FormData;
window.File = File;
window.Blob = Blob;

const appjs = readFileSync(`${root}/app.js`, 'utf8');
window.eval(appjs);

const $ = (id) => window.document.getElementById(id);
const sleep = (ms) => new Promise((r) => setTimeout(r, ms));
const waitFor = async (fn, timeout, label) => {
  const t0 = Date.now();
  while (Date.now() - t0 < timeout) {
    if (fn()) return true;
    await sleep(500);
  }
  check(false, label, `timed out after ${timeout / 1000}s`);
  return false;
};

// --- 1. boot ------------------------------------------------------------
await sleep(1200);
check(/pydock/.test($('env').textContent), 'health banner populated',
  $('env').textContent.trim());
check(!$('example').disabled, 'example button enabled when data is present');

// --- 2. site-mode switching --------------------------------------------
const modeBtns = [...$('siteMode').children];
modeBtns.find((b) => b.dataset.mode === 'centre').click();
const refBox = window.document.querySelector('.mode[data-mode=reference]');
const cenBox = window.document.querySelector('.mode[data-mode=centre]');
check(cenBox.classList.contains('hidden') === false, 'coordinate panel shown');
check(refBox.classList.contains('hidden'), 'reference panel hidden');
check([...refBox.querySelectorAll('input:not([type=file])')].every((i) => i.disabled),
  'hidden-mode value inputs disabled so they are not submitted');
check([...refBox.querySelectorAll('input[type=file]')].every((i) => !i.disabled),
  'reference file stays attached across modes, so RMSD survives');
check([...cenBox.querySelectorAll('input')].every((i) => !i.disabled),
  'active-mode inputs enabled');
modeBtns.find((b) => b.dataset.mode === 'reference').click();
check(!refBox.classList.contains('hidden'), 'switching back restores reference panel');

// --- 3. submit the bundled example and follow it ------------------------
$('example').click();
const started = await waitFor(
  () => !$('statusCard').classList.contains('hidden'), 15000,
  'status card appears after submitting');
if (started) {
  check($('intro').classList.contains('hidden'), 'intro panel hides on submit');
  await waitFor(() => Number($('fill').style.width.replace('%', '')) > 0,
    60000, 'progress bar advances');
  check($('log').textContent.includes('cavity'), 'live log streams engine output');

  const done = await waitFor(
    () => $('statusText').textContent === 'done', 420000, 'job reaches done');

  if (done) {
    check(!$('viewerCard').classList.contains('hidden'), 'viewer card revealed');
    await waitFor(
      () => $('poseTable').querySelectorAll('tbody tr').length > 0, 30000,
      'pose table populated');
    const rows = [...$('poseTable').querySelectorAll('tbody tr')];
    check(rows.length > 0, 'pose rows rendered', `${rows.length} rows`);
    check(rows[0].classList.contains('on'), 'first pose selected by default');

    const cells = [...rows[0].querySelectorAll('td')].map((c) => c.textContent.trim());
    check(/^#1$/.test(cells[0]), 'rank cell', cells[0]);
    check(!Number.isNaN(parseFloat(cells[1])), 'score cell numeric', cells[1]);
    check(/good|bad/.test(rows[0].innerHTML), 'RMSD coloured by the 2 A criterion');
    check($('summary').textContent.includes('binding mode'),
      'summary line written', $('summary').textContent);
    check($('dl').getAttribute('href').includes('/poses/0.sdf'),
      'download link points at the pose SDF');

    // Clicking another pose must swap the model and move the highlight.
    const before = calls.filter((c) => c.name === 'addModel').length;
    if (rows[1]) {
      rows[1].click();
      await sleep(400);
      check(rows[1].classList.contains('on') && !rows[0].classList.contains('on'),
        'clicking a row moves the selection');
      check(calls.filter((c) => c.name === 'addModel').length > before,
        'selecting a pose loads a new model into the viewer');
      check(calls.some((c) => c.name === 'removeModel'),
        'previous pose model removed');
    }
    const models = calls.filter((c) => c.name === 'addModel').length;
    check(models >= 3, 'receptor, reference and pose all loaded',
      `${models} addModel calls`);
    check(calls.some((c) => c.name === 'zoomTo'), 'view framed on the pose');
    check($('legend').textContent.includes('reference'),
      'legend lists the reference overlay');

    // Toggles must restyle without throwing.
    $('optSurface').checked = true;
    $('optSurface').dispatchEvent(new window.Event('change'));
    await sleep(200);
    check(calls.some((c) => c.name === 'addSurface'), 'surface toggle requests a surface');
    $('optPocket').checked = false;
    $('optPocket').dispatchEvent(new window.Event('change'));
    await sleep(200);
    check(true, 'pocket toggle handled');

    await sleep(600);
    check(!$('jobsCard').classList.contains('hidden'), 'recent jobs listed');
    check($('jobs').children.length > 0, 'recent job rows rendered',
      `${$('jobs').children.length} jobs`);
  }
}

// --- 4. server-side rejection surfaces in the UI ------------------------
const fd = new FormData();
fd.set('scoring', 'magic');
const r = await fetch(BASE + '/api/jobs/example', { method: 'POST', body: fd });
check(r.status === 400, 'bad scoring value rejected by the API', `HTTP ${r.status}`);

// --- report -------------------------------------------------------------
console.log('\nPASS');
for (const p of pass) console.log('  ok   ' + p);
if (fail.length) {
  console.log('\nFAIL');
  for (const f of fail) console.log('  FAIL ' + f);
}
if (pageErrors.length) {
  console.log('\nJS ERRORS IN PAGE');
  for (const e of [...new Set(pageErrors)].slice(0, 10)) console.log('  ! ' + e.slice(0, 220));
}
console.log(`\n${pass.length} passed, ${fail.length} failed, `
  + `${new Set(pageErrors).size} page errors`);
dom.window.close();
process.exit(fail.length || pageErrors.length ? 1 : 0);
