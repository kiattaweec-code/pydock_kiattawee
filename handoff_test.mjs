/* The docking page must adopt a structure prepared on the other page.
 *
 * This is the seam where the two pages meet, and a mistake here is quiet:
 * the page would simply dock an unprepared receptor, or send an empty file
 * input over the server's own copy.
 *
 *   node handoff_test.mjs <prep-id> [base-url] [static-dir]
 */
import { readFileSync } from 'node:fs';
import { JSDOM, VirtualConsole } from 'jsdom';

const PREP = process.argv[2];
const BASE = process.argv[3] || 'http://127.0.0.1:8079';
const root = process.argv[4] || '/home/claude/pydock_kiattawee/pydock_web/static';
if (!PREP) { console.error('give a prep id'); process.exit(2); }

const pass = [], fail = [];
const check = (ok, label, extra = '') =>
  (ok ? pass : fail).push(label + (extra && !ok ? `  (${extra})` : ''));

const calls = [];
const rec = (name) => (...a) => { calls.push({ name, args: a }); return { id: calls.length }; };
let modelId = 0;
const atoms = (n, chain) => Array.from({ length: n }, (_, i) => ({
  x: i * 0.4, y: 0, z: 0, chain, resi: 100 + i, serial: i + 1, elem: 'C',
}));
const model = (n, chain) => ({
  id: ++modelId, selectedAtoms: () => atoms(n, chain), setStyle: rec('m.setStyle'),
});
const viewer = {
  addModel: (...a) => { calls.push({ name: 'addModel', args: a }); return model(40, 'A'); },
  setStyle: rec('setStyle'), addStyle: rec('addStyle'), setClickable: rec('setClickable'),
  zoomTo: rec('zoomTo'), zoom: rec('zoom'), render: rec('render'), clear: rec('clear'),
  removeAllShapes: rec('removeAllShapes'), addLine: rec('addLine'),
  addSphere: rec('addSphere'), addSurface: rec('addSurface'),
  removeAllSurfaces: rec('removeAllSurfaces'), getModel: () => model(40, 'A'),
  removeModel: rec('removeModel'), center: rec('center'), rotate: rec('rotate'),
  resize: rec('resize'), setViewStyle: rec('setViewStyle'),
};

const vc = new VirtualConsole();
const noise = [];
vc.on('jsdomError', (e) => noise.push('jsdomError: ' + e.message));
vc.on('error', (...a) => noise.push('console.error: ' + a.join(' ')));

const dom = new JSDOM(readFileSync(`${root}/index.html`, 'utf8'), {
  url: `${BASE}/index.html?prep=${PREP}&ligand=LIG%20A%3A901`,
  runScripts: 'outside-only', pretendToBeVisual: true, virtualConsole: vc,
});
const { window } = dom;
window.$3Dmol = { createViewer: () => viewer, elementColors: { defaultColors: {} } };
window.fetch = (url, opts) =>
  fetch(typeof url === 'string' && url.startsWith('/') ? BASE + url : url, opts);
// Node's FormData cannot be built from a form element, which is exactly
// what the page does, so the fields are collected the way a browser would.
window.FormData = class extends FormData {
  constructor(form) {
    super();
    if (!form || !form.elements) return;
    for (const el of form.elements) {
      if (!el.name || el.disabled) continue;
      if (el.type === 'file') {
        for (const f of (el.files || [])) this.append(el.name, f);
      } else if (el.type === 'checkbox' || el.type === 'radio') {
        if (el.checked) this.append(el.name, el.value);
      } else {
        this.append(el.name, el.value);
      }
    }
  }
};
window.File = File;
window.Blob = Blob;
window.DataTransfer = class {
  constructor() { this.items = { _f: [], add: (f) => this.items._f.push(f) }; }
  get files() { return this.items._f; }
};

// jsdom implements neither DataTransfer nor an assignable FileList, so the
// two file inputs get plain writable properties before the page runs. The
// page's own technique (new DataTransfer -> input.files) is the standard one
// and works in browsers; only the stand-in needs this.
for (const name of ['receptor', 'reference', 'ligands']) {
  const el = dom.window.document.querySelector(`[name=${name}]`);
  if (el) Object.defineProperty(el, 'files', { value: [], writable: true, configurable: true });
}

window.eval(readFileSync(`${root}/app.js`, 'utf8'));
const doc = window.document;
const $ = (id) => doc.getElementById(id);
const form = $('form');

const sleep = (ms) => new Promise((r) => setTimeout(r, ms));
async function settle(test, ms = 90000) {
  const until = Date.now() + ms;
  while (Date.now() < until) { if (test()) return true; await sleep(150); }
  return false;
}

/* ---- the banner and the adopted files -------------------------------- */
const adopted = await settle(() => !$('preparedCard').classList.contains('hidden'));
check(adopted, 'the prepared-structure banner appears', $('preparedNote')?.textContent);

if (adopted) {
  check(/ligand[\s\S]*removed/i.test($('preparedNote').textContent),
        'the banner says the ligand was removed', $('preparedNote').textContent);
  check(/LIG A:901/.test($('preparedNote').textContent),
        'the banner names the extracted ligand', $('preparedNote').textContent);

  const rec0 = form.elements.receptor.files[0];
  check(!!rec0, 'the receptor input has been filled in');
  check(rec0 && /prepared/.test(rec0.name), 'with the prepared receptor', rec0 && rec0.name);
  check(form.elements.receptor.required === false,
        'the receptor is no longer a required upload');

  const ref0 = form.elements.reference.files[0];
  check(!!ref0, 'the extracted ligand is loaded as the reference');
  check(form.elements.ligands.required === false,
        'ligands are optional, since the prepared one can be redocked');
  check(!$('prepRedock').parentElement.classList.contains('hidden'),
        'the redock option is offered');

  /* ---- the site opens from the reference ----------------------------- */
  const sited = await settle(() => !$('siteCard').classList.contains('hidden'));
  check(sited, 'the binding-site card opens');
  // The card opens first and fills in once the server has scanned for
  // pockets, which takes a few seconds on a real protein.
  const read = await settle(() => !/reading structure/.test($('siteMeta').textContent));
  const meta = $('siteMeta').textContent;
  check(read && /atom/.test(meta), 'the receptor was read by the server', meta);
  check(calls.some((c) => c.name === 'addModel'), 'the structure is drawn',
        calls.map((c) => c.name).slice(0, 8).join(','));

  /* ---- submitting sends prep_id, not an empty receptor --------------- */
  let sent = null;
  const realFetch = window.fetch;
  window.fetch = async (url, opts) => {
    if (typeof url === 'string' && url.endsWith('/api/jobs') && opts && opts.body) {
      sent = opts.body;
      return new Response(JSON.stringify({ id: 'testjob0000', status: 'queued' }),
                          { status: 202, headers: { 'content-type': 'application/json' } });
    }
    return realFetch(url, opts);
  };
  $('prepRedock').checked = true;
  form.dispatchEvent(new window.Event('submit'));
  const posted = await settle(() => sent !== null, 15000);
  check(posted, 'the form submits');

  if (posted) {
    check(sent.get('prep_id') === PREP,
          'the preparation id is sent', String(sent.get('prep_id')));
    check(sent.get('dock_prepared_ligand') === 'true',
          'the redock choice is sent', String(sent.get('dock_prepared_ligand')));
    check(!sent.has('receptor'),
          'no receptor file is sent, so the server keeps its own copy');
  }

  /* ---- choosing their own receptor drops the prepared one ------------ */
  window.fetch = realFetch;
  const own = new window.File([readFileSync('/tmp/clean.pdb', 'utf8')],
                              'mine.pdb', { type: 'chemical/x-pdb' });
  form.elements.receptor.files = [own];
  form.elements.receptor.dispatchEvent(new window.Event('change'));
  await sleep(300);
  check($('preparedCard').classList.contains('hidden'),
        'uploading their own receptor hides the prepared banner');
  check(form.elements.ligands.required === true,
        'and makes a ligand file required again');
}

check(noise.length === 0, 'the page runs without script errors', noise.join(' | '));

for (const line of pass) console.log('  ok   ' + line);
for (const line of fail) console.log('  FAIL ' + line);
console.log(`\n${pass.length} passed, ${fail.length} failed\n`);
process.exit(fail.length ? 1 : 0);
