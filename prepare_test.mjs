/* Drive the preparation page against the live server, in jsdom.
 *
 * The backend is already tested through curl; what this checks is the part
 * curl cannot see: that the page wires itself up, that clicking a ligand
 * selects it, that Prepare sends the right form and renders the result, and
 * that the handoff link carries the preparation id.
 *
 *   node prepare_test.mjs [base-url] [static-dir]
 */
import { readFileSync } from 'node:fs';
import { JSDOM, VirtualConsole } from 'jsdom';

const BASE = process.argv[2] || 'http://127.0.0.1:8079';
const root = process.argv[3] || '/home/claude/pydock_kiattawee/pydock_web/static';
const pass = [], fail = [];
const check = (ok, label, extra = '') =>
  (ok ? pass : fail).push(label + (extra && !ok ? `  (${extra})` : ''));

/* ---- a viewer stub that records what the page draws ------------------- */
const calls = [];
const rec = (name) => (...a) => {
  calls.push({ name, args: a });
  return { id: calls.length };
};
const viewer = {
  addModel: rec('addModel'), setStyle: rec('setStyle'), addStyle: rec('addStyle'),
  setClickable: rec('setClickable'), zoomTo: rec('zoomTo'), zoom: rec('zoom'),
  render: rec('render'), clear: rec('clear'), removeAllShapes: rec('removeAllShapes'),
  addLine: rec('addLine'), addSphere: rec('addSphere'), getModel: () => null,
};

const vc = new VirtualConsole();
const noise = [];
vc.on('jsdomError', (e) => noise.push('jsdomError: ' + e.message));
vc.on('error', (...a) => noise.push('console.error: ' + a.join(' ')));

const dom = new JSDOM(readFileSync(`${root}/prepare.html`, 'utf8'), {
  url: `${BASE}/prepare.html`,
  runScripts: 'outside-only',
  pretendToBeVisual: true,
  virtualConsole: vc,
});
const { window } = dom;
window.$3Dmol = { createViewer: () => viewer };
// The browser resolves a relative URL against the page; Node's fetch does
// not, so the stub does it here.
window.fetch = (url, opts) =>
  fetch(typeof url === 'string' && url.startsWith('/') ? BASE + url : url, opts);
window.FormData = FormData;
window.File = File;
window.Blob = Blob;

window.eval(readFileSync(`${root}/prepare.js`, 'utf8'));
const doc = window.document;
const $ = (id) => doc.getElementById(id);

const sleep = (ms) => new Promise((r) => setTimeout(r, ms));
async function settle(test, ms = 60000) {
  const until = Date.now() + ms;
  while (Date.now() < until) {
    if (test()) return true;
    await sleep(120);
  }
  return false;
}

/* ---- 1. the page loads and is wired ---------------------------------- */
check(!!$('file'), 'the file input exists');
check($('pickSection').classList.contains('hidden'),
      'the ligand list is hidden before a file is loaded');
check(!!$('build') && !!$('toDock'), 'the Prepare and Send buttons exist');

/* ---- 2. loading a structure ------------------------------------------ */
const pdb = readFileSync('/tmp/complex.pdb', 'utf8');
// jsdom gives no real file picker, so the change event is fired by hand with
// a File the page can read -- the same object a picker would hand it.
const file = new window.File([pdb], 'complex.pdb', { type: 'chemical/x-pdb' });
Object.defineProperty($('file'), 'files', { value: [file], writable: true });
$('file').dispatchEvent(new window.Event('change'));

const loaded = await settle(() => !$('pickSection').classList.contains('hidden'));
check(loaded, 'the ligand list appears once the structure is read',
      $('loadStatus').textContent);

if (loaded) {
  const text = $('loadStatus').textContent;
  check(/2,498 atoms/.test(text), 'the atom count is reported', text);
  check(/2 chains/.test(text), 'the chain count is reported', text);

  const buttons = [...doc.querySelectorAll('#picks button')];
  check(buttons.length === 1, 'exactly one ligand candidate is offered',
        `${buttons.length}: ` + buttons.map((b) => b.textContent).join(' | '));
  check(/LIG A:901/.test(buttons[0]?.textContent || ''),
        'the ligand is the one in the file', buttons[0]?.textContent);
  check(/43 atoms/.test(buttons[0]?.textContent || ''),
        'its heavy-atom count is shown');
  check(buttons[0]?.getAttribute('aria-pressed') === 'true',
        'the only candidate is selected for them');
  check(!$('smilesSection').classList.contains('hidden'),
        'the SMILES box appears with a ligand selected');
  check(/43 heavy atoms/.test($('smilesHint').textContent),
        'the SMILES box states the expected atom count');

  // The viewer must be told to draw the hetero groups solid, or they cannot
  // be found on screen.
  const sphered = calls.filter((c) => c.name === 'setStyle'
    && JSON.stringify(c.args[1] || {}).includes('sphere'));
  check(sphered.length > 0, 'hetero groups are drawn as spheres');
  check(calls.some((c) => c.name === 'setClickable'),
        'the viewer is made clickable');
  check(/metal/.test($('legend').innerHTML), 'the legend names the metal');

  check($('optWater').checked && /6 present/.test($('wCount').textContent),
        'waters are counted and removal is on by default', $('wCount').textContent);
  check($('optH').disabled, 'the hydrogen option is disabled when there are none');
  check(/ZN/.test($('mCount').textContent), 'the metal is named', $('mCount').textContent);

  /* ---- 3. clicking a candidate toggles it ---------------------------- */
  buttons[0].dispatchEvent(new window.Event('click'));
  await sleep(50);
  check(doc.querySelector('#picks button').getAttribute('aria-pressed') === 'false',
        'clicking the selected ligand deselects it');
  check($('smilesSection').classList.contains('hidden'),
        'the SMILES box hides with nothing selected');
  doc.querySelector('#picks button').dispatchEvent(new window.Event('click'));
  await sleep(50);
  check(doc.querySelector('#picks button').getAttribute('aria-pressed') === 'true',
        'clicking again selects it');

  /* ---- 4. a wrong SMILES is reported, not swallowed ------------------ */
  $('smiles').value = 'CCO';
  $('build').dispatchEvent(new window.Event('click'));
  const refused = await settle(() => /not the same compound/.test($('buildStatus').textContent));
  check(refused, 'a mismatched SMILES is reported on the page',
        $('buildStatus').textContent);
  check($('resultSection').classList.contains('hidden'),
        'no results are shown after a refusal');

  /* ---- 5. preparing for real ---------------------------------------- */
  const smiles = readFileSync('/tmp/smiles.txt', 'utf8').trim();
  $('smiles').value = smiles;
  $('build').dispatchEvent(new window.Event('click'));
  const built = await settle(() => !$('resultSection').classList.contains('hidden'));
  check(built, 'the results section appears', $('buildStatus').textContent);

  if (built) {
    const summary = $('summary').textContent;
    check(/2,444 atoms/.test(summary), 'the receptor atom count is shown', summary);
    check(/43 ligands?/.test(summary) && /6 waters/.test(summary),
          'what was removed is itemised', summary);
    check(/from your SMILES/.test(summary),
          'the bond orders are credited to the SMILES', summary);

    const links = [...doc.querySelectorAll('#files a')].map((a) => a.getAttribute('href'));
    check(links.length === 3, 'three files are offered', links.join(' '));
    check(links.some((h) => h.endsWith('receptor.pdb')), 'the receptor is downloadable');
    check(links.some((h) => h.endsWith('ligand.mol2')), 'the MOL2 is downloadable');

    // Every link must actually serve something.
    for (const href of links) {
      const r = await fetch(BASE + href);
      check(r.ok, `${href.split('/').pop()} downloads`, `http ${r.status}`);
    }

    /* ---- saving onto this machine ------------------------------------ */
    check(!!$('saveBtn') && !!$('saveFolder'), 'the save controls exist');
    check($('dlZip').getAttribute('href').endsWith('/bundle.zip'),
          'the ZIP link is set', $('dlZip').getAttribute('href'));
    const zip = await fetch(BASE + $('dlZip').getAttribute('href'));
    check(zip.ok && (await zip.arrayBuffer()).byteLength > 1000,
          'the ZIP downloads');

    const dest = '/tmp/jsdom_save';
    $('saveFolder').value = dest;
    $('saveBtn').dispatchEvent(new window.Event('click'));
    const saved = await settle(() => /saved/.test($('saveStatus').textContent));
    check(saved, 'the page reports the files as saved', $('saveStatus').textContent);
    check(/jsdom_save/.test($('saveStatus').innerHTML),
          'and shows where they went');
    const { existsSync, readdirSync } = await import('node:fs');
    const onDisk = existsSync(dest) ? readdirSync(dest) : [];
    check(onDisk.length === 3, 'three files really are on disk', onDisk.join(', '));
    check(onDisk.some((f) => f.endsWith('.mol2')) &&
          onDisk.some((f) => f.endsWith('_receptor.pdb')),
          'the ligand MOL2 and the receptor are among them', onDisk.join(', '));

    // A second save must not quietly replace the first.
    $('saveBtn').dispatchEvent(new window.Event('click'));
    await settle(() => /numbered/.test($('saveStatus').textContent));
    check(/numbered/.test($('saveStatus').textContent),
          'saving twice numbers the copies instead of overwriting',
          $('saveStatus').textContent);
    check(readdirSync(dest).length === 6, 'so six files are there now',
          String(readdirSync(dest).length));

    const good = doc.querySelector('#resultSection .warn.good');
    check(!!good && /matches the SMILES/.test(good.textContent),
          'the page confirms the ligand matches the SMILES');

    /* ---- 6. the handoff ---------------------------------------------- */
    const went = $('toDock').getAttribute('href');
    check(typeof went === 'string' && went.startsWith('index.html?'),
          'Send to docking goes to the docking page', String(went));
    check(/prep=[0-9a-f]{32}/.test(went || ''),
          'the handoff carries the preparation id', String(went));
  }
}

check(noise.length === 0, 'the page runs without script errors', noise.join(' | '));

for (const line of pass) console.log('  ok   ' + line);
for (const line of fail) console.log('  FAIL ' + line);
console.log(`\n${pass.length} passed, ${fail.length} failed\n`);
process.exit(fail.length ? 1 : 0);
