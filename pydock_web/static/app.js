/* pydock web client: submit a job, follow it, then look at the poses. */
'use strict';

const $ = (id) => document.getElementById(id);
const api = (path, opts) => fetch(path, opts);

const state = {
  jobId: null,
  job: null,
  ligIndex: 0,
  poseRank: 1,
  poll: null,
  viewer: null,
  models: { receptor: null, pose: null, reference: null },
  poses: [],          // SDF text per pose, for the current ligand
  referenceSdf: null,
  framed: false,
  surfaceOn: false,
};

Object.assign(state, {
  prep: null,                 // the /api/prepare response for the loaded receptor
  siteViewer: null,
  siteModels: { receptor: null, reference: null },
  boxShape: [],
  preparePending: false,
  box: { centre: [0, 0, 0], size: [22, 22, 22] },
  pickMode: false,
  preparing: false,
});

/* -------------------------------------------------- the docking box -- */
const AXES = ['x', 'y', 'z'];
const SIZE_LABELS = ['w', 'h', 'd'];
const fmt = (v) => (Math.round(v * 10) / 10).toFixed(1);

function buildSliders() {
  for (const [host, labels, kind] of [
    [$('centreSliders'), AXES, 'centre'],
    [$('sizeSliders'), SIZE_LABELS, 'size'],
  ]) {
    host.innerHTML = '';
    labels.forEach((label, i) => {
      const cell = document.createElement('div');
      cell.className = 'sliderCell';
      const r = document.createElement('input');
      r.type = 'range';
      r.step = '0.2';
      if (kind === 'size') { r.min = '6'; r.max = '60'; } else { r.min = '-100'; r.max = '100'; }
      r.dataset.kind = kind;
      r.dataset.axis = String(i);
      const out = document.createElement('output');
      cell.innerHTML = `<span>${label}</span>`;
      cell.appendChild(r);
      cell.appendChild(out);
      host.appendChild(cell);
      r.addEventListener('input', () => {
        const next = { centre: [...state.box.centre], size: [...state.box.size] };
        next[kind][i] = parseFloat(r.value);
        setBox(next.centre, next.size, 'slider');
      });
    });
  }
}
buildSliders();

function sliderOf(kind, i) {
  const host = kind === 'centre' ? $('centreSliders') : $('sizeSliders');
  return host.querySelectorAll('.sliderCell')[i];
}

/** Single source of truth for the site box: updates sliders, the form
 *  fields the job is submitted with, the 3D box, and the readout. */
function setBox(centre, size, source) {
  state.box.centre = centre.map(Number);
  state.box.size = size.map((v) => Math.min(80, Math.max(6, Number(v))));

  const f = $('form').elements;
  AXES.forEach((a, i) => {
    f['centre_' + a].value = fmt(state.box.centre[i]);
    f['size_' + a].value = fmt(state.box.size[i]);
    for (const [kind, arr] of [['centre', state.box.centre], ['size', state.box.size]]) {
      const cell = sliderOf(kind, i);
      if (!cell) continue;
      const range = cell.querySelector('input');
      const out = cell.querySelector('output');
      if (kind === 'centre' && state.prep) {
        // Keep the slider range tied to the structure that is loaded.
        range.min = fmt(state.prep.bounds.lower[i] - 5);
        range.max = fmt(state.prep.bounds.upper[i] + 5);
      }
      if (source !== 'slider' || Number(range.value) !== arr[i]) range.value = arr[i];
      out.textContent = fmt(arr[i]);
    }
  });

  // Touching the box means the box is the site, so move off reference mode.
  if (source === 'slider' || source === 'pick' || source === 'pocket') {
    selectSiteMode('centre');
  }
  markActivePocket();
  drawBox();
  readout();
}

function readout() {
  const [w, h, d] = state.box.size;
  const vol = w * h * d;
  const big = vol > 27000;
  $('boxReadout').innerHTML =
    `box <b>${fmt(w)} × ${fmt(h)} × ${fmt(d)} Å</b> · volume <b>${Math.round(vol).toLocaleString()} Å³</b>`
    + ` · centre <b>${state.box.centre.map(fmt).join(', ')}</b>`
    + (big ? ' · <span class="warn">a box this large will dock slowly; '
           + 'narrow it to the pocket if you can</span>' : '');
}

function markActivePocket() {
  const near = (a, b) => Math.hypot(a[0] - b[0], a[1] - b[1], a[2] - b[2]) < 0.6;
  for (const btn of $('pockets').querySelectorAll('button[data-centre]')) {
    btn.classList.toggle('on', near(JSON.parse(btn.dataset.centre), state.box.centre));
  }
}

function selectSiteMode(mode) {
  const btn = [...$('siteMode').children].find((b) => b.dataset.mode === mode);
  if (btn && !btn.classList.contains('on')) btn.click();
}

/* ------------------------------------------- receptor preview fetch -- */
function fileOf(name) {
  const input = $('form').elements[name];
  return input && input.files && input.files[0] ? input.files[0] : null;
}

async function prepareReceptor() {
  const rec = fileOf('receptor');
  if (!rec) return;
  if (state.preparing) {
    // Picking a receptor and a reference in quick succession would otherwise
    // drop the second request and leave the reference unloaded.
    state.preparePending = true;
    return;
  }
  state.preparing = true;
  $('siteCard').classList.remove('hidden');
  $('intro').classList.add('hidden');
  $('siteMeta').textContent = 'reading structure…';
  $('pockets').innerHTML = '<span class="scanning">detecting pockets…</span>';
  showErr('');

  const fd = new FormData();
  fd.set('receptor', rec);
  const ref = fileOf('reference');
  if (ref) fd.set('reference', ref);
  fd.set('find_pockets', 'true');

  try {
    const r = await api('/api/prepare', { method: 'POST', body: fd });
    const body = await r.json();
    if (!r.ok) throw new Error(body.detail || `could not read the receptor (${r.status})`);
    state.prep = body;
    $('siteMeta').textContent =
      `${body.n_residues} residues · ${body.n_atoms} atoms`;
    renderPockets(body);
    await drawSiteStructure(body);

    // Opening box: the reference if there is one, else the best pocket.
    if (body.reference) {
      setBox(body.reference.centre, body.reference.size, 'init');
      $('siteFitRef').disabled = false;
    } else if (body.pockets.length) {
      setBox(body.pockets[0].centre, body.pockets[0].size, 'pocket');
    } else {
      setBox(body.bounds.centre, [22, 22, 22], 'init');
    }
    frameBox();
  } catch (err) {
    $('siteMeta').textContent = '';
    $('pockets').innerHTML = '';
    showErr(err.message || String(err));
    $('siteCard').classList.add('hidden');
  } finally {
    state.preparing = false;
    if (state.preparePending) {
      state.preparePending = false;
      prepareReceptor();
    }
  }
}

function renderPockets(body) {
  const host = $('pockets');
  host.innerHTML = '';
  if (!body.pockets.length) {
    host.innerHTML = '<span class="scanning">no enclosed pocket found — '
      + 'place the box by hand, or load a reference ligand</span>';
    return;
  }
  const lab = document.createElement('span');
  lab.className = 'plab';
  lab.textContent = 'detected pockets';
  host.appendChild(lab);
  for (const p of body.pockets) {
    const b = document.createElement('button');
    b.type = 'button';
    b.dataset.centre = JSON.stringify(p.centre);
    b.innerHTML = `#${p.rank}<span class="pvol">${p.volume.toLocaleString()} Å³</span>`;
    b.title = `centre ${p.centre.join(', ')}`;
    b.addEventListener('click', () => { setBox(p.centre, p.size, 'pocket'); frameBox(); });
    host.appendChild(b);
  }
}

async function drawSiteStructure(body) {
  if (!window.$3Dmol) {
    $('siteNote').textContent =
      '3Dmol.js could not be loaded, so the box cannot be shown here. '
      + 'The numeric centre and size fields still work.';
    $('siteNote').classList.remove('hidden');
    return;
  }
  if (!state.siteViewer) {
    state.siteViewer = $3Dmol.createViewer($('siteViewer'),
      { backgroundColor: '#090d13', antialias: true });
  }
  const v = state.siteViewer;
  v.clear();
  state.boxShape = [];
  state.siteModels.receptor = v.addModel(body.pdb, 'pdb');
  v.setStyle({ model: state.siteModels.receptor },
    { cartoon: { color: '#6d8099', opacity: 0.85, thickness: 0.3 } });
  if (body.reference) {
    state.siteModels.reference = v.addModel(body.reference.sdf, 'sdf');
    v.setStyle({ model: state.siteModels.reference },
      { stick: { radius: 0.11, colorscheme: 'greenCarbon' } });
  } else {
    state.siteModels.reference = null;
  }

  // Clicking an atom recentres the box on it.
  v.setClickable({ model: state.siteModels.receptor }, true, (atom) => {
    if (!state.pickMode || !atom) return;
    setBox([atom.x, atom.y, atom.z], state.box.size, 'pick');
  });

  $('siteLegend').innerHTML =
    '<span><i style="background:#45cfc4"></i>docking box</span>'
    + (body.reference ? '<span><i style="background:#6ec98a"></i>reference</span>' : '')
    + '<span><i style="background:#6d8099"></i>receptor</span>';
  v.zoomTo();
  v.render();
}

/* 3Dmol's wireframe box draws the triangulation, so every face gets a
 * diagonal through it.  Twelve explicit edges give a clean cage instead. */
const BOX_EDGES = [
  [0, 1], [1, 3], [3, 2], [2, 0],
  [4, 5], [5, 7], [7, 6], [6, 4],
  [0, 4], [1, 5], [2, 6], [3, 7],
];

function drawBox() {
  const v = state.siteViewer;
  if (!v) return;
  for (const s of state.boxShape || []) {
    try { v.removeShape(s); } catch { /* already gone */ }
  }
  const [cx, cy, cz] = state.box.centre;
  const [w, h, d] = state.box.size;
  const corner = (i) => ({
    x: cx + (i & 1 ? w : -w) / 2,
    y: cy + (i & 2 ? h : -h) / 2,
    z: cz + (i & 4 ? d : -d) / 2,
  });
  state.boxShape = BOX_EDGES.map(([a, b]) => v.addLine({
    start: corner(a), end: corner(b), color: '#45cfc4', linewidth: 2,
  }));
  v.render();
}

$('sitePick').addEventListener('change', (e) => {
  state.pickMode = e.target.checked;
  $('siteViewer').style.cursor = state.pickMode ? 'crosshair' : '';
});
$('siteFitRef').addEventListener('click', () => {
  if (state.prep && state.prep.reference) {
    setBox(state.prep.reference.centre, state.prep.reference.size, 'pocket');
    frameBox();
  }
});
$('siteWhole').addEventListener('click', () => {
  if (!state.prep) return;
  setBox(state.prep.bounds.centre, state.prep.bounds.size, 'pocket');
  if (state.siteViewer) { state.siteViewer.zoomTo(); state.siteViewer.render(); }
});
/** Frame the view on the box with some protein around it.
 *  zoomTo() on the whole structure leaves the site too small to aim at, and
 *  3Dmol's `within` selector is unreliable, so the nearby atoms are picked by
 *  serial number — which selects exactly. */
function frameBox(pad = 9) {
  const v = state.siteViewer;
  const m = state.siteModels.receptor;
  if (!v || !m) return;
  const [cx, cy, cz] = state.box.centre;
  const reach = Math.max(...state.box.size) / 2 + pad;
  let serials = [];
  try {
    for (const a of m.selectedAtoms({})) {
      if (Math.abs(a.x - cx) <= reach && Math.abs(a.y - cy) <= reach
          && Math.abs(a.z - cz) <= reach) serials.push(a.serial);
    }
  } catch { serials = []; }
  if (serials.length > 3) v.zoomTo({ model: m, serial: serials });
  else v.zoomTo();
  v.zoom(0.9);
  v.render();
}

$('siteRecentre').addEventListener('click', () => frameBox());

// Typing in the numeric fields drives the same box.
for (const name of ['centre_x', 'centre_y', 'centre_z', 'size_x', 'size_y', 'size_z']) {
  const el = $('form').elements[name];
  if (el) {
    el.addEventListener('change', () => {
      const f = $('form').elements;
      setBox(AXES.map((a) => parseFloat(f['centre_' + a].value) || 0),
             AXES.map((a) => parseFloat(f['size_' + a].value) || 22), 'field');
    });
  }
}

$('form').elements.receptor.addEventListener('change', prepareReceptor);
$('form').elements.reference.addEventListener('change', () => {
  if (fileOf('receptor')) prepareReceptor();
});

/* ----------------------------------------------------------- startup -- */
async function boot() {
  try {
    const r = await api('/api/health');
    const h = await r.json();
    $('env').innerHTML =
      `pydock <b>${h.pydock}</b> · rdkit <b>${h.rdkit}</b> · numba <b>${h.numba}</b>`;
    if (!h.example_available) {
      $('example').disabled = true;
      $('example').title = 'bundled example data not found on the server';
      $('example2').classList.add('hidden');
    }
  } catch {
    $('env').textContent = 'engine unreachable';
  }
  refreshJobs();
}

/* ------------------------------------------------------- site modes -- */
$('siteMode').addEventListener('click', (e) => {
  const btn = e.target.closest('button[data-mode]');
  if (!btn) return;
  const mode = btn.dataset.mode;
  for (const b of $('siteMode').children) b.classList.toggle('on', b === btn);
  for (const box of document.querySelectorAll('.mode')) {
    box.classList.toggle('hidden', box.dataset.mode !== mode);
    // Disabled inputs are not submitted, which is exactly what we want: the
    // server reads the site mode from which fields arrive.  File inputs are
    // the exception -- a reference ligand stays attached whatever the mode,
    // because it also supplies the RMSD the results are measured against.
    for (const inp of box.querySelectorAll('input:not([type=file])')) {
      inp.disabled = box.dataset.mode !== mode;
    }
  }
});
for (const box of document.querySelectorAll('.mode.hidden')) {
  for (const inp of box.querySelectorAll('input:not([type=file])')) inp.disabled = true;
}

/* ---------------------------------------------------------- submit --- */
$('form').addEventListener('submit', async (e) => {
  e.preventDefault();
  showErr('');
  const fd = new FormData($('form'));
  for (const key of ['detect_cavity', 'rotate_receptor_oh', 'keep_waters']) {
    fd.set(key, $('form').elements[key].checked ? 'true' : 'false');
  }
  await send('/api/jobs', fd);
});

/* A form that fails native validation never fires `submit`, so without this
 * the Dock button would look dead with no explanation. */
$('form').addEventListener('invalid', (e) => {
  const el = e.target;
  const where = el.type === 'file' ? 'Choose a file for' : 'Check';
  const name = (el.name || 'this field').replace(/_/g, ' ');
  showErr(`${where} ${name}: ${el.validationMessage}`);
}, true);

$('example').addEventListener('click', runExample);
$('example2').addEventListener('click', runExample);

async function runExample() {
  showErr('');
  const fd = new FormData();
  const scoring = $('form').querySelector('input[name=scoring]:checked');
  fd.set('scoring', scoring ? scoring.value : 'empirical');
  fd.set('n_runs', '4');
  fd.set('n_evals', '150000');
  await send('/api/jobs/example', fd);
}

async function send(url, fd) {
  $('go').disabled = true;
  $('example').disabled = true;
  try {
    const r = await api(url, { method: 'POST', body: fd });
    const body = await r.json().catch(() => ({}));
    if (!r.ok) throw new Error(body.detail || `request failed (${r.status})`);
    watch(body.id);
  } catch (err) {
    showErr(err.message || String(err));
    $('go').disabled = false;
    $('example').disabled = false;
  }
}

function showErr(msg) {
  $('formErr').textContent = msg;
  $('formErr').classList.toggle('hidden', !msg);
}

/* ----------------------------------------------------------- polling -- */
function watch(id) {
  state.jobId = id;
  state.framed = false;
  state.poses = [];
  state.referenceSdf = null;
  $('intro').classList.add('hidden');
  $('statusCard').classList.remove('hidden');
  $('viewerCard').classList.add('hidden');
  $('log').textContent = '';
  if (state.poll) clearInterval(state.poll);
  tick();
  state.poll = setInterval(tick, 1000);
}

async function tick() {
  if (!state.jobId) return;
  let job;
  try {
    const r = await api(`/api/jobs/${state.jobId}`);
    if (!r.ok) throw new Error('job is gone');
    job = await r.json();
  } catch {
    clearInterval(state.poll);
    return;
  }
  state.job = job;
  paintStatus(job);

  const finished = ['done', 'failed', 'cancelled'].includes(job.status);
  if (finished) {
    clearInterval(state.poll);
    $('go').disabled = false;
    $('example').disabled = false;
    $('cancel').classList.add('hidden');
    if (job.status === 'done' && job.results.length) await loadResults(job);
    if (job.status === 'failed') showErr(job.error || 'the job failed');
    refreshJobs();
  }
}

function paintStatus(job) {
  const running = ['queued', 'preparing', 'running'].includes(job.status);
  $('dot').className = 'dot ' + (running ? 'run'
    : job.status === 'done' ? 'done'
    : job.status === 'cancelled' ? '' : 'fail');
  $('statusText').textContent = job.status;
  const lig = job.ligands && job.ligands.length
    ? `${job.ligands.length} ligand${job.ligands.length > 1 ? 's' : ''} · ` : '';
  $('statusMeta').textContent =
    `${lig}${job.params.scoring} · ${job.params.n_runs} runs · ${job.elapsed}s`;
  $('fill').style.width = `${Math.round(job.fraction * 100)}%`;
  $('cancel').classList.toggle('hidden', !running);
  const log = (job.log || []).join('\n');
  if ($('log').textContent !== log) {
    const atBottom = $('log').scrollTop + $('log').clientHeight
      >= $('log').scrollHeight - 24;
    $('log').textContent = log;
    if (atBottom) $('log').scrollTop = $('log').scrollHeight;
  }
}

$('cancel').addEventListener('click', async () => {
  if (!state.jobId) return;
  $('cancel').disabled = true;
  await api(`/api/jobs/${state.jobId}/cancel`, { method: 'POST' });
  $('cancel').disabled = false;
});

/* ----------------------------------------------------------- results -- */
async function loadResults(job) {
  $('viewerCard').classList.remove('hidden');
  $('siteCard').classList.add('hidden');
  renderTabs(job);
  renderScreening(job);
  const best = job.results.length > 1 ? bestLigandIndex(job) : 0;
  await selectLigand(best);
}

function bestLigandIndex(job) {
  let best = 0;
  job.results.forEach((r, i) => {
    if (r.score < job.results[best].score) best = i;
  });
  return best;
}

function renderTabs(job) {
  const tabs = $('ligTabs');
  tabs.innerHTML = '';
  const span = document.createElement('span');
  span.className = 'muted';
  span.id = 'ligLabel';
  tabs.appendChild(span);
  // With a library the screening table is the navigator, so no tab strip.
  if (job.results.length > 1 && job.results.length <= 6) {
    job.results.forEach((r, i) => {
      const b = document.createElement('button');
      b.type = 'button';
      b.textContent = r.name;
      b.dataset.index = i;
      b.addEventListener('click', () => selectLigand(i));
      tabs.appendChild(b);
    });
  }
}

/* -------------------------------------------------------- screening -- */
function renderScreening(job) {
  const many = job.results.length > 1;
  $('screenCard').classList.toggle('hidden', !many);
  if (!many) return;

  $('dlRanking').href = `/api/jobs/${state.jobId}/ranking.csv`;
  const total = job.results.reduce((s, r) => s + r.seconds, 0);
  $('screenMeta').textContent =
    `${job.results.length} ligands · ${Math.round(total)}s total`;
  paintScreening(job);
}

function paintScreening(job) {
  const byLE = $('byEfficiency').checked;
  const rows = job.results.map((r) => ({ r, key: byLE ? r.efficiency : r.score }))
    .sort((a, b) => a.key - b.key);
  const tb = $('screenTable').querySelector('tbody');
  tb.innerHTML = '';
  const anyRmsd = job.results.some((r) => r.rmsd_to_reference != null);
  $('screenTable').querySelectorAll('th')[7].style.display = anyRmsd ? '' : 'none';

  rows.forEach(({ r }, i) => {
    const tr = document.createElement('tr');
    tr.dataset.index = r.index;
    const rmsdCell = anyRmsd
      ? `<td>${r.rmsd_to_reference == null ? '–'
        : `<span class="${r.rmsd_to_reference < 2 ? 'good' : 'bad'}">`
          + `${r.rmsd_to_reference.toFixed(2)}</span>`}</td>`
      : '';
    tr.innerHTML = `<td>${i + 1}</td>`
      + `<td class="name" title="${r.name}">${r.name}</td>`
      + `<td class="${i === 0 ? 'best' : ''}">${r.score.toFixed(3)}</td>`
      + `<td>${r.efficiency.toFixed(3)}</td>`
      + `<td>${r.n_heavy}</td><td>${r.n_rotatable}</td><td>${r.n_modes}</td>`
      + rmsdCell + `<td>${r.seconds}s</td>`;
    tr.addEventListener('click', () => selectLigand(r.index));
    tb.appendChild(tr);
  });
  highlightScreening();
}

function highlightScreening() {
  for (const tr of $('screenTable').querySelectorAll('tbody tr')) {
    tr.classList.toggle('on', Number(tr.dataset.index) === state.ligIndex);
  }
  for (const b of $('ligTabs').querySelectorAll('button')) {
    b.classList.toggle('on', Number(b.dataset.index) === state.ligIndex);
  }
}

$('byEfficiency').addEventListener('change', () => {
  if (state.job && state.job.results.length > 1) paintScreening(state.job);
});

async function selectLigand(index) {
  state.ligIndex = index;
  state.framed = false;
  const res = state.job.results[index];
  const label = $('ligLabel');
  if (label) {
    label.textContent = `${res.name} · ${res.n_heavy} heavy atoms · `
      + `${res.n_rotatable} rotors`;
  }
  highlightScreening();
  renderTable(res);
  $('dl').href = `/api/jobs/${state.jobId}/poses/${index}.sdf`;
  $('dl').classList.toggle('hidden', !res.sdf);

  const summary = [`${res.n_modes} distinct binding mode${res.n_modes > 1 ? 's' : ''}`,
    `${res.poses.length} poses kept`, `${res.seconds}s`];
  if (res.rmsd_to_reference != null) {
    summary.push(`top-ranked pose ${res.rmsd_to_reference} Å from the reference`);
  }
  $('summary').textContent = summary.join(' · ');

  await loadStructures(index);
  showPose(1);
}

function renderTable(res) {
  const tb = $('poseTable').querySelector('tbody');
  tb.innerHTML = '';
  const hasRmsd = res.poses.some((p) => p.rmsd != null);
  $('poseTable').querySelectorAll('th')[2].style.display = hasRmsd ? '' : 'none';

  for (const p of res.poses) {
    const tr = document.createElement('tr');
    tr.dataset.rank = p.rank;
    const terms = Object.entries(p.terms)
      .filter(([k]) => !['rotors', 'outside_box', 'metal', 'plp_steric'].includes(k))
      .map(([k, v]) => `${k} ${v}`).join('   ');
    const rmsdCell = hasRmsd
      ? `<td>${p.rmsd == null ? '–'
        : `<span class="${p.rmsd < 2 ? 'good' : 'bad'}">${p.rmsd.toFixed(2)}</span>`}</td>`
      : '';
    tr.innerHTML = `<td>#${p.rank}</td><td>${p.score.toFixed(3)}</td>`
      + rmsdCell + `<td>${p.run}</td><td class="terms">${terms}</td>`;
    tr.addEventListener('click', () => showPose(p.rank));
    tb.appendChild(tr);
  }
}

/* ------------------------------------------------------------ viewer -- */
function viewerReady() {
  if (!window.$3Dmol) {
    $('viewerNote').textContent =
      '3Dmol.js could not be loaded, so the structure cannot be drawn here. '
      + 'The scores and the pose files below are unaffected — download the SDF '
      + 'and open it in PyMOL or Chimera.';
    $('viewerNote').classList.remove('hidden');
    return false;
  }
  if (!state.viewer) {
    state.viewer = $3Dmol.createViewer($('viewer'),
      { backgroundColor: '#090d13', antialias: true });
  }
  return true;
}

async function loadStructures(index) {
  if (!viewerReady()) return;
  const v = state.viewer;
  v.clear();
  state.models = { receptor: null, pose: null, reference: null };

  const [pdb, sdf, ref] = await Promise.all([
    text(`/api/jobs/${state.jobId}/receptor.pdb`),
    text(`/api/jobs/${state.jobId}/poses/${index}.sdf`),
    text(`/api/jobs/${state.jobId}/reference.sdf`),
  ]);

  state.poses = splitSdf(sdf);
  state.referenceSdf = ref;

  if (pdb) state.models.receptor = v.addModel(pdb, 'pdb');
  if (ref) state.models.reference = v.addModel(ref, 'sdf');

  $('optRef').parentElement.classList.toggle('hidden', !ref);
  $('legend').innerHTML =
    `<span><i style="background:#e9a14b"></i>docked pose</span>`
    + (ref ? `<span><i style="background:#6ec98a"></i>reference</span>` : '')
    + `<span><i style="background:#56687e"></i>receptor</span>`;
}

async function text(url) {
  try {
    const r = await api(url);
    if (!r.ok) return null;
    const t = await r.text();
    return t && t.trim() ? t : null;
  } catch { return null; }
}

function splitSdf(sdf) {
  if (!sdf) return [];
  return sdf.split(/\$\$\$\$\s*\n?/).filter((s) => s.trim()).map((s) => s + '$$$$\n');
}

function showPose(rank) {
  state.poseRank = rank;
  for (const tr of $('poseTable').querySelectorAll('tbody tr')) {
    tr.classList.toggle('on', Number(tr.dataset.rank) === rank);
  }
  if (!state.viewer || !state.poses.length) return;
  const v = state.viewer;

  if (state.models.pose) { v.removeModel(state.models.pose); state.models.pose = null; }
  const record = state.poses[rank - 1];
  if (record) state.models.pose = v.addModel(record, 'sdf');

  restyle();
  if (!state.framed && state.models.pose) {
    v.zoomTo({ model: state.models.pose });
    v.zoom(0.75);
    state.framed = true;
  }
  v.render();
}

/* 3Dmol's `within` selector returns nothing once a `model` filter is combined
 * with it, so the pocket is worked out here from the coordinates instead:
 * every receptor residue with an atom inside `cut` of any ligand atom. */
function pocketResidues(cut) {
  const { receptor, pose } = state.models;
  if (!receptor || !pose) return [];
  let rec, lig;
  try {
    rec = receptor.selectedAtoms({});
    lig = pose.selectedAtoms({});
  } catch { return []; }
  if (!rec.length || !lig.length) return [];

  const c2 = cut * cut;
  const byChain = new Map();
  for (const a of rec) {
    for (const b of lig) {
      const dx = a.x - b.x, dy = a.y - b.y, dz = a.z - b.z;
      if (dx * dx + dy * dy + dz * dz <= c2) {
        const ch = a.chain || ' ';
        if (!byChain.has(ch)) byChain.set(ch, new Set());
        byChain.get(ch).add(a.resi);
        break;
      }
    }
  }
  return [...byChain.entries()].map(([chain, resis]) => ({
    chain, resi: [...resis],
  }));
}

function restyle() {
  const v = state.viewer;
  const { receptor, pose, reference } = state.models;
  v.removeAllSurfaces();
  v.setStyle({}, {});

  const pocket = ($('optPocket').checked || $('optSurface').checked)
    ? pocketResidues(5.0) : [];

  if (receptor) {
    v.setStyle({ model: receptor },
      { cartoon: { color: '#7d90a8', opacity: 0.6, thickness: 0.28 } });
    if ($('optPocket').checked) {
      for (const sel of pocket) {
        v.addStyle({ model: receptor, chain: sel.chain, resi: sel.resi },
          { stick: { radius: 0.1, colorscheme: 'lightGrayCarbon' } });
      }
    }
  }
  if (reference && $('optRef').checked) {
    v.setStyle({ model: reference },
      { stick: { radius: 0.08, colorscheme: 'greenCarbon' } });
  } else if (reference) {
    v.setStyle({ model: reference }, {});
  }
  if (pose) {
    v.setStyle({ model: pose },
      { stick: { radius: 0.17, colorscheme: 'orangeCarbon' } });
  }
  if ($('optSurface').checked && receptor && pocket.length) {
    try {
      for (const sel of pocket) {
        v.addSurface($3Dmol.SurfaceType.VDW,
          { opacity: 0.55, color: '#33455c' },
          { model: receptor, chain: sel.chain, resi: sel.resi });
      }
    } catch { /* surface generation can fail on odd residues; skip it */ }
  }
}

for (const id of ['optRef', 'optPocket', 'optSurface']) {
  $(id).addEventListener('change', () => {
    if (!state.viewer) return;
    restyle();
    state.viewer.render();
  });
}
$('reset').addEventListener('click', () => {
  if (!state.viewer) return;
  const sel = state.models.pose ? { model: state.models.pose } : {};
  state.viewer.zoomTo(sel);
  state.viewer.zoom(0.75);
  state.viewer.render();
});

/* ------------------------------------------------------- recent jobs -- */
async function refreshJobs() {
  let jobs = [];
  try {
    const r = await api('/api/jobs');
    jobs = (await r.json()).jobs || [];
  } catch { return; }
  $('jobsCard').classList.toggle('hidden', jobs.length === 0);
  const ul = $('jobs');
  ul.innerHTML = '';
  for (const j of jobs.slice(0, 12)) {
    const li = document.createElement('li');
    const score = j.results && j.results.length && j.results[0].poses.length
      ? j.results[0].poses[0].score.toFixed(2) : '';
    li.innerHTML = `<span class="jid">${j.id}</span>`
      + `<span class="jname">${(j.ligands || []).join(', ') || j.receptor}</span>`
      + `<span class="jst">${j.status}${score ? ' · ' + score : ''}</span>`;
    const open = document.createElement('button');
    open.type = 'button';
    open.textContent = 'open';
    open.addEventListener('click', () => watch(j.id));
    li.appendChild(open);
    const del = document.createElement('button');
    del.type = 'button';
    del.textContent = '×';
    del.title = 'delete this job and its files';
    del.addEventListener('click', async () => {
      await api(`/api/jobs/${j.id}`, { method: 'DELETE' });
      if (state.jobId === j.id) {
        state.jobId = null;
        $('statusCard').classList.add('hidden');
        $('viewerCard').classList.add('hidden');
        $('intro').classList.remove('hidden');
      }
      refreshJobs();
    });
    li.appendChild(del);
    ul.appendChild(li);
  }
}

boot();
