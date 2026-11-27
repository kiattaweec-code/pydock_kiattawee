/* The preparation page: load a complex, pick the ligand, get clean inputs.
 *
 * The viewer shows the structure exactly as uploaded -- ligand still in
 * place -- because choosing it is the task.  Hetero groups are drawn as
 * spheres so they stand out against the cartoon, and clicking one selects
 * it, which is quicker than reading a list of residue numbers.
 */
'use strict';

const $ = (id) => document.getElementById(id);

const state = {
  native: { folder_picker: false, default_folder: '' },
  prepId: null,
  meta: null,
  selected: null,
  viewer: null,
  built: null,
  busy: false,
};

const COLOURS = {
  ligand: '0x45cfc4',
  metal: '0xe9a14b',
  ion: '0x8c9ab0',
  buffer: '0x6b7889',
  water: '0x3a4656',
};

/* The app runs on the person's own machine, so it can write where they
 * ask. In its own window there is a real folder dialog; in a browser tab
 * there is not, and a typed path plus a ZIP is the next best thing. */
fetch('/api/prep/native').then((r) => r.json()).then((info) => {
  state.native = info;
  const box = $('saveFolder');
  if (info.folder_picker) {
    box.value = '';
    box.placeholder = 'choose a folder\u2026';
    box.readOnly = true;
    $('saveBtn').textContent = 'Choose folder\u2026';
  } else if (info.default_folder) {
    box.value = info.default_folder;
  }
}).catch(() => { /* saving is simply unavailable; downloads still work */ });

/* ------------------------------------------------------------- helpers -- */
function status(el, text, kind) {
  const node = $(el);
  node.textContent = text || '';
  node.className = 'status' + (kind ? ' ' + kind : '');
}

function bytes(n) {
  if (n < 1024) return n + ' B';
  if (n < 1024 * 1024) return (n / 1024).toFixed(0) + ' kB';
  return (n / (1024 * 1024)).toFixed(1) + ' MB';
}

async function post(url, form) {
  const reply = await fetch(url, { method: 'POST', body: form });
  let data = null;
  try { data = await reply.json(); } catch (e) { /* not json */ }
  if (!reply.ok) {
    throw new Error((data && data.detail) || `${reply.status} ${reply.statusText}`);
  }
  return data;
}

/* -------------------------------------------------------------- viewer -- */
function showStructure(pdbText) {
  const box = $('viewer');
  if (!window.$3Dmol) {
    status('loadStatus', 'the 3D viewer could not be loaded', 'err');
    return;
  }
  $('placeholder').classList.add('hidden');
  if (state.viewer) { state.viewer.clear(); }
  else {
    state.viewer = $3Dmol.createViewer(box, { backgroundColor: '#0d1117' });
  }
  const view = state.viewer;
  view.addModel(pdbText, 'pdb');

  view.setStyle({}, { cartoon: { color: '#3f4a5a', opacity: 0.85 } });
  // Anything that is not protein is drawn solid, so it can be found and
  // clicked; the protein stays a quiet backdrop.
  for (const group of (state.meta.ligands || []).concat(state.meta.others || [])) {
    view.setStyle(
      { chain: group.chain, resi: group.resseq },
      { stick: { radius: 0.17, colorscheme: 'default' },
        sphere: { radius: 0.34 } });
  }
  for (const group of state.meta.ligands || []) {
    view.addStyle({ chain: group.chain, resi: group.resseq },
                  { stick: { radius: 0.2, color: COLOURS.ligand } });
  }
  view.setClickable({}, true, (atom) => {
    const hit = pickFromAtom(atom);
    if (hit) select(hit.key);
  });
  view.zoomTo();
  view.render();
  highlight();
  legend();
}

function pickFromAtom(atom) {
  if (!atom) return null;
  const all = (state.meta.ligands || []).concat(state.meta.others || []);
  return all.find((g) => g.chain === atom.chain && g.resseq === atom.resi) || null;
}

function highlight() {
  const view = state.viewer;
  if (!view || !state.meta) return;
  for (const group of state.meta.ligands || []) {
    const on = group.key === state.selected;
    view.setStyle({ chain: group.chain, resi: group.resseq },
                  { stick: { radius: on ? 0.26 : 0.18,
                             color: on ? COLOURS.ligand : undefined,
                             colorscheme: on ? undefined : 'default' },
                    sphere: { radius: on ? 0.42 : 0.3 } });
  }
  const chosen = (state.meta.ligands || []).find((g) => g.key === state.selected);
  if (chosen) {
    view.zoomTo({ chain: chosen.chain, resi: chosen.resseq });
    view.zoom(0.35);
  }
  view.render();
}

function legend() {
  const parts = [];
  if ((state.meta.ligands || []).length) {
    parts.push(`<span><i style="background:#45cfc4"></i>ligand candidate</span>`);
  }
  const kinds = new Set((state.meta.others || []).map((o) => o.kind));
  if (kinds.has('metal')) parts.push(`<span><i style="background:#e9a14b"></i>metal</span>`);
  if (kinds.has('ion')) parts.push(`<span><i style="background:#8c9ab0"></i>ion</span>`);
  if (kinds.has('buffer')) parts.push(`<span><i style="background:#6b7889"></i>buffer</span>`);
  parts.push(`<span><i style="background:#3f4a5a"></i>protein</span>`);
  parts.push(`<span>click a molecule to choose it</span>`);
  $('legend').innerHTML = parts.join('');
}

/* ------------------------------------------------------------ the list -- */
function renderPicks() {
  const list = $('picks');
  const ligands = state.meta.ligands || [];
  list.innerHTML = '';
  for (const group of ligands) {
    const li = document.createElement('li');
    const button = document.createElement('button');
    button.type = 'button';
    button.setAttribute('aria-pressed', String(group.key === state.selected));
    button.innerHTML =
      `<span class="name">${group.label}</span>` +
      `<span class="meta">${group.n_heavy} atoms &middot; ${group.formula}</span>`;
    button.addEventListener('click', () => select(group.key));
    li.appendChild(button);
    list.appendChild(li);
  }
  const none = $('noLigand');
  if (!ligands.length) {
    none.textContent =
      'No bound ligand was found. You can still clean the structure and use ' +
      'it as a receptor, but you will have to place the docking box by hand.';
  } else {
    none.textContent = '';
  }
}

function select(key) {
  state.selected = (state.selected === key) ? null : key;
  renderPicks();
  highlight();
  const chosen = (state.meta.ligands || []).find((g) => g.key === state.selected);
  $('smilesSection').classList.toggle('hidden', !chosen);
  if (chosen) {
    $('smilesHint').textContent =
      `${chosen.label} has ${chosen.n_heavy} heavy atoms (${chosen.formula}). ` +
      `A SMILES with a different count is rejected.`;
  }
}

/* --------------------------------------------------------------- load --- */
$('file').addEventListener('change', async (event) => {
  const file = event.target.files && event.target.files[0];
  if (!file) return;
  status('loadStatus', `reading ${file.name}…`);
  for (const id of ['pickSection', 'smilesSection', 'optionSection', 'resultSection']) {
    $(id).classList.add('hidden');
  }
  state.selected = null;
  state.built = null;

  const form = new FormData();
  form.append('structure', file);
  try {
    const meta = await post('/api/prep/inspect', form);
    state.prepId = meta.prep_id;
    state.meta = meta;

    const reply = await fetch(`/api/prep/${meta.prep_id}/viewer.pdb`);
    showStructure(await reply.text());

    const chains = Object.keys(meta.chains).length;
    status('loadStatus',
      `${meta.name}: ${meta.n_atoms.toLocaleString()} atoms, ` +
      `${meta.polymer_residues} residues, ${chains} chain${chains === 1 ? '' : 's'}`,
      'good');

    renderPicks();
    $('pickSection').classList.remove('hidden');
    $('optionSection').classList.remove('hidden');

    $('hCount').textContent = meta.n_hydrogens ? `${meta.n_hydrogens} present` : 'none present';
    $('optH').disabled = !meta.n_hydrogens;
    $('wCount').textContent = meta.waters ? `${meta.waters} present` : 'none';
    $('optWater').disabled = !meta.waters;
    const buffers = (meta.others || []).filter((o) => o.kind === 'buffer');
    const metals = (meta.others || []).filter((o) => o.kind === 'metal' || o.kind === 'ion');
    $('bCount').textContent = buffers.length ? `${buffers.length} present` : 'none';
    $('optBuffer').disabled = !buffers.length;
    $('mCount').textContent = metals.length
      ? metals.map((m) => m.resname).join(', ') : 'none';
    $('optMetal').disabled = !metals.length;

    // The biggest hetero group is nearly always the ligand, so offer it
    // already selected rather than making the obvious click.
    if ((meta.ligands || []).length) select(meta.ligands[0].key);
  } catch (err) {
    status('loadStatus', String(err.message || err), 'err');
  }
});

/* -------------------------------------------------------------- build --- */
$('build').addEventListener('click', async () => {
  if (state.busy || !state.prepId) return;
  state.busy = true;
  $('build').disabled = true;
  status('buildStatus', 'preparing…');

  const form = new FormData();
  form.append('prep_id', state.prepId);
  if (state.selected) form.append('ligand_key', state.selected);
  form.append('smiles', $('smiles').value.trim());
  form.append('remove_hydrogens', $('optH').checked);
  form.append('remove_waters', $('optWater').checked);
  form.append('remove_buffers', $('optBuffer').checked);
  form.append('keep_metals', $('optMetal').checked);

  try {
    const out = await post('/api/prep/build', form);
    state.built = out;
    status('buildStatus', 'done', 'good');
    showResult(out);
  } catch (err) {
    status('buildStatus', String(err.message || err), 'err');
  } finally {
    state.busy = false;
    $('build').disabled = false;
  }
});

function showResult(out) {
  const rows = [];
  const dropped = out.receptor.dropped || {};
  rows.push(['Receptor', `${out.receptor.atoms.toLocaleString()} atoms`]);
  const removed = Object.entries(dropped)
    .map(([what, n]) => `${n} ${what}${n === 1 ? '' : 's'}`).join(', ');
  rows.push(['Removed', removed || 'nothing']);
  if (out.ligand) {
    rows.push(['Ligand', `${out.ligand.label} — ${out.ligand.n_atoms} atoms, ${out.ligand.formula}`]);
    rows.push(['Bond orders', out.ligand.bond_orders === 'template'
      ? 'from your SMILES'
      : out.ligand.bond_orders === 'perceived'
        ? 'worked out from the geometry'
        : 'not determined — all single']);
  }
  $('summary').innerHTML = rows
    .map(([k, v]) => `<dt>${k}</dt><dd>${v}</dd>`).join('');

  const links = [[`/api/prep/${state.prepId}/receptor.pdb`, 'receptor_clean.pdb', 'protein, ligand removed']];
  if (out.ligand) {
    links.push([out.ligand.mol2, 'ligand.mol2', 'extracted ligand, Sybyl types']);
    links.push([out.ligand.sdf, 'ligand.sdf', 'the same ligand as SDF']);
  }
  $('files').innerHTML = links.map(([href, name, what]) =>
    `<a href="${href}" download><span class="what">${name}</span>` +
    `<span class="size">${what}</span></a>`).join('');

  const warnings = (out.warnings || []).map((w) => `<p class="warn">${w}</p>`);
  if (out.ligand && out.ligand.trusted) {
    warnings.unshift('<p class="warn good">The ligand matches the SMILES you ' +
      'gave, and its coordinates are unchanged from the crystal structure.</p>');
  }
  $('files').insertAdjacentHTML('afterend', '');
  const old = document.querySelectorAll('#resultSection .warn');
  old.forEach((n) => n.remove());
  $('files').insertAdjacentHTML('beforebegin', warnings.join(''));

  // A real link, so it can be opened in a new tab and its target is visible
  // before clicking.
  const params = new URLSearchParams({ prep: state.prepId });
  if (out.ligand) params.set('ligand', out.ligand.label);
  $('toDock').href = 'index.html?' + params.toString();

  $('dlZip').href = `/api/prep/${state.prepId}/bundle.zip`;
  status('saveStatus', '');

  $('resultSection').classList.remove('hidden');
  // The panel is long enough that the results land below the fold, which
  // reads as nothing having happened.
  $('resultSection').scrollIntoView({ behavior: 'smooth', block: 'nearest' });
}


/* ------------------------------------------------------------- saving --- */
$('saveBtn').addEventListener('click', async () => {
  if (!state.prepId || !state.built) return;
  const button = $('saveBtn');
  button.disabled = true;
  status('saveStatus', 'saving…');

  const form = new FormData();
  if (state.native.folder_picker) form.append('browse', 'true');
  form.append('folder', $('saveFolder').value.trim());

  try {
    const out = await post(`/api/prep/${state.prepId}/save`, form);
    if (out.cancelled) { status('saveStatus', ''); return; }
    $('saveFolder').value = out.folder;
    const names = out.files.map((f) => f.name).join(', ');
    status('saveStatus',
      `${out.files.length} files saved` + (out.renamed
        ? ' — numbered, so nothing already there was replaced' : ''),
      'good');
    $('saveStatus').insertAdjacentHTML('beforeend',
      `<div class="savedPath">${out.folder}<br>${names}</div>`);
  } catch (err) {
    status('saveStatus', String(err.message || err), 'err');
  } finally {
    button.disabled = false;
  }
});
