"""Build a single self-contained HTML page that replays a real docking job.

The point is that the *actual* front end runs: index.html, style.css and
app.js are inlined unchanged, and only `fetch` is replaced by a stub that
serves recorded responses instead of a live server.  So what the page shows is
the real interface and real results, not a mock-up of them.

The box editor is genuinely interactive here, because once `/api/prepare` has
answered, placing the box is all client-side work.
"""
import base64
import json
from pathlib import Path

STATIC = Path("pydock_web/static")
ASSETS = Path("demo_assets")
OUT = Path("pydock_demo.html")

html = (STATIC / "index.html").read_text()
css = (STATIC / "style.css").read_text()
appjs = (STATIC / "app.js").read_text()
mol = (STATIC / "vendor" / "3Dmol-min.js").read_text()
ico = base64.b64encode((STATIC / "favicon.ico").read_bytes()).decode()

job = json.loads((ASSETS / "job.json").read_text())
prep = json.loads((ASSETS / "prepare.json").read_text())
poses = [ (ASSETS / f"poses_{i}.sdf").read_text()
          for i in range(len(job["results"])) ]

data = {
    "job": job,
    "prepare": prep,
    "receptor": (ASSETS / "receptor.pdb").read_text(),
    "reference": (ASSETS / "reference.sdf").read_text(),
    "poses": poses,
    "ranking": (ASSETS / "ranking.csv").read_text(),
    "health": {
        "status": "ok", "pydock": "0.1.0", "rdkit": "2026.09.1",
        "numba": "0.68.1", "queued": 0, "example_available": True,
        "scoring_functions": ["empirical", "chemplp", "forcefield"],
    },
}

BANNER = """
<div class="demoBanner">
  <strong>Static demo &mdash; but the box editor really works.</strong>
  An example receptor is loaded for you: drag the box, click a detected pocket,
  tick <em>click to centre</em> and click an atom. Pressing <em>Dock</em>
  replays a recorded 5-ligand screening run (ergotamine redocked to
  <b>0.72&nbsp;&Aring;</b> from its crystal pose). Nothing computes here;
  to dock your own structures, run the server from <code>pydock_web/</code>.
</div>
"""

BANNER_CSS = """
.demoBanner {
  margin: 0 16px; padding: 10px 14px; border-radius: 7px;
  background: #1b2430; border: 1px solid #2d4a52; border-left: 3px solid #45cfc4;
  color: #9aa7b8; font-size: 12.5px; line-height: 1.55; max-width: 1568px;
}
.demoBanner strong { color: #45cfc4; }
.demoBanner em { color: #dfe6f0; font-style: normal; }
.demoBanner b { color: #dfe6f0; font-weight: 600; font-family: var(--mono); }
.demoBanner code {
  font-family: var(--mono); font-size: 11.5px; color: #e9a14b;
  background: #0d1117; padding: 1px 5px; border-radius: 3px;
}
@media (min-width: 1601px) { .demoBanner { margin: 0 auto; } }
"""

# The fetch stub: the same URLs the client already calls, served from the
# recording.  Timings are compressed; everything else is verbatim.
STUB = r"""
(function () {
  const D = window.__PYDOCK_DEMO__;
  const REPLAY_MS = 19000;
  const JOB = D.job;
  const LOG = JOB.log.slice();
  const started = new Map();
  let counter = 0;

  const reply = (body) => Promise.resolve({
    ok: true, status: 200,
    json: () => Promise.resolve(typeof body === 'string' ? JSON.parse(body) : body),
    text: () => Promise.resolve(typeof body === 'string' ? body : JSON.stringify(body)),
  });
  const fail = (status, detail) => Promise.resolve({
    ok: false, status,
    json: () => Promise.resolve({ detail }),
    text: () => Promise.resolve(detail),
  });

  function snapshot(id) {
    const t0 = started.get(id);
    if (t0 == null) return null;
    const dt = Date.now() - t0;
    const f = Math.min(1, dt / REPLAY_MS);
    const j = JSON.parse(JSON.stringify(JOB));
    j.id = id;
    if (f < 1) {
      const n = Math.max(1, Math.floor(LOG.length * f));
      j.status = dt < 900 ? 'queued' : 'running';
      j.fraction = f;
      j.log = LOG.slice(0, n);
      j.message = (LOG[n - 1] || '').slice(11) || 'starting';
      j.results = [];
      j.elapsed = Math.round(dt / 100) / 10;
    } else {
      j.status = 'done';
      j.fraction = 1;
      j.log = LOG;
      j.message = 'finished';
    }
    return j;
  }

  window.fetch = function (url, opts) {
    url = String(url).replace(/^https?:\/\/[^/]+/, '');
    const method = (opts && opts.method) || 'GET';

    if (url === '/api/health') return reply(D.health);
    if (url === '/api/prepare') return reply(D.prepare);

    if (method === 'POST' && /^\/api\/jobs(\/example)?$/.test(url)) {
      const id = 'demo' + (++counter).toString().padStart(8, '0');
      started.set(id, Date.now());
      return reply({ id, status: 'queued' });
    }
    if (method === 'POST' && /\/cancel$/.test(url)) {
      return reply({ cancelled: false, status: 'running' });
    }
    if (method === 'DELETE') {
      started.delete(url.split('/').pop());
      return reply({ deleted: true });
    }
    if (url === '/api/jobs') {
      const jobs = [...started.keys()].reverse().map((id) => {
        const j = snapshot(id);
        return { id, status: j.status, receptor: j.receptor, ligands: j.ligands,
                 params: j.params, results: j.results, elapsed: j.elapsed };
      });
      return reply({ jobs, queued: 0 });
    }
    const m = url.match(/^\/api\/jobs\/([^/]+)$/);
    if (m) {
      const j = snapshot(m[1]);
      return j ? reply(j) : fail(404, 'no such job');
    }
    if (/\/receptor\.pdb/.test(url)) return reply(D.receptor);
    if (/\/reference\.sdf/.test(url)) return reply(D.reference);
    if (/\/ranking\.csv/.test(url)) return reply(D.ranking);
    const p = url.match(/\/poses\/(\d+)\.sdf/);
    if (p) return reply(D.poses[Number(p[1])] || D.poses[0]);
    return fail(404, 'not available in the static demo');
  };

  function attach(inputName, text, filename, type) {
    const input = document.querySelector(`input[name="${inputName}"]`);
    if (!input) return;
    const dt = new DataTransfer();
    dt.items.add(new File([text], filename, { type }));
    input.files = dt.files;
  }

  document.addEventListener('DOMContentLoaded', () => {
    for (const el of document.querySelectorAll('input[required]')) el.required = false;

    // Hand the page the example structures so the box editor is live without
    // the visitor having to find a PDB file.
    attach('receptor', D.prepare.pdb, 'receptor.pdb', 'chemical/x-pdb');
    attach('reference', D.reference, 'reference.sdf', 'chemical/x-mdl-sdfile');
    attach('ligands', D.poses.join(''), 'library.sdf', 'chemical/x-mdl-sdfile');
    setTimeout(() => {
      document.querySelector('input[name="receptor"]')
        .dispatchEvent(new Event('change', { bubbles: true }));
    }, 150);

    for (const id of ['dl', 'dlRanking']) {
      const a = document.getElementById(id);
      if (!a) continue;
      const label = a.textContent;
      a.addEventListener('click', (e) => {
        e.preventDefault();
        a.textContent = 'download needs the server';
        setTimeout(() => { a.textContent = label; }, 2200);
      });
    }
  });
})();
"""

page = html
page = page.replace('<link rel="icon" href="favicon.ico">',
                    f'<link rel="icon" href="data:image/x-icon;base64,{ico}">')
page = page.replace('<link rel="stylesheet" href="style.css">',
                    f"<style>\n{css}\n{BANNER_CSS}</style>")
page = page.replace("<title>pydock</title>",
                    "<title>pydock — flexible-ligand docking in the browser</title>")
page = page.replace('<main class="layout">', BANNER + '<main class="layout">')

start = page.index("<!-- 3Dmol.js is vendored")
end = page.index("</body>")
page = page[:start] + (
    f"<script>{mol}</script>\n"
    f"<script>window.__PYDOCK_DEMO__ = {json.dumps(data)};</script>\n"
    f"<script>{STUB}</script>\n"
    f"<script>{appjs}</script>\n"
) + page[end:]

OUT.write_text(page)
print(f"wrote {OUT}  {len(page) / 1024:.0f} KB")
assert "3Dmol.org" not in page[page.index("</style>"):], "no remote scripts remain"
print(f"inlined: 3Dmol, css, app.js, prepare ({len(prep['pockets'])} pockets), "
      f"{len(poses)} ligands' poses, receptor, reference, ranking")
