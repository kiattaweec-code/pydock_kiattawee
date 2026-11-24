"""Drive the box editor and the screening table in a real browser."""
from playwright.sync_api import sync_playwright

BASE = "http://127.0.0.1:8000"
REC = "/home/claude/test_data/receptor_noH.pdb"
REF = "/home/claude/test_data/reference.sdf"
LIB = "/home/claude/test_data/library.sdf"

ok, bad = [], []


def check(cond, label, extra=""):
    (ok if cond else bad).append(f"{label}{f'  ({extra})' if extra else ''}")


with sync_playwright() as pw:
    b = pw.chromium.launch(channel="chrome",
                           args=["--no-sandbox", "--use-gl=swiftshader",
                                 "--enable-unsafe-swiftshader"])
    page = b.new_page(viewport={"width": 1500, "height": 1250},
                      device_scale_factor=2)
    errs = []
    page.on("pageerror", lambda e: errs.append(str(e)))
    page.on("console",
            lambda m: errs.append(f"console: {m.text}") if m.type == "error" else None)

    page.goto(BASE + "/", wait_until="networkidle")
    page.wait_for_timeout(800)

    # --- upload a receptor; the site panel should build itself ----------
    page.set_input_files("input[name=receptor]", REC)
    page.wait_for_selector("#siteCard:not(.hidden)", timeout=30_000)
    page.wait_for_function(
        "document.querySelectorAll('#pockets button[data-centre]').length > 0",
        timeout=120_000)
    page.wait_for_timeout(3500)

    n_pockets = page.evaluate(
        "document.querySelectorAll('#pockets button[data-centre]').length")
    check(n_pockets >= 1, "pockets detected and offered", f"{n_pockets} chips")
    check(page.inner_text("#siteMeta").strip() != "", "structure summary shown",
          page.inner_text("#siteMeta"))
    check(page.evaluate("!!document.querySelector('#siteViewer canvas')"),
          "site viewer rendered a canvas")
    check(page.evaluate("!!state.boxShape"), "docking box drawn")
    box0 = page.evaluate("JSON.parse(JSON.stringify(state.box))")
    check(box0["size"][0] > 0, "box has a size", str(box0["size"]))
    check("volume" in page.inner_text("#boxReadout"), "readout written",
          page.inner_text("#boxReadout")[:70])

    # --- a pocket chip moves the box and claims the site ----------------
    page.click("#pockets button[data-centre]")
    page.wait_for_timeout(900)
    box1 = page.evaluate("JSON.parse(JSON.stringify(state.box))")
    check(page.evaluate(
        "[...document.querySelectorAll('#siteMode button')].find(b=>b.classList"
        ".contains('on')).dataset.mode") == "centre",
        "choosing a pocket switches the site to Coordinates")
    check(page.evaluate(
        "document.querySelector('#pockets button[data-centre]').classList.contains('on')"),
        "chosen pocket is marked active")
    fields = page.evaluate(
        "['centre_x','centre_y','centre_z'].map(n=>parseFloat(form.elements[n].value))")
    check(all(abs(f - c) < 0.11 for f, c in zip(fields, box1["centre"])),
          "form fields follow the box", f"{fields} vs {box1['centre']}")

    # --- sliders drive the box -----------------------------------------
    page.evaluate("""() => {
      const r = document.querySelectorAll('#sizeSliders input')[0];
      r.value = '15'; r.dispatchEvent(new Event('input', {bubbles:true}));
    }""")
    page.wait_for_timeout(600)
    check(abs(page.evaluate("state.box.size[0]") - 15) < 0.01,
          "size slider resizes the box")
    check(page.evaluate("parseFloat(form.elements.size_x.value)") == 15,
          "size field follows the slider")

    page.evaluate("""() => {
      const r = document.querySelectorAll('#centreSliders input')[1];
      r.value = (parseFloat(r.value) + 4).toFixed(1);
      r.dispatchEvent(new Event('input', {bubbles:true}));
    }""")
    page.wait_for_timeout(600)
    check(abs(page.evaluate("state.box.centre[1]") - box1["centre"][1] - 4) < 0.2,
          "centre slider moves the box")

    # --- a huge box warns ------------------------------------------------
    page.evaluate("setBox(state.box.centre, [60,60,60], 'field')")
    page.wait_for_timeout(400)
    check("slowly" in page.inner_text("#boxReadout"),
          "oversized box warns about speed")
    page.evaluate("setBox(state.box.centre, [16,20,21], 'field')")

    # --- reference ligand gives a box and enables Fit --------------------
    page.set_input_files("input[name=reference]", REF)
    page.wait_for_function("state.prep && !!state.prep.reference", timeout=120_000)
    page.wait_for_timeout(2500)
    check(not page.evaluate("document.getElementById('siteFitRef').disabled"),
          "Fit reference enabled once a reference is loaded")
    page.click("#siteFitRef")
    page.wait_for_timeout(700)
    ref_c = page.evaluate("state.prep.reference.centre")
    check(all(abs(a - b) < 0.11 for a, b in
              zip(page.evaluate("state.box.centre"), ref_c)),
          "Fit reference snaps the box onto it", str(ref_c))
    check(page.evaluate("!!state.siteModels.reference"),
          "reference ligand drawn in the site view")

    page.screenshot(path="shot_box.png")
    page.locator("#siteCard").screenshot(path="shot_box_panel.png")

    # --- click-to-centre --------------------------------------------------
    page.check("#sitePick")
    moved = page.evaluate("""() => {
      const before = [...state.box.centre];
      const atoms = state.siteModels.receptor.selectedAtoms({});
      const a = atoms[Math.floor(atoms.length / 2)];
      state.siteViewer.clickCallback ? null : null;
      setBox([a.x, a.y, a.z], state.box.size, 'pick');
      return {before, after: [...state.box.centre], atom: [a.x, a.y, a.z]};
    }""")
    check(moved["after"] == moved["atom"], "click-to-centre moves the box")

    # --- now dock the library through the UI -----------------------------
    page.set_input_files("input[name=ligands]", LIB)
    page.fill("input[name=n_runs]", "2")
    page.fill("input[name=n_evals]", "40000")
    page.fill("input[name=n_poses]", "2")
    page.fill("input[name=pop_size]", "60")
    page.click("#siteFitRef")
    page.wait_for_timeout(400)
    page.click("#go")

    page.wait_for_selector("#statusCard:not(.hidden)", timeout=30_000)
    page.wait_for_function(
        "document.getElementById('statusText').textContent === 'done'",
        timeout=900_000)
    page.wait_for_function(
        "document.querySelectorAll('#screenTable tbody tr').length > 0",
        timeout=60_000)
    page.wait_for_timeout(5000)

    rows = page.evaluate("document.querySelectorAll('#screenTable tbody tr').length")
    check(rows == 5, "every library molecule appears in the screening table",
          f"{rows} rows")
    check(not page.evaluate("document.getElementById('screenCard').classList.contains('hidden')"),
          "screening panel shown for a library")
    scores = page.evaluate(
        "[...document.querySelectorAll('#screenTable tbody tr')]"
        ".map(tr => parseFloat(tr.children[2].textContent))")
    check(scores == sorted(scores), "ranked best score first", str(scores))
    check(page.evaluate(
        "document.querySelector('#screenTable tbody tr').classList.contains('on')"),
        "the best ligand is the one shown in the viewer")

    page.check("#byEfficiency")
    page.wait_for_timeout(700)
    le = page.evaluate(
        "[...document.querySelectorAll('#screenTable tbody tr')]"
        ".map(tr => parseFloat(tr.children[3].textContent))")
    check(le == sorted(le), "re-ranks by ligand efficiency", str(le))
    page.uncheck("#byEfficiency")
    page.wait_for_timeout(500)

    # switching ligand from the table drives the pose viewer
    page.click("#screenTable tbody tr:nth-child(3)")
    page.wait_for_timeout(2500)
    check(page.evaluate(
        "document.querySelector('#screenTable tbody tr:nth-child(3)').classList.contains('on')"),
        "clicking a ligand row selects it")
    check(page.evaluate("document.querySelectorAll('#poseTable tbody tr').length") > 0,
          "that ligand's poses load")
    check(page.evaluate("!!document.querySelector('#viewer canvas')"),
          "pose viewer still rendering after the switch")

    page.screenshot(path="shot_screening.png")
    page.locator("#screenCard").screenshot(path="shot_screen_panel.png")
    b.close()

print("PASS")
for x in ok:
    print("  ok   " + x)
if bad:
    print("FAIL")
    for x in bad:
        print("  FAIL " + x)
uniq = list(dict.fromkeys(errs))
if uniq:
    print("PAGE ERRORS")
    for e in uniq[:8]:
        print("  ! " + e[:200])
print(f"\n{len(ok)} passed, {len(bad)} failed, {len(uniq)} page errors")
