"""Screenshot the preparation page in a real browser, with WebGL.

jsdom checks the wiring but has no graphics, so this is the only way to see
whether the structure actually renders and the panel reads properly.

    python prep_shot.py [base-url] [pdb-file] [smiles-file]
"""
import sys
from pathlib import Path

from playwright.sync_api import sync_playwright

BASE = sys.argv[1] if len(sys.argv) > 1 else "http://127.0.0.1:8079"
PDB = sys.argv[2] if len(sys.argv) > 2 else "/tmp/complex.pdb"
SMILES = Path(sys.argv[3] if len(sys.argv) > 3 else "/tmp/smiles.txt").read_text().strip()

errors = []

with sync_playwright() as pw:
    browser = pw.chromium.launch(channel="chrome",
                                 args=["--no-sandbox", "--use-gl=swiftshader",
                                       "--enable-unsafe-swiftshader"])
    page = browser.new_page(viewport={"width": 1440, "height": 900})
    page.on("pageerror", lambda e: errors.append(f"pageerror: {e}"))
    page.on("console", lambda m: errors.append(f"console.{m.type}: {m.text}")
            if m.type == "error" else None)

    page.goto(f"{BASE}/prepare.html")
    page.wait_for_timeout(1200)
    page.screenshot(path="/tmp/prep_1_empty.png")

    page.set_input_files("#file", PDB)
    page.wait_for_selector("#picks button", timeout=60_000)
    page.wait_for_timeout(2500)
    page.screenshot(path="/tmp/prep_2_loaded.png")
    print("loaded:", page.inner_text("#loadStatus"))
    print("candidate:", page.inner_text("#picks button"))

    page.fill("#smiles", SMILES)
    page.click("#build")
    page.wait_for_selector("#resultSection:not(.hidden)", timeout=120_000)
    page.wait_for_timeout(900)
    page.screenshot(path="/tmp/prep_3_built.png")
    print("summary:", " | ".join(page.inner_text("#summary").split("\n")))
    print("handoff:", page.get_attribute("#toDock", "href"))

    # And the docking page receiving it, which is where a broken handoff
    # would show up as an empty form.
    page.click("#toDock")
    page.wait_for_selector("#preparedCard:not(.hidden)", timeout=60_000)
    page.wait_for_selector("#siteCard:not(.hidden)", timeout=60_000)
    page.wait_for_function(
        "!document.getElementById('siteMeta').textContent.includes('reading')",
        timeout=180_000)
    page.wait_for_timeout(2500)
    page.screenshot(path="/tmp/prep_4_handoff.png")
    print("site:", page.inner_text("#siteMeta"))
    print("banner:", page.inner_text("#preparedNote"))

    browser.close()

print("\nbrowser errors:", errors if errors else "none")
