"""Screenshot the running app: idle, then with the finished example loaded."""
import sys

from playwright.sync_api import sync_playwright

BASE = "http://127.0.0.1:8000"
JOB = sys.argv[1] if len(sys.argv) > 1 else None

with sync_playwright() as pw:
    b = pw.chromium.launch(channel="chrome",
                           args=["--no-sandbox", "--use-gl=swiftshader",
                                 "--enable-unsafe-swiftshader"])
    page = b.new_page(viewport={"width": 1500, "height": 1180},
                      device_scale_factor=2)
    errs = []
    page.on("pageerror", lambda e: errs.append(str(e)))
    page.on("console", lambda m: errs.append(f"console.{m.type}: {m.text}")
            if m.type == "error" else None)

    page.goto(BASE + "/", wait_until="networkidle")
    page.wait_for_timeout(1500)
    page.screenshot(path="shot_idle.png")
    print("wrote shot_idle.png")
    print("3Dmol loaded:", page.evaluate("!!window.$3Dmol"))

    if JOB:
        # Open the finished job through the UI's own code path.
        page.evaluate(f"watch('{JOB}')")
        page.wait_for_selector("#viewerCard:not(.hidden)", timeout=120_000)
        page.wait_for_function(
            "document.querySelectorAll('#poseTable tbody tr').length > 0",
            timeout=60_000)
        page.wait_for_timeout(6000)          # let WebGL settle
        page.screenshot(path="shot_result.png")
        print("wrote shot_result.png")
        print("canvas present:",
              page.evaluate("!!document.querySelector('#viewer canvas')"))
        print("pose rows:",
              page.evaluate("document.querySelectorAll('#poseTable tbody tr').length"))
        page.locator("#viewerCard").screenshot(path="shot_viewer.png")
        print("wrote shot_viewer.png")

    if errs:
        print("\npage errors:")
        for e in dict.fromkeys(errs):
            print("  !", e[:200])
    else:
        print("\nno page errors")
    b.close()
