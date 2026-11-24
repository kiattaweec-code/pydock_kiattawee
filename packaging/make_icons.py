"""Generate the app icon in every format the three platforms want.

Writes icon.png (Linux), icon.ico (Windows) and icon.icns (macOS) from one
drawing, so the three stay in step.  Pure stdlib plus Pillow, which RDKit
already brings in.
"""
import math
import struct
from pathlib import Path

from PIL import Image, ImageDraw, ImageFont

HERE = Path(__file__).resolve().parent
TEAL = (69, 207, 196)
TEAL_DARK = (28, 125, 120)
AMBER = (233, 161, 75)
BG = (16, 22, 31)


def draw(size: int) -> Image.Image:
    """A ligand in a pocket: three bonded atoms inside a teal box."""
    s = size * 4                       # supersample, then downscale
    img = Image.new("RGBA", (s, s), (0, 0, 0, 0))
    d = ImageDraw.Draw(img)

    r = s * 0.17
    d.rounded_rectangle([0, 0, s - 1, s - 1], radius=r, fill=BG)

    # The docking box, drawn as an open cage in perspective.
    m, off = s * 0.19, s * 0.13
    box = [m, m + off, s - m - off, s - m]
    back = [m + off, m, s - m, s - m - off]
    w = max(2, int(s * 0.016))
    for rect in (back, box):
        d.rectangle(rect, outline=TEAL_DARK, width=w)
    for (x1, y1), (x2, y2) in (
            ((box[0], box[1]), (back[0], back[1])),
            ((box[2], box[1]), (back[2], back[1])),
            ((box[0], box[3]), (back[0], back[3])),
            ((box[2], box[3]), (back[2], back[3]))):
        d.line([x1, y1, x2, y2], fill=TEAL_DARK, width=w)

    # A small three-atom ligand sitting inside it.
    cx, cy, rad = s * 0.5, s * 0.52, s * 0.135
    pts = [(cx + rad * math.cos(a), cy + rad * math.sin(a))
           for a in (-math.pi / 2, math.pi / 6, math.pi * 5 / 6)]
    bond = max(2, int(s * 0.026))
    for i in range(3):
        d.line([*pts[i], *pts[(i + 1) % 3]], fill=TEAL, width=bond)
    for i, (x, y) in enumerate(pts):
        ar = s * (0.072 if i == 0 else 0.058)
        d.ellipse([x - ar, y - ar, x + ar, y + ar],
                  fill=AMBER if i == 0 else TEAL)

    return img.resize((size, size), Image.LANCZOS)


def write_icns(path: Path, sizes):
    """Minimal ICNS: PNG payloads under the standard type codes."""
    codes = {16: b"icp4", 32: b"icp5", 64: b"icp6", 128: b"ic07",
             256: b"ic08", 512: b"ic09", 1024: b"ic10"}
    chunks = b""
    for size in sizes:
        code = codes.get(size)
        if not code:
            continue
        tmp = path.with_suffix(f".{size}.png")
        draw(size).save(tmp, "PNG")
        payload = tmp.read_bytes()
        tmp.unlink()
        chunks += code + struct.pack(">I", len(payload) + 8) + payload
    path.write_bytes(b"icns" + struct.pack(">I", len(chunks) + 8) + chunks)


def dmg_background(w: int = 740, h: int = 460) -> Image.Image:
    """The window art behind the .dmg: title, an arrow, one instruction.

    The icon positions in build.py's AppleScript are matched to this: the app
    at x=180, Applications at x=560, both at y=210.
    """
    s = 2                                   # drawn at 2x for retina
    img = Image.new("RGB", (w * s, h * s), BG)
    d = ImageDraw.Draw(img)

    for i in range(h * s):                  # faint vertical gradient
        t = i / (h * s)
        d.line([(0, i), (w * s, i)],
               fill=(int(BG[0] + 10 * t), int(BG[1] + 12 * t), int(BG[2] + 16 * t)))

    def text(xy, msg, size, fill, anchor="mm"):
        try:
            f = ImageFont.truetype("DejaVuSans.ttf", size * s)
        except OSError:
            f = ImageFont.load_default(size * s)
        d.text((xy[0] * s, xy[1] * s), msg, font=f, fill=fill, anchor=anchor)

    text((w / 2, 52), "pydock", 34, (223, 230, 240))
    text((w / 2, 88), "flexible-ligand molecular docking", 14, (122, 136, 155))

    # Arrow from the app icon to the Applications folder.
    y = 210 * s
    x1, x2 = 268 * s, 472 * s
    d.line([(x1, y), (x2 - 16 * s, y)], fill=TEAL_DARK, width=3 * s)
    d.polygon([(x2, y), (x2 - 18 * s, y - 11 * s), (x2 - 18 * s, y + 11 * s)],
              fill=TEAL)

    text((w / 2, 300), "Drag pydock into Applications", 17, (223, 230, 240))
    text((w / 2, 412), "First time? Open READ ME FIRST.txt — "
                       "macOS needs a right-click → Open, once.",
         13, (233, 161, 75))

    return img.resize((w, h), Image.LANCZOS)


if __name__ == "__main__":
    draw(512).save(HERE / "icon.png", "PNG")
    base = draw(256)
    base.save(HERE / "icon.ico", "ICO",
              sizes=[(16, 16), (24, 24), (32, 32), (48, 48), (64, 64),
                     (128, 128), (256, 256)])
    write_icns(HERE / "icon.icns", [16, 32, 64, 128, 256, 512, 1024])
    dmg_background().save(HERE / "dmg_background.png", "PNG")
    for f in ("icon.png", "icon.ico", "icon.icns", "dmg_background.png"):
        p = HERE / f
        print(f"  {f:12s} {p.stat().st_size / 1024:7.1f} KB")
