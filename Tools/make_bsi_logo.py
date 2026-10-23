"""Builds the extension's bSI logo images from the buildingSMART logo (a PNG on a white background, e.g. the
200x200 bSI.png in the BIM for Bridges folder):

  IfcProjectImporter.ico  IDI_BSI: 32 bpp, transparent, 16 to 128 px. The (R) is left out below 48 px, where it's a smudge.
  IfcExtensions15.bmp     IDB_BSI_15: the 15x15 menu bitmap (the menu check-mark size on Windows 11 at 96 DPI; see
                          MenuLogo.h). Pixel art: each pixel takes the link color covering most of it, since smooth
                          resampling blends the overlapping links into mud at this size.

IfcExtensions.bmp (IDB_BSI, 13x13) is hand-drawn pixel art and isn't generated.

usage: python make_bsi_logo.py <bSI.png> [output folder, default: the extension's folder]
Needs Pillow (python -m pip install pillow).
"""
import os
import sys
from collections import deque
from PIL import Image

# the logo's link colors
LINK_COLORS = {"blue": (0x01, 0x50, 0x9f), "red": (0xe6, 0x01, 0x46), "purple": (0x9f, 0x36, 0x8b), "teal": (0x01, 0x9e, 0xb1)}

ICON_SIZES = [16, 20, 24, 32, 40, 48, 64, 96, 128]
ICON_REGISTERED_MARK_FROM = 48  # smallest icon size with the (R)


def color_to_alpha(img):
    """white -> transparent, keeping the anti-aliased edges (GIMP's color-to-alpha for white)"""
    w, h = img.size
    px = img.load()
    rgba = Image.new("RGBA", (w, h))
    out = rgba.load()
    for y in range(h):
        for x in range(w):
            r, g, b = px[x, y]
            a = max(255 - r, 255 - g, 255 - b)
            if a == 0:
                out[x, y] = (0, 0, 0, 0)
            else:
                f = 255.0 / a
                un = lambda c: max(0, min(255, round(255 - (255 - c) * f)))
                out[x, y] = (un(r), un(g), un(b), a)
    return rgba


def components(rgba):
    """the connected groups of visible pixels"""
    w, h = rgba.size
    px = rgba.load()
    seen = [[False] * w for _ in range(h)]
    found = []
    for y in range(h):
        for x in range(w):
            if seen[y][x] or px[x, y][3] < 16:
                continue
            comp, q = [], deque([(x, y)])
            seen[y][x] = True
            while q:
                cx, cy = q.popleft()
                comp.append((cx, cy))
                for nx, ny in ((cx + 1, cy), (cx - 1, cy), (cx, cy + 1), (cx, cy - 1)):
                    if 0 <= nx < w and 0 <= ny < h and not seen[ny][nx] and px[nx, ny][3] >= 16:
                        seen[ny][nx] = True
                        q.append((nx, ny))
            found.append(comp)
    return found


def without_registered_mark(rgba):
    """the logo without the (R): the pixel groups entirely in the lower right corner are removed, keeping the faint
    anti-aliased fringe of the links"""
    w, h = rgba.size

    def bbox(c):
        return min(p[0] for p in c), min(p[1] for p in c)

    groups = components(rgba)
    mark = [c for c in groups if bbox(c)[0] > 0.75 * w and bbox(c)[1] > 0.75 * h]
    if not mark:
        sys.exit("The (R) wasn't found in the lower right corner of the logo")
    logo = set(p for c in groups if c not in mark for p in c)

    out = rgba.copy()
    o = out.load()
    for y in range(h):
        for x in range(w):
            if (x, y) not in logo and not any((x + dx, y + dy) in logo for dx in range(-2, 3) for dy in range(-2, 3)):
                o[x, y] = (0, 0, 0, 0)
    return out


def square(rgba):
    """cropped to the logo and centered in a transparent square"""
    im = rgba.crop(rgba.getbbox())
    side = max(im.size)
    sq = Image.new("RGBA", (side, side), (0, 0, 0, 0))
    sq.paste(im, ((side - im.size[0]) // 2, (side - im.size[1]) // 2))
    return sq


def resized(rgba, s):
    """premultiplied resampling, so transparent pixels don't darken the edges"""
    return rgba.convert("RGBa").resize((s, s), Image.LANCZOS).convert("RGBA")


def pixel_art(rgba, n, logo, threshold=0.30):
    """n x n on white with the logo at logo x logo pixels, centered; each pixel takes the link color covering most of
    it if that color covers at least threshold of it, else it's white"""
    w, h = rgba.size
    px = rgba.load()
    out = Image.new("RGB", (n, n), (255, 255, 255))
    o = out.load()
    off = (n - logo) // 2
    for j in range(logo):
        for i in range(logo):
            x0, x1 = int(i * w / logo), int((i + 1) * w / logo)
            y0, y1 = int(j * h / logo), int((j + 1) * h / logo)
            cover = dict.fromkeys(LINK_COLORS, 0)
            for y in range(y0, y1):
                for x in range(x0, x1):
                    r, g, b, a = px[x, y]
                    if a >= 128:
                        cover[min(LINK_COLORS, key=lambda k: sum((c - p) ** 2 for c, p in zip((r, g, b), LINK_COLORS[k])))] += 1
            name = max(cover, key=cover.get)
            if cover[name] / ((x1 - x0) * (y1 - y0)) >= threshold:
                o[i + off, j + off] = LINK_COLORS[name]
    return out


def main():
    if len(sys.argv) < 2:
        sys.exit(__doc__)
    out_dir = sys.argv[2] if len(sys.argv) > 2 else os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

    rgba = color_to_alpha(Image.open(sys.argv[1]).convert("RGB"))
    full, small = square(rgba), square(without_registered_mark(rgba))

    frames = [resized(full if s >= ICON_REGISTERED_MARK_FROM else small, s) for s in ICON_SIZES]
    ico = os.path.join(out_dir, "IfcProjectImporter.ico")
    frames[-1].save(ico, format="ICO", sizes=[(s, s) for s in ICON_SIZES], append_images=frames[:-1])
    print("wrote", ico)

    bmp = os.path.join(out_dir, "IfcExtensions15.bmp")
    pixel_art(small, 15, 13).save(bmp, format="BMP")
    print("wrote", bmp)


if __name__ == "__main__":
    main()
