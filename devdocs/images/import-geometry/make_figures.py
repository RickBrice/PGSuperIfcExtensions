"""Writes the figures of devdocs/ImportGeometry.md (SVG, one file per figure, in this folder).

    python make_figures.py

The figures are schematic: proportions are chosen for legibility, not to scale. Each figure has a light background so it
reads the same in light and dark viewers.
"""

import math
from pathlib import Path

HERE = Path(__file__).resolve().parent

# colors
INK = "#222"
MUTED = "#666"
CONCRETE = "#e4e4e4"
CONCRETE_EDGE = "#6b6b6b"
GIRDER = "#cfd8e3"
GIRDER_EDGE = "#4a5a70"
ALIGN = "#1f6fb2"
MEAS = "#c0392b"
GREEN = "#2e8b57"
STEEL = "#8a5a2b"
ORANGE = "#d4801f"

STYLE = f"""
<style>
  text {{ font-family: 'Segoe UI', Helvetica, Arial, sans-serif; font-size: 13px; fill: {INK}; }}
  .muted {{ fill: {MUTED}; }}
  .small {{ font-size: 11.5px; }}
  .title {{ font-weight: 600; font-size: 14px; }}
  .meas {{ fill: {MEAS}; }}
  .align {{ fill: {ALIGN}; }}
  .green {{ fill: {GREEN}; }}
  .steel {{ fill: {STEEL}; }}
  .orange {{ fill: {ORANGE}; }}
  .mid {{ text-anchor: middle; }}
  .end {{ text-anchor: end; }}
</style>
<defs>
  <marker id="arrow" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="7" markerHeight="7" orient="auto-start-reverse">
    <path d="M0,0 L10,5 L0,10 z" fill="{MEAS}"/>
  </marker>
  <marker id="arrow-ink" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="7" markerHeight="7" orient="auto-start-reverse">
    <path d="M0,0 L10,5 L0,10 z" fill="{INK}"/>
  </marker>
  <marker id="arrow-align" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="7" markerHeight="7" orient="auto-start-reverse">
    <path d="M0,0 L10,5 L0,10 z" fill="{ALIGN}"/>
  </marker>
  <marker id="arrow-green" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="7" markerHeight="7" orient="auto-start-reverse">
    <path d="M0,0 L10,5 L0,10 z" fill="{GREEN}"/>
  </marker>
</defs>
"""


class Figure:
    def __init__(self, name, width, height, title):
        self.name, self.width, self.height, self.title = name, width, height, title
        self.items = []

    def add(self, s):
        self.items.append(s)

    def write(self):
        body = "\n".join(self.items)
        svg = (f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {self.width} {self.height}" width="{self.width}" height="{self.height}" role="img" aria-label="{esc(self.title)}">\n'
               f'<title>{esc(self.title)}</title>\n{STYLE}\n'
               f'<rect x="0" y="0" width="{self.width}" height="{self.height}" rx="8" fill="#ffffff" stroke="#d0d0d0"/>\n{body}\n</svg>\n')
        (HERE / f"{self.name}.svg").write_text(svg, encoding="utf-8")

    # drawing helpers
    def text(self, x, y, s, cls="", anchor=None):
        a = f' text-anchor="{anchor}"' if anchor else ""
        c = f' class="{cls}"' if cls else ""
        self.add(f'<text x="{x:.1f}" y="{y:.1f}"{c}{a}>{esc(s)}</text>')

    def lines(self, x, y, rows, cls="", anchor=None, dy=16):
        for i, row in enumerate(rows):
            self.text(x, y + i * dy, row, cls, anchor)

    def line(self, x1, y1, x2, y2, color=INK, width=1.2, dash=None, arrow=None, start_arrow=False):
        d = f' stroke-dasharray="{dash}"' if dash else ""
        m = ""
        if arrow:
            m += f' marker-end="url(#{arrow})"'
            if start_arrow:
                m += f' marker-start="url(#{arrow})"'
        self.add(f'<line x1="{x1:.1f}" y1="{y1:.1f}" x2="{x2:.1f}" y2="{y2:.1f}" stroke="{color}" stroke-width="{width}"{d}{m}/>')

    def poly(self, points, fill="none", stroke=INK, width=1.2, closed=True, dash=None, opacity=None):
        pts = " ".join(f"{x:.1f},{y:.1f}" for x, y in points)
        tag = "polygon" if closed else "polyline"
        d = f' stroke-dasharray="{dash}"' if dash else ""
        o = f' fill-opacity="{opacity}"' if opacity is not None else ""
        self.add(f'<{tag} points="{pts}" fill="{fill}" stroke="{stroke}" stroke-width="{width}"{d}{o} stroke-linejoin="round"/>')

    def path(self, d, fill="none", stroke=INK, width=1.2, dash=None):
        da = f' stroke-dasharray="{dash}"' if dash else ""
        self.add(f'<path d="{d}" fill="{fill}" stroke="{stroke}" stroke-width="{width}"{da} stroke-linejoin="round"/>')

    def circle(self, x, y, r, fill=INK, stroke="none", width=1):
        self.add(f'<circle cx="{x:.1f}" cy="{y:.1f}" r="{r}" fill="{fill}" stroke="{stroke}" stroke-width="{width}"/>')

    def rect(self, x, y, w, h, fill="none", stroke=INK, width=1.2, rx=0, dash=None):
        d = f' stroke-dasharray="{dash}"' if dash else ""
        self.add(f'<rect x="{x:.1f}" y="{y:.1f}" width="{w:.1f}" height="{h:.1f}" rx="{rx}" fill="{fill}" stroke="{stroke}" stroke-width="{width}"{d}/>')

    def dim(self, x1, y1, x2, y2, label=None, color=MEAS, lx=None, ly=None, anchor="middle", cls="meas"):
        """a dimension line with arrows at both ends and a label (at lx, ly, else the middle)"""
        self.line(x1, y1, x2, y2, color=color, width=1.3, arrow="arrow" if color == MEAS else "arrow-ink", start_arrow=True)
        if label:
            self.text(lx if lx is not None else (x1 + x2) / 2, ly if ly is not None else (y1 + y2) / 2 - 5, label, cls, anchor)

    def panel_title(self, x, y, s):
        self.text(x, y, s, "title")


def esc(s):
    return s.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;").replace('"', "&quot;")


def i_section(cx, top, scale, mirror_top_width=48.0):
    """the outline of an I-girder (WF-like), top center at (cx, top), inches * scale, y down on the page"""
    half = [(0, 0), (-mirror_top_width / 2, 0), (-mirror_top_width / 2, 3), (-12, 5), (-3, 8), (-3, 62), (-12.5, 67), (-12.5, 74), (0, 74)]
    left = [(cx + x * scale, top + y * scale) for x, y in half]
    right = [(cx - x * scale, top + y * scale) for x, y in reversed(half)]
    return left + right[1:-1]


def bezier(p0, p1, p2, t):
    x = (1 - t) ** 2 * p0[0] + 2 * (1 - t) * t * p1[0] + t * t * p2[0]
    y = (1 - t) ** 2 * p0[1] + 2 * (1 - t) * t * p1[1] + t * t * p2[1]
    dx = 2 * (1 - t) * (p1[0] - p0[0]) + 2 * t * (p2[0] - p1[0])
    dy = 2 * (1 - t) * (p1[1] - p0[1]) + 2 * t * (p2[1] - p1[1])
    n = math.hypot(dx, dy)
    return (x, y), (dx / n, dy / n)


# ---------------------------------------------------------------------------------------------------------------------
def fig_coordinates():
    f = Figure("coordinates", 760, 360, "Coordinate systems used by the import")

    # (a) plan
    f.panel_title(20, 28, "(a) Plan: station and offset, girder axis")
    p0, p1, p2 = (40, 120), (230, 70), (440, 120)
    pts = [bezier(p0, p1, p2, t / 40)[0] for t in range(41)]
    f.poly(pts, stroke=ALIGN, width=2, closed=False, dash="10 4 2 4")
    f.text(40, 108, "alignment", "align small")
    for t in (0.1, 0.35, 0.6, 0.85):
        (x, y), (dx, dy) = bezier(p0, p1, p2, t)
        f.line(x - dy * 6, y + dx * 6, x + dy * 6, y - dx * 6, ALIGN, 1.2)
    (x, y), (dx, dy) = bezier(p0, p1, p2, 0.6)
    f.text(x - 2, y - 14, "station", "align small", "middle")
    f.line(x + 14 * dx, y + 14 * dy, x + 50 * dx, y + 50 * dy, ALIGN, 1.5, arrow="arrow-align")
    f.text(x + 54 * dx + 2, y + 54 * dy - 4, "ahead", "align small")
    # offset to a point on the right
    nx, ny = -dy, dx  # right of the direction of travel (page y is down)
    px, py = x + 75 * nx, y + 75 * ny
    f.dim(x, y, px, py, None)
    f.text(px + 8, (y + py) / 2 + 4, "offset, + right", "meas small")
    f.circle(px, py, 3.2, MEAS)
    # girder in plan
    gx1, gy1, gx2, gy2 = 70, 205, 410, 230
    ang = math.atan2(gy2 - gy1, gx2 - gx1)
    w = 14
    ox, oy = -math.sin(ang) * w / 2, math.cos(ang) * w / 2
    f.poly([(gx1 - ox, gy1 - oy), (gx2 - ox, gy2 - oy), (gx2 + ox, gy2 + oy), (gx1 + ox, gy1 + oy)], GIRDER, GIRDER_EDGE)
    f.line(gx1, gy1, gx2, gy2, GIRDER_EDGE, 1, dash="8 3 2 3")
    f.circle(gx1, gy1, 3, INK)
    f.circle(gx2, gy2, 3, INK)
    f.text(gx1 - 4, gy1 - 12, "start (lesser station)", "small", "start")
    f.text(gx2, gy2 + 26, "end", "small", "end")
    # girder axes
    ux, uy = math.cos(ang), math.sin(ang)
    bx, by = gx1 + 40 * ux, gy1 + 40 * uy + 30
    f.line(bx, by, bx + 60 * ux, by + 60 * uy, GREEN, 1.6, arrow="arrow-green")
    f.text(bx + 64 * ux, by + 64 * uy + 4, "x (along the girder)", "green small")
    f.line(bx, by, bx - 40 * uy, by + 40 * ux, GREEN, 1.6, arrow="arrow-green")
    f.text(bx - 40 * uy - 6, by + 40 * ux + 16, "z (+ right)", "green small", "middle")
    f.lines(20, 315, ["Global plan coordinates are the model's, in PGSuper's units. A point's station and",
                      "offset are from the alignment, offsets positive to the right looking ahead on station."], "muted small", dy=15)

    # (b) elevation in girder coordinates
    f.panel_title(470, 28, "(b) Girder elevation")
    ex, ey = 490, 70
    f.rect(ex, ey, 230, 52, GIRDER, GIRDER_EDGE)
    f.line(ex, ey, ex + 245, ey, GREEN, 1.6, arrow="arrow-green")
    f.text(ex + 236, ey - 6, "x", "green small")
    f.line(ex, ey, ex, ey + 70, GREEN, 1.6, arrow="arrow-green")
    f.text(ex + 6, ey + 76, "y (depth below the top)", "green small")
    f.text(ex - 6, ey - 6, "top of girder at its start", "small", "start")

    # (c) section coordinates
    f.panel_title(470, 170, "(c) Section")
    cx, top, s = 600, 200, 1.35
    f.poly(i_section(cx, top, s), GIRDER, GIRDER_EDGE)
    f.line(cx, top, cx + 70, top, GREEN, 1.6, arrow="arrow-green")
    f.text(cx + 74, top + 4, "x", "green small")
    f.line(cx, top, cx, top - 25, GREEN, 1.6, arrow="arrow-green")
    f.text(cx + 6, top - 16, "y (up)", "green small")
    f.circle(cx, top, 3, GREEN)
    f.lines(470, 315, ["Sections are cut normal to the girder axis, origin at the",
                       "top center, x horizontal, y up."], "muted small", dy=15)
    f.write()


# ---------------------------------------------------------------------------------------------------------------------
def fig_mesh_steps():
    f = Figure("mesh-steps", 760, 300, "From an element's solid to its edges")
    titles = ["1. Triangulated solid", "2. Smooth components", "3. Boundary of the top", "4. Smooth segments"]
    for i, t in enumerate(titles):
        f.text(20 + i * 185, 28, t, "title")

    def box(ox, oy, fill_top, fill_side, edge, show_tris):
        # an oblique box: top face and front face
        L, W, H, sk = 140, 34, 34, 22
        top = [(ox, oy + sk), (ox + L, oy), (ox + L + 18, oy + W - 6), (ox + 18, oy + sk + W - 6)]
        front = [top[3], top[2], (top[2][0], top[2][1] + H), (top[3][0], top[3][1] + H)]
        end = [top[0], top[3], front[3], (top[0][0], top[0][1] + H)]
        f.poly(end, fill_side, edge)
        f.poly(front, fill_side, edge)
        f.poly(top, fill_top, edge)
        if show_tris:
            for k in range(1, 5):
                a = (top[0][0] + (top[1][0] - top[0][0]) * k / 5, top[0][1] + (top[1][1] - top[0][1]) * k / 5)
                b = (top[3][0] + (top[2][0] - top[3][0]) * k / 5, top[3][1] + (top[2][1] - top[3][1]) * k / 5)
                b0 = (top[3][0] + (top[2][0] - top[3][0]) * (k - 1) / 5, top[3][1] + (top[2][1] - top[3][1]) * (k - 1) / 5)
                f.line(*a, *b, edge, 0.7)
                f.line(*a, *b0, edge, 0.7)
            f.line(*top[1], *(top[2]), edge, 0.7)
            for k in range(1, 5):
                a = (front[0][0] + (front[1][0] - front[0][0]) * k / 5, front[0][1] + (front[1][1] - front[0][1]) * k / 5)
                f.line(a[0], a[1], a[0] - 28, a[1] + H, edge, 0.7)
        return top

    box(25, 70, GIRDER, GIRDER, GIRDER_EDGE, True)
    t2 = box(210, 70, "#9fc3e6", "#eeeeee", GIRDER_EDGE, False)
    f.text(285, 180, "top: highest when each", "small muted", "middle")
    f.text(285, 195, "component is pushed out", "small muted", "middle")
    f.text(285, 210, "along its average normal", "small muted", "middle")
    # 3. boundary wire of the top
    ox = 395
    L, W, sk = 140, 34, 22
    top = [(ox, 70 + sk), (ox + L, 70), (ox + L + 18, 70 + W - 6), (ox + 18, 70 + sk + W - 6)]
    f.poly(top, "#eef4fb", ALIGN, 2)
    for p in top:
        f.circle(*p, 2.5, ALIGN)
    f.text(470, 180, "the outline of the top face,", "small muted", "middle")
    f.text(470, 195, "a chain of edges", "small muted", "middle")
    # 4. segments
    ox = 580
    top = [(ox, 70 + sk), (ox + L, 70), (ox + L + 18, 70 + W - 6), (ox + 18, 70 + sk + W - 6)]
    f.line(*top[0], *top[1], GREEN, 3)
    f.line(*top[2], *top[3], GREEN, 3)
    f.line(*top[1], *top[2], MEAS, 3)
    f.line(*top[3], *top[0], MEAS, 3)
    for a, b in ((top[1], top[2]), (top[3], top[0])):
        f.circle((a[0] + b[0]) / 2, (a[1] + b[1]) / 2, 3.5, MEAS)
    f.text(655, 180, "shortest two: the ends", "small meas", "middle")
    f.text(655, 195, "(their midpoints are on the CL)", "small meas", "middle")
    f.text(655, 212, "longest two: the edges", "small green", "middle")
    f.lines(20, 250, ["The same steps find a girder's ends (segments split where the direction turns more than 1 degree) and the",
                      "deck's edges (10 degrees). Faces belong to one smooth component while adjacent faces turn less than 10 degrees."], "muted small", dy=15)
    f.write()


# ---------------------------------------------------------------------------------------------------------------------
def fig_pipeline():
    steps = [
        ("Alignment", "the model's IfcAlignment; a straight one through the deck ends if it has none"),
        ("Piers and abutments", "the number of spans, and the pier stations"),
        ("Girders", "ends and CL from each girder's top surface, span, order; bearings at the ends"),
        ("Pier orientation", "positioning referent, else the CL bearing line, else the girder ends"),
        ("Girder spacing", "where the girder CLs cross the pier lines"),
        ("Girder sections", "cut from the solids; a library entry matched, or a new one fitted"),
        ("End distance, bearings", "end distance and bearing offset from the CL bearings; bearing size and count"),
        ("Boundary conditions", "pier diaphragms, then support conditions"),
        ("Strands", "centerlines to the library's strand grid, else individual strands"),
        ("Girder bars", "longitudinal rows, stirrup zones, confinement bars"),
        ("Deck", "plan edges, profile, roadway section, traffic barriers, depths, haunch,"),
    ]
    f = Figure("pipeline", 760, 545, "Order of the import steps")
    f.text(20, 28, "Order of the steps: each one uses what the steps before it found", "title")
    y = 48
    for i, (name, what) in enumerate(steps):
        f.rect(30, y, 190, 30, "#eef4fb", ALIGN, 1.2, 5)
        f.text(125, y + 20, name, "", "middle")
        f.text(236, y + 20, what, "small")
        if i + 1 < len(steps):
            f.line(125, y + 30, 125, y + 41, ALIGN, 1.3, arrow="arrow-align")
        y += 42
    f.text(236, y - 42 + 36, "slab offset, deck reinforcement", "small")
    f.lines(20, 530, ["The deck is last: its cross section is measured over the girders, and the slab offset at the CL bearings."], "muted small")
    f.write()


# ---------------------------------------------------------------------------------------------------------------------
def fig_girder_layout():
    f = Figure("girder-layout", 760, 330, "Girder layout from the girder solids")
    f.text(20, 28, "Plan: each girder's ends and CL, its span, and its place across the span", "title")
    # alignment
    f.line(30, 70, 730, 70, ALIGN, 2, dash="10 4 2 4")
    f.text(735, 74, "", "align small")
    f.text(30, 62, "alignment", "align small")
    piers = [80, 380, 690]
    for i, px in enumerate(piers):
        f.line(px, 50, px, 290, INK, 1.2, dash="8 3 2 3")
        f.text(px, 305, ["Abutment 1", "Pier 2", "Abutment 3"][i], "small", "middle")
    # span 1 girders
    offsets = [110, 150, 190, 230]
    for k, y in enumerate(offsets):
        f.rect(88, y - 6, 284, 12, GIRDER, GIRDER_EDGE)
        f.line(88, y, 372, y, GIRDER_EDGE, 0.8, dash="6 3")
        f.circle(88, y, 2.8, MEAS)
        f.circle(372, y, 2.8, MEAS)
        f.text(70, y + 4, f"G{k + 1}", "small", "end")
    # span 2, one girder highlighted
    for k, y in enumerate(offsets):
        f.rect(388, y - 6, 294, 12, GIRDER, GIRDER_EDGE)
        f.circle(388, y, 2.8, MEAS)
        f.circle(682, y, 2.8, MEAS)
    f.circle(535, 150, 4, GREEN)
    f.line(535, 150, 560, 268, GREEN, 1, dash="3 3")
    f.lines(565, 268, ["midpoint: in span 2"], "green small")
    f.lines(240, 268, ["red: CL girder at each end,", "the midpoint of an end segment"], "meas small")
    f.lines(30, 322 - 4, [""], "small")
    # offset ordering arrow
    f.line(720, 100, 720, 240, MEAS, 1.3, arrow="arrow")
    f.text(730, 58, "girders ordered left to right by offset", "meas small", "end")
    f.write()


# ---------------------------------------------------------------------------------------------------------------------
def fig_spacing():
    f = Figure("spacing", 760, 340, "Girder spacing at a pier line")
    f.text(20, 28, "Plan at a pier: where each girder CL crosses the measurement line", "title")
    ax, cy = 380, 150  # the alignment runs down the page (ahead is down); the pier station is at (ax, cy)
    top, bottom = 50, 250
    f.line(ax, top, ax, bottom, ALIGN, 2, dash="10 4 2 4")
    f.line(ax + 14, bottom - 40, ax + 14, bottom - 10, ALIGN, 1.4, arrow="arrow-align")
    f.text(ax + 20, bottom - 18, "ahead", "align small")
    f.text(ax + 6, top + 8, "alignment", "align small")
    f.circle(ax, cy, 4, ALIGN)
    f.text(ax + 8, cy + 18, "pier station", "align small")
    # normal line
    f.line(60, cy, 700, cy, MEAS, 1.4)
    f.text(64, cy - 8, "normal to the alignment", "meas small")
    # skewed CL pier line
    sk = math.radians(15)
    dx, dy = math.cos(sk), math.sin(sk)
    f.line(ax - 330 * dx, cy - 330 * dy, ax + 330 * dx, cy + 330 * dy, GREEN, 1.4, dash="7 4")
    f.text(ax + 330 * dx, cy + 330 * dy + 16, "along the CL pier", "green small", "end")
    # girder CLs, nearly parallel to the alignment
    xs = [-240, -80, 80, 240]
    crossings = []
    for gx in xs:
        x = ax + gx
        f.line(x - 5, top, x + 5, bottom, GIRDER_EDGE, 1.2, dash="6 3")
        t = (cy - top) / (bottom - top)
        xn = x - 5 + 10 * t
        crossings.append(xn)
        f.circle(xn, cy, 3.5, MEAS)
        # crossing with the skewed line: x - 5 + 10 t = ax + u dx ; top + (bottom - top) t = cy + u dy
        a11, a12, b1 = 10, -dx, ax - (x - 5)
        a21, a22, b2 = bottom - top, -dy, cy - top
        det = a11 * a22 - a12 * a21
        t2 = (b1 * a22 - a12 * b2) / det
        f.circle(x - 5 + 10 * t2, top + (bottom - top) * t2, 3.5, GREEN)
    for i in range(3):
        f.dim(crossings[i], cy - 28, crossings[i + 1], cy - 28, None)
        f.text((crossings[i] + crossings[i + 1]) / 2 + (18 if i == 1 else 0), cy - 34, f"s{i + 1}", "meas")
    f.lines(20, 295, ["Spacing is the difference of the crossing offsets. It's uniform when every spacing at every pier line agrees within",
                      "1/16 in, tried normal to the alignment first, then along the CL piers; otherwise it's set at each pier face."], "muted small", dy=15)
    f.write()


# ---------------------------------------------------------------------------------------------------------------------
def fig_pier_orientation():
    f = Figure("pier-orientation", 760, 360, "CL pier direction")
    f.text(20, 28, "Plan at a pier: the CL pier direction, from the first source the model has", "title")
    ax, cy = 380, 175
    f.line(ax, 50, ax, 300, ALIGN, 2, dash="10 4 2 4")
    f.text(ax + 6, 62, "alignment", "align small")
    sk = math.radians(-12)
    dx, dy = math.cos(sk), math.sin(sk)
    for gx in (-240, -80, 80, 240):
        u = gx / dx
        px, py = ax + u * dx, cy + u * dy
        f.rect(px - 7, py - 70, 14, 58, GIRDER, GIRDER_EDGE)   # back span girder, its end at py - 12
        f.rect(px - 7, py + 12, 14, 58, GIRDER, GIRDER_EDGE)   # ahead span girder
        f.rect(px - 9, py - 28, 18, 9, "#f2d9c4", STEEL)
        f.rect(px - 9, py + 19, 18, 9, "#f2d9c4", STEEL)
        f.circle(px, py - 23.5, 2.6, MEAS)
        f.circle(px, py + 23.5, 2.6, MEAS)
        f.circle(px, py - 12, 2.6, GREEN)
        f.circle(px, py + 12, 2.6, GREEN)
    L = 320
    f.line(ax - L * dx, cy - L * dy - 23.5, ax + L * dx, cy + L * dy - 23.5, MEAS, 1.2, dash="5 3")
    f.line(ax - L * dx, cy - L * dy - 12, ax + L * dx, cy + L * dy - 12, GREEN, 1.2, dash="2 3")
    f.line(ax, cy, ax + 110 * dx, cy + 110 * dy, INK, 2, arrow="arrow-ink")
    f.circle(ax, cy, 3.5, INK)
    f.lines(482, 232, ["1. the positioning referent's", "explicit RefDirection (black)"], "small", dy=14)
    f.line(486, 220, ax + 60 * dx, cy + 60 * dy + 4, INK, 0.8)
    f.text(30, 70, "2. CL bearing line: best-fit line through the CL bearings (red)", "meas small")
    f.text(30, 86, "3. line through the girder ends (green)", "green small")
    f.lines(20, 315, ["The first available source sets the orientation, as a skew from the alignment normal; the others are checked",
                      "against it and a difference over 0.5 degree is logged. With none, the pier is normal to the alignment."], "muted small", dy=15)
    f.write()


# ---------------------------------------------------------------------------------------------------------------------
def fig_bearings():
    f = Figure("bearings", 760, 360, "Bearings, end distance, and bearing offset")
    f.text(20, 28, "(a) Plan at a pier face", "title")
    px, gy = 90, 150
    f.line(px, 60, px, 250, INK, 1.4, dash="10 4 2 4")
    f.text(px, 268, "CL pier", "small", "middle")
    x_end, x_cl = 150, 190
    f.rect(x_end, gy - 26, 300, 52, GIRDER, GIRDER_EDGE)
    f.line(x_end, gy, 460, gy, GIRDER_EDGE, 0.8, dash="6 3")
    f.text(462, gy + 4, "CL girder", "small")
    for by in (gy - 15, gy + 15):
        f.rect(x_cl - 14, by - 7, 28, 14, "#f2d9c4", STEEL)
        f.circle(x_cl, by, 2.4, STEEL)
    f.circle(x_cl, gy, 4, MEAS)
    f.circle(x_end, gy, 3, INK)
    # bearing offset: CL pier to CL bearing
    f.line(px, gy - 60, px, gy - 26, MEAS, 0.8, dash="2 2")
    f.line(x_cl, gy - 60, x_cl, gy - 22, MEAS, 0.8, dash="2 2")
    f.dim(px, gy - 52, x_cl, gy - 52, None)
    f.text((px + x_cl) / 2, gy - 60, "bearing offset", "meas small", "middle")
    # end distance: end of girder to CL bearing
    f.line(x_end, gy + 26, x_end, gy + 60, MEAS, 0.8, dash="2 2")
    f.line(x_cl, gy + 22, x_cl, gy + 60, MEAS, 0.8, dash="2 2")
    f.dim(x_end, gy + 52, x_cl, gy + 52, None)
    f.text(x_cl + 8, gy + 56, "end distance", "meas small")
    f.lines(250, gy - 52, ["CL bearing (red): on the CL girder, at the bearings'", "average distance from the girder end"], "meas small", dy=14)
    f.lines(20, 296, ["Both are measured normal to the CL pier; each pier face gets the median",
                      "of its girders (a spread over 1/4 in is logged). A bearing belongs to the",
                      "nearest girder end it's under: up to 6 ft in from the end, within 3 ft of the CL."], "muted small", dy=15)
    # (b) bearing sizes
    bx = 540
    f.text(bx, 28, "(b) Bearing size", "title")
    f.text(bx, 60, "rectangular: length along the girder,", "small")
    f.text(bx, 75, "width across; height from the z range", "small")
    f.rect(bx + 30, 92, 90, 44, "#f2d9c4", STEEL)
    f.dim(bx + 30, 150, bx + 120, 150, None)
    f.text(bx + 75, 166, "length", "meas small", "middle")
    f.dim(bx + 134, 92, bx + 134, 136, None)
    f.text(bx + 142, 118, "width", "meas small")
    f.text(bx, 198, "round: 8 or more distinct plan points", "small")
    f.text(bx, 213, "within 2% of a common radius", "small")
    f.add(f'<circle cx="{bx + 75}" cy="256" r="28" fill="#f2d9c4" stroke="{STEEL}" stroke-width="1.2"/>')
    f.dim(bx + 47, 256, bx + 103, 256, None)
    f.text(bx + 75, 300, "diameter", "meas small", "middle")
    f.write()


# ---------------------------------------------------------------------------------------------------------------------
def fig_continuity():
    f = Figure("continuity", 760, 250, "Boundary condition from a pier diaphragm")
    f.text(20, 28, "Elevation at a pier: continuous when a diaphragm engages the girder ends of both spans", "title")
    for ox, cont in ((40, True), (420, False)):
        # deck
        f.rect(ox, 60, 320, 16, CONCRETE, CONCRETE_EDGE)
        f.rect(ox, 80, 140, 50, GIRDER, GIRDER_EDGE)
        f.rect(ox + 180, 80, 140, 50, GIRDER, GIRDER_EDGE)
        f.rect(ox + 120, 140, 80, 50, CONCRETE, CONCRETE_EDGE)
        f.line(ox + 160, 50, ox + 160, 205, INK, 1, dash="8 3 2 3")
        f.rect(ox + 125, 130, 16, 10, "#f2d9c4", STEEL)
        f.rect(ox + 179, 130, 16, 10, "#f2d9c4", STEEL)
        if cont:
            f.rect(ox + 128, 76, 64, 64, "#b9dcc6", GREEN, 1.6)
            f.text(ox + 160, 225, "diaphragm spans both ends (within 1 in):", "green small", "middle")
            f.text(ox + 160, 240, "continuous after deck", "green small", "middle")
        else:
            f.text(ox + 160, 225, "no such diaphragm: the bearings' support", "small", "middle")
            f.text(ox + 160, 240, "condition (Fixed: hinge, Expansion: roller)", "small", "middle")
    f.write()


# ---------------------------------------------------------------------------------------------------------------------
def fig_girder_section():
    f = Figure("girder-section", 760, 380, "Girder section cut and measures")
    f.text(20, 28, "(a) Cut at mid-length, normal to the girder axis", "title")
    # girder elevation with a cut plane
    f.rect(30, 70, 300, 60, GIRDER, GIRDER_EDGE)
    f.line(180, 50, 180, 150, MEAS, 2)
    f.text(186, 162, "cut at 0.5 of the length", "meas small")
    f.text(30, 178, "(away from end blocks and tapers)", "muted small")
    # (b) rotated section then upright
    f.text(20, 215, "(b) Turned upright", "title")
    s = 0.95
    sec = i_section(0, 0, s)
    ang = math.radians(-6)
    rot = [(100 + x * math.cos(ang) - y * math.sin(ang), 240 + x * math.sin(ang) + y * math.cos(ang)) for x, y in sec]
    f.poly(rot, GIRDER, GIRDER_EDGE)
    f.text(100, 330 + 12, "rotated with the deck", "small", "middle")
    f.line(170, 280, 200, 280, INK, 1.4, arrow="arrow-ink")
    up = [(260 + x, 240 + y) for x, y in sec]
    f.poly(up, GIRDER, GIRDER_EDGE)
    f.text(260, 330 + 12, "least height: upright", "small", "middle")
    f.lines(20, 365, ["Turned within 15 degrees, in 0.05 degree steps; the largest loop is the outside, loops inside it are voids."], "muted small")
    # (c) measures
    f.text(400, 28, "(c) What's measured", "title")
    cx, top, sc = 520, 60, 3.1
    sec = i_section(cx, top, sc)
    f.poly(sec, GIRDER, GIRDER_EDGE)
    # width profile lines
    for k in range(0, 75, 6):
        y = top + k * sc
        # find the width at depth k from the half outline
        half = [(0, 0), (-24, 0), (-24, 3), (-12, 5), (-3, 8), (-3, 62), (-12.5, 67), (-12.5, 74), (0, 74)]
        w = None
        for (x1, y1), (x2, y2) in zip(half, half[1:]):
            if y1 <= k <= y2 and y1 != y2:
                w = -(x1 + (x2 - x1) * (k - y1) / (y2 - y1))
                break
        if w is None:
            continue
        f.line(cx - w * sc, y, cx + w * sc, y, MEAS, 0.9)
    f.dim(cx + 95, top, cx + 95, top + 74 * sc, "height", lx=cx + 101, ly=top + 37 * sc, anchor="start")
    f.circle(cx, top + 30 * sc, 4, GREEN)
    f.text(cx + 10, top + 30 * sc + 4, "centroid, I", "green small")
    f.lines(400, 315, ["height, area, centroid, I about the centroid, width;",
                       "the width at every 1/8 in of depth (red lines, drawn coarser)",
                       "and the area above every 1/2 in of depth"], "small", dy=15)
    f.write()


# ---------------------------------------------------------------------------------------------------------------------
def fig_section_fit():
    f = Figure("section-fit", 760, 420, "Choosing or fitting a girder library entry")
    f.text(20, 28, "(a) Half outlines: the model's (grey) and a beam factory's fitted section (blue)", "title")
    ox, oy, s = 290, 60, 3.4
    model = [(0, 0), (-24, 0), (-24, 3), (-12, 5), (-3, 8), (-3, 62), (-12.5, 67), (-12.5, 74), (0, 74)]
    fit = [(0, 0), (-23, 0), (-23, 3.6), (-11, 5.5), (-3.5, 9), (-3.5, 61), (-13.5, 67.5), (-13.5, 74), (0, 74)]
    f.poly([(ox + x * s, oy + y * s) for x, y in model] + [(ox, oy + 74 * s)], CONCRETE, CONCRETE_EDGE, closed=True)
    f.poly([(ox + x * s, oy + y * s) for x, y in fit], "none", ALIGN, 1.8, closed=False)
    for (x1, y1), (x2, y2) in zip(model[1:-1], fit[1:-1]):
        f.line(ox + x1 * s, oy + y1 * s, ox + x2 * s, oy + y2 * s, MEAS, 1)
        f.circle(ox + x1 * s, oy + y1 * s, 2.4, CONCRETE_EDGE)
        f.circle(ox + x2 * s, oy + y2 * s, 2.4, ALIGN)
    f.line(ox, oy - 10, ox, oy + 74 * s + 10, INK, 1, dash="8 3 2 3")
    f.text(ox + 6, oy + 74 * s + 4, "axis of symmetry", "small")
    f.lines(20, 380, ["The left halves only. Corners (red pairs) are paired in order, then the distance between the outlines (points spaced",
                      "along each outline to the nearest edge of the other) is minimized by least squares."], "muted small", dy=15)
    # (b) decision flow
    bx = 400
    f.text(bx, 58, "(b) The entry for a girder", "title")
    steps = [("Same section as a girder already done?", "use its entry"),
             ("A library entry's section matches?", "use it (a named entry is preferred)"),
             ("Fit every beam factory's dimensions", "best accepted fit, constant depth"),
             ("Make a new entry \"<type> (IFC)\"", "fitted dimensions, a 2 in strand grid")]
    y = 70
    for i, (q, a) in enumerate(steps):
        f.rect(bx, y, 330, 46, "#eef4fb", ALIGN, 1.2, 5)
        f.text(bx + 12, y + 19, q)
        f.text(bx + 12, y + 37, a, "small green")
        if i + 1 < len(steps):
            f.line(bx + 165, y + 46, bx + 165, y + 64, ALIGN, 1.3, arrow="arrow-align")
            f.text(bx + 172, y + 59, "no", "small muted")
        y += 64
    f.lines(bx, 330, ["Matching tolerances: height, centroid, and the width at", "each depth 1/8 in; area and I 1%."], "small", dy=15)
    f.write()


# ---------------------------------------------------------------------------------------------------------------------
def fig_strands():
    f = Figure("strands", 760, 330, "Strands in girder coordinates")
    f.text(20, 28, "Girder elevation: strand centerlines read from their solids, in girder coordinates", "title")
    x0, x1, top, depth = 60, 700, 60, 150
    f.rect(x0, top, x1 - x0, depth, GIRDER, GIRDER_EDGE)
    # temporary strand
    f.line(x0 + 10, top + 18, x1 - 10, top + 18, ORANGE, 2)
    f.text(x1 - 12, top + 14, "temporary (tendon status TEMPORARY)", "orange small", "end")
    # harped strand
    h1, h2 = x0 + 0.4 * (x1 - x0), x0 + 0.6 * (x1 - x0)
    f.poly([(x0, top + 40), (h1, top + 125), (h2, top + 125), (x1, top + 40)], "none", MEAS, 2, closed=False)
    f.circle(h1, top + 125, 3.5, MEAS)
    f.circle(h2, top + 125, 3.5, MEAS)
    f.text((h1 + h2) / 2, top + 118, "harped: its centerline bends", "meas small", "middle")
    f.text(h1, top + depth + 18, "harp point", "meas small", "middle")
    f.text(h2, top + depth + 18, "harp point", "meas small", "middle")
    # straight strands
    for y in (top + 135, top + 142):
        f.line(x0, y, x1, y, STEEL, 2)
    # debonded at start: straight strand starting later? draw a sheath
    f.line(x0, top + 128, x1, top + 128, STEEL, 2)
    f.rect(x0 + 2, top + 125, 70, 6, "none", INK, 1, dash="3 2")
    f.text(x0 + 80, top + 104, "debonded length (from properties)", "small", "start")
    f.line(x0 + 78, top + 108, x0 + 50, top + 124, INK, 0.8)
    # extended strand
    f.line(x1, top + 142, x1 + 30, top + 142, STEEL, 2)
    f.text(x1 + 34, top + 140, "", "small")
    f.lines(560, 248, ["projects beyond the end:", "extended there"], "small steel", dy=14)
    f.line(600, 250 - 18, x1 + 20, top + 146, STEEL, 0.8)
    f.text(x0 + 140, top + 162 + 10, "straight", "steel small")
    # axes
    f.line(x0, top, x0 + 60, top, GREEN, 1.6, arrow="arrow-green")
    f.text(x0 + 4, top - 6, "x", "green small")
    f.line(x0, top, x0, top + 40, GREEN, 1.6, arrow="arrow-green")
    f.text(x0 - 6, top + 36, "y", "green small", "end")
    f.lines(20, 285, ["Every strand on a grid position of the library entry (and harp points as the entry gives them): the grid is filled.",
                      "Otherwise the strands are input individually, with the harped strands' common harp points."], "muted small", dy=15)
    f.write()


# ---------------------------------------------------------------------------------------------------------------------
def fig_stirrups():
    f = Figure("stirrups", 760, 350, "Stirrup stations to shear zones")
    f.text(20, 28, "Girder elevation: stations of vertical legs, then zones of equal spacing", "title")
    x0, x1, top, depth = 40, 720, 80, 110
    f.rect(x0, top - 18, x1 - x0, 18, CONCRETE, CONCRETE_EDGE)  # deck
    f.rect(x0, top, x1 - x0, depth, GIRDER, GIRDER_EDGE)
    mid = (x0 + x1) / 2
    # stations: zone 1 spacing 10 px x 6, zone 2 spacing 22 px to mid, symmetric
    stations = []
    x = x0 + 6
    stations.append(x)
    for _ in range(6):
        x += 10
        stations.append(x)
    while x + 22 < mid:
        x += 22
        stations.append(x)
    allst = stations + [x0 + x1 - s for s in reversed(stations)]
    for s in allst:
        f.line(s, top - 12, s, top + depth - 8, STEEL, 1.6)
    # end tie at the start that doesn't engage the deck
    f.line(x0 + 2, top + 10, x0 + 2, top + depth - 8, MUTED, 1.6, dash="3 2")
    f.text(x0 + 30, top + depth + 58, "dashed: an end tie within the girder depth of an end that doesn't reach the deck, left out", "muted small")
    f.line(x0 + 2, top + depth - 4, x0 + 28, top + depth + 52, MUTED, 0.8)
    # zones
    z1 = (stations[0], stations[6])
    f.dim(x0, top + depth + 14, stations[6], top + depth + 14, None)
    f.text((x0 + stations[6]) / 2, top + depth + 30, "zone 1", "meas small", "middle")
    f.dim(stations[6], top + depth + 14, mid, top + depth + 14, None)
    f.text((stations[6] + mid) / 2, top + depth + 30, "zone 2, to midspan", "meas small", "middle")
    f.line(mid, top - 30, mid, top + depth + 20, INK, 1, dash="8 3 2 3")
    f.text(mid, top - 34, "midspan: the zones are symmetric", "small", "middle")
    f.lines(500, top - 40, ["legs reaching more than 1 in above the", "girder top engage the deck"], "small steel", dy=14)
    f.lines(20, 305, ["A leg is a straight piece within 15 degrees of vertical, at least half the girder depth. Bars within 1 in along the",
                      "girder are one station. Equal spaces (1/4 in) with the same bars make a zone; spacings round to 1/8 in (1 mm in SI)."], "muted small", dy=15)
    f.write()


# ---------------------------------------------------------------------------------------------------------------------
def fig_confinement():
    f = Figure("confinement", 760, 360, "Confinement bars from the bottom flange bars")
    f.text(20, 28, "Section: a bar is a confinement bar when it encloses the bottom row of strands", "title")
    cx, top, s = 230, 50, 3.6
    sec = i_section(cx, top, s)
    f.poly(sec, GIRDER, GIRDER_EDGE)
    # 40% line
    y40 = top + 0.6 * 74 * s
    f.line(cx - 120, y40, cx + 120, y40, INK, 1, dash="6 3")
    f.text(cx + 124, y40 + 4, "40% of the depth", "small")
    # bottom strands
    for k in range(-5, 6):
        f.circle(cx + k * 2 * s, top + 70 * s, 3.0, STEEL)
    for k in range(-4, 5):
        f.circle(cx + k * 2 * s, top + 66 * s, 3.0, STEEL)
    # confinement bar around the bottom flange
    b = [(-11, 64.5), (-11, 72.5), (11, 72.5), (11, 64.5)]
    f.poly([(cx + x * s, top + y * s) for x, y in b], "none", GREEN, 2.2)
    f.text(cx + 140, top + 70 * s + 4, "encloses the bottom row: confinement", "green small")
    f.line(cx + 11 * s, top + 70 * s, cx + 136, top + 70 * s, GREEN, 0.8)
    # a bar higher up that doesn't enclose
    f.line(cx - 2.2 * s, top + 40 * s, cx + 2.2 * s, top + 40 * s, MUTED, 2.2)
    f.lines(cx - 40, top + 40 * s - 4, ["above the 40% line:", "not confinement"], "muted small", "end", dy=14)
    f.line(cx - 3 * s, top + 40 * s, cx - 36, top + 40 * s, MUTED, 0.8)
    # stirrup legs
    f.line(cx - 1.2 * s, top + 2 * s, cx - 1.2 * s, top + 71 * s, STEEL, 1.4)
    f.line(cx + 1.2 * s, top + 2 * s, cx + 1.2 * s, top + 71 * s, STEEL, 1.4)
    f.lines(440, 80, ["The bar's path in the section, its ends joined,", "must enclose the strands within 1 in of the",
                      "lowest (with no strands: the CL girder 1 in", "above the bar's lowest point).", "",
                      "It's at the stirrup station nearest it (within", "3 in and half the spacing), whose zone gets its",
                      "size as the confinement bar size."], "small", dy=16)
    f.write()


# ---------------------------------------------------------------------------------------------------------------------
def fig_longitudinal_bars():
    f = Figure("longitudinal-bars", 760, 330, "Longitudinal girder bars")
    f.text(20, 28, "(a) Section: face and cover", "title")
    cx, top, s = 160, 50, 3.0
    f.poly(i_section(cx, top, s), GIRDER, GIRDER_EDGE)
    for x in (-16, -8, 8, 16):
        f.circle(cx + x * s, top + 2 * s, 3.4, STEEL)
    for x in (-8, 8):
        f.circle(cx + x * s, top + 71 * s, 3.4, STEEL)
    for y1, y2, x, label in ((top, top + 2 * s, cx + 28 * s, "cover from the top"), (top + 71 * s, top + 74 * s, cx + 18 * s, "cover from the bottom")):
        f.line(x - 16, y1, x + 4, y1, MEAS, 1)
        f.line(x - 16, y2, x + 4, y2, MEAS, 1)
        f.line(x, y1 - 12, x, y1, MEAS, 1.2, arrow="arrow")
        f.line(x, y2 + 12, x, y2, MEAS, 1.2, arrow="arrow")
        f.text(x + 8, (y1 + y2) / 2 + 4, label, "meas small")
    f.lines(20, 300, ["Bars pair at +/-z about the CL girder (a", "bar without a mirror bar is skipped)."], "muted small", dy=15)
    # (b) elevation layouts
    f.text(370, 28, "(b) Elevation: layout from the bar ends", "title")
    x0, x1, top2 = 380, 715, 60
    f.rect(x0, top2, x1 - x0, 90, GIRDER, GIRDER_EDGE)
    f.line(x0 + 4, top2 + 10, x1 - 4, top2 + 10, STEEL, 2.4)
    f.text(x1 - 6, top2 + 26, "full length (ends within 3 in of the end faces)", "small steel", "end")
    f.line(x0 + 4, top2 + 50, x0 + 150, top2 + 50, STEEL, 2.4)
    f.text(x0 + 156, top2 + 54, "from the left end", "small steel")
    f.line((x0 + x1) / 2 - 80, top2 + 78, (x0 + x1) / 2 + 80, top2 + 78, STEEL, 2.4)
    f.text((x0 + x1) / 2, top2 + 72, "centered", "small steel", "middle")
    f.line(x1 - 4, top2 + 10, x1 + 20, top2 + 10, STEEL, 2.4)
    f.lines(370, 190, ["A longitudinal bar's longest straight piece runs within 2 degrees",
                       "of the girder axis and is at least 2 ft long. Bars of one size that",
                       "lap (collinear within 1.5 bar diameters) are one bar. A bar that",
                       "projects beyond an end (right) is extended there. Bars with the",
                       "same face, cover, size, and layout make a row."], "small", dy=16)
    f.write()


# ---------------------------------------------------------------------------------------------------------------------
def fig_deck_plan():
    f = Figure("deck-plan", 760, 330, "Deck edges and deck points")
    f.text(20, 28, "Plan: the deck's long edges become deck points (station, left and right offsets)", "title")
    left = [(60, 90), (300, 90), (480, 70), (700, 70)]
    right = [(60, 230), (700, 230)]
    f.poly(left + list(reversed(right)), "#f3f3f3", CONCRETE_EDGE, 1.2)
    f.line(30, 160, 730, 160, ALIGN, 2, dash="10 4 2 4")
    f.text(70, 152, "alignment", "align small")
    f.text(684, 218, "the deck's two longest edges (green); its ends (red)", "small", "end")
    f.poly(left, "none", GREEN, 3, closed=False)
    f.poly(right, "none", GREEN, 3, closed=False)
    f.line(60, 90, 60, 230, MEAS, 3)
    f.line(700, 70, 700, 230, MEAS, 3)
    for st in (60, 300, 480):
        f.line(st, 50, st, 255, INK, 1, dash="3 3")
        f.circle(st, 160, 3, INK)
    f.text(60, 272, "deck point", "small", "middle")
    f.text(300, 272, "deck point", "small", "middle")
    f.text(480, 272, "deck point", "small", "middle")
    f.text(180, 84, "parallel", "small green", "middle")
    f.text(390, 70, "linear", "small green", "middle")
    f.text(590, 64, "parallel", "small green", "middle")
    f.lines(20, 295, ["Edges are simplified to within 1/8 in of offset (Douglas-Peucker); a stretch of constant offset parallels the alignment.",
                      "A top surface stopping short of the deck's outermost points (a chamfer or drip, 1/16 to 6 in) moves the edge out."], "muted small", dy=15)
    f.write()


# ---------------------------------------------------------------------------------------------------------------------
def fig_deck_section():
    f = Figure("deck-section", 760, 505, "Deck cross-section measurements")
    f.text(20, 28, "(a) Cross section normal to the alignment at the 1/4, 1/2, and 3/4 points of each span (left part)", "title")
    s = 3.0          # px per inch
    top = 90         # top of the deck
    gross, haunch, edge, fillet, tf = 8, 10.5, 7, 1.5, 48
    left_edge = 30
    g = [left_edge + (36 + tf / 2) * s + k * 72 * s for k in range(3)]   # girder CLs: 36 in overhang, 6 ft spacing

    def y(d):
        return top + d * s

    soffit = [(left_edge, y(edge)), (g[0] - tf / 2 * s, y(haunch))]
    for k in range(3):
        l, r = g[k] - tf / 2 * s, g[k] + tf / 2 * s
        if k > 0:
            soffit += [(l - fillet * s, y(gross)), (l, y(haunch))]
        soffit.append((r, y(haunch)))
        if k < 2:
            soffit.append((r + fillet * s, y(gross)))
    right = 730
    soffit.append((right, y(haunch)))
    deck = [(left_edge, top), (right, top)] + list(reversed(soffit))
    f.poly(deck, CONCRETE, CONCRETE_EDGE)

    # girders, top part only, with a break line
    for k in range(3):
        cx = g[k]
        if cx > right - 10:
            continue
        half = [(-24, 0), (-24, 3), (-12, 5), (-3, 8), (-3, 26)]
        pts = [(cx + x * s, y(haunch) + yy * s) for x, yy in half]
        pts += [(cx - 3 * s + 4, y(haunch) + 24 * s), (cx + 3 * s - 4, y(haunch) + 28 * s)]
        pts += [(cx - x * s, y(haunch) + yy * s) for x, yy in reversed(half)]
        f.poly(pts, GIRDER, GIRDER_EDGE)
    f.line(right, top, right, y(haunch), "#ffffff", 3)  # the cut edge on the right

    # gross depth: 2 in outside the flange tips, in the bays
    for k in range(2):
        for x in (g[k] + tf / 2 * s + 2 * s, g[k + 1] - tf / 2 * s - 2 * s):
            f.line(x, top, x, y(gross), MEAS, 1.6)
            f.circle(x, top, 2.2, MEAS)
            f.circle(x, y(gross), 2.2, MEAS)
    f.lines(g[0] + tf / 2 * s + 20, top - 30, ["gross depth: 2 in outside", "each flange tip, in the bays"], "meas small", dy=14)
    f.line(g[0] + tf / 2 * s + 6, top - 4, g[0] + tf / 2 * s + 18, top - 22, MEAS, 0.8)

    # haunch over the middle of a flange
    f.line(g[1], top, g[1], y(haunch), GREEN, 1.8)
    f.circle(g[1], y(haunch), 2.2, GREEN)
    f.lines(g[1] - 10, 236, ["over the middle of the flange: deeper", "than the gross depth, a haunch in the deck solid"], "green small", dy=14)
    f.line(g[1] + 1, y(haunch) + 2, g[1] - 4, 224, GREEN, 0.8)

    # fillet
    xr = g[1] + tf / 2 * s
    f.line(xr, y(haunch) + 6, xr, y(haunch) + 14, GREEN, 1)
    f.line(xr + fillet * s, y(gross) + 6, xr + fillet * s, y(haunch) + 14, GREEN, 1)
    f.lines(xr + 20, 176, ["fillet: the", "transition's width"], "green small", dy=14)
    f.line(xr + 3, y(haunch) + 16, xr + 18, 170, GREEN, 0.8)

    # edge depth and taper
    f.line(left_edge + 0.25 * s, top, left_edge + 0.25 * s, y(edge), MEAS, 1.8)
    f.lines(left_edge, y(edge) + 22, ["edge depth, 1/4 in", "inside the edge"], "meas small", dy=14)
    f.line(left_edge + 1 * s, y(edge) + 1, g[0] - tf / 2 * s, y(haunch), ORANGE, 1.6, dash="5 3")
    f.lines(30, 236, ["soffit line (dashed): its depth at the", "exterior flange tip sets the taper"], "orange small", dy=14)

    rows = [("Gross depth", "the median thickness 2 in outside each flange tip in the bays (clear of fillets and crowns)", MEAS),
            ("Haunch in the deck", "the thickness over the middle of the flanges exceeds the gross depth by 1/8 in", GREEN),
            ("Fillet", "the transition's width at the flange tips of interior bays (1/32 in steps); under 1/8 in: square", GREEN),
            ("Edge depth", "1/4 in inside the edge, or the soffit extrapolated to the edge where it's straight (within 1/32 in)", MEAS),
            ("Overhang taper", "the soffit's depth at the exterior flange tip compared with the edge and the haunch, within 1/4 in", ORANGE)]
    for i, (a, b, color) in enumerate(rows):
        cls = {MEAS: "small meas", GREEN: "small green", ORANGE: "small orange"}[color]
        f.text(30, 290 + i * 18, a, cls)
        f.text(160, 290 + i * 18, b, "small")
    f.text(30, 290 + 5 * 18 + 4, "Each is measured at every cut; the import uses the median (depths) or the most common (tapers).", "muted small")

    # (b) tapers
    f.text(20, 418, "(b) Overhang taper", "title")
    k = 1.8
    for i, (name, tip_depth) in enumerate((("none", None), ("to the top of the top flange", 10.5), ("to the bottom of the top flange", 13.5))):
        ox = 40 + i * 240
        tip = ox + 120
        yt = 436
        pts = [(ox, yt), (ox + 200, yt), (ox + 200, yt + 10.5 * k), (tip, yt + 10.5 * k)]
        if tip_depth is None:
            pts += [(tip, yt + 7 * k), (ox, yt + 7 * k)]
        elif tip_depth == 10.5:
            pts += [(ox, yt + 7 * k)]
        else:
            pts += [(tip, yt + 13.5 * k), (ox, yt + 7 * k)]
        f.poly(pts, CONCRETE, CONCRETE_EDGE)
        if tip_depth == 13.5:
            f.poly([(tip, yt + 10.5 * k), (ox + 200, yt + 10.5 * k), (ox + 200, yt + 13.5 * k + 2), (tip, yt + 13.5 * k + 2)], GIRDER, GIRDER_EDGE)
        else:
            f.rect(tip, yt + 10.5 * k, 80, 6, GIRDER, GIRDER_EDGE)
        f.text(ox + 100, yt + 52, name, "small", "middle")
    f.write()


# ---------------------------------------------------------------------------------------------------------------------
def fig_slab_offset():
    f = Figure("slab-offset", 760, 300, "Slab offset at the CL bearing")
    f.text(20, 28, "Elevation at a girder end: the slab offset (\"A\" dimension) is top of deck to top of girder at the CL bearing", "title")
    for ox, reaches in ((30, True), (400, False)):
        x_end = ox + 30
        deck_from = ox if reaches else ox + 90
        f.rect(deck_from, 70, ox + 330 - deck_from, 22, CONCRETE, CONCRETE_EDGE)
        f.rect(x_end, 100, 300, 60, GIRDER, GIRDER_EDGE)
        f.rect(x_end, 92, 300, 8, "#ececec", CONCRETE_EDGE)  # haunch
        bx = x_end + 30
        f.rect(bx - 12, 160, 24, 10, "#f2d9c4", STEEL)
        f.line(bx, 50, bx, 190, INK, 1, dash="8 3 2 3")
        f.text(bx, 205, "CL bearing", "small", "middle")
        if reaches:
            f.dim(bx + 14, 70, bx + 14, 100, None)
            f.text(bx + 20, 88, "slab offset", "meas small")
            f.text(ox + 165, 240, "measured at the CL bearing", "small", "middle")
        else:
            for k in range(4):
                x = deck_from + 10 + k * 30
                f.line(x, 70, x, 100, MEAS, 1.2)
                f.circle(x, 70, 2, MEAS)
            f.line(bx, 70 + 1, deck_from + 100, 70, MEAS, 1, dash="3 3")
            f.circle(bx, 71, 3, MEAS)
            f.text(ox + 165, 240, "no slab over the bearing (e.g. an integral end):", "small", "middle")
            f.text(ox + 165, 255, "a line through four samples 1 ft apart, extrapolated", "small", "middle")
    f.lines(20, 282, ["One slab offset for the bridge if all agree within 1/16 in, else one per bearing line, else one per girder end."], "muted small")
    f.write()


# ---------------------------------------------------------------------------------------------------------------------
def fig_roadway():
    f = Figure("roadway", 760, 320, "Roadway section from the top of the deck")
    f.text(20, 28, "Top of the deck at a cut: its corners are the ridge points, with one at the alignment", "title")
    base = 150
    pts = [(40, base + 30), (70, base + 18), (380, base - 22), (690, base + 18), (720, base + 30)]
    # deck body
    body = [(40, base + 30), (70, base + 18), (380, base - 22), (690, base + 18), (720, base + 30), (720, base + 60), (40, base + 60)]
    f.poly(body, CONCRETE, CONCRETE_EDGE)
    f.poly(pts[1:4], "none", MEAS, 2.4, closed=False)
    f.line(40, base + 30, 70, base + 18, MUTED, 2.4, dash="3 2")
    f.line(690, base + 18, 720, base + 30, MUTED, 2.4, dash="3 2")
    for p in pts[1:4]:
        f.circle(*p, 4, MEAS)
    f.line(380, 60, 380, 240, ALIGN, 2, dash="10 4 2 4")
    f.text(386, 72, "alignment: ridge point (a corner within 1/2 in moves here)", "align small")
    f.text(210, base - 12, "slope", "meas small", "middle")
    f.text(550, base - 12, "slope", "meas small", "middle")
    f.lines(560, 230, ["edges steeper than 25%:", "a chamfer or drip, left out"], "muted small", dy=14)
    f.lines(20, 270, ["Cuts at the quarter points of each span and as near the bridge ends as the deck is full width. One template if the",
                      "sections agree within 1/8 in; otherwise templates where the section changes (1/16 in), extended to the end piers."], "muted small", dy=15)
    f.write()


# ---------------------------------------------------------------------------------------------------------------------
def fig_profile():
    f = Figure("profile", 760, 300, "Profile fitted to the top of the deck")
    f.text(20, 28, "Elevation along the alignment: samples of the deck top, fitted with tangents and vertical curves", "title")
    import random
    random.seed(4)
    x0, x1 = 50, 710

    def vc(x):
        g1, g2, L = -0.25, 0.18, 180
        u = x - 300
        return 200 - 0.25 * (300 - x0) + g1 * u + (g2 - g1) / (2 * L) * u * u

    def y(x):
        if x < 300:
            return 200 - 0.25 * (x - x0)
        if x <= 480:
            return vc(x)
        return vc(480) + 0.18 * (x - 480)

    curve = [(x, 300 - y(x) - 20) for x in range(x0, x1 + 1, 5)]
    f.poly(curve, "none", ALIGN, 2, closed=False)
    for xs in range(x0, x1 + 1, 30):
        yy = 300 - y(xs) - 20 + random.uniform(-1.5, 1.5)
        if xs == 410:
            f.circle(xs, yy - 16, 3.5, "#ffffff", MEAS, 1.6)
            f.text(xs, yy - 40, "outlier (1/2 in off the local fit): left out", "meas small", "middle")
        else:
            f.circle(xs, yy, 3, INK)
    f.circle(300, 300 - y(300) - 20, 4, GREEN)
    f.circle(480, 300 - y(480) - 20, 4, GREEN)
    f.text(390, 300 - y(390) - 20 + 30, "symmetric parabolic vertical curve (green: its ends)", "green small", "middle")
    f.lines(20, 262, ["Samples about every 10 ft (at least 8 per span). The simplest profile within 1/8 in becomes the profile when the",
                      "alignment has none; when it has one, the deck is compared with it and a difference is logged."], "muted small", dy=15)
    f.write()


# ---------------------------------------------------------------------------------------------------------------------
def fig_barriers():
    f = Figure("barriers", 760, 330, "Traffic barrier section")
    f.text(20, 28, "(a) Cut at mid-bridge, normal to the alignment", "title")
    # deck with a cross slope and a barrier at the left edge, rotated with the slope
    slope = math.radians(-5)
    ox, oy = 40, 200

    def r(x, y):
        return (ox + x * math.cos(slope) - y * math.sin(slope), oy + x * math.sin(slope) + y * math.cos(slope))

    deck = [r(0, 0), r(300, 0), r(300, 20), r(0, 20)]
    f.poly(deck, CONCRETE, CONCRETE_EDGE)
    barrier = [(0, 0), (40, 0), (34, -20), (18, -32), (16, -90), (0, -90)]
    f.poly([r(x, y) for x, y in barrier], "#d8cfc4", "#7c6a55", 1.4)
    f.circle(*r(0, 0), 3.5, MEAS)
    f.text(r(0, 0)[0] - 4, r(0, 0)[1] + 34, "deck edge at the top of the deck", "meas small")
    f.text(170, 120, "the outermost barrier on each side", "small", "middle")
    f.text(170, 135, "is the exterior barrier", "small", "middle")
    # (b) barrier coordinates
    f.text(400, 28, "(b) In traffic barrier coordinates", "title")
    bx, by, s = 470, 230, 1.5
    f.poly([(bx + x * s, by + y * s) for x, y in barrier], "#d8cfc4", "#7c6a55", 1.4)
    f.poly([(bx + x * s, by + y * s) for x, y in [(0, 0), (41, 0), (35, -19), (18, -32.5), (16, -90), (0, -90)]], "none", ALIGN, 1.4, dash="5 3")
    f.line(bx, by, bx + 100, by, GREEN, 1.6, arrow="arrow-green")
    f.text(bx + 104, by + 4, "toward the roadway", "green small")
    f.line(bx, by, bx, by - 150, GREEN, 1.6, arrow="arrow-green")
    f.circle(bx, by, 3.5, MEAS)
    f.text(bx + 6, by + 18, "origin: the deck edge", "meas small")
    f.lines(560, 110, ["turned back by the cross slope", "at the edge (PGSuper turns it to", "the slope again); the right barrier", "is mirrored"], "small", dy=15)
    f.text(560, 190, "dashed: a library entry", "align small")
    f.lines(20, 285, ["A library entry whose outline is within 1/8 in (Hausdorff distance) is used, the one named in the model if it's that close;",
                      "otherwise a new entry \"<name> (IFC)\" is made from the section, its weight computed."], "muted small", dy=15)
    f.write()


# ---------------------------------------------------------------------------------------------------------------------
def fig_deck_rebar():
    f = Figure("deck-rebar", 760, 350, "Deck reinforcement as lump sums")
    f.text(20, 28, "(a) Deck section: primary lines in two mats", "title")
    f.rect(40, 70, 300, 70, CONCRETE, CONCRETE_EDGE)
    for x in range(60, 330, 30):
        f.circle(x, 88, 3.6, STEEL)
        f.circle(x, 124, 3.6, STEEL)
    f.dim(356, 70, 356, 88, None)
    f.text(362, 84, "top cover", "meas small")
    f.dim(356, 124, 356, 140, None)
    f.text(362, 136, "bottom cover", "meas small")
    f.dim(60, 160, 330 + 30, 160, None)
    f.text(210, 180, "width covered: first line to last, plus one spacing", "meas small", "middle")
    f.lines(20, 210, ["lump sum = area of one bar of each primary line / width covered",
                      "cover: the depth of the primary bars' centroid (top mat), or its",
                      "height above the bottom (bottom mat)",
                      "mats: the primary lines' depths split at their largest gap (over 1 in)"], "small", dy=16)
    f.text(500, 28, "(b) Plan: supplemental bars at a pier", "title")
    f.line(500, 70, 740, 70, CONCRETE_EDGE, 1)
    f.line(500, 170, 740, 170, CONCRETE_EDGE, 1)
    f.line(620, 55, 620, 185, INK, 1.2, dash="8 3 2 3")
    f.text(620, 200, "CL pier", "small", "middle")
    for y in (85, 100, 115, 130, 145, 160):
        f.line(500, y, 740, y, STEEL, 1.2)
    for y in (92, 122, 152):
        f.line(560, y, 690, y, ORANGE, 2.2)
    f.dim(560, 214, 620, 214, None)
    f.dim(620, 214, 690, 214, None)
    f.text(625, 232, "cutoffs: median distances to the bar ends", "meas small", "middle")
    f.lines(500, 256, ["primary lines cover at least 75% of the bridge",
                       "length; the others are supplemental, at the",
                       "pier they cross nearest their middle"], "small", dy=15)
    f.lines(20, 318, ["Deck bars: centerlines 90% inside the deck slab, more along the alignment than across it. Bars within 2 in of offset",
                      "and 1 in of depth are a line (e.g. bars lapped end to end)."], "muted small", dy=15)
    f.write()


if __name__ == "__main__":
    for fn in [fig_coordinates, fig_mesh_steps, fig_pipeline, fig_girder_layout, fig_spacing, fig_pier_orientation, fig_bearings,
               fig_continuity, fig_girder_section, fig_section_fit, fig_strands, fig_stirrups, fig_confinement, fig_longitudinal_bars,
               fig_deck_plan, fig_deck_section, fig_slab_offset, fig_roadway, fig_profile, fig_barriers, fig_deck_rebar]:
        fn()
    print("figures written to", HERE)
