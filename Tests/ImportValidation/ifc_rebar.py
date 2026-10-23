"""Reduces the longitudinal girder bars and stirrups of an IFC model to a flat summary, in the form of a pgs_extract.py summary, so
compare_bridge.py can compare them (see ifc_strands.py, which does the same for strands).

The round trip compares the bars of the export of the source project with the bars of an export of the imported
project. Both exports are written by PGSuper's exporter from PGSuper's own rebar layout, so the comparison is
physical: it doesn't depend on how the bars are arranged in rows.

The bars of a girder are the IfcReinforcingBar parts of the IfcElementAssembly parts of its IfcBeam; an element with
several items in its body representation is a group of bars, one for each item. A longitudinal bar is a bar whose
centerline is one straight line along the girder (within 2 degrees of the beam's x axis). A stirrup is a bar with a
vertical leg: a straight piece of its centerline within 15 degrees of vertical and at least 1 ft long (enough for the
exporter's bars: its confinement bars' legs are shorter). Its position is the middle of its legs; stirrups within 1 in
of each other are at one place.

Values for each girder (group[i].girder[j].longitudinal_bars, from the beam name "Span n, Girder m", else the beam name):
  count       number of longitudinal bars
  sizes       the number of bars of each size (e.g. "6 #4; 4 #5"), from the bar type's name
  material    bar material names
  paths       the bar centerlines (start and end), in girder coordinates (x along the girder from the beam's
              placement, y to the side, z up from the top of the girder)

Values for the stirrups of each girder (group[i].girder[j].stirrups):
  count       number of stirrups
  sizes       the number of stirrups of each size
  places      where the stirrups are along the girder (x of each place)
  legs        the number of places with each number of legs (e.g. "64 x 2"), and of legs above the top of the girder

A confinement bar is a bar across the girder (not longitudinal, not a stirrup) whose centerline is all below 40% of the
girder's depth (taken from the lowest bar of the girder), e.g. the exporter's G9 and G10 bars. Values for each girder
(group[i].girder[j].confinement):
  count       number of confinement bars
  sizes       the number of confinement bars of each size
  places      where they are along the girder (x of each place, bars within 1 in at one place)

Values are SI.

Usage: python ifc_rebar.py <file.ifc> [--out summary.json]
"""
import argparse
import collections
import json
import math
import re

import numpy as np
import ifcopenshell
import ifcopenshell.util.element as element
import ifcopenshell.util.placement as placement
import ifcopenshell.util.unit as unit

from ifc_strands import girder_key

MAX_ANGLE = math.radians(2.0)
MAX_LEG_SLOPE = math.tan(math.radians(15.0))
MIN_LEG = 0.3048
PLACE_TOLERANCE = 0.0254
PROJECTION = 0.0254
BAR_DIAMETERS = {3: 0.375, 4: 0.5, 5: 0.625, 6: 0.75, 7: 0.875, 8: 1.0, 9: 1.128, 10: 1.27, 11: 1.41, 14: 1.693, 18: 2.257}  # in


def bars_of(beam):
    return [b for part in element.get_parts(beam) if part.is_a("IfcElementAssembly")
            for b in element.get_parts(part) if b.is_a("IfcReinforcingBar")]


def directrix_points(curve, scale):
    """the points of a polyline or an indexed poly curve, in its coordinates (SI), or None"""
    if curve.is_a("IfcPolyline"):
        coordinates = [p.Coordinates for p in curve.Points]
    elif curve.is_a("IfcIndexedPolyCurve"):
        coordinates = curve.Points.CoordList
    else:
        return None
    return [[v * scale for v in c] + [0.0] * (3 - len(c)) for c in coordinates]


def directrix_lines(curve, scale):
    """the straight pieces of a polyline or an indexed poly curve (its arcs left out), in its coordinates (SI)"""
    points = directrix_points(curve, scale)
    if points is None:
        return []
    if curve.is_a("IfcIndexedPolyCurve") and curve.Segments:
        lines = []
        for segment in curve.Segments:
            if segment.is_a("IfcLineIndex"):
                indices = segment.wrappedValue
                lines += [(points[a - 1], points[b - 1]) for a, b in zip(indices, indices[1:])]
        return lines
    return list(zip(points, points[1:]))


def centerlines(bar, to_girder, scale, pieces=False):
    """the centerlines of the bar's swept disks in girder coordinates (SI): one for each item of its body representation (a
    group of bars of one type is one element with a mapped item for each bar). With pieces, each centerline is its
    straight pieces (pairs of points) instead of its points"""
    def swept_disks(representation_items, m):
        for item in representation_items:
            if item.is_a("IfcMappedItem"):
                yield from swept_disks(item.MappingSource.MappedRepresentation.Items, m @ placement.get_mappeditem_transformation(item))
            elif item.is_a("IfcSweptDiskSolid"):
                yield item, m

    m0 = to_girder @ placement.get_placement(bar.ObjectPlacement, should_return_si=True)
    lines = []
    for representation in bar.Representation.Representations if bar.Representation else []:
        if representation.RepresentationIdentifier not in (None, "Body"):
            continue
        for top_item in representation.Items:
            for item, m in swept_disks([top_item], m0):
                to = lambda p: [float(v) for v in (m @ np.array(list(p) + [1.0]))[:3]]
                if pieces:
                    segments = directrix_lines(item.Directrix, scale)
                    if segments:
                        lines.append([(to(a), to(b)) for a, b in segments])
                        break
                else:
                    points = directrix_points(item.Directrix, scale)
                    if points:
                        lines.append([to(p) for p in points])
                        break
    return lines


def is_longitudinal(points):
    if not points or len(points) != 2:
        return False
    d = np.array(points[1]) - np.array(points[0])
    length = np.linalg.norm(d)
    return 0 < length and math.acos(min(1.0, abs(d[0]) / length)) <= MAX_ANGLE


def stirrup_of(segments):
    """(x, number of legs, number of legs above the top of the girder) of a bar's straight pieces, or None"""
    legs = []
    for a, b in segments:
        dz = abs(b[2] - a[2])
        if MIN_LEG <= dz and math.hypot(b[0] - a[0], b[1] - a[1]) <= MAX_LEG_SLOPE * dz:
            legs.append(((a[0] + b[0]) / 2, PROJECTION < max(a[2], b[2])))
    if not legs:
        return None
    return sum(x for x, _ in legs) / len(legs), len(legs), sum(1 for _, above in legs if above)


def places_of(stirrups):
    """[(x, legs, legs above the top)] of the stirrups (x, legs, above), those within PLACE_TOLERANCE of the first at a place"""
    places = []
    for x, legs, above in sorted(stirrups):
        if places and x - places[-1][3] <= PLACE_TOLERANCE:
            place = places[-1]
            place[0].append(x)
            place[1] += legs
            place[2] += above
        else:
            places.append([[x], legs, above, x])
    return [(sum(xs) / len(xs), legs, above) for xs, legs, above, _ in places]


def is_confinement(points, depth):
    """the bar (its centerline points) is across the girder and all below 40% of its depth"""
    xs = [p[0] for p in points]
    ys = [p[1] for p in points]
    return 3 <= len(points) and max(xs) - min(xs) <= 2 * (max(ys) - min(ys)) and all(p[2] <= -0.4 * depth for p in points)


def bar_size(bar, scale):
    """the bar size from the bar's type name (e.g. "Girder_Longitudinal_Bar_#5"), else from its nominal diameter"""
    bar_type = element.get_type(bar)
    m = re.search(r"#\d+", (bar_type.Name if bar_type else None) or bar.Name or "")
    if m:
        return m[0]
    diameter = getattr(bar_type, "NominalDiameter", None) or getattr(bar, "NominalDiameter", None)
    if not diameter:
        return "?"
    diameter *= scale / 0.0254  # in
    size, d = min(BAR_DIAMETERS.items(), key=lambda item: abs(item[1] - diameter))
    return f"#{size}" if abs(d - diameter) <= 0.01 else f"{diameter:.3f} in"


def extract(path_to_ifc):
    f = ifcopenshell.open(str(path_to_ifc))
    scale = unit.calculate_unit_scale(f)
    values = {}

    def put(key, value, kind):
        values[key] = {"value": value, "kind": kind}

    for beam in f.by_type("IfcBeam"):
        bars = bars_of(beam)
        if not bars and element.get_predefined_type(beam) not in ("BEAM", "GIRDER_SEGMENT", "T_BEAM"):  # e.g. diaphragms
            continue
        key = f"{girder_key(beam)}.longitudinal_bars"
        to_girder = np.linalg.inv(placement.get_placement(beam.ObjectPlacement, should_return_si=True))

        longitudinal = []
        stirrups = []
        others = []  # neither longitudinal nor stirrups
        for bar in bars:
            for points in centerlines(bar, to_girder, scale):
                if is_longitudinal(points):
                    longitudinal.append((bar, sorted(points)))  # start first
            for segments, points in zip(centerlines(bar, to_girder, scale, pieces=True), centerlines(bar, to_girder, scale)):
                stirrup = stirrup_of(segments)
                if stirrup:
                    stirrups.append((bar, stirrup))
                elif not is_longitudinal(points):
                    others.append((bar, points))

        depth = -min((p[2] for _, points in others for p in points), default=0.0)
        confinement = [(b, points) for b, points in others if is_confinement(points, depth)]
        put(f"{girder_key(beam)}.confinement.count", len(confinement), "count")
        if confinement:
            ckey = f"{girder_key(beam)}.confinement"
            sizes = collections.Counter(bar_size(b, scale) for b, _ in confinement)
            put(f"{ckey}.sizes", "; ".join(f"{n} {size}" for size, n in sorted(sizes.items())), "text")
            places = places_of([(sum(p[0] for p in points) / len(points), 1, 0) for _, points in confinement])
            put(f"{ckey}.places", [x for x, _, _ in places], "lengths")

        put(f"{girder_key(beam)}.stirrups.count", len(stirrups), "count")
        if stirrups:
            skey = f"{girder_key(beam)}.stirrups"
            sizes = collections.Counter(bar_size(b, scale) for b, _ in stirrups)
            put(f"{skey}.sizes", "; ".join(f"{n} {size}" for size, n in sorted(sizes.items())), "text")
            places = places_of([s for _, s in stirrups])
            put(f"{skey}.places", [x for x, _, _ in places], "lengths")
            legs = collections.Counter((n, above) for _, n, above in places)
            put(f"{skey}.legs", "; ".join(f"{count} x {n}" + (f" ({above} above the top)" if above else "") for (n, above), count in sorted(legs.items())), "text")

        put(f"{key}.count", len(longitudinal), "count")
        if not longitudinal:
            continue
        sizes = collections.Counter(bar_size(b, scale) for b, _ in longitudinal)
        put(f"{key}.sizes", "; ".join(f"{n} {size}" for size, n in sorted(sizes.items())), "text")
        materials = [m.Name for m in (element.get_material(b) for b, _ in longitudinal) if m is not None and m.is_a("IfcMaterial")]
        if materials:
            put(f"{key}.material", "; ".join(sorted(set(materials))), "text")
        put(f"{key}.paths", [p for _, p in longitudinal], "bar_paths")

    return {"meta": {"file": str(path_to_ifc)}, "values": values}


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("file")
    parser.add_argument("--out")
    args = parser.parse_args()
    text = json.dumps(extract(args.file), indent=1)
    if args.out:
        with open(args.out, "w", encoding="utf-8") as f:
            f.write(text)
    else:
        print(text)


if __name__ == "__main__":
    main()
