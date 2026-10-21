"""Reduces the longitudinal girder bars of an IFC model to a flat summary, in the form of a pgs_extract.py summary, so
compare_bridge.py can compare them (see ifc_strands.py, which does the same for strands).

The round trip compares the bars of the export of the source project with the bars of an export of the imported
project. Both exports are written by PGSuper's exporter from PGSuper's own rebar layout, so the comparison is
physical: it doesn't depend on how the bars are arranged in rows.

The bars of a girder are the IfcReinforcingBar parts of the IfcElementAssembly parts of its IfcBeam; an element with
several items in its body representation is a group of bars, one for each item. A longitudinal bar is a bar whose
centerline is one straight line along the girder (within 2 degrees of the beam's x axis).

Values for each girder (group[i].girder[j].longitudinal_bars, from the beam name "Span n, Girder m", else the beam name):
  count       number of longitudinal bars
  sizes       the number of bars of each size (e.g. "6 #4; 4 #5"), from the bar type's name
  material    bar material names
  paths       the bar centerlines (start and end), in girder coordinates (x along the girder from the beam's
              placement, y to the side, z up from the top of the girder)

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


def centerlines(bar, to_girder, scale):
    """the centerlines of the bar's swept disks in girder coordinates (SI): one for each item of its body representation (a
    group of bars of one type is one element with a mapped item for each bar)"""
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
                points = directrix_points(item.Directrix, scale)
                if points:
                    lines.append([[float(v) for v in (m @ np.array(p + [1.0]))[:3]] for p in points])
                    break
    return lines


def is_longitudinal(points):
    if not points or len(points) != 2:
        return False
    d = np.array(points[1]) - np.array(points[0])
    length = np.linalg.norm(d)
    return 0 < length and math.acos(min(1.0, abs(d[0]) / length)) <= MAX_ANGLE


def bar_size(bar):
    """the bar size from the bar's type name (e.g. "Girder_Longitudinal_Bar_#5"), else its nominal diameter"""
    bar_type = element.get_type(bar)
    m = re.search(r"#\d+", (bar_type.Name if bar_type else None) or bar.Name or "")
    if m:
        return m[0]
    diameter = getattr(bar_type, "NominalDiameter", None) or getattr(bar, "NominalDiameter", None)
    return f"{diameter:.4f}" if diameter else "?"


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
        for bar in bars:
            for points in centerlines(bar, to_girder, scale):
                if is_longitudinal(points):
                    longitudinal.append((bar, sorted(points)))  # start first

        put(f"{key}.count", len(longitudinal), "count")
        if not longitudinal:
            continue
        sizes = collections.Counter(bar_size(b) for b, _ in longitudinal)
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
