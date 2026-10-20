"""Reduces the prestressing strands of an IFC model to a flat summary, in the form of a pgs_extract.py summary, so
compare_bridge.py can compare them.

The round trip compares the strands of the export of the source project with the strands of an export of the
imported project. Both exports are written by PGSuper's exporter from PGSuper's own strand geometry, so the
comparison is physical: it doesn't depend on how the strands are defined in PGSuper (total number, straight and
harped numbers, grid fill, rows, or individual strands), and harped strand offsets, extended strands, and the
strand grid of the library entry are already applied.

The strands of a girder are the IfcTendon parts of the IfcElementAssembly parts of its IfcBeam. Each strand is
classified by its path and status:
  harped     the path has interior points (the harp points)
  temporary  Status TEMPORARY (Pset_ElementComponentCommon or usBrPset_Common, on the strand or its type)
  straight   otherwise

Values for each girder (group[i].girder[j].strands, from the beam name "Span n, Girder m", else the beam name):
  count, <type>.count                   number of strands (every type, 0 if none)
  <type>.material                       strand material names
  <type>.prestress                      IfcTendon.PreStress (jacking stress)
  <type>.tension_force                  IfcTendon.TensionForce (jacking force of one strand)
  <type>.paths                          the strand centerlines, in girder coordinates (x along the girder from
                                        the beam's placement, y to the side, z up from the top of the girder)
  <type>.debonded.count                 strands with a debond length (usBrPset_TendonDebonding)
  <type>.debond_lengths                 (start, end) debond lengths of the debonded strands
  harp_points                           the x of the harp points of the harped strands

Values are SI.

Usage: python ifc_strands.py <file.ifc> [--out summary.json]
"""
import argparse
import json
import re

import numpy as np
import ifcopenshell
import ifcopenshell.util.element as element
import ifcopenshell.util.placement as placement
import ifcopenshell.util.unit as unit

STRAND_TYPES = ("straight", "harped", "temporary")


def girder_key(beam):
    m = re.fullmatch(r"Span (\d+), Girder (\d+)", beam.Name or "")
    return f"group[{int(m[1]) - 1}].girder[{int(m[2]) - 1}]" if m else f"girder[{beam.Name or beam.GlobalId}]"


def strands_of(beam):
    return [t for part in element.get_parts(beam) if part.is_a("IfcElementAssembly")
            for t in element.get_parts(part) if t.is_a("IfcTendon")]


def status(tendon):
    psets = element.get_psets(tendon)  # includes the type's property sets
    for pset in ("Pset_ElementComponentCommon", "usBrPset_Common"):
        value = psets.get(pset, {}).get("Status")
        if value:
            return (value[0] if isinstance(value, (list, tuple)) else value).upper()
    return None


def path(tendon, to_girder, scale):
    """the points of the strand's centerline in girder coordinates (SI), or None if it isn't a polyline swept disk"""
    for representation in tendon.Representation.Representations if tendon.Representation else []:
        for item in representation.Items:
            if item.is_a("IfcSweptDiskSolid") and item.Directrix.is_a("IfcPolyline"):
                m = to_girder @ placement.get_placement(tendon.ObjectPlacement, should_return_si=True)
                points = []
                for p in item.Directrix.Points:
                    c = [v * scale for v in p.Coordinates] + [0.0] * (3 - len(p.Coordinates))
                    points.append([float(v) for v in (m @ np.array(c + [1.0]))[:3]])
                return points
    return None


def length_property(obj, pset_name, name):
    """a length property of obj, else of its type, in metres (with the property's own unit, else the project's), or None"""
    psets = [rel.RelatingPropertyDefinition for rel in getattr(obj, "IsDefinedBy", None) or [] if rel.is_a("IfcRelDefinesByProperties")]
    obj_type = element.get_type(obj)
    psets += list(getattr(obj_type, "HasPropertySets", None) or [])
    for pset in psets:
        if not pset.is_a("IfcPropertySet") or pset.Name != pset_name:
            continue
        for p in pset.HasProperties:
            if p.Name == name and p.is_a("IfcPropertySingleValue") and p.NominalValue is not None:
                u = unit.get_property_unit(p, obj.file)
                value = p.NominalValue.wrappedValue
                return unit.convert(value, getattr(u, "Prefix", None), u.Name, None, "METRE") if u else value
    return None


def debond_lengths(tendon):
    start = length_property(tendon, "usBrPset_TendonDebonding", "DebondLengthStart") or 0.0
    end = length_property(tendon, "usBrPset_TendonDebonding", "DebondLengthEnd") or 0.0
    return (start, end) if 0 < start or 0 < end else None


def single_or_list(values):
    """the value if all values are the same, else the sorted distinct values"""
    distinct = sorted(set(values))
    return distinct[0] if len(distinct) == 1 else distinct


def extract(path_to_ifc):
    f = ifcopenshell.open(str(path_to_ifc))
    scale = unit.calculate_unit_scale(f)
    values = {}

    def put(key, value, kind):
        values[key] = {"value": value, "kind": kind}

    for beam in f.by_type("IfcBeam"):
        tendons = strands_of(beam)
        if not tendons and element.get_predefined_type(beam) not in ("BEAM", "GIRDER_SEGMENT", "T_BEAM"):  # e.g. diaphragms
            continue
        key = f"{girder_key(beam)}.strands"
        to_girder = np.linalg.inv(placement.get_placement(beam.ObjectPlacement, should_return_si=True))

        by_type = {t: [] for t in STRAND_TYPES}
        for tendon in tendons:
            points = path(tendon, to_girder, scale)
            if status(tendon) == "TEMPORARY":
                strand_type = "temporary"
            elif points and len(points) > 2:
                strand_type = "harped"
            else:
                strand_type = "straight"
            by_type[strand_type].append((tendon, points))

        put(f"{key}.count", len(tendons), "count")
        for strand_type, strands in by_type.items():
            put(f"{key}.{strand_type}.count", len(strands), "count")
            if not strands:
                continue
            materials = [m.Name for m in (element.get_material(t) for t, _ in strands) if m is not None and m.is_a("IfcMaterial")]
            if materials:
                put(f"{key}.{strand_type}.material", "; ".join(sorted(set(materials))), "text")
            prestress = [t.PreStress for t, _ in strands if t.PreStress is not None]
            if prestress:
                put(f"{key}.{strand_type}.prestress", single_or_list(prestress), "stress")
            tension = [t.TensionForce for t, _ in strands if t.TensionForce is not None]
            if tension:
                put(f"{key}.{strand_type}.tension_force", single_or_list(tension), "force")
            put(f"{key}.{strand_type}.paths", [p for _, p in strands if p], "paths")
            debonds = [d for d in (debond_lengths(t) for t, _ in strands) if d]
            put(f"{key}.{strand_type}.debonded.count", len(debonds), "count")
            if debonds:
                put(f"{key}.{strand_type}.debond_lengths", sorted([s, e] for s, e in debonds), "length_pairs")

        harp_points = sorted({round(p[0], 3) for _, points in by_type["harped"] if points for p in points[1:-1]})
        if harp_points:
            put(f"{key}.harp_points", harp_points, "lengths")

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
