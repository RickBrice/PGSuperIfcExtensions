"""Reduces a PGSuper project (.pgs) or template (.pgt) to a flat summary of the values the
IFC importer is responsible for.

PGSuper files store values in system units (kg, m, s, Pa, rad), so every value in the summary
is SI. Conversion to US units happens only when a report is written (see compare_bridge.py).

The summary is {"meta": {...}, "values": {key: {"value": v, "kind": k}}}. Keys are paths like
"pier[1].station" or "group[0].girder[2].fc". The kind tells compare_bridge.py how to compare
and display the value (see KINDS in compare_bridge.py).

Alignment and profile are compared by evaluating them at the pier stations, not by comparing
their defining data, because an imported alignment may be defined from a different reference
station or point (or with different elements) and still be the same alignment.

Usage: python pgs_extract.py <file.pgs> [--out summary.json]
"""
import argparse
import json
import math
import xml.etree.ElementTree as ET


def _f(e, tag, default=None):
    """float value of child tag, or default"""
    if e is None:
        return default
    c = e.find(tag)
    if c is None or c.text is None:
        return default
    return float(c.text)


def _t(e, tag, default=None):
    """text value of child tag, or default"""
    if e is None:
        return default
    c = e.find(tag)
    return default if c is None or c.text is None else c.text


class HorizontalAlignment:
    """Evaluates a PGSuper alignment made of tangents and circular curves.

    Mirrors CBridgeAgentImp's alignment construction: PI stations are measured along the back
    tangent from the previous curve, FwdTangent is either the forward tangent direction or the
    curve delta (forward = back + delta), and the alignment is positioned so the reference
    station is at the reference point. Spirals are not supported (supported == False).
    """

    def __init__(self, alignment):
        self.supported = True
        self.ref_station = _f(alignment, "RefStation", 0.0)
        self.ref_point = (_f(alignment, "RefPointEasting", 0.0), _f(alignment, "RefPointNorthing", 0.0))
        direction = _f(alignment, "Direction", 0.0)

        curves = alignment.findall("HorzCurveData")
        start_station = min([self.ref_station] + [_f(c, "PIStation") for c in curves]) - 1.0e5

        # elements in a local coordinate system, starting at start_station at (0,0)
        self.elements = []  # (start station, end station, kind, data)
        station, point, back = start_station, (0.0, 0.0), direction
        for curve in curves:
            if 0 < _f(curve, "EntrySpiral", 0.0) or 0 < _f(curve, "ExitSpiral", 0.0):
                self.supported = False
                return

            pi_station = _f(curve, "PIStation")
            fwd = _f(curve, "FwdTangent")
            if _t(curve, "FwdTangentIsBearing", "0") == "0":
                fwd += back  # FwdTangent is the curve delta
            delta = math.atan2(math.sin(fwd - back), math.cos(fwd - back))
            radius = _f(curve, "Radius", 0.0)
            pi = (point[0] + (pi_station - station) * math.cos(back), point[1] + (pi_station - station) * math.sin(back))

            if radius <= 0 or abs(delta) < 1.0e-12:
                # angle point
                self.elements.append((station, pi_station, "line", (point, back)))
                station, point, back = pi_station, pi, fwd
                continue

            T = radius * math.tan(abs(delta) / 2)
            ts_station = pi_station - T
            ts = (pi[0] - T * math.cos(back), pi[1] - T * math.sin(back))
            self.elements.append((station, ts_station, "line", (point, back)))

            arc_length = radius * abs(delta)
            side = 1 if 0 < delta else -1  # 1 = curve to the left
            center = (ts[0] - side * radius * math.sin(back), ts[1] + side * radius * math.cos(back))
            self.elements.append((ts_station, ts_station + arc_length, "arc", (ts, back, radius, side, center)))

            station = ts_station + arc_length
            point = (pi[0] + T * math.cos(fwd), pi[1] + T * math.sin(fwd))
            back = fwd

        self.elements.append((station, math.inf, "line", (point, back)))

        # shift so the reference station is at the reference point
        x, y = self._local(self.ref_station)[0]
        self.shift = (self.ref_point[0] - x, self.ref_point[1] - y)

    def _local(self, station):
        for start, end, kind, data in self.elements:
            if station <= end or end == math.inf:
                if kind == "line":
                    point, direction = data
                    d = station - start
                    return (point[0] + d * math.cos(direction), point[1] + d * math.sin(direction)), direction
                ts, back, radius, side, center = data
                angle = (station - start) / radius
                direction = back + side * angle
                return (center[0] + side * radius * math.sin(direction), center[1] - side * radius * math.cos(direction)), direction
        raise ValueError(station)

    def evaluate(self, station):
        """(easting, northing), direction at station"""
        (x, y), direction = self._local(station)
        return (x + self.shift[0], y + self.shift[1]), direction


class Profile:
    """Evaluates a PGSuper profile (grades and vertical curves). See CBridgeAgentImp.
    If L2 is zero, the curve is symmetric and L1 is the full curve length."""

    def __init__(self, profile):
        self.station = _f(profile, "Station", 0.0)
        self.elevation = _f(profile, "Elevation", 0.0)
        self.grade = _f(profile, "Grade", 0.0)
        self.curves = []  # (pvi station, pvi elevation, g1, g2, L1, L2)
        station, elevation, grade = self.station, self.elevation, self.grade
        for curve in profile.findall("VertCurveData"):
            pvi = _f(curve, "PVIStation")
            L1, L2 = _f(curve, "L1", 0.0), _f(curve, "L2", 0.0)
            if L2 == 0:
                L1 = L2 = L1 / 2
            exit_grade = _f(curve, "ExitGrade", grade)
            pvi_elevation = elevation + grade * (pvi - station)
            self.curves.append((pvi, pvi_elevation, grade, exit_grade, L1, L2))
            station, elevation, grade = pvi, pvi_elevation, exit_grade

    def evaluate(self, s):
        """elevation, grade at station s"""
        if not self.curves:
            return self.elevation + self.grade * (s - self.station), self.grade

        for pvi, elev, g1, g2, L1, L2 in self.curves:
            if s <= pvi + L2:
                bvc, evc = pvi - L1, pvi + L2
                if s <= bvc:
                    return elev + g1 * (s - pvi), g1
                if 0 < L1 + L2:
                    e = L1 * L2 * (g2 - g1) / (2 * (L1 + L2))  # offset at the PVI
                    if s <= pvi:
                        u = (s - bvc) / L1
                        return elev + g1 * (s - pvi) + e * u * u, g1 + 2 * e * u / L1
                    u = (evc - s) / L2
                    return elev + g2 * (s - pvi) + e * u * u, g2 - 2 * e * u / L2
        pvi, elev, g1, g2, L1, L2 = self.curves[-1]
        return elev + g2 * (s - pvi), g2


def _dms(parts):
    """degrees from [ddd[, mm[, ss.s]]] strings"""
    return sum(float(v) / 60 ** i for i, v in enumerate(parts))


def pier_skew(orientation, alignment_direction):
    """skew of a pier (radians, positive is left/counterclockwise) from its PGSuper orientation string:
    NORMAL, a skew angle ("20 L", "15 30 R", "-10"), or a bearing of the pier line ("N 35 E")"""
    text = orientation.strip().upper()
    if text in ("N", "NORMAL"):
        return 0.0
    parts = text.split()
    if parts[0] in ("N", "S") and parts[-1] in ("E", "W"):
        angle = math.radians(_dms(parts[1:-1]))
        direction = {("N", "E"): math.pi / 2 - angle, ("N", "W"): math.pi / 2 + angle,
                     ("S", "E"): -math.pi / 2 + angle, ("S", "W"): -math.pi / 2 - angle}[(parts[0], parts[-1])]
        # PGSuper's normal pier line points to the left of the alignment
        skew = direction - (alignment_direction + math.pi / 2)
        skew = math.atan2(math.sin(skew), math.cos(skew))
        if skew <= -math.pi / 2:
            skew += math.pi
        elif math.pi / 2 < skew:
            skew -= math.pi
        return skew
    sign = 1.0
    if parts[-1] in ("L", "R"):
        sign = -1.0 if parts[-1] == "R" else 1.0
        parts = parts[:-1]
    elif parts[-1][-1:] in ("L", "R"):
        sign = -1.0 if parts[-1][-1] == "R" else 1.0
        parts[-1] = parts[-1][:-1]
    return sign * math.radians(_dms(parts))


def deck_edges_at(deck_points, station):
    """(left, right) deck edge offsets at station, measurement type of the governing point.
    Parallel transitions hold the edge until the next point, linear and spline transitions
    are interpolated linearly (spline edges are approximate)."""
    if not deck_points:
        return None
    if station <= deck_points[0]["station"]:
        p = deck_points[0]
        return p["left"], p["right"], p["measurement"]
    for p0, p1 in zip(deck_points, deck_points[1:]):
        if station <= p1["station"]:
            u = (station - p0["station"]) / (p1["station"] - p0["station"]) if p1["station"] != p0["station"] else 0
            left = p0["left"] if p0["left_transition"] == 2 else p0["left"] + u * (p1["left"] - p0["left"])
            right = p0["right"] if p0["right_transition"] == 2 else p0["right"] + u * (p1["right"] - p0["right"])
            return left, right, p0["measurement"]
    p = deck_points[-1]
    return p["left"], p["right"], p["measurement"]


def expand_spacing(spacing):
    """list of spacings between adjacent girders from a GirderSpacing element"""
    gaps = []
    for group in spacing.findall("SpacingGroup"):
        first, last = int(_f(group, "FirstGirderIndex")), int(_f(group, "LastGirderIndex"))
        gaps.extend([_f(group, "Spacing")] * max(0, last - first))
    return gaps


def extract(path, stations=None):
    """Summary of the project at path. Alignment, profile, and deck edges are evaluated at
    stations; if stations is None, the pier stations and mid-span stations of this project are used."""
    root = ET.parse(path).getroot()
    bridge = root.find(".//BridgeDescription")
    values = {}

    def put(key, value, kind):
        values[key] = {"value": value, "kind": kind}

    # bridge
    put("bridge.girder_family", _t(bridge, "GirderFamilyName"), "text")
    put("bridge.girder_orientation", _t(bridge, "GirderOrientation"), "text")
    put("bridge.same_girder_for_entire_bridge", _t(bridge, "UseSameGirderForEntireBridge"), "text")
    put("bridge.girder_spacing_type", _t(bridge, "GirderSpacingType"), "text")
    if _t(bridge, "GirderSpacingType") == "0":  # uniform spacing, a single value for the bridge
        datum = f"type {_t(bridge, 'MeasurementType')}, location {_t(bridge, 'MeasurementLocation')}"
        put("bridge.girder_spacing", {"value": _f(bridge, "GirderSpacing"), "datum": datum}, "spacing")
    put("bridge.slab_offset_type", _t(bridge, "SlabOffsetType"), "text")
    put("bridge.slab_offset", _f(bridge, "SlabOffset"), "length_in")
    put("bridge.fillet", _f(bridge, "Fillet"), "length_in")

    # bearings: the bearing that applies at each pier face (first girder), whatever level it is defined at
    bearing_type = _t(bridge, "BearingType")
    bridge_bearing = bridge.find("BearingData2")

    def put_bearing(prefix, bd):
        if bd is None:
            return
        put(f"{prefix}.shape", _t(bd, "Shape"), "text")
        put(f"{prefix}.length", _f(bd, "Length"), "length_in")
        if _t(bd, "Shape") != "1":  # width isn't used for round bearings
            put(f"{prefix}.width", _f(bd, "Width"), "length_in")
        put(f"{prefix}.height", _f(bd, "Height"), "length_in")
        put(f"{prefix}.count", int(_f(bd, "BearingCount", 1)), "count")
        if int(_f(bd, "BearingCount", 1)) > 1:
            put(f"{prefix}.spacing", _f(bd, "Spacing"), "length_in")
        put(f"{prefix}.fixed_x", _t(bd, "FixedX"), "text")
        put(f"{prefix}.fixed_y", _t(bd, "FixedY"), "text")

    # piers
    piers = bridge.find("Piers").findall("PierDataDetails")
    put("bridge.pier_count", len(piers), "count")
    pier_stations = []
    for i, pier in enumerate(piers):
        station = _f(pier, "Station")
        pier_stations.append(station)
        put(f"pier[{i}].station", station, "station")
        put(f"pier[{i}].connection_type", _t(pier, "PierConnectionType"), "text")
        for face_tag, face in (("Back", "back"), ("Ahead", "ahead")):
            f = pier.find(face_tag)
            outer_abutment_face = (i == 0 and face == "back") or (i == len(piers) - 1 and face == "ahead")  # no girders
            if outer_abutment_face:
                pass
            elif bearing_type == "0":
                put_bearing(f"pier[{i}].{face}.bearing", bridge_bearing)
            elif f is not None and f.find("BearingDataVec") is not None:
                put_bearing(f"pier[{i}].{face}.bearing", f.find("BearingDataVec").find("BearingData2"))
            if f is None:
                continue
            put(f"pier[{i}].{face}.end_distance", _f(f, "GirderEndDistance"), "length_in")
            put(f"pier[{i}].{face}.end_distance_measure", _t(f, "EndDistanceMeasurementType"), "text")
            put(f"pier[{i}].{face}.bearing_offset", _f(f, "GirderBearingOffset"), "length_in")
            put(f"pier[{i}].{face}.bearing_offset_measure", _t(f, "BearingOffsetMeasurementType"), "text")
            spacing = f.find("GirderSpacing")
            if spacing is not None:
                gaps = expand_spacing(spacing)
                datum = f"type {_t(spacing, 'MeasurementType')}, location {_t(spacing, 'MeasurementLocation')}"
                for j, gap in enumerate(gaps):
                    put(f"pier[{i}].{face}.spacing[{j}]", {"value": gap, "datum": datum}, "spacing")

    # girders
    groups = bridge.find("GirderGroups").findall("GirderGroup")
    put("bridge.group_count", len(groups), "count")
    for i, group in enumerate(groups):
        put(f"group[{i}].piers", f"{_t(group, 'StartPier')}-{_t(group, 'EndPier')}", "text")
        girders = group.find("Girders").findall("Girder")
        put(f"group[{i}].girder_count", len(girders), "count")
        for j, girder in enumerate(girders):
            # with UseSameGirderForEntireBridge, the girder type is only stored for the bridge
            put(f"group[{i}].girder[{j}].type", _t(girder, "GirderType", _t(bridge, "Girder")), "text")
            segment = girder.find("PrecastSegment")
            concrete = segment.find(".//Concrete") if segment is not None else None
            # effective slab offset at the ends of the girder, wherever it is stored for the slab offset type
            for end, pier_idx, face, seg_tag in (("start", int(_t(group, "StartPier")), "Ahead", "StartSlabOffset"), ("end", int(_t(group, "EndPier")), "Back", "EndSlabOffset")):
                slab_offset_type = _t(bridge, "SlabOffsetType")
                if slab_offset_type == "0":
                    value = _f(bridge, "SlabOffset")
                elif slab_offset_type == "1":
                    value = _f(piers[pier_idx], f"{face}SlabOffset")
                else:
                    value = _f(segment, seg_tag)
                put(f"group[{i}].girder[{j}].slab_offset_{end}", value, "length_in")
            put(f"group[{i}].girder[{j}].fc", _f(concrete, "Fc"), "stress")
            put(f"group[{i}].girder[{j}].fci", _f(concrete, "Fci"), "stress")

    # deck
    deck = bridge.find("Deck")
    put("deck.type", _t(deck, "SlabType"), "text")
    put("deck.gross_depth", _f(deck, "GrossDepth"), "length_in")
    put("deck.left_overhang_edge_depth", _f(deck, "LeftOverhangEdgeDepth"), "length_in")
    put("deck.right_overhang_edge_depth", _f(deck, "RightOverhangEdgeDepth"), "length_in")
    put("deck.left_overhang_taper", _t(deck, "LeftOverhangTaperType"), "text")
    put("deck.right_overhang_taper", _t(deck, "RightOverhangTaperType"), "text")
    put("deck.haunch_shape", _t(deck, "HaunchShape"), "text")
    put("deck.fc", _f(deck.find("Concrete"), "Fc"), "stress")
    deck_points = [{
        "station": _f(p, "Station"), "left": _f(p, "LeftEdge"), "right": _f(p, "RightEdge"),
        "measurement": int(_f(p, "MeasurementType", 0)),
        "left_transition": int(_f(p, "LeftTransitionType", 0)), "right_transition": int(_f(p, "RightTransitionType", 0)),
    } for p in deck.find("DeckEdgePoints").findall("DeckPoint")]
    put("deck.edge_point_count", len(deck_points), "count")

    # railings
    for side in ("Left", "Right"):
        railing = bridge.find(f"{side}RailingSystem/RailingSystem")
        put(f"railing.{side.lower()}.exterior", _t(railing, "ExteriorRailingName"), "text")

    # evaluation stations
    if stations is None:
        stations = list(pier_stations)
        stations += [(a + b) / 2 for a, b in zip(pier_stations, pier_stations[1:])]
        stations = sorted(stations)

    alignment = HorizontalAlignment(root.find(".//AlignmentData"))

    # pier skew, so skews entered as a bearing compare with skews entered as an angle
    for i, pier in enumerate(piers):
        orientation = _t(pier, "Orientation")
        if not alignment.supported:
            put(f"pier[{i}].skew", "alignment has spirals", "unsupported")
        elif orientation:
            put(f"pier[{i}].skew", pier_skew(orientation, alignment.evaluate(pier_stations[i])[1]), "angle")
    profile = Profile(root.find(".//ProfileData"))
    alignment_offset = _f(bridge, "AlignmentOffset", 0.0)
    for k, station in enumerate(stations):
        tag = f"at[{k}]"  # the station is in meta["stations"][k]
        if alignment.supported:
            (x, y), direction = alignment.evaluate(station)
            put(f"{tag}.alignment.point", [x, y], "point")
            put(f"{tag}.alignment.direction", direction, "angle")
        else:
            put(f"{tag}.alignment.point", "alignment has spirals", "unsupported")
        elevation, grade = profile.evaluate(station)
        put(f"{tag}.profile.elevation", elevation, "length_ft")
        put(f"{tag}.profile.grade", grade, "grade")
        edges = deck_edges_at(deck_points, station)
        if edges:
            left, right, measurement = edges
            if measurement == 1:  # measured from CL bridge - convert to alignment
                left, right = left + alignment_offset, right - alignment_offset
            put(f"{tag}.deck.left_edge", left, "length_ft")
            put(f"{tag}.deck.right_edge", right, "length_ft")

    meta = {"file": str(path), "stations": stations, "pier_stations": pier_stations}
    return {"meta": meta, "values": values}


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("file")
    parser.add_argument("--out")
    args = parser.parse_args()
    summary = extract(args.file)
    text = json.dumps(summary, indent=1)
    if args.out:
        with open(args.out, "w", encoding="utf-8") as f:
            f.write(text)
    else:
        print(text)


if __name__ == "__main__":
    main()
