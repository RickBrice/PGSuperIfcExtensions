"""Compares an imported PGSuper project against expected values and writes a Markdown report.

Expected values come from the project the IFC model was exported from (round trip), or from a
hand-made expected file (e.g. values taken from the drawings). Values are compared in SI with
a tolerance for each kind of value and reported in US units.

Each value gets a status:
  match            within tolerance
  default match    within tolerance, but the imported value is also the template default, so the
                   match doesn't show that the importer set it (not counted as a match in the score)
  mismatch         outside tolerance
  not imported     the imported value is still the template default, but a different value is expected
  not comparable   the values are measured differently (e.g. girder spacing datum) or can't be evaluated
  missing          expected, but not in the imported project
  extra            in the imported project, but not expected (not counted in the score)

An expected file has the same form as a pgs_extract.py summary, or just {"values": {key: value}}
with values in SI. Only the keys it lists are compared.

Usage: python compare_bridge.py <imported.pgs> <expected.pgs|expected.json> [--template template.pgt] [--out report.md]
"""
import argparse
import json
import math
import re

import pgs_extract

FT = 0.3048
IN = 0.0254
KSI = 6894757.293168361

# kind: (tolerance in SI, formatter)
KINDS = {
    "station": (0.003, lambda v: _station(v)),
    "length_ft": (0.003, lambda v: f"{v / FT:.3f} ft"),
    "length_in": (0.003, lambda v: f"{v / IN:.3f} in"),
    "point": (0.003, lambda v: f"({v[0] / FT:.3f}, {v[1] / FT:.3f}) ft"),
    "angle": (1.0e-4, lambda v: f"{math.degrees(v):.4f} deg"),
    "grade": (1.0e-5, lambda v: f"{v * 100:.4f} %"),
    "stress": (0.05 * KSI, lambda v: f"{v / KSI:.3f} ksi"),
    "spacing": (0.003, lambda v: f"{v['value'] / FT:.3f} ft ({v['datum']})"),
    "count": (0, str),
    "text": (None, str),
    "unsupported": (None, str),
}


def _station(v):
    ft = v / FT
    sign = "-" if ft < 0 else ""
    ft = abs(ft)
    return f"{sign}{int(ft // 100)}+{ft % 100:05.2f}"


def fmt(value, kind):
    if value is None:
        return ""
    try:
        return KINDS.get(kind, (None, str))[1](value)
    except (TypeError, ValueError, KeyError, IndexError):
        return str(value)


def same(a, b, kind):
    """True if a and b are the same value within the tolerance for kind"""
    if a is None or b is None:
        return a is None and b is None
    tol = KINDS.get(kind, (None, None))[0]
    if kind == "point":
        return math.dist(a, b) <= tol
    if kind == "spacing":
        return a["datum"] == b["datum"] and abs(a["value"] - b["value"]) <= tol
    if kind == "angle":
        return abs(math.atan2(math.sin(a - b), math.cos(a - b))) <= tol
    if tol is None or isinstance(a, str) or isinstance(b, str):
        return str(a) == str(b)
    return abs(a - b) <= tol


def load(path, stations=None):
    """summary of a .pgs/.pgt file or an expected .json file"""
    if path.lower().endswith((".pgs", ".pgt")):
        return pgs_extract.extract(path, stations)
    with open(path, encoding="utf-8") as f:
        data = json.load(f)
    values = {}
    for key, v in data.get("values", {}).items():
        values[key] = v if isinstance(v, dict) and "value" in v else {"value": v, "kind": None}
    return {"meta": data.get("meta", {"file": path}), "values": values}


def generic(key):
    """key with its indices replaced by [*], so a value can be compared with the template's value for
    any pier, girder, etc. Values at stations (at[n]) keep their index - the template is evaluated at the same stations"""
    return re.sub(r"(?<!^at)\[\d+\]", "[*]", key)


def compare(imported, expected, template=None):
    """list of (key, status, expected value, imported value, kind, note)"""
    template_values = {}
    if template:
        for key, v in template["values"].items():
            template_values.setdefault(generic(key), []).append(v["value"])

    rows = []
    for key, exp in expected["values"].items():
        imp = imported["values"].get(key)
        kind = exp["kind"] or (imp["kind"] if imp else None)
        e = exp["value"]
        if imp is None:
            rows.append((key, "missing", e, None, kind, ""))
            continue
        i = imp["value"]
        if kind == "unsupported" or imp["kind"] == "unsupported":
            rows.append((key, "not comparable", e, i, kind, "not supported by pgs_extract"))
        elif kind == "spacing" and e["datum"] != i["datum"]:
            rows.append((key, "not comparable", e, i, kind, "measured differently"))
        elif same(e, i, kind):
            is_default = any(same(t, i, kind) for t in template_values.get(generic(key), []))
            rows.append((key, "default match" if is_default else "match", e, i, kind, "template default" if is_default else ""))
        elif any(same(t, i, kind) for t in template_values.get(generic(key), [])):
            rows.append((key, "not imported", e, i, kind, "template default"))
        else:
            rows.append((key, "mismatch", e, i, kind, ""))

    for key, imp in imported["values"].items():
        if key not in expected["values"]:
            rows.append((key, "extra", None, imp["value"], imp["kind"], ""))
    return rows


def inventory(imported, template):
    """For a model without expected values: rows of (key, status, template value, imported value, kind, note)
    where status is "imported" if the import changed the value from the template default, else "default".
    Values with a count ("[n]") the template doesn't have are compared with the template's other instances."""
    template_values = {}
    for key, v in template["values"].items():
        template_values.setdefault(generic(key), []).append(v["value"])

    rows = []
    for key, imp in imported["values"].items():
        defaults = template_values.get(generic(key), [])
        is_default = any(same(t, imp["value"], imp["kind"]) for t in defaults)
        rows.append((key, "default" if is_default else "imported", defaults[0] if defaults else None, imp["value"], imp["kind"], ""))
    return rows


def inventory_report(rows, title, imported_file, stations=None):
    imported_count = sum(1 for r in rows if r[1] == "imported")
    lines = [f"# {title}", "", f"- Imported: `{imported_file}`", "- No expected values - this shows which values the import set", ""]
    if stations:
        lines += ["Alignment, profile, and deck edges are evaluated at " + ", ".join(f"at[{k}] = {_station(s)}" for k, s in enumerate(stations)), ""]
    lines += [f"**{imported_count} of {len(rows)} values set by the import**, {len(rows) - imported_count} are template defaults", ""]
    for status, heading in (("default", "template default"), ("imported", "set by the import")):
        group = [r for r in rows if r[1] == status]
        if not group:
            continue
        lines += [f"## {heading} ({len(group)})", "", "| Value | Template | Imported |", "|---|---|---|"]
        for key, _, t, i, kind, _ in group:
            lines.append(f"| `{key}` | {fmt(t, kind)} | {fmt(i, kind)} |")
        lines.append("")
    return "\n".join(lines)


STATUS_ORDER = ["mismatch", "not imported", "missing", "not comparable", "default match", "match", "extra"]


def score(rows):
    counts = {s: 0 for s in STATUS_ORDER}
    for row in rows:
        counts[row[1]] += 1
    compared = sum(counts[s] for s in STATUS_ORDER if s != "extra")
    return counts, compared


def report(rows, title, imported_file, expected_file, show_extra=False, stations=None):
    counts, compared = score(rows)
    lines = [f"# {title}", "",
             f"- Imported: `{imported_file}`",
             f"- Expected: `{expected_file}`", ""]
    if stations:
        lines += ["Alignment, profile, and deck edges are compared at " + ", ".join(f"at[{k}] = {_station(s)}" for k, s in enumerate(stations)), ""]
    lines += [
             f"**{counts['match']} of {compared} match** - " + ", ".join(f"{s}: {counts[s]}" for s in STATUS_ORDER if s != "match"), "",
             "A *default match* equals the expected value only because the template has the same value.", ""]
    for status in STATUS_ORDER:
        group = [r for r in rows if r[1] == status]
        if not group or (status == "extra" and not show_extra):
            continue
        lines += [f"## {status} ({len(group)})", "", "| Value | Expected | Imported | Note |", "|---|---|---|---|"]
        for key, _, e, i, kind, note in group:
            lines.append(f"| `{key}` | {fmt(e, kind)} | {fmt(i, kind)} | {note} |")
        lines.append("")
    return "\n".join(lines)


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("imported")
    parser.add_argument("expected")
    parser.add_argument("--template", help="template the project was imported into, to detect values that weren't imported")
    parser.add_argument("--title", default="IFC import comparison")
    parser.add_argument("--show-extra", action="store_true")
    parser.add_argument("--out")
    args = parser.parse_args()

    expected = load(args.expected)
    stations = expected["meta"].get("stations")
    imported = load(args.imported, stations)
    template = load(args.template, stations) if args.template else None
    text = report(compare(imported, expected, template), args.title, args.imported, args.expected, args.show_extra, stations or imported["meta"]["stations"])
    if args.out:
        with open(args.out, "w", encoding="utf-8") as f:
            f.write(text)
    else:
        print(text)


if __name__ == "__main__":
    main()
