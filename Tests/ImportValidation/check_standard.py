"""Checks the exports and the standard mapping table against the standards in standards.json.

IDS check: an export is checked with each of the standard's IDS files (ifctester). The failures are
written to results/<run>.<standard>-ids.txt, grouped by IDS file, specification, requirement, entity,
and reason, with the number of elements and a few of their names. The file is sorted, so the change from
one run to the next shows up as a difference in the repository.

Workbook check: the standard mapping table is compared with the standard's data dictionary workbook,
role by role. A role is matched to a workbook class by its classification (usBridge_GirderPrecastConcrete
is the class GirderPrecastConcrete). The findings go to results/Standard-vs-<standard>-workbook.txt:
  CLASS       the role's entity, predefined type, or ObjectType differs from the class, or the class isn't in the workbook
  REQUIRED    a required property of a property set the role declares isn't declared (an IDS failure)
  NO-VALUE    a required property is declared without a target or a value, so it's exported without a value
  VALUE       a constant isn't one of the workbook's enumeration values (IDS values are case sensitive)
  TYPE        the data type differs
  NAME        a property name differs only in case
  ENUM        the enumeration values differ
  RESTRICTED  the workbook limits the property set to some classes ("Applies only to ...") and the role doesn't
              declare it: check whether it applies
  PSET        a property set of the class isn't declared for the role. IDS specifications only check the
              property sets an element has, and the workbook doesn't say which ones an element must have
  OPTIONAL    an optional property of a property set the role declares isn't declared
  EXTRA       the table declares a property or property set the class doesn't have
  NO-CLASS    the role has no classification of the standard, so it isn't compared

Usage:
  python check_standard.py --ifc results/PGSuper-Display.ifc [--standard usBridge]
  python check_standard.py --workbook [--standard usBridge]

ifctester: the one in %IFCOPENSHELL_DIR%\\src\\ifctester is used if it's there (the ifctester installed
with the local IfcOpenShell build's Python is older than its ifcopenshell), otherwise the installed one.
"""
import argparse
import json
import os
import re
import sys
import time
from collections import defaultdict
from pathlib import Path

HERE = Path(__file__).resolve().parent
RESULTS = HERE / "results"
STANDARDS = HERE / "standards.json"
STANDARD_TABLE = HERE.parent.parent / "MappingTables" / "Standard.ifcmap"

CATEGORIES = ["CLASS", "REQUIRED", "NO-VALUE", "VALUE", "TYPE", "NAME", "ENUM", "PSET", "RESTRICTED", "OPTIONAL", "EXTRA", "NO-CLASS"]
INFORMATIONAL = ("PSET", "RESTRICTED", "OPTIONAL", "EXTRA", "NO-CLASS")


def load_standards(names=None):
    """the standards of standards.json (all, or the named ones), with "path" set to their folder"""
    standards = json.loads(STANDARDS.read_text(encoding="utf-8"))["standards"]
    for s in standards:
        s["path"] = Path(os.path.expandvars(s["folder"]))
    return [s for s in standards if not names or s["name"] in names]


_ifctester = None


def ifctester_ids():
    """the ifctester.ids module (see the module notes). Raises ImportError if there isn't one"""
    global _ifctester
    if _ifctester is None:
        source = Path(os.environ.get("IFCOPENSHELL_DIR", "")) / "src" / "ifctester"
        if os.environ.get("IFCOPENSHELL_DIR") and (source / "ifctester").is_dir():
            sys.path.insert(0, str(source))
        import ifctester.ids
        _ifctester = ifctester.ids
    return _ifctester


def unavailable(standard):
    """why the standard's files can't be used, or an empty string"""
    if not standard["path"].is_dir():
        return f"folder not found: {standard['path']}"
    return ""


# ---------------------------------------------------------------------------------------------------
# IDS check

def _element_name(element):
    name = getattr(element, "Name", None)
    return name if name else f"#{element.id()}"


def _requirement_text(requirement, specification):
    try:
        text = requirement.to_string("requirement", specification, requirement)
    except Exception:
        text = type(requirement).__name__
    return re.sub(r" at 0x[0-9A-Fa-f]+", "", text)  # no object addresses, so the file is the same from run to run


def check_ids(standard, ifc_path, out_path, samples=3):
    """checks an export with the standard's IDS files and writes the failures. Returns (failed, total) specifications"""
    import ifcopenshell
    ids = ifctester_ids()
    ifc = ifcopenshell.open(str(ifc_path))
    lines = []
    header = []
    failed_total = specs_total = 0
    for ids_file in sorted(standard["path"].glob(standard["ids"])):
        specs = ids.open(str(ids_file))
        specs.validate(ifc)
        failed = 0
        for spec in specs.specifications:
            groups = defaultdict(list)
            for requirement in spec.requirements:
                for failure in getattr(requirement, "failures", None) or []:
                    element = failure["element"]
                    key = (_requirement_text(requirement, spec), element.is_a(), str(failure.get("reason") or ""))
                    groups[key].append(_element_name(element))
            if not groups and spec.status is False:
                applicable = len(getattr(spec, "applicable_entities", []) or [])
                groups[("the specification's cardinality", "", f"{applicable} applicable elements")] = []
            if groups:
                failed += 1
            for (requirement, entity, reason), names in sorted(groups.items()):
                shown = ", ".join(f'"{n}"' for n in sorted(set(names))[:samples])
                more = f" (and {len(set(names)) - samples} more)" if len(set(names)) > samples else ""
                count = f"{len(names)} elements: {shown}{more}" if names else ""
                lines.append(f"{ids_file.name} | {spec.name} | {requirement} | {entity} | {reason} | {count}")
        header.append(f"{ids_file.name}: {failed} of {len(specs.specifications)} specifications fail")
        failed_total += failed
        specs_total += len(specs.specifications)
    text = [f"{standard['name']} IDS check of {Path(ifc_path).name}", ""] + header + [""] + sorted(lines)
    out_path.write_text("\n".join(text) + "\n", encoding="utf-8")
    return failed_total, specs_total


# ---------------------------------------------------------------------------------------------------
# Workbook check

def _text(cell):
    return cell.strip() if isinstance(cell, str) else ("" if cell is None else str(cell))


def _parse_ifc_name(text):
    """("IfcBeam", {"BEAM"}, None) from "IfcBeam.Beam"; a USERDEFINED type gives its ObjectType"""
    text = " ".join(text.split())
    m = re.match(r"(Ifc\w+)\s*\.?\s*(.*)", text)
    if not m:
        return text, set(), None
    entity, rest = m.group(1), m.group(2).strip()
    if rest.upper().startswith("USERDEFINED"):
        return entity, {"USERDEFINED"}, rest.split(":", 1)[1].strip() if ":" in rest else None
    words = {w.upper() for w in re.findall(r"[A-Za-z_]+", rest) if w.lower() != "or"}
    if words & {"PREDEFINEDTYPE", "PREDEFINED", "TYPE"}:
        words = set()  # any predefined type
    return entity, words, None


def read_workbook(path):
    """(classes, sheets): classes[dd name] = list of {sheet, industry, entity, predefined, object_type};
    sheets[sheet] = list of property sets {name, note, properties: [{name, cardinality, type, property_type, enumeration}]}"""
    import openpyxl
    wb = openpyxl.load_workbook(str(path), read_only=True, data_only=True)
    classes = defaultdict(list)
    sheets = {}
    for ws in wb.worksheets:
        section = None
        cols = None
        psets = []
        current = None
        for row in ws.iter_rows(values_only=True):
            cells = [_text(c) for c in row]
            first = next((c for c in cells if c), "")
            if not first:
                continue
            if first == "Entity":
                section, cols = "entity", None
                continue
            if first == "IFC Attributes":
                section, cols = "attributes", None
                continue
            if first.endswith("Property Sets and Properties"):
                section, cols, current = "psets", None, None
                continue
            if section == "entity":
                if cols is None:
                    if "IFC Name" in cells:
                        cols = {c: i for i, c in enumerate(cells) if c}
                    continue
                dd = cells[cols["usBridge DD Name"]] if "usBridge DD Name" in cols else ""
                ifc_name = cells[cols["IFC Name"]]
                if not dd or not ifc_name or dd == "Not Included":
                    continue
                entity, predefined, object_type = _parse_ifc_name(ifc_name)
                classes[dd].append({"sheet": ws.title, "industry": cells[cols["Industry Name"]], "entity": entity,
                                    "predefined": predefined, "object_type": object_type, "ifc_name": " ".join(ifc_name.split())})
            elif section == "psets":
                if cols is None:
                    if "Property Set Name" in cells:
                        cols = {c: i for i, c in enumerate(cells) if c}
                    continue
                pset_cell = cells[cols["Property Set Name"]]
                if pset_cell:
                    if pset_cell == "None":
                        current = None
                        continue
                    # text after the name limits the property set to some of the sheet's classes, e.g. "(Applies only to Pier Cap)"
                    name_line, _, note = pset_cell.partition("\n")
                    note = " ".join(note.split())
                    if note.startswith("(") and note.endswith(")"):
                        note = note[1:-1]
                    current = {"name": name_line.strip(), "note": note, "properties": []}
                    psets.append(current)
                name = cells[cols["Property Name"]]
                if current is None or not name:
                    continue
                property_type = cells[cols["PropertyType"]] if "PropertyType" in cols else ""
                enumeration = cells[cols["Enumeration"]] if "Enumeration" in cols else ""
                values = [v.strip() for v in re.split(r";|\n", enumeration) if v.strip()] if "Enumerated" in property_type else []
                current["properties"].append({"name": name, "cardinality": cells[cols["IDS Cardinality"]], "type": cells[cols["Data Type"]],
                                              "property_type": property_type, "enumeration": values})
        sheets[ws.title] = psets
    return classes, sheets


def _roles(applies_to):
    return applies_to if isinstance(applies_to, list) else [applies_to]


def _declared(table):
    """declared[role][(pset, property)] = the table's declaration"""
    declared = defaultdict(dict)
    for key, items in (("property_sets", "properties"), ("quantity_sets", "quantities")):
        for p in table.get(key, []):
            if p.get("remove"):
                continue
            for role in _roles(p["applies_to"]):
                for prop in p[items]:
                    enumeration = prop.get("enumeration") or {}
                    declared[role][(p["name"], prop["name"])] = {
                        "type": prop["type"], "target": prop.get("target"), "value": prop.get("value"),
                        "enumeration": enumeration.get("values"), "attach": p.get("attach", "occurrence")}
    return declared


def compare_workbook(standard, table):
    """findings (category, role, item, detail) of the table compared with the standard's workbook"""
    classes, sheets = read_workbook(standard["path"] / standard["workbook"])
    material_psets = sheets.get("IfcMaterial", [])
    prefix = standard.get("classification_prefix", "")
    role_class = {}
    for c in table.get("classifications", []):
        if c.get("system") == standard["classification_system"]:
            for role in _roles(c["applies_to"]):
                role_class[role] = c["identification"].removeprefix(prefix)
    declared = _declared(table)
    findings = []

    for role, selector in table.get("elements", {}).items():
        dd = role_class.get(role)
        if not dd:
            findings.append(("NO-CLASS", role, "", f"{selector.get('entity', '')}: no {standard['classification_system']} classification in the table"))
            continue
        label = f"{role} ({dd})"
        candidates = classes.get(dd, [])
        if not candidates:
            findings.append(("CLASS", label, "", f"the class {dd} isn't in the workbook"))
            continue
        cls = next((c for c in candidates if c["entity"].lower() == selector.get("entity", "").lower()), candidates[0])

        # the IFC class
        entity = selector.get("entity", "")
        predefined = (selector.get("predefined_type") or "").upper()
        object_type = (selector.get("attributes") or {}).get("ObjectType")
        if entity.lower() != cls["entity"].lower():
            findings.append(("CLASS", label, "entity", f"table {entity}, workbook {cls['ifc_name']}"))
        if cls["predefined"] and predefined not in cls["predefined"]:
            findings.append(("CLASS", label, "predefined type", f"table {predefined or '(none)'}, workbook {cls['ifc_name']}"))
        if cls["object_type"] and object_type != cls["object_type"]:
            findings.append(("CLASS", label, "ObjectType", f"table {object_type or '(none)'}, workbook {cls['ifc_name']}"))

        # the properties
        mine = declared.get(role, {})
        by_lower = {(p.lower(), n.lower()): (p, n) for p, n in mine}
        compared = set()
        workbook_psets = [(p, False) for p in sheets.get(cls["sheet"], [])]
        workbook_psets += [(p, True) for p in material_psets if any(k[0] == p["name"] and d["attach"] == "material" for k, d in mine.items())]
        known_psets = {p["name"].lower() for p, _ in workbook_psets}
        for pset, _ in workbook_psets:
            if not any(k[0].lower() == pset["name"].lower() for k in mine):
                if pset["note"]:
                    findings.append(("RESTRICTED", label, pset["name"], f"not declared; the workbook says \"{pset['note']}\""))
                else:
                    n = sum(1 for wp in pset["properties"] if wp["cardinality"].lower() == "required")
                    findings.append(("PSET", label, pset["name"], f"not declared ({len(pset['properties'])} properties, {n} required)"))
                continue
            for wp in pset["properties"]:
                key = by_lower.get((pset["name"].lower(), wp["name"].lower()))
                item = f"{pset['name']}.{wp['name']}"
                required = wp["cardinality"].lower() == "required"
                if not key:
                    findings.append(("REQUIRED" if required else "OPTIONAL", label, item, f"{wp['type']} not declared"))
                    continue
                compared.add(key)
                d = mine[key]
                if key != (pset["name"], wp["name"]):
                    findings.append(("NAME", label, item, f"the table has {key[0]}.{key[1]}"))
                if required and not d["target"] and d["value"] is None:
                    findings.append(("NO-VALUE", label, item, "required, but declared without a target or a value"))
                enumerated = "Enumerated" in wp["property_type"] or wp["type"].lower().startswith("penum")
                if not enumerated and wp["type"] and d["type"].lower() != wp["type"].lower():
                    findings.append(("TYPE", label, item, f"table {d['type']}, workbook {wp['type']}"))
                if wp["enumeration"]:
                    allowed = wp["enumeration"]
                    if d["value"] is not None and str(d["value"]) not in allowed:
                        case = [v for v in allowed if v.lower() == str(d["value"]).lower()]
                        hint = f" (the workbook has \"{case[0]}\")" if case else ""
                        findings.append(("VALUE", label, item, f"\"{d['value']}\" isn't one of {allowed}{hint}"))
                    if d["enumeration"] and d["enumeration"] != allowed:
                        findings.append(("ENUM", label, item, f"table {d['enumeration']}, workbook {allowed}"))
        for key, d in sorted(mine.items()):
            if key in compared:
                continue
            if key[0].lower() in known_psets:
                findings.append(("EXTRA", label, f"{key[0]}.{key[1]}", "not in the workbook's property set for this class"))
        for pset_name in sorted({k[0] for k in mine if k[0].lower() not in known_psets}):
            findings.append(("EXTRA", label, pset_name, "property set not in the workbook for this class"))
    return findings


def write_workbook_findings(standard, findings, out_path):
    counts = {c: sum(1 for f in findings if f[0] == c) for c in CATEGORIES}
    order = {c: i for i, c in enumerate(CATEGORIES)}
    lines = [f"Standard mapping table compared with the {standard['name']} workbook {standard['workbook']} (categories: see check_standard.py)", "",
             "Counts: " + ", ".join(f"{c} {n}" for c, n in counts.items() if n), ""]
    for category, role, item, detail in sorted(findings, key=lambda f: (order[f[0]], f[1], f[2], f[3])):
        lines.append(f"{category} | {role} | {item} | {detail}")
    out_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return counts


def workbook_message(counts, out_name):
    shown = [f"{n} {c}" for c, n in counts.items() if n and c not in INFORMATIONAL]
    other = [f"{n} {c}" for c, n in counts.items() if n and c in INFORMATIONAL]
    return (", ".join(shown) or "no differences") + (f" ({', '.join(other)})" if other else "") + f" (see {out_name})"


def check_workbook(standard, table_path=STANDARD_TABLE):
    """compares the table with the standard's workbook. Returns (ok, message)"""
    problem = unavailable(standard)
    if problem:
        return False, f"skipped: {problem}"
    table = json.loads(Path(table_path).read_text(encoding="utf-8"))
    out = RESULTS / f"Standard-vs-{standard['name']}-workbook.txt"
    counts = write_workbook_findings(standard, compare_workbook(standard, table), out)
    ok = not any(counts[c] for c in ("CLASS", "REQUIRED", "NO-VALUE", "VALUE", "TYPE", "NAME"))
    return ok, workbook_message(counts, out.name)


def ids_message(standard, ifc_path, run):
    """checks an export with the standard's IDS files. Returns a message for the summary"""
    problem = unavailable(standard)
    if problem:
        return f", {standard['name']} IDS skipped: {problem}"
    try:
        ifctester_ids()
    except ImportError as e:
        return f", {standard['name']} IDS skipped: ifctester not available ({e})"
    out = RESULTS / f"{run}.{standard['name']}-ids.txt"
    failed, total = check_ids(standard, ifc_path, out)
    if failed == 0:
        return f", {standard['name']} IDS: all {total} specifications pass"
    return f", {standard['name']} IDS: {failed} of {total} specifications fail (see {out.name})"


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--ifc", nargs="*", default=[], help="exports to check with the standard's IDS files")
    parser.add_argument("--workbook", action="store_true", help="compare the standard mapping table with the standard's workbook")
    parser.add_argument("--standard", nargs="*", help="names of the standards (default: all in standards.json)")
    args = parser.parse_args()
    if not args.ifc and not args.workbook:
        parser.error("give --ifc or --workbook")
    RESULTS.mkdir(exist_ok=True)
    for standard in load_standards(args.standard):
        if args.workbook:
            ok, message = check_workbook(standard)
            print(f"Standard table vs {standard['name']} workbook: {message}")
        for ifc in args.ifc:
            start = time.time()
            run = Path(ifc).stem
            print(f"{run}{ids_message(standard, Path(ifc), run)} ({time.time() - start:.0f} s)")


if __name__ == "__main__":
    main()
