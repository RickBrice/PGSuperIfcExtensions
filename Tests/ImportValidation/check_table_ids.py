"""Round trip of a mapping table through a general IDS (devdocs/MappingTablesDesign.md, M5).

The standard table is written as a general IDS (/IfcTableToIds), a standalone table is generated back
from that IDS (/IfcIdsToTable with /IfcTableExtends=none), and the declarations of the two tables are compared: property sets, quantity
sets, and classifications of every element role. Material property sets aren't compared - a general
IDS can't hold them (IDS property facets check an element's own and type property sets).

Writes results/Standard.ids, results/Standard-from-ids.json, and their logs.

Usage: python check_table_ids.py [--bridgelink BridgeLink.exe]  (run_validation.py runs it too)
"""
import argparse
import json
import os
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
RESULTS = HERE / "results"
STANDARD = HERE.parent.parent / "MappingTables" / "Standard.json"
TEMPLATE = HERE.parent.parent / "IfcImportTemplate.pgt"
BSDD = "https://identifier.buildingsmart.org/uri/usTransportation/usBridge/1/"


def resolve(uri, name, kind):
    if not uri:
        return ""
    if uri == "bsdd":
        return f"{BSDD}{kind}/{name}"
    if uri.startswith("bsdd:"):
        return f"{BSDD}{kind}/{uri[5:]}"
    return uri


def roles(applies_to):
    return applies_to if isinstance(applies_to, list) else [applies_to]


def canonical(table):
    """set of lines describing what the table exports, per element role"""
    lines = set()
    for key, items_key in (("property_sets", "properties"), ("quantity_sets", "quantities")):
        for p in table.get(key, []):
            attach = p.get("attach", "occurrence")
            if attach == "material" or p.get("remove"):
                continue
            quantities = key == "quantity_sets"
            head = f'{p["name"]} attach={attach} condition={p.get("condition", "classify")} uri={resolve(p.get("uri", ""), p["name"], "class")}'
            if quantities:
                head += f' method={p.get("method", "")}'
            for role in roles(p["applies_to"]):
                for n, prop in enumerate(p[items_key]):
                    enumeration = prop.get("enumeration")
                    enum_text = f'{enumeration["name"]}{enumeration["values"]}' if enumeration else ""
                    value = prop.get("value", "")
                    lines.add(f'{role} | {head} | {n}: {prop["name"]} type={prop["type"].upper()} target={prop.get("target", "")} '
                              f'import={prop.get("import", not quantities) if prop.get("target") else ""} value={json.dumps(value)} '
                              f'enumeration={enum_text} uri={resolve(prop.get("uri", ""), prop["name"], "prop") if not quantities else ""}')
    for c in table.get("classifications", []):
        for role in roles(c["applies_to"]):
            location = c.get("location", "")
            if location == "bsdd":
                location = f'{BSDD}class/{c["identification"]}'
            lines.add(f'{role} | classification {c["system"]} {c["identification"]} name={c.get("name", "")} location={location}')
    return lines


def run(exe, arguments):
    command = f'"{exe}" ' + " ".join(arguments)
    return subprocess.run(command, timeout=600).returncode


def check(exe):
    """runs the round trip. Returns (ok, message)"""
    ids = RESULTS / "Standard.ids"
    generated = RESULTS / "Standard-from-ids.json"
    for f in (ids, generated):
        f.unlink(missing_ok=True)

    run(exe, [f'"/IfcTableToIds={ids}"', f'"{TEMPLATE}"'])
    if not ids.exists():
        return False, "the IDS wasn't written (see results/Standard.ids.log)"
    run(exe, [f'"/IfcIdsToTable={ids}"', f'"{TEMPLATE}"', f'"/IfcTable={generated}"', "/IfcTableExtends=none"])
    if not generated.exists():
        return False, "the table wasn't generated (see results/Standard-from-ids.json.log)"

    log = (RESULTS / "Standard-from-ids.json.log").read_text(errors="replace")
    if "The generated table is valid." not in log:
        return False, "the generated table isn't valid (see results/Standard-from-ids.json.log)"

    standard = canonical(json.loads(STANDARD.read_text(encoding="utf-8")))
    back = canonical(json.loads(generated.read_text(encoding="utf-8")))
    missing = sorted(standard - back)
    extra = sorted(back - standard)
    diff_file = RESULTS / "Standard-from-ids.diff.txt"
    diff_file.unlink(missing_ok=True)
    if missing or extra:
        diff_file.write_text("\n".join([f"- {l}" for l in missing] + [f"+ {l}" for l in extra]) + "\n", encoding="utf-8")
        return False, f"{len(missing)} declarations lost, {len(extra)} changed or added (see {diff_file.name})"
    return True, f"same: {len(standard)} declarations (material property sets not compared)"


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--bridgelink", default=str(Path(os.environ.get("ARPDIR", ".")) / "BridgeLink" / "RegFreeCOM" / "x64" / "Release" / "BridgeLink.exe"))
    args = parser.parse_args()
    RESULTS.mkdir(exist_ok=True)
    ok, message = check(args.bridgelink)
    print(("OK: " if ok else "FAILED: ") + message)
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
