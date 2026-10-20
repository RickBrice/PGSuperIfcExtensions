"""Runs the IFC import validation for the models in models.json and writes reports.

For a model with a "pgs" file (round trip):
   export the project to IFC (once for each "property_units" mode), import the IFC into a new
   project, and compare the imported project to the original project. "remove_from_export" lists IFC
   classes removed from the export before the import, with what they nest (e.g. "IfcAlignmentVertical",
   so the import measures the profile from the deck); "tolerances" replaces the comparison tolerance
   (SI) for the keys matching each pattern (e.g. "at[*].profile.grade", for values measured from geometry).
   The strands are compared physically (ifc_strands.py): the imported project is exported again
   (results/reexport/, not tracked), and its strands are compared with the strands of the first export.
   The template's girders have no strands, so a strand count of 0 is reported as the template default.
For a model with an "ifc" file (.ifc, or .zip containing one .ifc):
   import the IFC into a new project and compare to expected/<name>.json if it exists,
   otherwise report which values the import set (inventory).

A model can name a "mapping" table (/IfcMapping); without one, the configured mapping table (File >
Configure BridgeLink) is used, or the standard table installed with the extension if none is configured.
The committed results are made without a setting; the summary says when one is set. "expect_log" lists text the import log must contain (e.g. mapping table hints);
the summary reports any that is missing.

The table-driven content of each round-trip export (property sets, quantity sets, classifications,
element names; see compare_export.py) is compared with export_baseline/<run>.txt.gz, and the summary
reports any difference, with the details in results/<run>.export-diff.txt. The design-value IDS
written with each round-trip export (/IfcIds) is compared with export_baseline/<run>.ids.txt the same
way (see compare_ids.py). --update-export-baseline saves the current exports and IDS files as the
baseline instead (after an intended change to the export).

The standard table's round trip through a general IDS (/IfcTableToIds, /IfcIdsToTable; see
check_table_ids.py) runs after the models, as the "Standard table <-> IDS" row, unless --models is given.

Standards (standards.json, check_standard.py): each round-trip export is checked with the IDS files of
each standard (results/<run>.<standard>-ids.txt), and, unless --models is given, the standard table is
compared with each standard's data dictionary workbook (results/Standard-vs-<standard>-workbook.txt, one
row per standard). --skip-standards leaves these checks out.

Reports, logs, and imported projects are written to results/, replacing the previous results, so
the change from one run to the next shows up as a difference in the repository. summary.md has the
score of every run. Exported IFC files and unzipped models are written to results/ but not tracked.

Paths in models.json are relative to models.json. BridgeLink is the Release build in %ARPDIR%
unless --bridgelink is given.

--configuration switches the PGSuper configuration (e.g. "WSDOT (Local)" for exact girder
name matches) for the duration of the run and switches back to the current configuration
when done.

The girder label format (GirderLabelFormat in PGSuper's Settings key, "Alpha" by default) is set to
Numeric for the run, as the baselines have, and restored when done.

Usage: python run_validation.py [--models PGSuper Iowa] [--configuration "WSDOT (Local)"] [--bridgelink BridgeLink.exe] [--skip-standards]
"""
import argparse
import datetime
import json
import os
import subprocess
import time
import winreg
import zipfile
from pathlib import Path
from xml.etree import ElementTree

import compare_bridge
import compare_export
import compare_ids
import check_standard
import check_table_ids
import ifc_strands

HERE = Path(__file__).resolve().parent
RESULTS = HERE / "results"
EXPORT_BASELINE = HERE / "export_baseline"
TIMEOUT = 600  # seconds for one BridgeLink run
PGSUPER_OPTIONS_KEY = r"Software\Washington State Department of Transportation\PGSuper\Options"
PGSUPER_SETTINGS_KEY = r"Software\Washington State Department of Transportation\PGSuper\Settings"
GIRDER_LABEL_FORMAT = "Numeric"  # girder names in the exports and baselines are "Girder 1", ...
IFC_EXTENSIONS_KEY = r"Software\Washington State Department of Transportation\BridgeLink\IfcExtensions"


def rel(path):
    """path relative to this folder, for reports"""
    try:
        return Path(path).resolve().relative_to(HERE).as_posix()
    except ValueError:
        return str(path)


def quote(text):
    return f'"{text}"'


def run_bridgelink(exe, arguments):
    """Runs BridgeLink with a command line (list of already quoted arguments). Returns (ok, message)."""
    command = f'"{exe}" ' + " ".join(arguments)
    start = time.time()
    try:
        # BridgeLink is a GUI application - the command line is passed as is so quoting is under our control
        result = subprocess.run(command, timeout=TIMEOUT)
    except subprocess.TimeoutExpired:
        return False, f"timed out after {TIMEOUT} s (a dialog may be waiting for input)"
    elapsed = time.time() - start
    if result.returncode != 0:
        return False, f"exit code {result.returncode:#x} after {elapsed:.0f} s"
    return True, f"{elapsed:.0f} s"


def get_configuration():
    with winreg.OpenKey(winreg.HKEY_CURRENT_USER, PGSUPER_OPTIONS_KEY) as key:
        server = winreg.QueryValueEx(key, "CatalogServer2")[0]
        publisher = winreg.QueryValueEx(key, "Publisher2")[0]
    return server, publisher


def get_mapping_table_setting():
    """the configured mapping table (File > Configure BridgeLink), or an empty string"""
    try:
        with winreg.OpenKey(winreg.HKEY_CURRENT_USER, IFC_EXTENSIONS_KEY) as key:
            return winreg.QueryValueEx(key, "MappingTable")[0].strip()
    except OSError:
        return ""


def get_girder_label_format():
    """PGSuper's girder label format ("Alpha" or "Numeric"), or None when it isn't set (PGSuper's default is Alpha)"""
    try:
        with winreg.OpenKey(winreg.HKEY_CURRENT_USER, PGSUPER_SETTINGS_KEY) as key:
            return winreg.QueryValueEx(key, "GirderLabelFormat")[0]
    except OSError:
        return None


def set_girder_label_format(label_format):
    """sets the girder label format, or removes the setting when label_format is None"""
    with winreg.CreateKey(winreg.HKEY_CURRENT_USER, PGSUPER_SETTINGS_KEY) as key:
        if label_format is None:
            try:
                winreg.DeleteValue(key, "GirderLabelFormat")
            except OSError:
                pass
        else:
            winreg.SetValueEx(key, "GirderLabelFormat", 0, winreg.REG_SZ, label_format)


def set_configuration(exe, server, publisher):
    ok, message = run_bridgelink(exe, ["/App=PGSuper", f'/Configuration="{server}":"{publisher}"'])
    if not ok or get_configuration() != (server, publisher):
        raise RuntimeError(f'Unable to set configuration "{server}":"{publisher}" ({message})')


def read_log(path):
    try:
        return Path(path).read_text(errors="replace")
    except OSError:
        return ""


def unpack(path):
    """the .ifc file for a model - zipped models are extracted to results/models"""
    path = Path(path)
    if path.suffix.lower() != ".zip":
        return path
    with zipfile.ZipFile(path) as z:
        name = next(n for n in z.namelist() if n.lower().endswith(".ifc"))
        target = RESULTS / "models" / name
        if not target.exists() or target.stat().st_mtime < path.stat().st_mtime:
            z.extract(name, RESULTS / "models")
    return target


def import_ifc(exe, ifc, template, pgs, mapping=None):
    """imports ifc into a new project pgs, with a mapping table if given. Returns (ok, message)"""
    arguments = [quote(f"/IfcImport={ifc}"), quote(template), quote(f"/IfcOut={pgs}")]
    if mapping:
        arguments.append(quote(f"/IfcMapping={mapping}"))
    ok, message = run_bridgelink(exe, arguments)
    if not ok or not pgs.exists():
        return False, f"import failed: {message}"
    if "IFC import succeeded" not in read_log(pgs.with_suffix(".log")):
        message += ", import reported errors - see log"
    return True, message


def check_log(model, pgs):
    """message about the text of the model's "expect_log" that is missing from the import log"""
    expected = model.get("expect_log", [])
    if not expected:
        return ""
    log = read_log(pgs.with_suffix(".log"))
    missing = [text for text in expected if text not in log]
    if missing:
        return ", log check failed: missing " + "; ".join(f'"{text}"' for text in missing)
    return f", log check: {len(expected)} of {len(expected)} found"


def check_export(run, ifc, update_baseline):
    """compares the export's listing with its baseline, or saves it as the baseline. Returns a message for the summary"""
    current = compare_export.listing(ifc)
    baseline_file = EXPORT_BASELINE / f"{run}.txt.gz"
    diff_file = RESULTS / f"{run}.export-diff.txt"
    diff_file.unlink(missing_ok=True)
    if update_baseline or not baseline_file.exists():
        EXPORT_BASELINE.mkdir(exist_ok=True)
        compare_export.save(current, baseline_file)
        return ", export baseline saved"
    baseline = compare_export.listing(baseline_file)
    n = compare_export.count_differences(baseline, current)
    if n == 0:
        return ", export same as baseline"
    diff_file.write_text("\n".join(compare_export.compare(baseline, current, limit=5000)) + "\n", encoding="utf-8")
    return f", export differs from baseline: {n} lines (see {diff_file.name})"


def check_ids(run, ids, update_baseline):
    """compares the design-value IDS with its baseline, or saves it as the baseline. Returns a message for the summary"""
    if not ids.exists():
        return ", IDS not written"
    current = compare_ids.listing(ids)
    baseline_file = EXPORT_BASELINE / f"{run}.ids.txt"
    diff_file = RESULTS / f"{run}.ids-diff.txt"
    diff_file.unlink(missing_ok=True)
    if update_baseline or not baseline_file.exists():
        EXPORT_BASELINE.mkdir(exist_ok=True)
        compare_ids.save(current, baseline_file)
        return ", IDS baseline saved"
    baseline = compare_ids.listing(baseline_file)
    n = compare_ids.count_differences(baseline, current)
    if n == 0:
        return ", IDS same as baseline"
    diff_file.write_text("\n".join(compare_ids.compare(baseline, current, limit=5000)) + "\n", encoding="utf-8")
    return f", IDS differs from baseline: {n} lines (see {diff_file.name})"


def remove_from_ifc(ifc, classes):
    """removes the entities of the IFC classes from the file, with what they nest (e.g. an IfcAlignmentVertical and
    its segments), to test the import of a model without them"""
    import ifcopenshell
    f = ifcopenshell.open(str(ifc))
    doomed = []
    pending = [e for c in classes for e in f.by_type(c)]
    while pending:
        e = pending.pop()
        doomed.append(e)
        for rel in getattr(e, "IsNestedBy", None) or []:
            pending += list(rel.RelatedObjects)
            doomed.append(rel)
    ids = {e.id() for e in doomed}
    for rel in f.by_type("IfcRelNests"):
        if rel.id() not in ids and any(o.id() in ids for o in rel.RelatedObjects):
            rel.RelatedObjects = [o for o in rel.RelatedObjects if o.id() not in ids]
    for e in doomed:
        if e.id() in ids:
            ids.discard(e.id())
            f.remove(e)
    f.write(str(ifc))


def template_has_strands(template):
    """True if any girder of the template (.pgt) has strands"""
    root = ElementTree.parse(template).getroot()
    for prestress in root.iter("PrestressData"):
        counts = [prestress.findtext(tag) for tag in ("NumStraightStrands", "NumHarpedStrands", "NumTempStrands", "NumPermanentStrands")]
        if any(c and int(c) for c in counts) or any(int(rows.findtext("Count") or 0) for rows in prestress.iter("StrandRows")):
            return True
    return False


def compare_strands(exe, run, pgs, units, expected_strands, config, tolerances):
    """Compares the strands of the source export (expected_strands, an ifc_strands summary) with the strands of an export
    of the imported project pgs. Returns (comparison rows, message for the summary)."""
    reexport = RESULTS / "reexport" / f"{run}.ifc"
    reexport.parent.mkdir(exist_ok=True)
    reexport.unlink(missing_ok=True)
    ok, message = run_bridgelink(exe, [quote(f"/IfcExport={reexport}"), quote(pgs), f"/IfcPropertyUnits={units}"])
    if not ok or not reexport.exists():
        return [], f", strands not compared: the export of the imported project failed ({message})"
    # the template's girders have no strands, so a strand count of 0 is the template default
    template = None
    if not config["template_has_strands"]:
        template = {"values": {key: {"value": 0, "kind": "count"} for key, v in expected_strands["values"].items() if v["kind"] == "count"}}
    return compare_bridge.compare(ifc_strands.extract(reexport), expected_strands, template, tolerances), ""


def round_trip(model, config, exe):
    """export -> import -> compare for a PGSuper project. Returns summary rows."""
    rows = []
    source = config["base"] / model["pgs"]
    expected = compare_bridge.load(str(source))
    stations = expected["meta"]["stations"]
    template = compare_bridge.load(str(config["template"]), stations)

    for units in model.get("property_units", ["Display"]):
        run = f"{model['name']}-{units}"
        ifc = RESULTS / f"{run}.ifc"
        pgs = RESULTS / f"{run}.pgs"

        ids = RESULTS / f"{run}.ids"
        ids.unlink(missing_ok=True)
        ok, message = run_bridgelink(exe, [quote(f"/IfcExport={ifc}"), quote(source), f"/IfcPropertyUnits={units}", quote(f"/IfcIds={ids}")])
        if not ok or not ifc.exists():
            rows.append((run, f"export failed: {message}", None))
            continue
        expected_strands = ifc_strands.extract(ifc)  # before anything is removed from the export
        export_message = check_export(run, ifc, config["update_export_baseline"]) + check_ids(run, ids, config["update_export_baseline"])
        for standard in config["standards"]:
            export_message += check_standard.ids_message(standard, ifc, run)
        if model.get("remove_from_export"):
            remove_from_ifc(ifc, model["remove_from_export"])
            export_message += ", removed " + ", ".join(model["remove_from_export"]) + " before the import"

        ok, message = import_ifc(exe, ifc, config["template"], pgs)
        message += export_message
        if not ok:
            rows.append((run, message, None))
            continue

        imported = compare_bridge.load(str(pgs), stations)
        comparison = compare_bridge.compare(imported, expected, template, model.get("tolerances"))
        strand_rows, strand_message = compare_strands(exe, run, pgs, units, expected_strands, config, model.get("tolerances"))
        comparison += strand_rows
        message += strand_message
        report = compare_bridge.report(comparison, f"{run}: round trip", rel(pgs), rel(source), stations=stations)
        (RESULTS / f"{run}.md").write_text(report, encoding="utf-8")
        rows.append((run, message, compare_bridge.score(comparison)))
    return rows


def ifc_model(model, config, exe):
    """import -> compare (or inventory) for an IFC model. Returns summary rows."""
    run = model["name"]
    pgs = RESULTS / f"{run}.pgs"
    mapping = (config["base"] / model["mapping"]).resolve() if "mapping" in model else None
    ok, message = import_ifc(exe, unpack(config["base"] / model["ifc"]), config["template"], pgs, mapping)
    if not ok:
        return [(run, message, None)]
    message += check_log(model, pgs)

    expected_file = HERE / "expected" / f"{run}.json"
    if expected_file.exists():
        expected = compare_bridge.load(str(expected_file))
        imported = compare_bridge.load(str(pgs), expected["meta"].get("stations"))
        stations = imported["meta"]["stations"]
        template = compare_bridge.load(str(config["template"]), stations)
        comparison = compare_bridge.compare(imported, expected, template)
        report = compare_bridge.report(comparison, f"{run}: IFC import", rel(pgs), rel(expected_file), stations=stations)
        (RESULTS / f"{run}.md").write_text(report, encoding="utf-8")
        return [(run, message, compare_bridge.score(comparison))]

    imported = compare_bridge.load(str(pgs))
    stations = imported["meta"]["stations"]
    template = compare_bridge.load(str(config["template"]), stations)
    rows = compare_bridge.inventory(imported, template)
    report = compare_bridge.inventory_report(rows, f"{run}: IFC import", rel(pgs), stations)
    (RESULTS / f"{run}.md").write_text(report, encoding="utf-8")
    set_by_import = sum(1 for r in rows if r[1] == "imported")
    return [(run, message + f", no expected values: {set_by_import} of {len(rows)} values set by the import", None)]


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--models", nargs="*", help="names of the models to run (default: all)")
    parser.add_argument("--configuration", help='PGSuper configuration for the run, "server" or "server:publisher"')
    parser.add_argument("--bridgelink", help=r"BridgeLink.exe (default: %%ARPDIR%%\BridgeLink\RegFreeCOM\x64\Release\BridgeLink.exe)")
    parser.add_argument("--config-file", default=str(HERE / "models.json"))
    parser.add_argument("--update-export-baseline", action="store_true", help="save the round-trip exports as the export baseline")
    parser.add_argument("--skip-standards", action="store_true", help="leave out the checks against the standards in standards.json")
    args = parser.parse_args()

    config_file = Path(args.config_file).resolve()
    config = json.loads(config_file.read_text(encoding="utf-8"))
    config["base"] = config_file.parent
    config["template"] = (config["base"] / config["template"]).resolve()
    config["template_has_strands"] = template_has_strands(config["template"])
    config["update_export_baseline"] = args.update_export_baseline
    config["standards"] = [] if args.skip_standards else check_standard.load_standards()
    models = [m for m in config["models"] if not args.models or m["name"] in args.models]

    exe = args.bridgelink
    if not exe:
        if "ARPDIR" not in os.environ:
            parser.error("ARPDIR is not set - use --bridgelink")
        exe = str(Path(os.environ["ARPDIR"]) / "BridgeLink" / "RegFreeCOM" / "x64" / "Release" / "BridgeLink.exe")
    if not Path(exe).exists():
        parser.error(f"{exe} not found")

    RESULTS.mkdir(exist_ok=True)

    mapping_setting = get_mapping_table_setting()
    if mapping_setting:
        print(f"WARNING: the configured IFC mapping table is {mapping_setting}. Runs without a \"mapping\" use it instead of the standard table.")

    original_configuration = None
    if args.configuration:
        server, _, publisher = args.configuration.partition(":")
        original_configuration = get_configuration()
        print(f'Configuration "{server}":"{publisher or server}" (was "{original_configuration[0]}":"{original_configuration[1]}")')
        set_configuration(exe, server, publisher or server)

    original_label_format = get_girder_label_format()
    if original_label_format != GIRDER_LABEL_FORMAT:
        print(f'Girder label format {GIRDER_LABEL_FORMAT} (was {original_label_format or "not set"})')
        set_girder_label_format(GIRDER_LABEL_FORMAT)

    summary = []
    try:
        for model in models:
            print(f"Running {model['name']}...", flush=True)
            summary += round_trip(model, config, exe) if "pgs" in model else ifc_model(model, config, exe)
        if not args.models:
            print("Running the standard table <-> IDS round trip...", flush=True)
            ok, message = check_table_ids.check(exe)
            summary.append(("Standard table <-> IDS", ("OK: " if ok else "FAILED: ") + message, None))
            for standard in config["standards"]:
                print(f"Comparing the standard table with the {standard['name']} workbook...", flush=True)
                ok, message = check_standard.check_workbook(standard)
                summary.append((f"Standard table vs {standard['name']} workbook", ("OK: " if ok else "DIFFERS: ") + message, None))
    finally:
        if original_label_format != GIRDER_LABEL_FORMAT:
            set_girder_label_format(original_label_format)
            print(f'Girder label format restored to {original_label_format or "not set"}')
        if original_configuration:
            set_configuration(exe, *original_configuration)
            print(f'Configuration restored to "{original_configuration[0]}":"{original_configuration[1]}"')

    configuration = args.configuration or "{}:{}".format(*get_configuration())
    lines = [f"# IFC import validation", "", f"- Run: {datetime.datetime.now():%Y-%m-%d %H:%M}", f"- Configuration: {configuration}"]
    if mapping_setting:
        lines.append(f"- Configured mapping table: {mapping_setting} (used by the runs without a \"mapping\")")
    lines += ["",
             "| Run | Result | Match | Default match | Mismatch | Not imported | Not comparable | Missing |",
             "|---|---|---|---|---|---|---|---|"]
    for run, message, score in summary:
        if score:
            counts, compared = score
            lines.append(f"| [{run}]({run}.md) | {message} | {counts['match']} of {compared} | {counts['default match']} | {counts['mismatch']} | "
                         f"{counts['not imported']} | {counts['not comparable']} | {counts['missing']} |")
        else:
            link = f"[{run}]({run}.md)" if (RESULTS / f"{run}.md").exists() else run
            lines.append(f"| {link} | {message} | | | | | | |")
    text = "\n".join(lines) + "\n"
    (RESULTS / "summary.md").write_text(text, encoding="utf-8")
    print(text)
    print(f"Reports in {RESULTS}")


if __name__ == "__main__":
    main()
