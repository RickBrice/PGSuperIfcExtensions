"""Checks the export extension API with the test export extension (TestExportExtension.h, devdocs/ExportProfilesPlan.md).

The main PGSuper model is exported with /IfcTestExtension and mappings/TestExtension.ifcmap, a table that extends the
standard table and "requires" the extension:

- hooks: the export succeeds; the hooks are called in the documented order (devdocs/ExtensionsDesign.md), once for each
  pier, bearing, and girder segment of the file; the extension's marker element has its property set, and every girder
  and pier has the property set of the extension's targets on them; the rest of the export is the standard export (the
  export baseline of PGSuper-System), and the file is valid
- the refusals: a hook that throws (throw), an extension built for another API version (api) or with other build
  settings (fingerprint), and the table in a standard export (without the extension) each stop the export with a
  message, and no file is written
- the general IDS of the table with the extension's targets (/IfcTableToIds with /IfcTestExtension)

Results in results/ExtensionApi.md. run_validation.py runs it (the "Export extension API" row).

Usage: python check_extension.py [--bridgelink BridgeLink.exe]
"""
import argparse
import itertools
import os
import subprocess
import sys
from pathlib import Path

import ifcopenshell
import ifcopenshell.util.element
import ifcopenshell.validate

import compare_export

HERE = Path(__file__).resolve().parent
RESULTS = HERE / "results"
OUT = RESULTS / "extension"
MODEL = HERE / "models" / "PGSuper_Import_Model.pgs"
TABLE = HERE / "mappings" / "TestExtension.ifcmap"
TEMPLATE = HERE.parent.parent / "IfcImportTemplate.pgt"
BASELINE = HERE / "export_baseline" / "PGSuper-System.txt.gz"
EXTENSION_ID = "PGSuperIfc.TestExtension"
USBRIDGE = "https://identifier.buildingsmart.org/uri/usTransportation/usBridge/1/"  # the standard table's default bSDD dictionary
TEST_DICTIONARY = "https://identifier.buildingsmart.org/uri/example/test/1/"  # the test table's second one

# the hooks in the order IfcExport::IHooks documents (OnExportEnd isn't recorded: it's called after the model is written)
HOOK_ORDER = ["OnExportBegin", "OnDeck", "OnBridgeParts", "OnPier", "OnBearing", "OnBarrier", "OnSegment",
              "OnSuperstructureComplete", "OnBridgeComplete", "OnModelComplete"]
REPEATED = {"OnPier", "OnBearing", "OnBarrier", "OnSegment"}


def quote(text):
    return f'"{text}"'


def export(run_bridgelink, name, arguments):
    """exports the main model with the arguments. Returns (ifc, log text)"""
    OUT.mkdir(parents=True, exist_ok=True)
    ifc = OUT / f"{name}.ifc"
    log = ifc.with_suffix(".export.log")
    ifc.unlink(missing_ok=True)
    log.unlink(missing_ok=True)
    run_bridgelink([quote(f"/IfcExport={ifc}"), quote(MODEL), "/IfcPropertyUnits=System"] + arguments)
    return ifc, log.read_text(errors="replace") if log.exists() else ""


def check_hooks(run_bridgelink):
    """the export with the extension's hooks. Returns a list of problems"""
    ifc, log = export(run_bridgelink, "Hooks", [quote(f"/IfcMapping={TABLE}"), "/IfcTestExtension"])
    if "IFC export succeeded" not in log or not ifc.exists():
        return ["the export failed: " + (log.strip().splitlines()[-1] if log.strip() else "no log")]

    problems = [line for line in log.splitlines() if line.startswith("Test extension problem")]
    if f"Export extension {EXTENSION_ID}" not in log:
        problems.append("the log doesn't name the extension")

    model = ifcopenshell.open(str(ifc))
    psets = ifcopenshell.util.element.get_psets
    markers = [e for e in model.by_type("IfcBuildingElementProxy") if e.ObjectType == "TEST_MARKER"]
    if len(markers) != 1:
        return problems + [f"{len(markers)} marker elements (expected 1)"]
    marker = psets(markers[0]).get("TestPset_Marker", {})
    calls = (marker.get("HookCalls") or "").split(";")
    if marker.get("HookCallCount") != len(calls):
        problems.append(f"the marker's HookCallCount is {marker.get('HookCallCount')}, but it has {len(calls)} calls")

    # the order: each hook once (or once for each element), in the documented order
    groups = [(hook, len(list(g))) for hook, g in itertools.groupby(calls)]
    if [hook for hook, n in groups] != HOOK_ORDER:
        problems.append("the hooks were called in this order: " + ", ".join(f"{hook} x{n}" for hook, n in groups))
    counts = dict(groups)
    for hook, n in groups:
        if hook not in REPEATED and n != 1:
            problems.append(f"{hook} was called {n} times")

    girders = [b for b in model.by_type("IfcBeam") if "TestPset_Girder" in psets(b)]
    keys = [psets(b)["TestPset_Girder"]["SegmentKey"] for b in girders]
    piers = [p for p in model.by_type("IfcBridgePart") if "TestPset_Pier" in psets(p)]
    bearings = model.by_type("IfcBearing")
    for hook, elements, what in (("OnSegment", girders, "girder segments with TestPset_Girder"), ("OnPier", piers, "piers and abutments with TestPset_Pier"),
                                 ("OnBearing", bearings, "bearings")):
        if counts.get(hook) != len(elements):
            problems.append(f"{hook} was called {counts.get(hook)} times, for {len(elements)} {what}")
    if len(set(keys)) != len(keys) or not keys:
        problems.append(f"the girders' SegmentKey values aren't one per segment: {keys}")
    if sorted(psets(p)["TestPset_Pier"]["PierIndex"] for p in piers) != list(range(len(piers))):
        problems.append("the piers' PierIndex values aren't 0, 1, ...")

    # bSDD URIs in two dictionaries: the test table's own ("bsdd@test", "bsdd@test:Calls") and the standard table's default ("bsdd")
    def specification(element, pset, name):
        for rel in element.IsDefinedBy or []:
            definition = getattr(rel, "RelatingPropertyDefinition", None)
            if definition and definition.is_a("IfcPropertySet") and definition.Name == pset:
                for p in definition.HasProperties:
                    if p.Name == name:
                        return p.Specification
        return None
    for element, pset, name, expected in ((markers[0], "TestPset_Marker", "HookCalls", TEST_DICTIONARY + "prop/Calls"),
                                          (girders[0] if girders else None, "TestPset_Girder", "SegmentKey", TEST_DICTIONARY + "prop/SegmentKey"),
                                          (piers[0] if piers else None, "TestPset_Pier", "PierIndex", USBRIDGE + "prop/PierIndex")):
        if element is not None and specification(element, pset, name) != expected:
            problems.append(f"{pset}.{name} has the URI {specification(element, pset, name)!r}, not {expected!r}")

    # the rest is the standard export
    if BASELINE.exists():
        current = [line for line in compare_export.listing(ifc) if "TestPset_" not in line and '"Test extension marker"' not in line]
        n = compare_export.count_differences(compare_export.listing(BASELINE), current)
        if n:
            diff = RESULTS / "ExtensionApi.export-diff.txt"
            diff.write_text("\n".join(compare_export.compare(compare_export.listing(BASELINE), current, limit=2000)) + "\n", encoding="utf-8")
            problems.append(f"without the extension's property sets and marker, the export differs from the PGSuper-System baseline: {n} lines (see {diff.name})")

    logger = ifcopenshell.validate.json_logger()
    ifcopenshell.validate.validate(model, logger, express_rules=True)
    if logger.statements:
        problems.append(f"{len(logger.statements)} schema errors, e.g. {logger.statements[0].get('message', '')}")
    return problems


def check_refused(run_bridgelink, name, arguments, expected):
    """an export that must fail with the expected text in its log, without writing a file. Returns a list of problems"""
    ifc, log = export(run_bridgelink, name, arguments)
    problems = []
    if "IFC export failed" not in log:
        problems.append("the export didn't fail")
    if expected not in log:
        problems.append(f"the log doesn't say \"{expected}\"")
    if ifc.exists():
        problems.append("a file was written")
    return problems


def check_table_ids(run_bridgelink):
    """the general IDS of the table, with the extension's targets and roles. Returns a list of problems"""
    OUT.mkdir(parents=True, exist_ok=True)
    ids = OUT / "TestExtension.ids"
    ids.unlink(missing_ok=True)
    run_bridgelink([quote(f"/IfcTableToIds={ids}"), quote(TEMPLATE), quote(f"/IfcMapping={TABLE}"), "/IfcTestExtension"])
    if not ids.exists():
        log = Path(str(ids) + ".log")
        return ["the IDS wasn't written: " + (log.read_text(errors="replace").strip().splitlines()[-1] if log.exists() else "no log")]
    text = ids.read_text(encoding="utf-8", errors="replace")
    return [f"the IDS has no {name}" for name in ("TestPset_Marker", "TestPset_Girder", "TestPset_Pier", "usBrPset_ACI_BarShape") if name not in text]


def check_unknown_dictionary(run_bridgelink):
    """a table whose URI names a dictionary no table declares can't be loaded. Returns a list of problems"""
    OUT.mkdir(parents=True, exist_ok=True)
    table = OUT / "UnknownDictionary.ifcmap"
    table.write_text('{"format": "PGSuperIfcMapping", "version": 1, "name": "Unknown dictionary", "extends": "standard", '
                     '"property_sets": [{"name": "TestPset_Bridge", "applies_to": "bridge", "uri": "bsdd@missing", '
                     '"properties": [{"name": "Name", "type": "IfcLabel"}]}]}', encoding="utf-8")
    ids = OUT / "UnknownDictionary.ids"
    ids.unlink(missing_ok=True)
    run_bridgelink([quote(f"/IfcTableToIds={ids}"), quote(TEMPLATE), quote(f"/IfcMapping={table}")])
    log = Path(str(ids) + ".log")
    text = log.read_text(errors="replace") if log.exists() else ""
    expected = "The bSDD dictionary \"missing\" isn't in the \"bsdd_dictionaries\""
    problems = []
    if ids.exists():
        problems.append("the IDS was written")
    if expected not in text:
        problems.append(f"the log doesn't say \"{expected}\"")
    return problems


def check(run_bridgelink):
    """runs the checks. Returns (ok, message)"""
    cases = [
        ("hooks", lambda: check_hooks(run_bridgelink)),
        ("a hook that throws", lambda: check_refused(run_bridgelink, "Throw", [quote(f"/IfcMapping={TABLE}"), "/IfcTestExtension=throw"],
                                                    f'The export extension "{EXTENSION_ID}" failed in OnPier: the test extension fails here on purpose')),
        ("another API version", lambda: check_refused(run_bridgelink, "ApiVersion", [quote(f"/IfcMapping={TABLE}"), "/IfcTestExtension=api"],
                                                     f'The export extension "{EXTENSION_ID}" was built for version')),
        ("another build", lambda: check_refused(run_bridgelink, "Fingerprint", [quote(f"/IfcMapping={TABLE}"), "/IfcTestExtension=fingerprint"],
                                               f'The export extension "{EXTENSION_ID}" wasn\'t built with the same tools')),
        ("the table without the extension", lambda: check_refused(run_bridgelink, "Requires", [quote(f"/IfcMapping={TABLE}")],
                                                                 f'the table needs the export extension "{EXTENSION_ID}"')),
        ("the general IDS of the table", lambda: check_table_ids(run_bridgelink)),
        ("a bSDD dictionary that isn't declared", lambda: check_unknown_dictionary(run_bridgelink)),
    ]
    lines = ["# Export extension API", "", "The test export extension (TestExportExtension.h) with mappings/TestExtension.ifcmap (check_extension.py).", "",
             "| Check | Result |", "|---|---|"]
    failed = []
    for name, case in cases:
        problems = case()
        lines.append(f"| {name} | {'OK' if not problems else 'FAILED: ' + '; '.join(problems)} |")
        if problems:
            failed.append(name)
    (RESULTS / "ExtensionApi.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    if failed:
        return False, f"{len(failed)} of {len(cases)} checks failed: {', '.join(failed)}"
    return True, f"{len(cases)} of {len(cases)} checks passed (hooks, the four refusals, the general IDS, and bSDD dictionaries)"


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--bridgelink", help=r"BridgeLink.exe (default: %%ARPDIR%%\BridgeLink\RegFreeCOM\x64\Release\BridgeLink.exe)")
    args = parser.parse_args()
    exe = args.bridgelink or str(Path(os.environ["ARPDIR"]) / "BridgeLink" / "RegFreeCOM" / "x64" / "Release" / "BridgeLink.exe")

    def run_bridgelink(arguments):
        subprocess.run(f'"{exe}" ' + " ".join(arguments), timeout=600)

    ok, message = check(run_bridgelink)
    print(("OK: " if ok else "FAILED: ") + message)
    print(f"Details in {RESULTS / 'ExtensionApi.md'}")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
