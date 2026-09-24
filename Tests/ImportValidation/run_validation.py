"""Runs the IFC import validation for the models in models.json and writes reports.

For a model with a "pgs" file (round trip):
   export the project to IFC (once for each "property_units" mode), import the IFC into a new
   project, and compare the imported project to the original project.
For a model with an "ifc" file (.ifc, or .zip containing one .ifc):
   import the IFC into a new project and compare to expected/<name>.json if it exists,
   otherwise report which values the import set (inventory).

Reports, logs, and imported projects are written to results/, replacing the previous results, so
the change from one run to the next shows up as a difference in the repository. summary.md has the
score of every run. Exported IFC files and unzipped models are written to results/ but not tracked.

Paths in models.json are relative to models.json. BridgeLink is the Release build in %ARPDIR%
unless --bridgelink is given.

--configuration switches the PGSuper configuration (e.g. "WSDOT (Local)" for exact girder
name matches) for the duration of the run and switches back to the current configuration
when done.

Usage: python run_validation.py [--models PGSuper Iowa] [--configuration "WSDOT (Local)"] [--bridgelink BridgeLink.exe]
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

import compare_bridge

HERE = Path(__file__).resolve().parent
RESULTS = HERE / "results"
TIMEOUT = 600  # seconds for one BridgeLink run
PGSUPER_OPTIONS_KEY = r"Software\Washington State Department of Transportation\PGSuper\Options"


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


def import_ifc(exe, ifc, template, pgs):
    """imports ifc into a new project pgs. Returns (ok, message)"""
    ok, message = run_bridgelink(exe, [quote(f"/IfcImport={ifc}"), quote(template), quote(f"/IfcOut={pgs}")])
    if not ok or not pgs.exists():
        return False, f"import failed: {message}"
    if "IFC import succeeded" not in read_log(pgs.with_suffix(".log")):
        message += ", import reported errors - see log"
    return True, message


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

        ok, message = run_bridgelink(exe, [quote(f"/IfcExport={ifc}"), quote(source), f"/IfcPropertyUnits={units}"])
        if not ok or not ifc.exists():
            rows.append((run, f"export failed: {message}", None))
            continue

        ok, message = import_ifc(exe, ifc, config["template"], pgs)
        if not ok:
            rows.append((run, message, None))
            continue

        imported = compare_bridge.load(str(pgs), stations)
        comparison = compare_bridge.compare(imported, expected, template)
        report = compare_bridge.report(comparison, f"{run}: round trip", rel(pgs), rel(source), stations=stations)
        (RESULTS / f"{run}.md").write_text(report, encoding="utf-8")
        rows.append((run, message, compare_bridge.score(comparison)))
    return rows


def ifc_model(model, config, exe):
    """import -> compare (or inventory) for an IFC model. Returns summary rows."""
    run = model["name"]
    pgs = RESULTS / f"{run}.pgs"
    ok, message = import_ifc(exe, unpack(config["base"] / model["ifc"]), config["template"], pgs)
    if not ok:
        return [(run, message, None)]

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
    args = parser.parse_args()

    config_file = Path(args.config_file).resolve()
    config = json.loads(config_file.read_text(encoding="utf-8"))
    config["base"] = config_file.parent
    config["template"] = (config["base"] / config["template"]).resolve()
    models = [m for m in config["models"] if not args.models or m["name"] in args.models]

    exe = args.bridgelink
    if not exe:
        if "ARPDIR" not in os.environ:
            parser.error("ARPDIR is not set - use --bridgelink")
        exe = str(Path(os.environ["ARPDIR"]) / "BridgeLink" / "RegFreeCOM" / "x64" / "Release" / "BridgeLink.exe")
    if not Path(exe).exists():
        parser.error(f"{exe} not found")

    RESULTS.mkdir(exist_ok=True)

    original_configuration = None
    if args.configuration:
        server, _, publisher = args.configuration.partition(":")
        original_configuration = get_configuration()
        print(f'Configuration "{server}":"{publisher or server}" (was "{original_configuration[0]}":"{original_configuration[1]}")')
        set_configuration(exe, server, publisher or server)

    summary = []
    try:
        for model in models:
            print(f"Running {model['name']}...", flush=True)
            summary += round_trip(model, config, exe) if "pgs" in model else ifc_model(model, config, exe)
    finally:
        if original_configuration:
            set_configuration(exe, *original_configuration)
            print(f'Configuration restored to "{original_configuration[0]}":"{original_configuration[1]}"')

    configuration = args.configuration or "{}:{}".format(*get_configuration())
    lines = [f"# IFC import validation", "", f"- Run: {datetime.datetime.now():%Y-%m-%d %H:%M}", f"- Configuration: {configuration}", "",
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
