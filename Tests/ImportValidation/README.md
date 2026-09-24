# IFC import validation

Checks how well the IFC importer maps an IFC model into a PGSuper project. It runs BridgeLink from the command line (`/IfcExport`, `/IfcImport`), so no UI interaction is needed. Everything it needs is in this folder: the test models, the import template (`IfcImportTemplate.pgt` at the repository root), and the latest results.

Requirements: Python 3 (standard library only) and a Release build of BridgeLink with this extension. BridgeLink is found through `%ARPDIR%` (`%ARPDIR%\BridgeLink\RegFreeCOM\x64\Release\BridgeLink.exe`) unless `--bridgelink` is given.

```
python run_validation.py                                   # all models
python run_validation.py --models PGSuper                  # selected models
python run_validation.py --configuration "WSDOT (Local)"   # run with another PGSuper configuration
```

## Results

`results/` has the results of the last run: `summary.md` with a score for every run, a report (`.md`), the import and export logs, and the imported project (`.pgs`) for each model. Each run replaces the previous results, so committing them after an importer change records how the change affected the imports. Exported IFC files and unzipped models are also written to `results/` but are not tracked.

## Models

`models.json` lists the models, with paths relative to it. The model files are in `models/`. Large IFC files are zipped, and the runner extracts them to `results/models/`.

- **`pgs` models (round trip):** the project is exported to IFC once for each `property_units` mode (`Display`, `System`), the IFC is imported into a new project, and the result is compared to the original project.
- **`ifc` models:** the IFC is imported. If `expected/<name>.json` exists, the result is compared to it. Otherwise the report is an inventory of which values the import set and which are still template defaults.

An expected file lists values in SI (PGSuper system units) by the keys that `pgs_extract.py` produces, e.g. `{"values": {"pier[1].station": 16.4592, "group[0].girder_count": 5}}`. Only the listed values are compared.

| Model | Source |
|---|---|
| `PGSuper_Import_Model.pgs` | PGSuper project: 3 spans, 4/5/4 WF girders, curved alignment |
| `PennDOT_Rearick_Rd_323010451_STR1.ifc` | PennDOT Rearick Road Bridge #1, from OpenBridge (Quadri): single span, no alignment |
| `Iowa_ADCMS_Pilot_US59_over_IA92_2025_07_09.ifc.zip` | Iowa DOT ADCMS pilot, US 59 over IA 92: 3 spans, no alignment |

## Scripts

- `pgs_extract.py`: reduces a `.pgs`/`.pgt` to a flat summary of the values the importer is responsible for. Alignment, profile and deck edges are evaluated at the pier and mid-span stations, so an alignment defined from a different reference point still compares equal. Alignments with spirals are not evaluated yet.
- `compare_bridge.py`: compares a project to expected values with a tolerance for each kind of value, and writes a Markdown report in US units.
  - A *default match* equals the expected value only because the import template has the same value, so it doesn't show the importer set it.
  - Values measured differently (e.g. girder spacing along vs normal to the pier) are *not comparable*.
- `run_validation.py`: runs BridgeLink and the comparison for each model. `--configuration` switches the PGSuper configuration for the run, then switches back to the current one.

Exact girder name matches need a configuration whose library has the girders, e.g. `WSDOT (Local)` for the WF girders in the PGSuper model. The committed results use the `Regression` configuration.
