# PGSuperIfcExtensions: notes for Claude

Project knowledge that isn't in the code or git history. The design notes are in `devdocs/` (start with `IfcImportPlan.md` and `MappingTablesDesign.md`).

## Working rules
- **Commits:** commit code and docs first, then the validation results (`Tests/ImportValidation/results`, `export_baseline`) in a separate "Updates ... validation results" commit. Never push; the user pushes after the overnight PGSuper regression tests. Changes in WBFL and PGSuper stay uncommitted until the user says so. Don't commit the stray `PGSuper_Import_Model.dbp/.dbr/.ifc` files.
- **PGSuper configuration:** it's fine to switch the configuration for tests with `BridgeLink.exe /App=PGSuper /Configuration="<server>":"<publisher>"`. Record the current names first (`HKCU\Software\Washington State Department of Transportation\PGSuper\Options`, `CatalogServer2`/`Publisher2`) and switch back afterwards. The dev machine normally runs `"Regression":"Regression"`, which the regression tests depend on. `"WSDOT (Local)":"WSDOT (Local)"` gives exact matches for the WSDOT WF girders of the PGSuper export.
- **Units:** PGSuper's system units are SI (kg, m, Pa). In-memory data and `.pgs` files are SI. Import code converts IFC values with `WBFL::Units::ConvertToSysUnits`.
- **IfcOpenShell:** the upstream branch that matters is `v0.9.0` (not v0.8.0 or master). Compare against it and target PRs to it. The user's working branch is `rab_infrastructure`.

## Direction (2026-10-15)
- **The standard is data, not code.** The current target is usBridge: the IDS files and `usBridge_Data_Dictionary_Workbook_Version_2.0.xlsx` (IFC classes, predefined types, usBridge classification names, property sets with cardinality and enumerations, and the spatial hierarchy) in the user's OneDrive (BIM for Bridges\usBridge). They will change, and other standards must be usable too.
- Strive first to minimize, then to eliminate, C++ changes and DLL rebuilds when a standard changes. Expect this to need PGSuper API improvements and documentation (data PGSuper exposes, with names, units, and meaning).
- Prefer solutions that move standard-specific content (property sets, names, classes, structure) into mapping tables or other data files over adding it to the C++.

## Design decisions to keep
- **Importer template:** the importer starts from `IfcImportTemplate.pgt` on purpose (interim, fault tolerant across PGSuper versions). The goal is building `CBridgeDescription2` from scratch, so keep new import code independent of template values. Don't propose removing the template as a side fix.
- **Design-value IDS vs general IDS:** the IDS from `IdsExporter` is model specific: each specification is pinned to one element and requires PGSuper's exact design values. It is not a general (agency) IDS and is never a source for mapping tables. Mapping tables are generated from, and written as, general IDS files (`/IfcIdsToTable`, `/IfcTableToIds`).
- **Mapping tables:** PGSuper data items are targets defined in code (`IfcTargets.cpp`); tables only say where they live in IFC. `MappingTables/Standard.ifcmap` was written by hand from what the exporter wrote before the tables existed (TPF/usBridge property sets); it wasn't generated from the usBridge IDS files. The Iowa and PennDOT tables in `Tests/ImportValidation/mappings` extend it with names found in their example models.
- **Girder type StrengthClass:** the `IfcBeamType`'s `Pset_ConcreteElementGeneral.StrengthClass` is "Override on entity instance", and each `IfcBeam` overrides it with its value. usBridge IDS specs that apply to anything with the property set also check type objects, and a required property needs a value.
- **Tendon type/material names:** `IfcTendonType.Name` equals the strand `IfcMaterial.Name` on purpose. Bonsai's IDS auditor wrongly matches the type with a material specification (applicability facets are ANDed in IDS 1.0; the reference ifctester gets it right). The collision stays so the auditor bug stays visible; don't rename to work around it. Check other same-name symptoms against the reference ifctester first, and ask.
- **GRF005 warning:** `IfcMapConversion.Scale` includes the grid scale factor (`create_map_conversion` in `MapConversion.h`), as bSI's geo-referencing guide recommends. The bSI validator's GRF005 checks only the unit ratio (open issue buildingSMART/ifc-gherkin-rules#522). Don't drop the grid scale factor.

## Pending
- **CoatingSpecification (2026-10-15):** check whether the usBridge bSDD defines a list of values for `CoatingSpecification` (https://identifier.buildingsmart.org/uri/usTransportation/usBridge/1/prop/CoatingSpecification); bSDD was down. The reinforcing material's is "NONE" (the usBridge IDS says to use NONE when there's no coating). The tendon material's has no value: its `CoatingSpecificationVersion` is required and there's no version for NONE. Remind the user at the start of mapping table or usBridge property work until it's resolved.

## Checking an export against an IDS
`Tests/ImportValidation/check_ids.py <file.ifc> <file.ids> [--spec NAME]` runs ifctester and lists the failures. The usBridge property set IDS is `usBrSpec_PropertySets.ids` in the user's OneDrive (BIM for Bridges\usBridge). The Python that comes with the local IfcOpenShell build has an ifctester (0.8.5) that is older than its ifcopenshell; use the ifctester from the IfcOpenShell source tree: `PYTHONPATH=<IfcOpenShell>/src/ifctester`.
