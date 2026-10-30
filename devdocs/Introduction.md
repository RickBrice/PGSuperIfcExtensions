# PGSuper IFC Extensions: Developer Documentation

The design notes, plans, and explanations for the PGSuper IFC extensions: importing IFC bridge models into PGSuper, exporting PGSuper projects to IFC, mapping tables, and round trips.

Start with the import plan and the mapping tables design.

| Document | What it covers |
|---|---|
| [IFC import plan](IfcImportPlan.md) | The phased import plan, its status, decisions, and what's next |
| [Import geometry](ImportGeometry.md) | How PGSuper input is derived from the model's geometry, with figures |
| [Girder import plan](GirderImportPlan.md) | Girder sections, library entries, strands, and girder reinforcement |
| [Mapping tables design](MappingTablesDesign.md) | The element registry, targets, the table format, and how the import and export use them |
| [Mapping table targets](Targets.md) | Each target and the PGSuper data and code behind it |
| [Round trip design](RoundTripDesign.md) | Exporting a project as an updated copy of the model it was imported from |
| [PGSplice plan](PGSplicePlan.md) | Import and export for PGSplice |
| [3D view plan](3DViewPlan.md) | Exploratory: a 3D view of the bridge model, native or from IFC |
| [Dependencies](Dependencies.md) | Building IfcOpenShell and the other dependencies |

## Building this documentation

Run Doxygen in this folder:

    doxygen Doxyfile.dox

The HTML is written to `API/html` (open `API/html/index.html`); it isn't tracked by git.
