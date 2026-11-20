# Export profiles: one export mechanism for the standard export and export extensions

Evaluated 2026-11-19 (after the export extensions were merged) and implemented 2026-11-20. The question was whether the standard (usBridge) export should use the same hook mechanism as an export extension's export ([ExtensionsDesign.md](ExtensionsDesign.md)), so that both are specializations of one generic export.

## The answer

Yes in substance, but the generalization runs the other way: the standard export wasn't turned into a hook-based extension.

- **There was already one build.** `CIfcExporter::BuildModel` builds both kinds of export. What makes an export usBridge is data: the mapping table's property sets, constants, classifications, and element selectors.
- **What was usBridge-specific in C++ wasn't hook-shaped.** Hooks add elements and property sets at fixed points of the build. The usBridge code was either data that belonged in the table (it's there now) or decisions about how the model is built, which every export shares: the exporter's IFC conventions ([MappingTablesDesign.md](MappingTablesDesign.md), "The exporter's IFC conventions").
- **So every export is a profile:** a mapping table, hooks, and option overrides. usBridge is the built-in profile with hooks that do nothing; an export extension's export is a profile with hooks. Both run the same build, session, and hook code, so the hook path is used by every export and tested by every validation run, not only by the extensions' own tests.

## The profile

`CExportProfile` (`ExportProfile.h`):

- `CExportProfile::Standard()`: no extension, the table of the export options.
- `CExportProfile::Extension(extension)`: the extension's targets, element roles, and hooks, and its table.
- `CreateHooks()`: the extension's hooks, or a do-nothing `IfcExport::IHooks`. Never null.
- `Apply(options)`: the profile's option overrides. An extension's table is reported as the extension's (`MappingTableSource::ExportExtension`), and an extension's export can't be a round trip.

`CIfcExporter::BuildModel` takes a profile. The three entry points make one: File > Export (`PGSuperDataExporter`, the standard profile), `/IfcExport` (the standard profile, or the test extension's with `/IfcTestExtension`), and `IIfcModelBuilder::BuildModel` (the extension's, or the standard profile without one). The session always has hooks, so the build has no special cases for "no extension".

## What moved out of the C++ (2026-11-20)

- `usBrPset_ACI_BarShape` is declared in `Standard.ifcmap`, with targets for the bend shape's standard and version ("ACI 315-99", the standard of PGSuper's shape codes), the shape code, the bend radius, and one target per ACI dimension (`rebar.dimension_a` ... `rebar.dimension_v`, absent when the bar's shape doesn't have it). usBridge lists each dimension as its own property, so no new table mechanism was needed.
- The bSDD dictionaries are table data (`bsdd_dictionaries`, with a default `bsdd_dictionary`); the standard table has usBridge's.
- The design-value IDS requires every constant the table declares for a role, instead of naming property sets.
- The empty usBridge bar type classification placeholder is gone; the classification helpers are in `Classifications.h`.

Each change was checked with the export baselines: the exports are unchanged.

## The smaller extension API (2026-11-20)

Extensions pass `IfcExport::BuildOptions` and `IfcExport::ElementKeys`, and their targets read an element through `IfcExport::Element` (accessors). `ExportContext`, `TargetDef`, and `CIfcExportOptions` aren't in the SDK. Version 1 stays provisional until it's released; extensions built against the earlier SDK must be changed to match. See [ExtensionsDesign.md](ExtensionsDesign.md), "The SDK".

## The test profile

`CTestExportExtension` (`TestExportExtension.h`) and `check_extension.py`: see [ExtensionsDesign.md](ExtensionsDesign.md), "Testing the API".

## What it doesn't do

- **It doesn't remove an extension's coupling to this DLL.** An extension still shares C++ objects with it (IfcOpenShell types in the hooks), so the two must be built together. The remedies are a per-release API version, a pinned standard table, and a contract test (the test extension is this repository's side of it).
- **It doesn't make usBridge an external plug-in.** Nothing usBridge-specific needs to live outside this DLL; what's left in code is the exporter's conventions.
- **The model's shape isn't a profile option.** Make a convention a profile option only when a standard needs something different.
