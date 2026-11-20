# Export extensions

Decision G7 in [MappingTablesDesign.md](MappingTablesDesign.md). Proposed 2026-10-05, merged 2026-11-19, and the SDK made smaller 2026-11-20 ([ExportProfilesPlan.md](ExportProfilesPlan.md)). Its users are plug-ins with their own export commands that write an IFC file for another application (an exchange file).

## Why

The IFC export writes standard and usBridge property sets only. Other applications need more: an exchange with an analysis program, for example, may need its own property sets and elements the standard export doesn't make (stirrup zones, deck reinforcement, user loads). Putting them in this project would mix application-specific data into a model that is meant to be "pure".

An export extension lets another plug-in build the model exactly as the IFC export does and add to it, without this project knowing what is added. The standard export command never loads an extension.

Both are export profiles (`CExportProfile`, `ExportProfile.h`): a mapping table, hooks, and option overrides. The standard export is the profile without an extension, whose hooks do nothing; an extension's export is the profile with its targets, roles, and hooks. Every export runs the same build, session, and hook code (`CIfcExporter::BuildModel`, `CIfcExportSession`), so the hook path is used by every export, and the test extension (below) tests it in every validation run.

## How it works (template method)

`CIfcExporter::BuildModel` builds the model in a fixed order. At fixed points it calls the profile's hooks. The extension supplies three things:

| What | Where | How the IFC extension uses it |
|---|---|---|
| A mapping table | its own file, `"extends": "standard"`, `"requires": ["<id>"]` | loaded like any table, so its property sets are written by the same code as the standard ones |
| Targets and element roles | `IIfcExportExtension::RegisterContent` | the table binds properties to the targets. Targets can belong to this project's roles (e.g. `girder`) or to the extension's own roles |
| Hooks | `IIfcExportExtension::CreateHooks` | called during the build, to add elements and write their property sets through `IServices` |

Another plug-in reaches the builder through the broker: `IIfcModelBuilder`, registered by `CIfcExtensionAgent`. `BuildModel(options, file, extension)` builds and writes the model; `WriteTableIds(table, ids, extension)` writes the general IDS of the extension's table (as `/IfcTableToIds`), with the extension's targets and roles. Neither throws: problems are in the `IfcBuildResult`, with the table files and what the hooks logged.

## The SDK

The header-only SDK is in `Include\IfcExport\`:

- `IfcExtensionTypes.h`: the value types (`ElementKind`, `ValueKind`, `ExportUnit`, `TargetValue`, `ExportValue`, `IntermediateDiaphragmDetails`), and what an extension's targets see of an element: `IfcExport::Element`, with an accessor for each key (`Broker()`, `Segment()`, `Pier()`, `Span()`, `Index()`, `Strand()`, `StrandType()`, `StrandIndex()`, `Rebar()`, `BarrierSide()`, `Bearing()`, `BearingElevation()`, `IntermediateDiaphragm()`, `Extension()`), and `IfcExport::ElementKeys`, the keys an extension gives its own elements (`segment`, `pier`, `span`, `index`, `extension`).
- `IfcModelBuilder.h`: `IIfcModelBuilder`, `IIfcExportExtension`, `IfcExport::BuildOptions` (what an extension sets: the mapping table, rebar, quantities, and display units; everything else is the IFC export's default: the alignment and bridge, classified), `IfcExport::IHooks`, `IServices`, `IRegistry`, `RoleSpec`, `TargetSpec`, `HookContext`, `ExportedElements`, `IfcBuildResult`, and `BuildFingerprint`.

`IfcExtensionSdk.props`, at the repository root (where `PGSuper.props` and `WBFL.props` are in theirs), puts `Include` on the include path, so the headers are `<IfcExport/...>`, and has the IfcOpenShell settings both projects must share (see "ABI"). Both projects import it.

**Small on purpose** (decided 2026-11-20, before version 1 was released): the IFC extension's own export context (`ExportContext`, `TargetDef` in `IfcTargets.h`) and options (`CIfcExportOptions` in `IfcExportOptions.h`) aren't part of it, so they can change without changing the SDK. An extension target's getter (`IfcExport::Getter`) reads `IfcExport::Element`, which the IFC extension implements over its context (`CElementView`, `IfcExtensionRegistry.cpp`); a new key is a new accessor at the end of the interface. The agent turns `BuildOptions` into `CIfcExportOptions`.

## Hooks

In order: `OnExportBegin`, `OnDeck`, `OnBridgeParts`, `OnPier` (each), `OnBearing` (each), `OnBarrier` (each), `OnSegment` (each), `OnSuperstructureComplete`, `OnBridgeComplete`, `OnModelComplete`, `OnExportEnd`. The bridge hooks aren't called when the bridge isn't exported; a girder-only export calls `OnSegment` for its girder. `check_extension.py` checks the order, and that `OnPier`, `OnBearing`, and `OnSegment` are called once for each pier, bearing, and segment in the file.

`HookContext` has the file, the broker, the `BuildOptions`, the services, and `elements`, the elements made so far (project, site, bridge, its parts, the deck slab, piers, barriers, and the girder segments by key). The diaphragms of the standard export (`ElementKind::Diaphragm`: the pier diaphragms, and the intermediate diaphragms, one per bay with `DiaphragmType` "Intermediate", whose targets see `Span()`, `Index()`, `Segment()` (the bay's left girder), and `IntermediateDiaphragm()`) exist by `OnSuperstructureComplete`.

An extension's own elements get their classification and property sets through `IServices::Classify` and `IServices::WritePropertySets(role, objects, keys)`: the `ElementKeys` are what the role's targets read; the IFC extension adds the broker. `IServices::RebarMaterial` gives the `IfcMaterial` the export uses for a bar type and grade (with the material property sets the table declares for the `rebar` role), so an extension's elements can be associated with the same material as the standard export's bars.

An export with an extension can't be a round trip (`CIfcExportOptions::source_model`): `CExportProfile::Apply` refuses it with a `CIfcExtensionException`.

A hook that throws stops the export. The exception becomes a `CIfcExtensionException` that names the extension and the hook, and no file is written (the file is written only at the end). What the hooks logged before is still in `IfcBuildResult::log`.

## Visibility

`IfcExtensionRegistry.h`: an extension's content is registered for each export and owned by it (`CIfcExportSession`, or `WriteTableIds`), declared before the mapping table so it outlives the table's `const TargetDef*` pointers. Its getters are code in the extension's DLL, which is loaded throughout, because the extension is the caller. Role ids are `ElementKind::FirstExtension` (1000) and up, for each export.

`FindTargetDef`, `GetElementKind`, and `GetElementRoleName` find an extension's targets and roles only while a `CIfcExtensionScope` for it exists, that is, during `IIfcModelBuilder` calls that name the extension. The standard export, the import, `/IfcTableToIds`, and the mapping table editor never see them, so a table that `"requires"` an extension can't be used there. `GetTargetDefs()` lists only this project's targets.

Rules, checked when the extension registers:

- names start with the extension's prefix and a period (`targetApp.girder.k1`, `targetApp.diaphragm`)
- no duplicates
- every target has a getter (export only), and its role exists
- a table can't give an extension target import locations (`"targets"`)

## ABI

The extension and this project share C++ objects: the IfcOpenShell file and entity handles, `std::function`, `std::string`. Both must be built with the same compiler, run time library, IfcOpenShell build, and defines: `IfcExtensionSdk.props` has the IfcOpenShell defines (`HAS_SCHEMA_4x3_add2`, `IFOPSH_SAFE_INSTANCE`, `IFOPSH_WITH_ROCKSDB`, `BOOST_ALL_NO_LIB`, as IfcOpenShell's own build uses them) and the parse libraries, so the two can't differ. Before an extension is used, the IFC extension compares:

- `IFC_EXPORT_EXTENSION_API_VERSION`: bump it when anything in `Include\IfcExport` changes. The first released version is 1; it's provisional until then, and changes made before that release don't bump it.
- `IfcExport::BuildFingerprint()`: compiler version, iterator debug level, debug or release, the IfcOpenShell defines, and the sizes of the shared types (`File`, `ExportValue`, `ElementKeys`, `TargetSpec`, `BuildOptions`, `ExportedElements`, `HookContext`, `IfcBuildResult`), computed in each module.

A mismatch gives a message that says to install versions released together. The IfcOpenShell version isn't in the fingerprint: IfcOpenShell gives it at run time (`IFCOPENSHELL_VERSION`, from its DLL), so both modules would see the same one.

A change to the SDK breaks the extensions built against it until they're changed too: build them after this project.

## Testing the API

`TestExportExtension.h` is a test-only extension in this DLL, used with `/IfcTestExtension[=hooks|throw|api|fingerprint]` on the `/IfcExport` and `/IfcTableToIds` command lines and `Tests\ImportValidation\mappings\TestExtension.ifcmap` (which `"requires"` it). It goes through the SDK as another DLL's extension does. `check_extension.py` (the ExtensionApi row of `run_validation.py`, details in `results/ExtensionApi.md`) checks:

- the hooks: their order and counts; a property set on the extension's own element (`ifctest.marker`, an `IfcBuildingElementProxy`) and on the IFC extension's girders and piers; the rest of the export is the standard export (the `PGSuper-System` export baseline); the file is valid
- the refusals, each with its message and without a file: a hook that throws, another API version, another build fingerprint, and the table in a standard export
- the general IDS of the table with the extension's targets (`IIfcModelBuilder::WriteTableIds`)
- bSDD URIs in two dictionaries (the test table adds a `test` dictionary to the standard table's `usBridge`), and that a table naming an undeclared dictionary is refused

## Value kinds

`Time`, `LinearForce`, `MassDensity`, `Real`, and `Moment` were added for extension targets. They are generic physical quantity kinds. `AddUsedProjectUnits` declares SI project units for them: s, and the derived units N/m, kg/m³, and N·m (IFC 4.3 writes moments as `IfcTorqueMeasure`). Mapping tables can also use `IfcModulusOfElasticityMeasure` (N/m²) as a property type, e.g. for `Pset_MaterialMechanical.YoungModulus`.

## Not done yet

- The design-value IDS (`IdsExporter`) doesn't know about extension targets. Use `IIfcModelBuilder::WriteTableIds` for a general IDS of the extension's table.
- PGSplice: the builder is registered by the agent in both, but the extension hooks have only been designed with PGSuper bridges in mind.
- General PGSuper data that extensions register as their own targets (e.g. Ec, densities, stress limits, distribution factors) belongs in this project's target catalog, where every table and the import can use it. Elements extensions build directly with IfcOpenShell (e.g. an overlay, user loads, stirrup zones) are candidates for the standard export or a declarative `IServices::CreateElement`. Once the hooks no longer need IfcOpenShell types, an extension no longer has to match the IfcOpenShell build.
