# Export extensions

Drafted 2026-10-02, proposed for review 2026-10-05. Decision G7 in [MappingTablesDesign.md](MappingTablesDesign.md). Its intended users are plug-ins with their own export commands that write an IFC file for another application (an exchange file).

## Why

The IFC export writes standard and usBridge property sets only. Other applications need more: an exchange with an analysis program, for example, may need its own property sets and elements the standard export doesn't make (diaphragms, stirrup zones, deck reinforcement). Putting them in this project would mix application-specific data into a model that is meant to be "pure".

An export extension lets another plug-in build the model exactly as the IFC export does and add to it, without this project knowing what is added. The standard export command never loads an extension.

## How it works (template method)

`CIfcExporter::BuildModel` builds the model in a fixed order. At fixed points it calls the extension's hooks. The extension supplies three things:

| What | Where | How the IFC extension uses it |
|---|---|---|
| A mapping table | its own JSON file, `"extends": "standard"`, `"requires": ["<id>"]` | loaded like any table, so its property sets are written by the same code as the standard ones |
| Targets and element roles | `IIfcExportExtension::RegisterContent` | the table binds properties to the targets. Targets can belong to this project's roles (e.g. `girder`) or to the extension's own roles |
| Hooks | `IIfcExportExtension::CreateHooks` | called during the build, to add elements and write their property sets through `IServices` |

Another plug-in reaches the builder through the broker: `IIfcModelBuilder`, registered by `CIfcExtensionAgent`. It's declared in the header-only SDK, `Include\IfcExport\`:

- `IfcExtensionTypes.h`: the target types (`ElementKind`, `ValueKind`, `ExportContext`, `TargetDef`, ...), moved out of `IfcTargets.h`
- `IfcExportOptions.h`: `CIfcExportOptions`, moved out of `IfcExporter.h`
- `IfcModelBuilder.h`: `IIfcModelBuilder`, `IIfcExportExtension`, `IfcExport::IHooks`, `IServices`, `IRegistry`, `HookContext`, `ExportedElements`
- `IfcExtensionSdk.props`: the IfcOpenShell settings both projects must share

## Hooks

In order: `OnExportBegin`, `OnDeck`, `OnBridgeParts`, `OnPier` (each), `OnBearing` (each), `OnBarrier` (each), `OnSegment` (each), `OnSuperstructureComplete`, `OnBridgeComplete`, `OnModelComplete`, `OnExportEnd`. The bridge hooks aren't called when the bridge isn't exported. `HookContext.elements` has the elements made so far. The pier diaphragms of the standard export (`ElementKind::Diaphragm`) exist by `OnSuperstructureComplete`.

An export with an extension can't be a round trip (`CIfcExportOptions::source_model`): `BuildModel` refuses it with a `CIfcExtensionException`.

A hook that throws stops the export. The exception becomes a `CIfcExtensionException` that names the extension and the hook, and no file is written (the file is written only at the end).

## Visibility

`IfcExtensionRegistry.h`: an extension's content is registered for each export and owned by it (`CIfcExportSession`, or `WriteTableIds`), declared before the mapping table so it outlives the table's `const TargetDef*` pointers. Its getters are code in the extension's DLL, which is loaded throughout, because the extension is the caller. Role ids are `ElementKind::FirstExtension` (1000) and up, for each export.

`FindTargetDef`, `GetElementKind`, and `GetElementRoleName` find an extension's targets and roles only while a `CIfcExtensionScope` for it exists, that is, during `IIfcModelBuilder` calls that name the extension. The standard export, the import, `/IfcTableToIds`, and the mapping table editor never see them. `GetTargetDefs()` still lists only this project's targets.

Rules, checked when the extension registers:

- names start with the extension's prefix and a period (`targetApp.girder.k1`, `targetApp.diaphragm`)
- no duplicates
- every target has a getter (export only), and its role exists
- a table can't give an extension target import locations (`"targets"`)

## ABI

The extension and this project share C++ objects: the IfcOpenShell file and entity handles, `std::function`, `std::string`. Both must be built with the same compiler, run time library, IfcOpenShell build, and defines. Both projects import `IfcExtensionSdk.props`. Before an extension is used, the IFC extension compares:

- `IFC_EXPORT_EXTENSION_API_VERSION`: bump it when anything in `Include\IfcExport` changes. The first released version is 1; changes made before that release don't bump it
- `IfcExport::BuildFingerprint()`: compiler version, iterator debug level, debug or release, the IfcOpenShell defines (`IFOPSH_SAFE_INSTANCE`, `IFOPSH_WITH_ROCKSDB`), and the sizes of shared types, computed in each module

A mismatch gives a message that says to install versions released together.

`IfcExtensionSdk.props` has the IfcOpenShell defines (`HAS_SCHEMA_4x3_add2`, `IFOPSH_SAFE_INSTANCE`, `IFOPSH_WITH_ROCKSDB`, `BOOST_ALL_NO_LIB`, as IfcOpenShell's own build uses them) and the parse libraries, so the IFC extension and its extensions can't differ.

## Value kinds

`Time`, `LinearForce`, `MassDensity`, `Real`, and `Moment` were added for extension targets. They are generic physical quantity kinds. `AddUsedProjectUnits` declares SI project units for them: s, and the derived units N/m, kg/m³, and N·m (IFC 4.3 writes moments as `IfcTorqueMeasure`). Mapping tables can also use `IfcModulusOfElasticityMeasure` (N/m²) as a property type, e.g. for `Pset_MaterialMechanical.YoungModulus`.

`IServices::RebarMaterial` gives an extension the `IfcMaterial` the export uses for a bar type and grade (with the material property sets the table declares for the `rebar` role), so the extension's own elements can be associated with the same material as the standard export's bars.

## Not done yet

- The design-value IDS (`IdsExporter`) doesn't know about extension targets. Use `IIfcModelBuilder::WriteTableIds` for a general IDS of the extension's table.
- PGSplice: the builder is registered by the agent in both, but the extension hooks have only been designed with PGSuper bridges in mind.
