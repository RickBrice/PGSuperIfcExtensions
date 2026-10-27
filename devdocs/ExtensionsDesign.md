# Export extensions

Drafted 2026-10-02, proposed for review 2026-10-05. Decision G7 in [MappingTablesDesign.md](MappingTablesDesign.md). Its first user is BridgeSight's BrDr IFC exporter (File > Export > Bridge Model to BrDr IFC Model), which writes the PGSuper to AASHTOWare BrDR exchange file.

## Why

The IFC export writes standard and usBridge property sets only. Other applications need more: the BridgeSight BrDR exchange, for example, needs about 35 `pgsPset_*` property sets and elements the standard export doesn't make (diaphragms, stirrup zones, deck reinforcement). Putting them in this project would mix application-specific data into a model that is meant to be "pure".

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

In order: `OnExportBegin`, `OnDeck`, `OnBridgeParts`, `OnPier` (each), `OnBearing` (each), `OnBarrier` (each), `OnSegment` (each), `OnSuperstructureComplete`, `OnBridgeComplete`, `OnModelComplete`, `OnExportEnd`. The bridge hooks aren't called when the bridge isn't exported. `HookContext.elements` has the elements made so far.

A hook that throws stops the export. The exception becomes a `CIfcExtensionException` that names the extension and the hook, and no file is written (the file is written only at the end).

## Visibility

`IfcExtensionRegistry.h`: an extension's content is registered for each export and owned by it (`CIfcExportSession`, or `WriteTableIds`), declared before the mapping table so it outlives the table's `const TargetDef*` pointers. Its getters are code in the extension's DLL, which is loaded throughout, because the extension is the caller. Role ids are `ElementKind::FirstExtension` (1000) and up, for each export.

`FindTargetDef`, `GetElementKind`, and `GetElementRoleName` find an extension's targets and roles only while a `CIfcExtensionScope` for it exists, that is, during `IIfcModelBuilder` calls that name the extension. The standard export, the import, `/IfcTableToIds`, and the mapping table editor never see them. `GetTargetDefs()` still lists only this project's targets.

Rules, checked when the extension registers:

- names start with the extension's prefix and a period (`brdr.girder.k1`, `brdr.diaphragm`)
- no duplicates
- every target has a getter (export only), and its role exists
- a table can't give an extension target import locations (`"targets"`)

## ABI

The extension and this project share C++ objects: the IfcOpenShell file and entity handles, `std::function`, `std::string`. Both must be built with the same compiler, run time library, IfcOpenShell build, and defines. Both projects import `IfcExtensionSdk.props`. Before an extension is used, the IFC extension compares:

- `IFC_EXPORT_EXTENSION_API_VERSION`: bump it when anything in `Include\IfcExport` changes
- `IfcExport::BuildFingerprint()`: compiler version, iterator debug level, debug or release, and the sizes of shared types, computed in each module

A mismatch gives a message that says to install versions released together.

## Value kinds

`Time`, `LinearForce`, `MassDensity`, `Real`, and `Moment` were added for extension targets. They are generic physical quantity kinds. `AddUsedProjectUnits` declares SI project units for them: s, and the derived units N/m, kg/m³, and N·m (IFC 4.3 writes moments as `IfcTorqueMeasure`).

## Not done yet

- The design-value IDS (`IdsExporter`) doesn't know about extension targets. Use `IIfcModelBuilder::WriteTableIds` for a general IDS of the extension's table.
- PGSplice: the builder is registered by the agent in both, but the extension hooks have only been designed with PGSuper bridges in mind.
