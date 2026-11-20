# Export profiles: one export mechanism for usBridge and BrDR

An evaluation (2026-11-19, after PR #3 was merged), implemented 2026-11-20 (see [Status](#status-2026-11-20)). The question:
should the standard (usBridge) export use the same hook mechanism as BridgeSight's BrDR export, so that both are
specializations of one generic export?

## Short answer

Yes in substance, but the generalization runs the other way. Most of it already exists, and turning usBridge into a
hook-based extension would move little.

- **There is already one mechanism.** `CIfcExporter::BuildModel(broker, options, file, extension)` builds both:
  - **the standard export:** `extension` is null, with `Standard.ifcmap`;
  - **BrDR:** `extension` is BridgeSight's, with `BrDrExchange.ifcmap`, which extends the standard table.

  What makes an export usBridge is mostly data: the mapping table's property sets, constants, classification, and element
  selectors. The export dialog's Mapping Table tab already chooses among such data-only profiles.
- **What's usBridge-specific in C++ is small, and most of it isn't hook-shaped.** Hooks add elements and property sets at
  fixed points of the build. The usBridge code that's left is either data that belongs in the table, or decisions about how
  the base model is shaped, which every profile shares (see the inventory below).
- **So the useful generalization is an *export profile*: a mapping table, optional hooks, and options.**
  - **usBridge:** the built-in profile, with no hooks.
  - **BrDR:** a profile with hooks.
  - **Both:** run the same build, through the same session and hook code.

  This makes the standard export a client of the hook path, which is the main gain. The extension path is exercised by
  every validation run instead of only by BrDR's builds.

## Inventory: usBridge in the C++

| Where | What | Kind | Where it should go |
|---|---|---|---|
| `PropertySets.h`, `Create_usBrPset_ACI_BarShape` | `usBrPset_ACI_BarShape`, with the standard name ("ACI 315-99"), version, bend shape, inside bend radius, and one property per bar dimension | data | the standard table: targets for the bend shape and radius, and a list or per-dimension targets for the dimensions |
| `bSDD.h`, `BSDD_URI` | the usBridge bSDD dictionary URI behind `"uri": "bsdd"` | data | a table header property (the dictionary URI), so another standard's table names its own |
| `IdsExporter.cpp` | `usBrPset_ACI_TendonMaterial` and other usBridge names in the design-value IDS | data | read through the table's locations, as most of the IDS already is |
| `USBridge_Classifications.h` | `Classify_usBridge_ReinforcingBarType`, a placeholder (usBridge has no class for bar types) | data | the table's classification list, which already handles the other roles |
| `Rebar.h`, `IfcExporter.cpp` (strands) | strand and bar data on their types, with types by path and by debonding | model shape | the builder's convention (good IFC practice anyway); a profile option only if a standard needs data on occurrences |
| `IfcExporter.cpp` | the spatial structure (bridge parts, containment of diaphragms in the superstructure), `TENDON_BUNDLE` assembly types, the bridge-part rules from the usBridge IDM | model shape | the builder's convention |

None of these is an element that a hook would add. The two "model shape" rows are policies, and the profile would carry
them as options when a second standard needs something different, not before.

## Plan

1. **An internal profile object.**
   - `CExportProfile` carries the mapping table, an `IfcExport::IHooks` (a null object for the standard export), and the
     profile's option overrides.
   - The three entry points build one: File > Export (`PGSuperDataExporter`), `/IfcExport`, and `IIfcModelBuilder`.
   - `BuildModel` takes the profile. `CIfcExportSession` always has hooks, so the `extension == nullptr` special cases go
     away.
   - The ABI-facing API doesn't change.
2. **A test profile in the validation.**
   - A test-only extension in this repository, enabled with a command-line flag such as `/IfcTestExtension`, that:
     - implements every hook and records the calls and their order;
     - registers a role and a target, and writes a marker element and property set;
     - can be made to throw from a hook, to report another API version, and to use a table that `"requires"` it.
   - `run_validation.py` runs it and checks the results. This also covers the three error paths PR #3 left untested.
3. **Data out of the C++.**
   - The bar shape property set goes into `Standard.ifcmap`.
   - The bSDD dictionary URI becomes a table property.
   - The design-value IDS reads its usBridge names through the table.

   Each change is checked with the export baseline: the export must be unchanged.
4. **Write down the builder's conventions.** List the model-shape decisions as the builder's IFC conventions, in
   `MappingTablesDesign.md`, rather than as usBridge's. Make one a profile option only when a standard needs it.
5. **Keep the extension API small**, as recommended for PR #3: a small options struct, context keys through accessors, and
   API version 1 provisional until then. With the standard export as a client, each change to the API shows up in this
   repository's own validation first, which is the point.

## What it doesn't do

- **It doesn't reduce BridgeSight's coupling.** BrDR still shares C++ objects with this DLL (see the PR #3 review: a
  per-release API version, a pinned standard table, and a contract test). It makes the hook path tested, and the standard
  export one profile among others.
- **It doesn't make usBridge an external plug-in.** Nothing in the inventory needs to live outside this DLL.

## Effort

| Step | Effort |
|---|---|
| 1. The profile object | Small: a refactoring of the session and the entry points |
| 2. The test profile | Moderate, and the most valuable step |
| 3. Data out of the C++ | Small per item. The bar shape set needs a way to write one property per bar dimension from a table |

## Status (2026-11-20)

All five steps are done. The standard export is unchanged: the export baselines and the design-value IDS checks are as before (one IDS change, under step 3).

1. **The profile object:** `CExportProfile` (`ExportProfile.h`): `Standard()` or `Extension(extension)`, its hooks (`CreateHooks`, a do-nothing `IfcExport::IHooks` for the standard export), and its option overrides (`Apply`: an extension's table is reported as the extension's, and an extension's export can't be a round trip). `CIfcExporter::BuildModel` takes a profile, and the three entry points make one. The session always has hooks, so `RecordExportedElement` and `CallExportHook` no longer check for an extension. The hooks' log is kept when an export fails.
2. **The test profile:** `CTestExportExtension` (`TestExportExtension.h`), `/IfcTestExtension[=hooks|throw|api|fingerprint]` on `/IfcExport` and `/IfcTableToIds`, `Tests/ImportValidation/mappings/TestExtension.ifcmap`, and `check_extension.py`, the "Export extension API" row of `run_validation.py`: the hooks' order and counts, property sets on the extension's own element and on the IFC extension's girders and piers, the rest of the export against the `PGSuper-System` baseline, and the four refusals (a hook that throws, another API version, another build fingerprint, the table without the extension). The general IDS of a table with an extension's targets (`WriteTableIds`) is checked too.
3. **Data out of the C++:**
   - `usBrPset_ACI_BarShape` is in `Standard.ifcmap`, with targets (`rebar.shape_standard`, `rebar.shape_standard_version`, `rebar.bend_shape`, `rebar.bend_radius`, and `rebar.dimension_a` ... `rebar.dimension_v`, one per ACI dimension, absent when the bar's shape doesn't have it). No new table mechanism was needed: usBridge lists each dimension as its own property. ACI 315-99 is the standard of PGSuper's shape codes, so it's a target value, not a table constant. `PropertySets.h` is gone.
   - The bSDD dictionary is the table's `"bsdd_dictionary"` (the standard table's is usBridge's; an extending table uses its base's unless it names its own). `bSDD.h` is gone; the mapping table editor shortens URIs with the edited table's dictionary.
   - The design-value IDS requires every constant the table declares for the girder, the girder type, and the strand material, instead of naming them by property set. It names no usBridge property set any more. The IDS now also requires the girders' `usBrPset_Common.Status` and `Pset_PrecastConcreteElementGeneral.SupportDuringTransportDescription`, and the girder type's `StrengthClass`, which the export already wrote; the strand material's `usBrPset_ACI_TendonMaterial.Specification`, which the table attaches to the tendon type, was never required.
   - `Classify_usBridge_ReinforcingBarType` is gone; `USBridge_Classifications.h` is `Classifications.h`.
4. **The builder's conventions** are in [MappingTablesDesign.md](MappingTablesDesign.md), "The exporter's IFC conventions".
5. **A smaller API** (still version 1, provisional until it's released):
   - Extensions pass `IfcExport::BuildOptions` (the table, rebar, quantities, display units) to `IIfcModelBuilder::BuildModel` and see it in `HookContext::options`, instead of `CIfcExportOptions`, which is no longer in the SDK (`IfcExportOptions.h` is in the project).
   - Extension targets read an element's keys through `IfcExport::Element` (accessors), and hooks give their own elements' keys in `IfcExport::ElementKeys`; `ExportContext` and `TargetDef` are no longer in the SDK.
   - The build fingerprint has the sizes of the types that are still shared.
   - BridgeSight's BrDR exporter was changed to match (the getters' parameter and accessors, its hooks' keys, and its build options) and builds; the change is in `BSDev\BrDrIfcExporter`, for BridgeSight to review.

What's left from "What it doesn't do" is unchanged: BrDR still shares C++ objects with this DLL (IfcOpenShell types in the hooks).
